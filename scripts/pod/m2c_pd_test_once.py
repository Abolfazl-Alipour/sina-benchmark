#!/usr/bin/env python
"""m2c_pd_test_once — Task A: score the hidden PD test set ONCE per checkpoint.

Loss-best vs rank-best pair for the PD Arm E run (mondo__0005180):
  - epoch_001.pth (epoch 1 — best val MRR 0.0195)
  - best.pth as it stands after epoch 4 (epoch 3 weights, best val loss 0.000658;
    epoch 4 did not improve so best.pth was not overwritten — tensor-identical to
    epoch_003.pth, verified in reports/m2c/e2_pd_test_once.md)

Filtered evaluator over the pinned 17,448 genes, runner semantics: candidates =
pinned minus (train2 + valid positives of the unit head); rank = optimistic-tie bisect
(identical to m2c_arm_e_train.py validation and to Probe 1, which reproduced the
trainer's logged val MRR/H@10 exactly). All test queries of a unit share one disease
head + relation DIS-GEN (339 PD / 621 AD), so ONE forward per checkpoint over the
pinned tails scores every query exactly.

Usage (bpn120 env; CPU is the validated path — zero GPU needed):
    CUDA_VISIBLE_DEVICES="" ~/miniconda3/envs/bpn120/bin/python -u \
        scripts/pod/m2c_pd_test_once.py --out /tmp/pd_test_once.json

    no-ANA variant (KEEP-graph exports + noana checkpoints):
    CUDA_VISIBLE_DEVICES="" ~/miniconda3/envs/bpn120/bin/python -u \
        scripts/pod/m2c_pd_test_once.py \
        --exports-dir data/processed/m2c/completion_exports \
        --ckpt-root reports/m2c/checkpoints_completion_noana \
        --checkpoints epoch_001.pth --out /tmp/pd_test_once_noana.json
"""
import argparse
import bisect
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "third_party" / "biopathnet"))
sys.path.insert(0, str(REPO_ROOT / "src"))

EXPORTS_DIR = REPO_ROOT / "data/processed/m2c/completion_exports_ana"
PINNED_PATH = REPO_ROOT / "data/processed/m2c/pinned_genes_17448.txt"
CKPT_ROOT = REPO_ROOT / "reports/m2c/checkpoints_completion_ana"


def utcnow_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bootstrap_ci(values: np.ndarray, n_boot: int, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(values)
    reps = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        reps[b] = float(np.mean(values[idx]))
    lo, hi = np.percentile(reps, [2.5, 97.5])
    return float(lo), float(hi)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit", default="mondo__0005180")
    parser.add_argument("--config", type=Path, default=REPO_ROOT / "configs/m2c_arm_e.yaml")
    parser.add_argument("--checkpoints", nargs="+", default=["epoch_001.pth", "best.pth"])
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1337)
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--out", type=Path, default=Path("/tmp/pd_test_once.json"))
    parser.add_argument("--exports-dir", type=Path, default=EXPORTS_DIR,
                        help="export root containing <unit>/ (default: ANA exports)")
    parser.add_argument("--ckpt-root", type=Path, default=CKPT_ROOT,
                        help="checkpoint root containing <unit>/<name> (default: ANA wave)")
    args = parser.parse_args()

    import torch
    from biopathnet import dataset as bpn_dataset  # noqa: F401
    from biopathnet import model as bpn_model  # noqa: F401
    from biopathnet import task as bpn_task  # noqa: F401
    from biopathnet import util
    from torchdrug import core

    torch.set_num_threads(args.threads)
    print(f"[{utcnow_iso()}] pd-test-once start unit={args.unit} device=cpu", flush=True)

    pinned = PINNED_PATH.read_text(encoding="utf-8").strip().split()
    assert len(pinned) == 17448

    unit_dir = args.exports_dir / args.unit
    cfg = util.load_config(str(args.config), context={"gpus": []})
    cfg.dataset.path = str(unit_dir)
    cfg.engine.gpus = []
    dataset_obj = core.Configurable.load_config_dict(cfg.dataset)
    task_cfg = cfg.task
    task_cfg.model.num_relation = dataset_obj.num_relation
    task_obj = core.Configurable.load_config_dict(task_cfg)
    train_set, valid_set, test_set = dataset_obj.split()
    task_obj.preprocess(train_set, valid_set, test_set)
    model = task_obj.model
    model.eval()
    fact_graph = task_obj.fact_graph
    print(f"fact_graph edges {fact_graph.num_edge}", flush=True)

    inv_entity = dataset_obj.inv_entity_vocab

    def edges_of(subset):
        return [tuple(dataset_obj.graph.edge_list[i].tolist()) for i in subset.indices]

    train2_edges = edges_of(train_set)
    valid_edges = edges_of(valid_set)
    test_edges = edges_of(test_set)
    rel_ids = {r for _, _, r in train2_edges}
    assert len(rel_ids) == 1
    dis_gen_rel = next(iter(rel_ids))
    disease_idx = inv_entity[args.unit.replace("__", ":")]

    test_queries = [(h, t) for h, t, r in test_edges if r == dis_gen_rel and h == disease_idx]
    assert len(test_queries) == len(test_edges)  # 339 PD / 621 AD; all share the unit head
    train_pos = {t for h, t, r in train2_edges if r == dis_gen_rel and h == disease_idx}
    valid_pos = {t for h, t, r in valid_edges if r == dis_gen_rel and h == disease_idx}
    excluded = train_pos | valid_pos  # runner filtered-evaluator semantics
    print(f"test queries {len(test_queries)}; filtered-out positives {len(excluded)}",
          flush=True)

    pinned_indices = [inv_entity[g] for g in pinned]
    pos_of_gene = {g: i for i, g in enumerate(pinned_indices)}
    gw_filt_pos = [i for i, g in enumerate(pinned_indices) if g not in excluded]
    print(f"filtered candidate universe {len(gw_filt_pos)} of {len(pinned_indices)}",
          flush=True)

    num_pinned = len(pinned_indices)
    h_index = torch.full((1, num_pinned), disease_idx, dtype=torch.long)
    r_index = torch.full((1, num_pinned), dis_gen_rel, dtype=torch.long)
    t_index = torch.tensor(pinned_indices, dtype=torch.long).unsqueeze(0)

    results = {}
    for name in args.checkpoints:
        path = args.ckpt_root / args.unit / name
        ckpt = torch.load(path, map_location="cpu")
        model.load_state_dict(ckpt["model"])
        model.eval()
        t0 = time.perf_counter()
        with torch.no_grad():
            scores_t = model(fact_graph, h_index, t_index, r_index,
                             all_loss=None, metric=None)
        fwd_sec = time.perf_counter() - t0
        scores = scores_t[0].detach().cpu().numpy().astype(np.float64)
        del scores_t

        gw_sorted = sorted((scores[i] for i in gw_filt_pos), reverse=True)
        neg_sorted = [-s for s in gw_sorted]
        ranks, rr, h10 = [], [], []
        for _, t in test_queries:
            s_t = float(scores[pos_of_gene[t]])
            rank = bisect.bisect_left(neg_sorted, -s_t) + 1
            ranks.append(rank)
            rr.append(1.0 / rank)
            h10.append(1.0 if rank <= 10 else 0.0)
        rr, h10 = np.array(rr), np.array(h10)
        mrr_lo, mrr_hi = bootstrap_ci(rr, args.bootstrap, args.seed)
        h10_lo, h10_hi = bootstrap_ci(h10, args.bootstrap, args.seed)
        results[name] = {
            "epoch": ckpt.get("epoch"),
            "sha256": sha256_file(path),
            "forward_seconds": fwd_sec,
            "test_mrr": float(np.mean(rr)),
            "test_mrr_ci": [mrr_lo, mrr_hi],
            "test_hits10": float(np.mean(h10)),
            "test_hits10_ci": [h10_lo, h10_hi],
            "median_rank": float(np.median(ranks)),
            "n_queries": len(test_queries),
        }
        print(f"  {name} (epoch {ckpt.get('epoch')}): test MRR {np.mean(rr):.5f} "
              f"[{mrr_lo:.5f}, {mrr_hi:.5f}] H@10 {np.mean(h10):.5f} "
              f"[{h10_lo:.5f}, {h10_hi:.5f}] median rank {np.median(ranks):.0f} "
              f"({fwd_sec:.1f}s)", flush=True)

    out = {
        "finished": utcnow_iso(),
        "unit": args.unit,
        "n_test": len(test_queries),
        "n_filtered_candidates": len(gw_filt_pos),
        "bootstrap": args.bootstrap,
        "seed": args.seed,
        "checkpoints": results,
    }
    args.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
