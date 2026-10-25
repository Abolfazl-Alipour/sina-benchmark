#!/usr/bin/env python
"""m2c_rerank_probe — 100-candidate rerank diagnostic on PD checkpoints (zero training).

Question: does the current Arm E PD model (mondo__0005180) already carry a
product-shaped signal in a realistic 100-candidate rerank setting?

Design notes:
- All valid/test queries share one head (mondo:0005180) and one relation
  (DIS-GEN = ASSOCIATED_WITH). NBFNet scoring with all_loss=None depends only on
  (h, r, graph, params); the tail index is just a gather. So ONE forward per
  checkpoint over the pinned 17,448 tails yields every query's scores, both
  genome-wide and inside any shortlist.
- Shortlists (both deterministic, top-degree pinned genes from train1.txt,
  tiebreak gene token):
    - "primary" (pre-registered construction): top-100 pinned genes by
      undirected train1 degree, known train2+valid positives INCLUDED if they
      rank there; they are masked out of the ranking exactly like the runner's
      filtered evaluator (excluded = train2 + valid positives of the query
      disease, minus the query answer) — filtered non-hits, never countable.
    - "excl_known" (post-hoc sensitivity): top-100 pinned genes by train1
      degree among genes NOT already train2/valid positives of the disease —
      the realistic "100 new candidates" panel. Labeled sensitivity only.
- A query positive outside the shortlist is an automatic miss (intended
  realism); coverage is reported for both shortlists.
- Null: degree ranking of the same shortlist through the same filtered
  evaluator. Bootstrap CIs: 1000 resamples, seed 1337, over queries.
- Calibration: pooled logistic fit of hit indicator on score over
  queries x filtered shortlist candidates (Newton-Raphson, numpy only).

Usage (CPU, bpn120 env — preferred; never touches the GPU):
    CUDA_VISIBLE_DEVICES="" ~/miniconda3/envs/bpn120/bin/python -u \
        scripts/pod/m2c_rerank_probe.py --smoke 3 --out /tmp/rerank_probe_smoke.json

Outputs: JSON results blob (default /tmp/rerank_probe_pd.json); the committed
artifact is the hand-written reports/m2c/rerank_probe_pd.md.
"""
import argparse
import bisect
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "third_party" / "biopathnet"))
sys.path.insert(0, str(REPO_ROOT / "src"))

from biopathnet import dataset as bpn_dataset  # noqa: E402, F401
from biopathnet import model as bpn_model  # noqa: E402, F401
from biopathnet import task as bpn_task  # noqa: E402, F401
from biopathnet import util  # noqa: E402
from torchdrug import core  # noqa: E402

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


def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def load_pinned() -> list[str]:
    txt = PINNED_PATH.read_text(encoding="utf-8").strip().split()
    assert len(txt) == 17448, f"pinned len {len(txt)} != 17448"
    return txt


def gene_degrees_from_train1(train1_path: Path, pinned_set: set[str]) -> dict[str, int]:
    """Undirected degree of each pinned gene in train1 (count head+tail occurrences)."""
    deg = dict.fromkeys(pinned_set, 0)
    with train1_path.open(encoding="utf-8") as f:
        for line in f:
            h, _, t = line.rstrip("\n").split("\t")
            if h in pinned_set:
                deg[h] += 1
            if t in pinned_set:
                deg[t] += 1
    return deg


def filtered_rank(sorted_desc_scores: list[float], true_score: float) -> int:
    """1-based rank with optimistic ties — identical to the runner's bisect_left scheme."""
    neg_sorted = [-s for s in sorted_desc_scores]
    return bisect.bisect_left(neg_sorted, -true_score) + 1


def logistic_fit(x: np.ndarray, y: np.ndarray, max_iter: int = 100, tol: float = 1e-10):
    """Newton-Raphson logistic regression of binary y on scalar x. Returns (intercept, slope)."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if y.sum() < 2 or y.sum() > len(y) - 2:
        return None, None  # degenerate: (nearly) no hits or no non-hits
    beta = np.zeros(2)
    X = np.column_stack([np.ones_like(x), x])
    for _ in range(max_iter):
        eta = X @ beta
        p = 1.0 / (1.0 + np.exp(-eta))
        w = np.clip(p * (1 - p), 1e-12, None)
        grad = X.T @ (y - p)
        hess = (X.T * w) @ X
        try:
            step = np.linalg.solve(hess, grad)
        except np.linalg.LinAlgError:
            return None, None
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break
    return float(beta[0]), float(beta[1])


def bootstrap_ci(values: np.ndarray, n_boot: int, seed: int) -> tuple[float, float]:
    """Percentile 95% CI of the mean over resampled queries."""
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
    parser.add_argument("--shortlist-size", type=int, default=100)
    parser.add_argument("--hits-k", type=int, default=10)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1337)
    parser.add_argument("--smoke", type=int, default=None,
                        help="limit valid/test queries per checkpoint (forward cost unchanged)")
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu")
    parser.add_argument("--checkpoints", nargs="+",
                        default=["best.pth", "epoch_002.pth", "epoch_003.pth", "last.pth"])
    parser.add_argument("--out", type=Path, default=Path("/tmp/rerank_probe_pd.json"))
    parser.add_argument("--exports-dir", type=Path, default=EXPORTS_DIR,
                        help="export root containing <unit>/ (default: ANA exports)")
    parser.add_argument("--ckpt-root", type=Path, default=CKPT_ROOT,
                        help="checkpoint root containing <unit>/<name> (default: ANA wave)")
    args = parser.parse_args()

    torch.set_num_threads(args.threads)
    if args.device == "cuda":
        assert torch.cuda.is_available(), "cuda requested but not available"
        torch.cuda.set_per_process_memory_fraction(0.02, 0)
        torch.cuda.set_device(0)
        device = torch.device("cuda", 0)
    else:
        device = torch.device("cpu")

    started = utcnow_iso()
    print(f"[{started}] rerank probe start unit={args.unit} device={device}", flush=True)

    pinned_genes = load_pinned()
    pinned_set = set(pinned_genes)

    unit_dir = args.exports_dir / args.unit
    hashes = {
        name: sha256_file(unit_dir / name)
        for name in ("train1.txt", "train2.txt", "valid.txt", "test.txt")
    }
    hashes["pinned_genes_17448.txt"] = sha256_file(PINNED_PATH)

    # --- dataset & task setup: identical machinery to m2c_arm_e_train.py ---
    cfg = util.load_config(str(args.config), context={"gpus": [0] if args.device == "cuda" else []})
    cfg.dataset.path = str(unit_dir)
    cfg.engine.gpus = [0] if args.device == "cuda" else []
    dataset_obj = core.Configurable.load_config_dict(cfg.dataset)
    print(f"dataset loaded: {dataset_obj.num_entity} entities, "
          f"{dataset_obj.num_relation} relations", flush=True)
    task_cfg = cfg.task
    task_cfg.model.num_relation = dataset_obj.num_relation
    task_obj = core.Configurable.load_config_dict(task_cfg)
    train_set, valid_set, test_set = dataset_obj.split()
    print(f"#train2 {len(train_set)} #valid {len(valid_set)} #test {len(test_set)}", flush=True)
    task_obj.preprocess(train_set, valid_set, test_set)
    if args.device == "cuda":
        task_obj = task_obj.cuda(0)
    model = task_obj.model
    model.eval()
    fact_graph = task_obj.fact_graph
    print(f"fact_graph edges {fact_graph.num_edge}", flush=True)

    inv_entity = dataset_obj.inv_entity_vocab
    entity_vocab = dataset_obj.entity_vocab

    def edges_of(subset):
        return [tuple(dataset_obj.graph.edge_list[i].tolist()) for i in subset.indices]

    train2_edges = edges_of(train_set)
    valid_edges = edges_of(valid_set)
    test_edges = edges_of(test_set)

    # DIS-GEN relation id: train2 carries exactly one relation (runner semantics).
    rel_ids = {r for _, _, r in train2_edges}
    assert len(rel_ids) == 1, f"train2 relations not unique: {rel_ids}"
    dis_gen_rel = next(iter(rel_ids))

    disease_token = args.unit.replace("__", ":")
    assert disease_token in inv_entity, f"disease {disease_token} not in vocab"
    disease_idx = inv_entity[disease_token]

    # Queries: all valid/test rows are (disease, DIS-GEN, gene) for this unit.
    valid_queries = [(h, t) for h, t, r in valid_edges if r == dis_gen_rel and h == disease_idx]
    test_queries = [(h, t) for h, t, r in test_edges if r == dis_gen_rel and h == disease_idx]
    assert len(valid_queries) == len(valid_edges), "non-PD valid edges present"
    assert len(test_queries) == len(test_edges), "non-PD test edges present"

    # Filtered-evaluator exclusion set (runner semantics): train2 + valid positives.
    train_pos = {t for h, t, r in train2_edges if r == dis_gen_rel and h == disease_idx}
    valid_pos = {t for _, t in valid_queries}
    excluded_base = train_pos | valid_pos
    print(f"queries: valid {len(valid_queries)} test {len(test_queries)}; "
          f"filtered-out positives train2 {len(train_pos)} valid {len(valid_pos)}", flush=True)

    # Pinned gene indices.
    missing = [g for g in pinned_genes if g not in inv_entity]
    assert not missing, f"{len(missing)} pinned genes missing from vocab e.g. {missing[:5]}"
    pinned_indices = [inv_entity[g] for g in pinned_genes]
    assert len(pinned_indices) == 17448
    pinned_idx_set = set(pinned_indices)
    pos_of_gene = {g: i for i, g in enumerate(pinned_indices)}  # entity idx -> pinned position
    assert train_pos.issubset(pinned_idx_set) and valid_pos.issubset(pinned_idx_set)
    for _, t in valid_queries + test_queries:
        assert t in pinned_idx_set, "query tail outside pinned universe"

    # --- shortlists: top-degree pinned genes from train1 ---
    t0 = time.perf_counter()
    deg = gene_degrees_from_train1(unit_dir / "train1.txt", pinned_set)
    deg_sorted_tokens = sorted(deg, key=lambda g: (-deg[g], g))  # deterministic
    print(f"degrees computed ({time.perf_counter() - t0:.1f}s)", flush=True)

    token_of_idx = {g: pinned_genes[i] for i, g in enumerate(pinned_indices)}
    excluded_tokens = {token_of_idx[g] for g in excluded_base}

    def build_shortlist(tokens: list[str]) -> dict:
        idx = [inv_entity[g] for g in tokens]
        # Filtered candidate list: drop train2+valid positives (they can never be
        # hits and occupy no rank position, exactly like the runner's evaluator).
        filt = [(pos_of_gene[g], g) for g in idx if g not in excluded_base]
        # Degree-null ranking of the same shortlist through the same filter;
        # a query answer ranks by key insertion (runner adds the answer back).
        deg_order = sorted(filt, key=lambda pg: (-deg[token_of_idx[pg[1]]], token_of_idx[pg[1]]))
        deg_keys = [(-deg[token_of_idx[g]], token_of_idx[g]) for _, g in deg_order]

        def degree_rank(t: int):
            if t not in set(idx):
                return None
            tok = token_of_idx[t]
            return bisect.bisect_left(deg_keys, (-deg[tok], tok)) + 1

        return {
            "tokens": tokens,
            "idx_set": set(idx),
            "filt": filt,  # list of (pinned_position, entity_idx)
            "sha256": sha256_text("\n".join(tokens)),
            "n_known_pos_inside": len(tokens) - len(filt),
            "degree_rank": degree_rank,
        }

    shortlists = {
        "primary": build_shortlist(deg_sorted_tokens[: args.shortlist_size]),
        "excl_known": build_shortlist(
            [g for g in deg_sorted_tokens if g not in excluded_tokens][: args.shortlist_size]
        ),
    }
    for name, sl in shortlists.items():
        print(f"shortlist[{name}] k={len(sl['tokens'])} sha256 {sl['sha256']} "
              f"known-positives-inside {sl['n_known_pos_inside']}", flush=True)

    queries_by_split = {"valid": valid_queries, "test": test_queries}
    if args.smoke:
        queries_by_split = {k: v[: args.smoke] for k, v in queries_by_split.items()}
        print(f"SMOKE mode: {len(queries_by_split['valid'])} valid + "
              f"{len(queries_by_split['test'])} test queries", flush=True)

    # --- coverage + degree null (checkpoint-independent) ---
    variants = {}
    for name, sl in shortlists.items():
        entry = {"sha256": sl["sha256"], "n_known_pos_inside": sl["n_known_pos_inside"],
                 "n_scorable": len(sl["filt"]), "coverage": {}, "degree_null": {}}
        for split, queries in queries_by_split.items():
            hits, rrs, covered = [], [], []
            for _, t in queries:
                r = sl["degree_rank"](t)
                covered.append(1.0 if t in sl["idx_set"] else 0.0)
                hits.append(1.0 if (r is not None and r <= args.hits_k) else 0.0)
                rrs.append(1.0 / r if r is not None else 0.0)
            hits, rrs = np.array(hits), np.array(rrs)
            lo, hi = bootstrap_ci(hits, args.bootstrap, args.seed)
            entry["coverage"][split] = float(np.mean(covered))
            entry["degree_null"][split] = {"recall_at_k": float(np.mean(hits)),
                                           "mrr": float(np.mean(rrs)),
                                           "ci_lo": lo, "ci_hi": hi}
        variants[name] = entry
        print(f"variant[{name}] test coverage {entry['coverage']['test']:.4f} "
              f"null recall@{args.hits_k} {entry['degree_null']['test']['recall_at_k']:.4f} "
              f"CI [{entry['degree_null']['test']['ci_lo']:.4f}, "
              f"{entry['degree_null']['test']['ci_hi']:.4f}]", flush=True)

    # Genome-wide filtered candidate positions (shared by all checkpoints).
    gw_filt_pos = [i for i, g in enumerate(pinned_indices) if g not in excluded_base]

    # --- per-checkpoint scoring ---
    num_pinned = len(pinned_indices)
    h_index = torch.full((1, num_pinned), disease_idx, dtype=torch.long, device=device)
    r_index = torch.full((1, num_pinned), dis_gen_rel, dtype=torch.long, device=device)
    t_index = torch.tensor(pinned_indices, dtype=torch.long, device=device).unsqueeze(0)

    results = {}
    for name in args.checkpoints:
        path = args.ckpt_root / args.unit / name
        ckpt = torch.load(path, map_location="cpu")
        model.load_state_dict(ckpt["model"])
        model.eval()
        epoch = ckpt.get("epoch")
        t0 = time.perf_counter()
        with torch.no_grad():
            scores_t = model(fact_graph, h_index, t_index, r_index,
                             all_loss=None, metric=None)
        fwd_sec = time.perf_counter() - t0
        scores = scores_t[0].detach().cpu().numpy().astype(np.float64)
        del scores_t
        assert scores.shape == (num_pinned,)
        print(f"[{utcnow_iso()}] {name}: epoch {epoch}, forward {fwd_sec:.1f}s", flush=True)

        # Genome-wide filtered ranking (shared across variants).
        gw_sorted = sorted((scores[i] for i in gw_filt_pos), reverse=True)
        sl_sorted = {v: sorted((scores[i] for i, _ in shortlists[v]["filt"]), reverse=True)
                     for v in shortlists}

        ckpt_variants = {}
        for vname, sl in shortlists.items():
            per_split = {}
            for split, queries in queries_by_split.items():
                rows = []
                for _, t in queries:
                    s_t = float(scores[pos_of_gene[t]])
                    rank_gw = filtered_rank(gw_sorted, s_t)
                    # rank_sl None = automatic miss: positive outside shortlist
                    rank_sl = (filtered_rank(sl_sorted[vname], s_t)
                               if t in sl["idx_set"] else None)
                    hit = 1.0 if (rank_sl is not None and rank_sl <= args.hits_k) else 0.0
                    rr = 1.0 / rank_sl if rank_sl is not None else 0.0
                    rows.append({"gene": entity_vocab[t], "score": s_t, "rank_gw": rank_gw,
                                 "rank_in100": rank_sl, "hit": hit, "rr": rr,
                                 "shift": (rank_gw - rank_sl) if rank_sl is not None else None})
                hits = np.array([r["hit"] for r in rows])
                rrs = np.array([r["rr"] for r in rows])
                lo, hi = bootstrap_ci(hits, args.bootstrap, args.seed)
                covered = [r for r in rows if r["rank_in100"] is not None]
                per_split[split] = {
                    "n": len(rows),
                    "recall_at_k": float(np.mean(hits)),
                    "recall_ci": [lo, hi],
                    "mrr_in100": float(np.mean(rrs)),
                    "median_rank_gw": float(np.median([r["rank_gw"] for r in rows])),
                    "median_rank_in100": (float(np.median([r["rank_in100"] for r in covered]))
                                          if covered else None),
                    "median_shift": (float(np.median([r["shift"] for r in covered]))
                                     if covered else None),
                    "rows": rows,
                }
            # Calibration: pooled (score, hit) over test queries x filtered candidates.
            if args.smoke:
                b0, b1, n_pairs, n_hits = None, None, 0, 0
            else:
                cand = [(g, scores[i]) for i, g in sl["filt"]]
                xs, ys = [], []
                for _, t in queries_by_split["test"]:
                    for g, s in cand:
                        xs.append(float(s))
                        ys.append(1.0 if g == t else 0.0)
                b0, b1 = logistic_fit(np.array(xs), np.array(ys))
                n_pairs, n_hits = len(xs), int(sum(ys))
            ckpt_variants[vname] = {"splits": per_split,
                                    "calibration_intercept": b0,
                                    "calibration_slope": b1,
                                    "calibration_n_pairs": n_pairs,
                                    "calibration_n_hits": n_hits}
            trec = per_split["test"]
            print(f"  {name} [{vname}] test recall@{args.hits_k}-in-{args.shortlist_size} "
                  f"{trec['recall_at_k']:.4f} CI [{lo:.4f}, {hi:.4f}] "
                  f"mrr {trec['mrr_in100']:.4f} slope {b1}", flush=True)

        results[name] = {"epoch": epoch, "sha256": sha256_file(path),
                         "forward_seconds": fwd_sec, "variants": ckpt_variants}

    out = {
        "started": started,
        "finished": utcnow_iso(),
        "unit": args.unit,
        "device": str(device),
        "smoke": args.smoke,
        "shortlist_size": args.shortlist_size,
        "hits_k": args.hits_k,
        "bootstrap": args.bootstrap,
        "seed": args.seed,
        "input_hashes": hashes,
        "variants": variants,
        "shortlist_tokens": {v: shortlists[v]["tokens"] for v in shortlists},
        "n_filtered_gw_candidates": len(gw_filt_pos),
        "checkpoints": results,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"[{utcnow_iso()}] wrote {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
