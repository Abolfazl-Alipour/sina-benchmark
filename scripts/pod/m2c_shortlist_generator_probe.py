#!/usr/bin/env python
"""m2c_shortlist_generator_probe — Probe 1.5: shortlist-GENERATOR comparison (PD, zero training).

Probe 1 (reports/m2c/rerank_probe_pd.md) showed the rerank bottleneck is candidate
GENERATION: degree-top-100 shortlists contain held-out PD test positives only 2/339
of the time. This probe asks whether better generators (disease-personalized PPR)
produce shortlists that actually contain held-out genes, and whether the current
checkpoints still rank them top-10 when present.

Two subcommands (two envs, JSON handoff via /tmp):

  generate  (sina env, CPU: PYTHONPATH=src ~/miniconda3/envs/sina/bin/python ...)
      Generators (all restricted to pinned 17,448; membership EXCLUDES known
      train2(PD-head)+valid positives — the realistic not-yet-associated panel):
        - degree_train1 : undirected train1 degree (Probe 1 baseline)
        - ppr_train1    : PPR personalized on mondo:0005180, train1 graph
        - ppr_train1train2 : same, graph = train1+train2 (matches C1.4 baseline)
      PPR: sina.baselines.ppr_baseline, alpha 0.5, max_iter 100, tol 1e-6 (C1.4 params).
      K in {100, 200, 500}. Coverage over the 339 test positives (held out) and the
      967 valid positives (reported under train2-only exclusion, since valid genes are
      excluded by construction under the primary rule).

  score  (bpn120 env, CPU: CUDA_VISIBLE_DEVICES="" ...)
      Best-coverage generator at K=100 (+ degree reference): one NBFNet forward per
      checkpoint (best.pth, epoch_002.pth) over the pinned 17,448 tails — all test
      queries share (head, relation), so this scores every query exactly (Probe 1
      machinery). Two-stage recovery = P(test gene in shortlist AND model rank <= 10);
      rerank null = generator's own score ordering on the SAME shortlist. Bootstrap
      1000 resamples, seed 1337, over the 339 test queries.
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
DISEASE_TOKEN = "mondo:0005180"
GENERATORS = ("degree_train1", "ppr_train1", "ppr_train1train2")
KS = (100, 200, 500)


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


def read_tsv(path: Path) -> list[tuple[str, str, str]]:
    rows = []
    with path.open(encoding="utf-8") as f:
        for line in f:
            h, r, t = line.rstrip("\n").split("\t")
            rows.append((h, r, t))
    return rows


def top_k(score: dict[str, float], eligible: set[str], k: int) -> list[str]:
    """Top-k eligible tokens by (-score, token); missing scores count as 0.0."""
    return sorted(eligible, key=lambda g: (-score.get(g, 0.0), g))[:k]


def bootstrap_ci(values: np.ndarray, n_boot: int, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(values)
    reps = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        reps[b] = float(np.mean(values[idx]))
    lo, hi = np.percentile(reps, [2.5, 97.5])
    return float(lo), float(hi)


def cmd_generate(args) -> int:
    import polars as pl

    from sina.baselines import degree_baseline, ppr_baseline

    unit_dir = EXPORTS_DIR / args.unit
    pinned = load_pinned()
    pinned_set = set(pinned)
    train1 = read_tsv(unit_dir / "train1.txt")
    train2 = read_tsv(unit_dir / "train2.txt")
    valid = read_tsv(unit_dir / "valid.txt")
    test = read_tsv(unit_dir / "test.txt")
    print(f"edges: train1 {len(train1)} train2 {len(train2)} valid {len(valid)} "
          f"test {len(test)}", flush=True)

    test_genes = [t for _, _, t in test]
    valid_genes = [t for _, _, t in valid]
    assert all(h == DISEASE_TOKEN for h, _, _ in valid + test)
    train2_pos = {t for h, _, t in train2 if h == DISEASE_TOKEN}
    known = train2_pos | set(valid_genes)  # primary membership exclusion
    print(f"known positives: train2(PD) {len(train2_pos)} + valid {len(valid_genes)} "
          f"-> excluded {len(known)}", flush=True)

    train1_df = pl.DataFrame(
        {"src_id": [h for h, _, _ in train1], "relation": [r for _, r, _ in train1],
         "dst_id": [t for _, _, t in train1]})
    train2_df = pl.DataFrame(
        {"src_id": [h for h, _, _ in train2], "relation": [r for _, r, _ in train2],
         "dst_id": [t for _, _, t in train2]})
    train12_df = pl.concat([train1_df, train2_df])

    scores: dict[str, dict[str, float]] = {}
    t0 = time.perf_counter()
    deg_df = degree_baseline(train1_df, pinned)
    scores["degree_train1"] = dict(zip(deg_df["node_id"].to_list(),
                                       deg_df["degree"].to_list(), strict=True))
    print(f"degree_train1 done ({time.perf_counter() - t0:.1f}s)", flush=True)
    for name, df in (("ppr_train1", train1_df), ("ppr_train1train2", train12_df)):
        t0 = time.perf_counter()
        ppr_df = ppr_baseline(df, [DISEASE_TOKEN], pinned, alpha=args.ppr_alpha,
                              max_iter=args.ppr_max_iter, tol=args.ppr_tol)
        scores[name] = dict(zip(ppr_df["candidate_id"].to_list(),
                                ppr_df["ppr_score"].to_list(), strict=True))
        n_edges = df.height
        print(f"{name} done on {n_edges} edges ({time.perf_counter() - t0:.1f}s)",
              flush=True)

    eligible = pinned_set - known            # primary: exclude train2+valid positives
    eligible_v = pinned_set - train2_pos     # supplementary: valid treated as held-out
    out_shortlists: dict[str, dict[str, list[str]]] = {g: {} for g in GENERATORS}
    coverage: dict[str, dict[str, dict[str, float]]] = {g: {} for g in GENERATORS}
    for g in GENERATORS:
        for k in KS:
            sl = top_k(scores[g], eligible, k)
            sl_v = top_k(scores[g], eligible_v, k)
            out_shortlists[g][str(k)] = sl
            coverage[g][str(k)] = {
                "test": float(np.mean([t in set(sl) for t in test_genes])),
                "valid_train2only_excl": float(np.mean([t in set(sl_v) for t in valid_genes])),
            }
            print(f"{g} K={k}: test coverage {coverage[g][str(k)]['test']:.4f} "
                  f"valid(t2-only excl) "
                  f"{coverage[g][str(k)]['valid_train2only_excl']:.4f}", flush=True)

    result = {
        "started": utcnow_iso(),
        "unit": args.unit,
        "ppr": {"alpha": args.ppr_alpha, "max_iter": args.ppr_max_iter, "tol": args.ppr_tol},
        "input_hashes": {n: sha256_file(unit_dir / n)
                         for n in ("train1.txt", "train2.txt", "valid.txt", "test.txt")}
                        | {"pinned_genes_17448.txt": sha256_file(PINNED_PATH)},
        "n_known_excluded": len(known),
        "shortlists": out_shortlists,
        "shortlist_sha256": {g: {k: sha256_text("\n".join(sl))
                                 for k, sl in per_k.items()}
                             for g, per_k in out_shortlists.items()},
        "coverage": coverage,
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.out}", flush=True)
    return 0


def cmd_score(args) -> int:
    import torch

    if args.out.resolve() == args.shortlists.resolve():
        raise ValueError("--out must differ from --shortlists in score mode")

    from biopathnet import dataset as bpn_dataset  # noqa: F401
    from biopathnet import model as bpn_model  # noqa: F401
    from biopathnet import task as bpn_task  # noqa: F401
    from biopathnet import util
    from torchdrug import core

    gen_blob = json.loads(args.shortlists.read_text())
    unit_dir = EXPORTS_DIR / args.unit
    pinned = load_pinned()

    torch.set_num_threads(args.threads)

    # --- dataset & task setup: identical machinery to m2c_arm_e_train.py ---
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
    disease_idx = inv_entity[DISEASE_TOKEN]
    valid_queries = [(h, t) for h, t, r in valid_edges if r == dis_gen_rel and h == disease_idx]
    test_queries = [(h, t) for h, t, r in test_edges if r == dis_gen_rel and h == disease_idx]
    assert len(valid_queries) == len(valid_edges) and len(test_queries) == len(test_edges)

    pinned_indices = [inv_entity[g] for g in pinned]
    pos_of_gene = {g: i for i, g in enumerate(pinned_indices)}

    # Shortlists under evaluation: best-coverage generator + degree reference, K=100.
    cov = gen_blob["coverage"]
    # tiebreak prefers later-listed generator (ppr_train1train2 > ppr_train1 > degree)
    best_gen = max(GENERATORS, key=lambda g: (cov[g]["100"]["test"], GENERATORS.index(g)))
    # Always score: best-coverage generator, degree reference, and the PPR
    # representative (ppr_train1train2) so model-vs-null exists for a PPR panel too.
    eval_gens = sorted({best_gen, "degree_train1", "ppr_train1train2"})
    print(f"best-coverage generator: {best_gen}; scoring {eval_gens} at K=100", flush=True)

    shortlists = {g: gen_blob["shortlists"][g]["100"] for g in eval_gens}
    sl_pos = {g: [pos_of_gene[inv_entity[tok]] for tok in toks]
              for g, toks in shortlists.items()}  # pinned positions (scores-array indices)
    sl_set = {g: {inv_entity[tok] for tok in toks} for g, toks in shortlists.items()}

    # Generator-score null ordering of each shortlist (ascending key for bisect):
    # null score = generator's own score; we recompute it cheaply from the blob? No —
    # scores are not stored; rerank null = membership order in the shortlist itself,
    # which IS the generator's (-score, token) order. Answer inserts by its own score.
    # We therefore store the generator score per shortlist gene from generate-time order:
    # position in shortlist = generator rank for members; for a member answer the null
    # rank is simply its position (deterministic, tie-consistent).
    null_rank_of_tok = {g: {tok: i + 1 for i, tok in enumerate(toks)}
                        for g, toks in shortlists.items()}

    num_pinned = len(pinned_indices)
    h_index = torch.full((1, num_pinned), disease_idx, dtype=torch.long)
    r_index = torch.full((1, num_pinned), dis_gen_rel, dtype=torch.long)
    t_index = torch.tensor(pinned_indices, dtype=torch.long).unsqueeze(0)

    results = {}
    for ckpt_name in args.checkpoints:
        path = CKPT_ROOT / args.unit / ckpt_name
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
        print(f"[{utcnow_iso()}] {ckpt_name}: epoch {ckpt.get('epoch')}, "
              f"forward {fwd_sec:.1f}s", flush=True)

        per_gen = {}
        for g in eval_gens:
            sl_sorted = sorted((scores[i] for i in sl_pos[g]), reverse=True)
            rows = []
            for _, t in test_queries:
                in_sl = t in sl_set[g]
                if in_sl:
                    s_t = float(scores[pos_of_gene[t]])
                    neg = [-s for s in sl_sorted]
                    rank = bisect.bisect_left(neg, -s_t) + 1
                else:
                    rank = None
                tok = dataset_obj.entity_vocab[t]
                nr = null_rank_of_tok[g].get(tok)  # None if not in shortlist
                rows.append({
                    "gene": tok, "in_shortlist": bool(in_sl), "model_rank": rank,
                    "hit": 1.0 if (rank is not None and rank <= args.hits_k) else 0.0,
                    "null_rank": nr,
                    "null_hit": 1.0 if (nr is not None and nr <= args.hits_k) else 0.0,
                })
            hits = np.array([r["hit"] for r in rows])
            null_hits = np.array([r["null_hit"] for r in rows])
            cov_g = float(np.mean([r["in_shortlist"] for r in rows]))
            lo, hi = bootstrap_ci(hits, args.bootstrap, args.seed)
            nlo, nhi = bootstrap_ci(null_hits, args.bootstrap, args.seed + 1)
            covered = [r for r in rows if r["in_shortlist"]]
            per_gen[g] = {
                "n": len(rows),
                "coverage": cov_g,
                "two_stage_recall": float(np.mean(hits)),
                "two_stage_ci": [lo, hi],
                "recall_given_coverage": (float(np.mean([r["hit"] for r in covered]))
                                          if covered else None),
                "null_two_stage": float(np.mean(null_hits)),
                "null_ci": [nlo, nhi],
                "rows": rows,
            }
            print(f"  {ckpt_name} [{g}] coverage {cov_g:.4f} two-stage "
                  f"{per_gen[g]['two_stage_recall']:.4f} CI [{lo:.4f}, {hi:.4f}] "
                  f"| null {per_gen[g]['null_two_stage']:.4f} CI [{nlo:.4f}, {nhi:.4f}] "
                  f"| recall|cov {per_gen[g]['recall_given_coverage']}", flush=True)
        results[ckpt_name] = {"epoch": ckpt.get("epoch"), "sha256": sha256_file(path),
                              "forward_seconds": fwd_sec, "generators": per_gen}

    out = {
        "finished": utcnow_iso(),
        "unit": args.unit,
        "best_generator": best_gen,
        "eval_generators": eval_gens,
        "k": 100,
        "hits_k": args.hits_k,
        "bootstrap": args.bootstrap,
        "seed": args.seed,
        "checkpoints": results,
    }
    args.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.out}", flush=True)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["generate", "score"])
    parser.add_argument("--unit", default="mondo__0005180")
    parser.add_argument("--config", type=Path, default=REPO_ROOT / "configs/m2c_arm_e.yaml")
    parser.add_argument("--hits-k", type=int, default=10)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1337)
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--ppr-alpha", type=float, default=0.5)
    parser.add_argument("--ppr-max-iter", type=int, default=100)
    parser.add_argument("--ppr-tol", type=float, default=1e-6)
    parser.add_argument("--shortlists", type=Path, default=Path("/tmp/shortlist_gen.json"))
    parser.add_argument("--checkpoints", nargs="+", default=["best.pth", "epoch_002.pth"])
    parser.add_argument("--out", type=Path, default=Path("/tmp/shortlist_gen.json"))
    args = parser.parse_args()
    if args.command == "generate":
        return cmd_generate(args)
    return cmd_score(args)


if __name__ == "__main__":
    raise SystemExit(main())
