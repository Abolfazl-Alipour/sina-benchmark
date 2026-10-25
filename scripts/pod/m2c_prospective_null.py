#!/usr/bin/env python
"""m2c_prospective_null — temporal-split gate on the supervision-only null.

Owner directive 2026-09-06 (temporal track). Random holdout rewards popularity
interpolation; a prospective split (train on pre-cutoff knowledge, test genes
first reported post-cutoff) measures discovery. The null is a pure function of
the labels, so it goes first: if the counting null collapses on future genes
while holding on random holdout, popularity was its fuel. If it HOLDS on
future genes too, prospective discovery from co-occurrence is real and a
retrained model epoch is justified.

First-reported years: GWAS Catalog bulk
(gwas-catalog-download-associations-alt.*.tsv, genome-wide only) via
scripts/pod/m2c_gwas_loci_build.py -> reports/m2c/gwas_loci_<unit>.json.
Cutoff default 2020 (PD direct: 177 pre / 59 post; AD direct: 91 pre /
180 post — see build output).

Semantics per variant (candidate side identical to every completion number:
pinned minus train2/valid positives of d*; post-cutoff test genes stay
candidates — at cutoff time they were unknown):
  FULL       G = all train2 positives of d*, test = post-cutoff hidden genes
             (control: same null as 0.01484/0.02950, restricted test set)
  PRE_ALL    G = train2 positives minus genes first-reported >= cutoff (any
             trait incl. proxies), same prospective test
  PRE_DIRECT G = train2 positives minus genes first-reported >= cutoff on
             direct-trait rows only, same prospective test

Eval: filtered pinned 17,448, optimistic-tie rank, hidden test scored,
bootstrap 1000 / seed 1337 (same convention as m2c_reweighted_null.py).

Usage (sina env, CPU):
    ~/miniconda3/envs/sina/bin/python scripts/pod/m2c_prospective_null.py \
        --unit mondo__0005180 --loci-json reports/m2c/gwas_loci_mondo__0005180.json \
        --out reports/m2c/prospective_null_pd.json
"""

import argparse
import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
PINNED_PATH = REPO_ROOT / "data" / "processed" / "m2c" / "pinned_genes_17448.txt"


def utcnow_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def bootstrap_ci(values: np.ndarray, n_boot: int, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(values)
    means = [float(np.mean(rng.choice(values, size=n, replace=True)))
             for _ in range(n_boot)]
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def load_splits(exports_dir: Path) -> tuple[list, list, list]:
    def rows(fname):
        out = []
        with (exports_dir / fname).open(encoding="utf-8") as f:
            for line in f:
                h, r, t = line.rstrip("\n").split("\t")
                out.append((h, r, t))
        return out
    return rows("train2.txt"), rows("valid.txt"), rows("test.txt")


def null_scores(train2_rows: list, g_set: set[str],
                 drop_genes: set[str] | None = None) -> dict[str, float]:
    drop = drop_genes or set()
    genes_of: dict[str, set[str]] = defaultdict(set)
    for h, _, t in train2_rows:
        if t in drop:
            continue  # cutoff semantics: post-cutoff genes were unknown
        genes_of[h].add(t)
    inter = {d: len(genes & g_set) for d, genes in genes_of.items()}
    scores: dict[str, float] = defaultdict(float)
    for d, genes in genes_of.items():
        w = inter[d]
        if w == 0:
            continue
        for g in genes:
            scores[g] += w
    return scores


def evaluate(scores: dict[str, float], candidates: list[str],
             test_genes: list[str], n_boot: int, seed: int) -> dict:
    cand_scores = np.array([scores.get(g, 0.0) for g in candidates])
    rr = np.zeros(len(test_genes))
    hits10 = np.zeros(len(test_genes))
    for i, g in enumerate(test_genes):
        s = scores.get(g, 0.0)
        rank = 1 + int(np.sum(cand_scores > s))  # optimistic tie
        rr[i] = 1.0 / rank
        hits10[i] = 1.0 if rank <= 10 else 0.0
    mrr_lo, mrr_hi = bootstrap_ci(rr, n_boot, seed)
    h_lo, h_hi = bootstrap_ci(hits10, n_boot, seed + 1)
    return {"mrr": float(np.mean(rr)), "mrr_ci": [mrr_lo, mrr_hi],
            "h10": float(np.mean(hits10)), "h10_ci": [h_lo, h_hi],
            "n": len(test_genes)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--unit", required=True)
    ap.add_argument("--exports-dir", type=Path,
                    default=REPO_ROOT / "data" / "processed" / "m2c" / "completion_exports")
    ap.add_argument("--loci-json", type=Path, required=True)
    ap.add_argument("--cutoff", type=int, default=2020)
    ap.add_argument("--test-mode", choices=["hidden_post", "all_post"],
                    default="all_post")
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1337)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    unit = args.unit
    disease = unit.replace("__", ":")
    exports = args.exports_dir / unit
    pinned = PINNED_PATH.read_text(encoding="utf-8").strip().split()
    assert len(pinned) == 17448

    loci = json.loads(args.loci_json.read_text(encoding="utf-8"))
    fy_all = loci["gene_first_year_all"]
    fy_direct = loci["gene_first_year_direct"]

    train2_rows, valid_rows, test_rows = load_splits(exports)
    train_pos = {t for h, _, t in train2_rows if h == disease}
    valid_pos = {t for h, _, t in valid_rows if h == disease}
    test_genes = [t for h, _, t in test_rows if h == disease]
    assert all(h == disease for h, _, _ in test_rows)
    excluded = train_pos | valid_pos
    candidates = [g for g in pinned if g not in excluded]
    print(f"train positives {len(train_pos)}, candidates {len(candidates)}, "
          f"test {len(test_genes)}", flush=True)

    # Prospective test: hidden genes first reported >= cutoff (genetics clock).
    test_post = [g for g in test_genes if fy_all.get(g, -1) >= args.cutoff]
    test_post_direct = [g for g in test_genes if fy_direct.get(g, -1) >= args.cutoff]
    print(f"prospective test genes (all-traits clock): {len(test_post)}/{len(test_genes)}",
          flush=True)
    print(f"prospective test genes (direct-trait clock): "
          f"{len(test_post_direct)}/{len(test_genes)}", flush=True)
    if not test_post:
        raise SystemExit("empty prospective test set — lower --cutoff")

    post_all = {g for g, y in fy_all.items() if y >= args.cutoff}
    post_direct = {g for g, y in fy_direct.items() if y >= args.cutoff}
    if args.test_mode == "all_post":
        # True prospective semantics: every post-cutoff pinned gene is a test
        # query; pre-cutoff knowns are the only exclusions (post-cutoff genes
        # were unknown at cutoff, so they stay candidates).
        test_eval = sorted(post_all)
        pre_known = (set(train_pos) | set(valid_pos)) - post_all
        candidates = [g for g in pinned if g not in pre_known]
        drop = set(post_all)
        print(f"all_post mode: test {len(test_eval)} post-cutoff genes, "
              f"candidates {len(candidates)}", flush=True)
    else:
        test_eval = test_post
        drop = set()
    # FULL keeps full knowledge (no drop): "are future genes countable at all".
    # PRE_* enforce cutoff semantics: post-cutoff genes were unknown.
    variants = [("FULL", set(train_pos), None),
                ("PRE_ALL", set(train_pos) - post_all, drop),
                ("PRE_DIRECT", set(train_pos) - post_direct, drop)]
    results = {}
    for v, g_set, vdrop in variants:
        scores = null_scores(train2_rows, g_set, drop_genes=vdrop)
        rep = evaluate(scores, candidates, test_eval, args.bootstrap, args.seed)
        rep["g_size"] = len(g_set)
        results[v] = rep
        print(f"{v:10s} |G|={len(g_set):5d}  MRR {rep['mrr']:.5f} "
              f"[{rep['mrr_ci'][0]:.5f}, {rep['mrr_ci'][1]:.5f}]  "
              f"H@10 {rep['h10']:.5f} [{rep['h10_ci'][0]:.5f}, {rep['h10_ci'][1]:.5f}] "
              f"(n={rep['n']})", flush=True)

    out = {"unit": unit, "generated_at": utcnow_iso(), "seed": args.seed,
           "bootstrap": args.bootstrap, "cutoff": args.cutoff,
           "test_mode": args.test_mode,
           "loci_json": str(args.loci_json),
           "test_n": len(test_genes), "test_post_n": len(test_post),
           "test_post_direct_n": len(test_post_direct),
           "test_eval_n": len(test_eval),
           "train_positives": len(train_pos), "candidates": len(candidates),
           "variants": results,
           "eval": "filtered pinned 17448, candidates = pinned minus train2/valid "
                   "positives (post-cutoff genes stay candidates: unknown at cutoff), "
                   "optimistic-tie rank",
           "limitation": "first-reported clock is genetics-only (GWAS Catalog); "
                         "literature-only genes carry no year and stay in G_pre — "
                         "generous to the null, so a collapse is decisive"}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(f"wrote {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
