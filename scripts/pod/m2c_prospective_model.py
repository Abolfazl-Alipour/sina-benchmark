#!/usr/bin/env python
"""m2c_prospective_model — model-side prospective + per-locus probes (CPU).

Owner directive 2026-09-06 (temporal + per-locus tracks). Reads a pinned-order
model score vector (m2c_ot_panel_eval.py score --dump-scores-dir) and the GWAS
loci JSON (m2c_gwas_loci_build.py). Zero training; the running GPU job is
never touched.

Probe A — prospective model: the null scores post-cutoff genes ~0 by
construction (prospective_null_*). But the MODEL sees graph features (GO,
phenotypes, expression), not just co-occurrence. Clean prospective queries =
post-cutoff pinned genes in NEITHER the train2 NOR valid splits of the unit
(the checkpoint never saw them as positives). Standard filtered genome-wide
eval on that query set: MRR / H@10 / median rank + bootstrap CIs. If the model
ranks them above rock-bottom while the null gives zeros, graph-side needle
signal exists. If the model also gives ~0, discovery is dead on this stack.

Probe B — per-locus choice: for each GWAS locus with >=2 pinned genes, rank
the in-locus genes by model score. No causal gold standard exists in-repo, so
the measurable question is known-gene recovery: in loci containing >=1 KNOWN
unit positive (train2 head positives), is the model's in-locus top-1 a known
gene? Baselines: nearest-gene (min SNP distance from Catalog up/downstream +
snp_gene = 0) and chance (1/n_locus_genes). Report top-1-known rate, MRR of
known genes in-locus, and model-vs-nearest agreement. Product-shaped: "of
these 5, which one" — validated where the answer is known.

Usage (sina env, CPU):
    ~/miniconda3/envs/sina/bin/python scripts/pod/m2c_prospective_model.py \
        --unit mondo__0005180 \
        --scores /tmp/pd_scores/mondo__0005180__epoch_001.scores.json \
        --loci-json reports/m2c/gwas_loci_mondo__0005180.json \
        --cutoff 2020 --out reports/m2c/prospective_model_pd.json
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


def load_splits(exports_dir: Path, unit: str, disease: str):
    def tails(fname):
        out = []
        with (exports_dir / unit / fname).open(encoding="utf-8") as f:
            for line in f:
                h, _, t = line.rstrip("\n").split("\t")
                if h == disease:
                    out.append(t)
        return out
    return set(tails("train2.txt")), set(tails("valid.txt")), tails("test.txt")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--unit", required=True)
    ap.add_argument("--scores", type=Path, required=True)
    ap.add_argument("--loci-json", type=Path, required=True)
    ap.add_argument("--exports-dir", type=Path,
                    default=REPO_ROOT / "data" / "processed" / "m2c" / "completion_exports")
    ap.add_argument("--cutoff", type=int, default=2020)
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1337)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    unit, disease = args.unit, args.unit.replace("__", ":")
    pinned = PINNED_PATH.read_text(encoding="utf-8").strip().split()
    assert len(pinned) == 17448
    pos_of = {g: i for i, g in enumerate(pinned)}

    dump = json.loads(args.scores.read_text(encoding="utf-8"))
    skeys = ["scores" in dump, "pinned_scores" in dump, "score_vector" in dump]
    skey = (["scores", "pinned_scores", "score_vector"])[skeys.index(True)] \
        if any(skeys) else None
    if skey is None:
        raise SystemExit(f"unknown score-dump keys: {sorted(dump.keys())}")
    assert dump.get("pinned_order") == "pinned_genes_17448.txt line order", \
        dump.get("pinned_order")
    scores = np.array(dump[skey], dtype=float)
    assert len(scores) == 17448, len(scores)
    print(f"score vector: {skey}, n={len(scores)}, ckpt={dump.get('checkpoint')}",
          flush=True)

    loci = json.loads(args.loci_json.read_text(encoding="utf-8"))
    fy_all = loci["gene_first_year_all"]
    post = {g for g, y in fy_all.items() if y >= args.cutoff}

    train_pos, valid_pos, _test_genes = load_splits(args.exports_dir, unit, disease)
    excluded = train_pos | valid_pos
    cand_scores = np.array([scores[pos_of[g]] for g in pinned if g not in excluded])

    def filtered_rank(g: str) -> int:
        s = scores[pos_of[g]]
        return 1 + int(np.sum(cand_scores > s))  # optimistic tie

    # --- Probe A: clean prospective queries (never seen as positives) ---
    clean = sorted(g for g in post if g in pos_of and g not in excluded)
    seen_post = sorted(g for g in post if g in excluded)
    ranks = np.array([filtered_rank(g) for g in clean]) if clean else np.zeros(0)
    rr = 1.0 / np.maximum(ranks, 1)
    if len(clean):
        mrr_lo, mrr_hi = bootstrap_ci(rr, args.bootstrap, args.seed)
        h_lo, h_hi = bootstrap_ci((ranks <= 10).astype(float), args.bootstrap,
                                  args.seed + 1)
        probe_a = {"n_clean": len(clean), "n_seen_post": len(seen_post),
                   "mrr": float(np.mean(rr)), "mrr_ci": [mrr_lo, mrr_hi],
                   "h10": float(np.mean(ranks <= 10)), "h10_ci": [h_lo, h_hi],
                   "median_rank": float(np.median(ranks)),
                   "clean_genes": clean}
    else:
        probe_a = {"n_clean": 0, "n_seen_post": len(seen_post)}
    print(f"Probe A: clean post-cutoff queries {len(clean)} "
          + (f"MRR {probe_a['mrr']:.5f} H@10 {probe_a['h10']:.5f} "
             f"median-rank {probe_a['median_rank']:.0f}" if len(clean) else "(none)"),
          flush=True)

    # Degree baseline for Probe B (popularity confound check): train1 degree.
    deg = defaultdict(int)
    with (args.exports_dir / unit / "train1.txt").open(encoding="utf-8") as f:
        for line in f:
            h, _, t = line.rstrip("\n").split("\t")
            if h in pos_of:
                deg[h] += 1
            if t in pos_of:
                deg[t] += 1

    # --- Probe B: per-locus known-gene recovery ---
    def locus_gene_dist(loc: dict) -> dict:
        dist = {}
        for a in loc["members"]:
            for g, roles in a["genes"].items():
                if "snp_gene" in roles and (g not in dist or dist[g] > 0):
                    dist[g] = 0
            for col, dcol in (("upstream_id", "upstream_dist"),
                              ("downstream_id", "downstream_dist")):
                gid, dval = a.get(col), a.get(dcol)
                if (gid and gid.startswith("ensembl:") and dval is not None
                        and (gid not in dist or dval < dist[gid])):
                    dist[gid] = dval
        return dist

    b_rows = []
    for loc in loci["loci"]:
        genes = loc["genes_pinned"]
        if len(genes) < 2:
            continue
        known = [g for g in genes if g in train_pos]
        if not known:
            continue
        by_model = sorted(genes, key=lambda g: (-scores[pos_of[g]], g))
        dist = locus_gene_dist(loc)
        with_d = [g for g in genes if g in dist]
        by_near = sorted(with_d, key=lambda g: (dist[g], g)) if with_d else []
        ranks_in = {}
        for i, g in enumerate(by_model):
            ranks_in[g] = i + 1
        top1 = by_model[0]
        by_deg = sorted(genes, key=lambda g: (-deg.get(g, 0), g))
        b_rows.append({
            "locus": loc["id"], "chr": loc["chr"], "n_genes": len(genes),
            "n_known": len(known), "model_top1": top1,
            "model_top1_known": top1 in train_pos,
            "degree_top1": by_deg[0],
            "degree_top1_known": by_deg[0] in train_pos,
            "model_degree_agree": by_deg[0] == top1,
            "nearest_top1": by_near[0] if by_near else None,
            "nearest_top1_known": (by_near[0] in train_pos) if by_near else None,
            "agree": (by_near[0] == top1) if by_near else None,
            "known_ranks_model": sorted(ranks_in[g] for g in known),
            "chance_top1": len(known) / len(genes),
        })
    n_b = len(b_rows)
    if n_b:
        m_top1 = float(np.mean([r["model_top1_known"] for r in b_rows]))
        near_rows = [r for r in b_rows if r["nearest_top1"] is not None]
        n_top1 = float(np.mean([r["nearest_top1_known"] for r in near_rows])) \
            if near_rows else None
        chance = float(np.mean([r["chance_top1"] for r in b_rows]))
        agree = float(np.mean([r["agree"] for r in near_rows])) if near_rows else None
        mrr_known = float(np.mean([1.0 / r for row in b_rows
                                   for r in row["known_ranks_model"]]))
        d_top1 = float(np.mean([r["degree_top1_known"] for r in b_rows]))
        md_agree = float(np.mean([r["model_degree_agree"] for r in b_rows]))
        probe_b = {"n_loci_eval": n_b, "n_loci_nearest": len(near_rows),
                   "model_top1_known_rate": m_top1,
                   "degree_top1_known_rate": d_top1,
                   "model_degree_agreement": md_agree,
                   "nearest_top1_known_rate": n_top1,
                   "chance_top1_known_rate": chance,
                   "model_nearest_agreement": agree,
                   "mrr_known_in_locus": mrr_known, "rows": b_rows}
    else:
        probe_b = {"n_loci_eval": 0}
    print(f"Probe B: {n_b} multi-gene loci with known genes"
          + (f" | model top1-known {probe_b['model_top1_known_rate']:.3f} vs nearest "
             f"{probe_b['nearest_top1_known_rate']} vs chance "
             f"{probe_b['chance_top1_known_rate']:.3f}" if n_b else ""), flush=True)

    out = {"unit": unit, "generated_at": utcnow_iso(), "seed": args.seed,
           "bootstrap": args.bootstrap, "cutoff": args.cutoff,
           "scores": str(args.scores), "ckpt": dump.get("checkpoint"),
           "probe_a_prospective_model": probe_a,
           "probe_b_per_locus": probe_b}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(f"wrote {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
