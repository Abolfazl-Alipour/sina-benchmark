#!/usr/bin/env python
"""m2c_reweighted_null — label-reweighting pre-test on the supervision-only null.

Owner directive 2026-09-05: before spending GPU on label reweighting, test its
premise on CPU. The supervision-only null (reports/m2c/cpu_nulls_and_coverage.md
section 4) is a pure function of the train2 labels:

    Score(d*, g) = sum_{d' in D(g)} |G_{d*} intersect genes(d')|     (train2 only)

    D(g)      = diseases with a train2 DIS-GEN edge to g
    genes(d') = train2 gene tails of d'
    G_{d*}    = train2 positives of the unit disease d*

If reweighting/down-weighting literature-derived edges is worth a GPU epoch, then
the same reweighting applied to the null must NOT collapse its test score. Two
outcomes, both decisive:
  - reweighted null DROPS  -> the top-10 signal is literature-popularity-shaped;
    reweighting the model's labels would hurt. Lever dies on CPU, zero GPU.
  - reweighted null HOLDS  -> non-literature labels carry recoverable signal;
    a reweighted K=all epoch is justified.

Variants (G-set construction; candidate side D(g) stays full because per-disease
OT evidence exists in-repo only for PD/AD — stated limitation):
  FULL        all train2 positives of d*             (control: PD must reproduce
                                                     MRR 0.01484 / H@10 0.02950)
  OT_MATCHED  positives with an OT association row   (missing-data control)
  NONLIT      OT row with any non-literature class score > 0
  GENETICS    OT row with genetics class score > 0
  STRICT_GF   OT row with genetics or functional_other > 0 (lit + curated out)
  LITONLY     OT row with literature > 0 and every other class == 0

Eval: identical to every completion-track number this week — filtered over the
pinned 17,448 (candidates = pinned minus train2/valid positives of d*),
optimistic-tie rank, hidden test scored. Bootstrap CIs (1000, seed 1337), same
convention as m2c_shortlist_generator_probe.py.

Usage (sina env, CPU):
    ~/miniconda3/envs/sina/bin/python scripts/pod/m2c_reweighted_null.py \
        --unit mondo__0005180 --out reports/m2c/reweighted_null_pd.json
"""
import argparse
import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
PINNED_PATH = REPO_ROOT / "data" / "processed" / "m2c" / "pinned_genes_17448.txt"
OT_RAW_DIR = REPO_ROOT / "data" / "processed" / "m2c" / "ot_panels"

# Same datasource -> evidence-class map as the committed circularity audit
# (scripts/pod/m2c_ot_panel_eval.py:61-75 / m2c_generator_wave2.py:54-61).
CLASS_OF_DATASOURCE = {
    "europepmc": "literature", "uniprot_literature": "literature",
    "genomics_england": "curated_db", "clingen": "curated_db",
    "orphanet": "curated_db", "gene2phenotype": "curated_db",
    "cancer_gene_census": "curated_db", "intogen": "curated_db",
    "gwas_credible_sets": "genetics", "ot_genetics_portal": "genetics",
    "eva": "genetics", "gene_burden": "genetics", "uniprot_variants": "genetics",
}
# Anything unlisted (impc, crispr_screen, expression_atlas, ...) = functional_other.


def utcnow_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def bootstrap_ci(values: np.ndarray, n_boot: int, seed: int) -> tuple[float, float]:
    """Percentile bootstrap over per-query scores (matches shortlist probe)."""
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


def ot_class_set(unit: str) -> dict[str, dict[str, float]]:
    """ensembl token -> {evidence class: max datasource score} from the OT raw pull."""
    path = OT_RAW_DIR / f"raw_{unit}.jsonl"
    out: dict[str, dict[str, float]] = {}
    with path.open(encoding="utf-8") as f:
        for line in f:
            row = json.loads(line)
            tok = "ensembl:" + row["target"]["id"]
            classes: dict[str, float] = defaultdict(float)
            for ds in row["datasourceScores"]:
                cls = CLASS_OF_DATASOURCE.get(ds["id"], "functional_other")
                classes[cls] = max(classes[cls], float(ds["score"]))
            out[tok] = dict(classes)
    return out


def variant_mask(g_set: set[str], ot: dict[str, dict[str, float]], variant: str) -> set[str]:
    if variant == "FULL":
        return set(g_set)
    if variant == "OT_MATCHED":
        return {g for g in g_set if g in ot}
    if variant == "NONLIT":
        return {g for g in g_set
                if any(s > 0 for c, s in ot.get(g, {}).items() if c != "literature")}
    if variant == "GENETICS":
        return {g for g in g_set if ot.get(g, {}).get("genetics", 0.0) > 0}
    if variant == "STRICT_GF":
        return {g for g in g_set
                if ot.get(g, {}).get("genetics", 0.0) > 0
                or ot.get(g, {}).get("functional_other", 0.0) > 0}
    if variant == "LITONLY":
        return {g for g in g_set
                if ot.get(g, {}).get("literature", 0.0) > 0
                and all(s == 0 for c, s in ot[g].items() if c != "literature")}
    raise ValueError(variant)


def null_scores(train2_rows: list, disease: str, g_set: set[str]) -> dict[str, float]:
    """Score(g) = sum_{d' in D(g)} |G intersect genes(d')|, train2 only."""
    genes_of: dict[str, set[str]] = defaultdict(set)
    for h, _, t in train2_rows:
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
    pos_of = {g: i for i, g in enumerate(candidates)}
    rr = np.zeros(len(test_genes))
    hits10 = np.zeros(len(test_genes))
    for i, g in enumerate(test_genes):
        s = scores.get(g, 0.0)
        rank = 1 + int(np.sum(cand_scores > s))  # optimistic tie, runner semantics
        if g in pos_of and s == 0.0:
            rank = 1 + int(np.sum(cand_scores > 0.0))
        rr[i] = 1.0 / rank
        hits10[i] = 1.0 if rank <= 10 else 0.0
    mrr_lo, mrr_hi = bootstrap_ci(rr, n_boot, seed)
    h_lo, h_hi = bootstrap_ci(hits10, n_boot, seed + 1)
    return {"mrr": float(np.mean(rr)), "mrr_ci": [mrr_lo, mrr_hi],
            "h10": float(np.mean(hits10)), "h10_ci": [h_lo, h_hi]}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--unit", required=True)
    ap.add_argument("--exports-dir", type=Path,
                    default=REPO_ROOT / "data" / "processed" / "m2c" / "completion_exports")
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1337)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    unit = args.unit
    disease = unit.replace("__", ":")
    exports = args.exports_dir / unit
    pinned = PINNED_PATH.read_text(encoding="utf-8").strip().split()
    assert len(pinned) == 17448

    train2_rows, valid_rows, test_rows = load_splits(exports)
    print(f"train2 {len(train2_rows)} valid {len(valid_rows)} test {len(test_rows)}",
          flush=True)

    train_pos = {t for h, _, t in train2_rows if h == disease}
    valid_pos = {t for h, _, t in valid_rows if h == disease}
    test_genes = [t for h, _, t in test_rows if h == disease]
    assert all(h == disease for h, _, _ in test_rows)
    excluded = train_pos | valid_pos
    candidates = [g for g in pinned if g not in excluded]
    print(f"train positives {len(train_pos)}, candidates {len(candidates)}, "
          f"test {len(test_genes)}", flush=True)

    ot = ot_class_set(unit)
    matched = sum(1 for g in train_pos if g in ot)
    print(f"OT rows matched to train positives: {matched}/{len(train_pos)}", flush=True)

    variants = ["FULL", "OT_MATCHED", "NONLIT", "GENETICS", "STRICT_GF", "LITONLY"]
    results = {}
    for v in variants:
        g_variant = variant_mask(train_pos, ot, v)
        scores = null_scores(train2_rows, disease, g_variant)
        rep = evaluate(scores, candidates, test_genes, args.bootstrap, args.seed)
        rep["g_size"] = len(g_variant)
        results[v] = rep
        print(f"{v:10s} |G|={len(g_variant):5d}  MRR {rep['mrr']:.5f} "
              f"[{rep['mrr_ci'][0]:.5f}, {rep['mrr_ci'][1]:.5f}]  "
              f"H@10 {rep['h10']:.5f} [{rep['h10_ci'][0]:.5f}, {rep['h10_ci'][1]:.5f}]",
              flush=True)

    # Context: evidence-class profile of the hidden test positives themselves.
    test_profile = {"ot_row": 0, "lit_only": 0, "has_genetics": 0,
                    "has_nonlit": 0, "no_ot_row": 0}
    for g in test_genes:
        if g not in ot:
            test_profile["no_ot_row"] += 1
            continue
        test_profile["ot_row"] += 1
        classes = ot[g]
        if classes.get("genetics", 0.0) > 0:
            test_profile["has_genetics"] += 1
        if any(s > 0 for c, s in classes.items() if c != "literature"):
            test_profile["has_nonlit"] += 1
        else:
            test_profile["lit_only"] += 1

    out = {"unit": unit, "generated_at": utcnow_iso(), "seed": args.seed,
           "bootstrap": args.bootstrap, "test_n": len(test_genes),
           "train_positives": len(train_pos), "ot_matched": matched,
           "candidates": len(candidates), "variants": results,
           "test_positive_profile": test_profile,
           "eval": "filtered pinned 17448, candidates = pinned minus train2/valid "
                   "positives, optimistic-tie rank",
           "limitation": "candidate side D(g) uses full train2; OT evidence classes "
                         "available in-repo only for the PD/AD pulls, so reweighting "
                         "applies to the disease-side gene set G only"}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(f"wrote {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
