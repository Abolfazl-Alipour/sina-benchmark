#!/usr/bin/env python
"""m2c_generator_wave2 — candidate-GENERATION wave 2 (CPU-only, zero training).

Probe 1/1.5 + cpu_nulls_and_coverage.md sections 5-6 proved generation is the bottleneck:
degree-top-100 covers 2/339 hidden PD positives (8.85% excl-knowns), OT overall
covers 87/339, and the reranker's conversion (~3.3x null) is stack-invariant.
Owner directive 2026-09-05: attack candidate generation on CPU/API (decision 2c).

Generators (all: pinned 17,448 universe, train2+valid positives of the disease
excluded, tiebreak gene token, K in {100,200,500}):
  ot_overall     — reference, re-derived from the committed raw pull (§5 champion)
  ot_genetics    — OT datasource harmonic-share restricted to genetics class
  ot_literature  — same, literature class (europepmc, uniprot_literature)
  ot_functional  — same, functional_other class (impc, crispr_screen, ...)
  pheno_shared   — disease PHENOTYPE_PRESENT phenotypes (KEEP train1,
                   disease->hp, 11,850 edges) x gene via hp-ASSOCIATED_WITH->gene
                   (768,844 edges); score = #shared phenotypes
  drug_indication— drugs -INDICATION-> disease (1,741) then drug -TARGET-> gene
                   (1,886); score = #indicated drugs targeting the gene
  union_best     — rank fusion across generators: gene's panel rank = best
                   (minimum) rank in any single generator; ties broken by token.

Outputs: reports/m2c/panels_wave2/{unit}_{gen}_top{K}.txt + manifest.json
(coverage@K over the hidden test positives per generator per K + union table).

  stage: two-stage recall@10 per panel using a DUMPED SCORE VECTOR
  (reports/m2c/ot_panels/score_vectors/*.scores.json — no model load, no GPU):
      python m2c_generator_wave2.py stage --scores <dump.json> [--gen all]

sina env: ~/miniconda3/envs/sina/bin/python scripts/pod/m2c_generator_wave2.py build
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

PINNED_PATH = REPO_ROOT / "data/processed/m2c/pinned_genes_17448.txt"
DEFAULT_EXPORTS = REPO_ROOT / "data/processed/m2c/completion_exports"
RAW_DIR = REPO_ROOT / "data/processed/m2c/ot_panels"
OUT_DIR = REPO_ROOT / "reports/m2c/panels_wave2"
OT_PANELS_DIR = REPO_ROOT / "reports/m2c/ot_panels"

UNITS = {
    "mondo__0005180": {"label": "PD", "n_test": 339},
    "mondo__0004975": {"label": "AD", "n_test": 621},
}
KS = (100, 200, 500)

CLASS_OF_DATASOURCE = {
    "europepmc": "literature", "uniprot_literature": "literature",
    "genomics_england": "curated_db", "clingen": "curated_db",
    "orphanet": "curated_db", "gene2phenotype": "curated_db",
    "cancer_gene_census": "curated_db", "intogen": "curated_db",
    "gwas_credible_sets": "genetics", "ot_genetics_portal": "genetics",
    "eva": "genetics", "gene_burden": "genetics", "uniprot_variants": "genetics",
}


def utcnow_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_pinned() -> list[str]:
    pinned = PINNED_PATH.read_text(encoding="utf-8").strip().split()
    assert len(pinned) == 17448
    return pinned


def split_tails(exports_dir: Path, unit: str, fname: str) -> set[str]:
    head = unit.replace("__", ":")
    out = set()
    for line in (exports_dir / unit / fname).read_text(encoding="utf-8").splitlines():
        h, _r, t = line.split("\t")
        if h == head:
            out.add(t)
    return out


def harmonic_total(datasource_scores: list[dict], want_class: str | None) -> float:
    """OT harmonic-sum overall score restricted to one evidence class (or all)."""
    ss = sorted((float(d["score"]) for d in datasource_scores
                 if want_class is None
                 or CLASS_OF_DATASOURCE.get(d["id"], "functional_other") == want_class),
                reverse=True)
    return sum(s / (i + 1) ** 2 for i, s in enumerate(ss))


def build_ot_class_panels(unit: str, known: set[str], pinned: set[str]):
    """Return {gen_name: [(token, score)] sorted} for ot_overall + 3 classes.
    ot_overall ranks by the API's own overall score (identical construction to
    m2c_ot_panel_eval.py pull — reproduces the committed §5 panels); class panels
    rank by harmonic sum restricted to the class's datasources."""
    rows = {}
    with (RAW_DIR / f"raw_{unit}.jsonl").open(encoding="utf-8") as f:
        for line in f:
            r = json.loads(line)
            tok = f"ensembl:{r['target']['id']}"
            if tok in pinned and tok not in known:
                rows[tok] = (float(r["score"]), r["datasourceScores"])
    out = {}
    for gen, cls in (("ot_overall", None), ("ot_genetics", "genetics"),
                     ("ot_literature", "literature"), ("ot_functional", "functional_other")):
        if cls is None:
            scored = [(t, s0) for t, (s0, _ds) in rows.items()]
        else:
            scored = [(t, harmonic_total(ds, cls)) for t, (_s0, ds) in rows.items()]
        scored = [x for x in scored if x[1] > 0.0]
        scored.sort(key=lambda x: (-x[1], x[0]))
        out[gen] = scored
    return out


def build_graph_panels(unit: str, exports_dir: Path, known: set[str], pinned: set[str]):
    """pheno_shared + drug_indication from KEEP train1 edge types."""
    head = unit.replace("__", ":")
    pheno_of_disease: set[str] = set()
    gene_phenos: dict[str, set[str]] = {}
    drug_indicated: set[str] = set()
    gene_drugs: dict[str, set[str]] = {}
    for line in (exports_dir / unit / "train1.txt").read_text(encoding="utf-8").splitlines():
        h, r, t = line.split("\t")
        if r == "PHENOTYPE_PRESENT" and h == head:
            pheno_of_disease.add(t)
        elif r == "ASSOCIATED_WITH":
            # phenotype -> gene (HPO associations; 768,844 edges)
            if h.startswith("hp:") and t.startswith("ensembl:"):
                gene_phenos.setdefault(t, set()).add(h)
        elif r == "INDICATION" and t == head:
            drug_indicated.add(h)
        elif r == "TARGET":
            gene_drugs.setdefault(t, set()).add(h)
    pheno_scored = []
    for g, ph in gene_phenos.items():
        if g in pinned and g not in known:
            n = len(ph & pheno_of_disease)
            if n > 0:
                pheno_scored.append((g, float(n)))
    pheno_scored.sort(key=lambda x: (-x[1], x[0]))
    drug_scored = []
    for g, dr in gene_drugs.items():
        if g in pinned and g not in known:
            n = len(dr & drug_indicated)
            if n > 0:
                drug_scored.append((g, float(n)))
    drug_scored.sort(key=lambda x: (-x[1], x[0]))
    return {"pheno_shared": pheno_scored, "drug_indication": drug_scored}, {
        "disease_phenotypes": len(pheno_of_disease),
        "genes_with_any_phenotype": len(gene_phenos),
        "genes_sharing_ge1": len(pheno_scored),
        "indicated_drugs": len(drug_indicated),
        "genes_targeted_by_indicated": len(drug_scored),
    }


def cmd_build(args) -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pinned_list = load_pinned()
    pinned = set(pinned_list)
    manifest = {"built_at": utcnow_iso(), "units": {}}
    for unit, spec in UNITS.items():
        known = split_tails(args.exports_dir, unit, "train2.txt") | \
            split_tails(args.exports_dir, unit, "valid.txt")
        test_pos = split_tails(args.exports_dir, unit, "test.txt")
        assert len(test_pos) == spec["n_test"]
        gens = build_ot_class_panels(unit, known, pinned)
        graph_gens, graph_stats = build_graph_panels(unit, args.exports_dir, known, pinned)
        gens.update(graph_gens)

        # union_best: min rank across all single generators
        best_rank: dict[str, int] = {}
        for scored in gens.values():
            for rank0, (tok, _s) in enumerate(scored):
                if tok not in best_rank or rank0 < best_rank[tok]:
                    best_rank[tok] = rank0
        union = sorted(best_rank.items(), key=lambda x: (x[1], x[0]))
        gens["union_best"] = [(t, float(r)) for t, r in union]

        unit_man = {"label": spec["label"], "knowns_excluded": len(known),
                    "graph_stats": graph_stats, "generators": {}}
        for gen, scored in gens.items():
            gen_man = {"n_scored": len(scored), "panels": {}}
            for k in KS:
                panel = [t for t, _s in scored[:k]]
                path = OUT_DIR / f"{unit}_{gen}_top{k}.txt"
                path.write_text("".join(f"{t}\n" for t in panel), encoding="utf-8")
                cov = sum(1 for t in panel if t in test_pos)
                gen_man["panels"][str(k)] = {
                    "sha256": sha256_file(path),
                    "test_coverage": cov / len(test_pos),
                    "test_hits_in_panel": cov}
            unit_man["generators"][gen] = gen_man
            covs = " ".join(
                f"K{k}={gen_man['panels'][str(k)]['test_hits_in_panel']}"
                for k in KS)
            print(f"[{utcnow_iso()}] {unit} {gen:16s} scored={len(scored):6d} {covs}",
                  flush=True)
        manifest["units"][unit] = unit_man
    mpath = OUT_DIR / "manifest.json"
    mpath.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                     encoding="utf-8")
    print(f"[{utcnow_iso()}] manifest -> {mpath}", flush=True)
    return 0


def cmd_stage(args) -> int:
    """Two-stage recall@10 per panel from a dumped score vector (no model load)."""
    import numpy as np

    dump = json.loads(Path(args.scores).read_text(encoding="utf-8"))
    unit = dump["unit"]
    n_test = UNITS[unit]["n_test"]
    scores = dump["scores"]
    pinned = load_pinned()
    pos_of = {t: i for i, t in enumerate(pinned)}
    test_pos = sorted(split_tails(args.exports_dir, unit, "test.txt"))
    test_positions = [pos_of[t] for t in test_pos]
    assert len(test_positions) == n_test

    def bootstrap_ci(values, n_boot=1000, seed=1337):
        rng = np.random.default_rng(seed)
        reps = np.empty(n_boot)
        n = len(values)
        for b in range(n_boot):
            reps[b] = float(np.mean(values[rng.integers(0, n, n)]))
        lo, hi = np.percentile(reps, [2.5, 97.5])
        return float(lo), float(hi)

    results = {"unit": unit, "scores_input": args.scores,
               "scores_ckpt_sha256": dump["ckpt_sha256"], "panels": {}}
    for k in KS:
        for pf in sorted(OUT_DIR.glob(f"{unit}_*_top{k}.txt")):
            gen = pf.name.replace(f"{unit}_", "").replace(f"_top{k}.txt", "")
            toks = pf.read_text(encoding="utf-8").split()
            panel_idx = [pos_of[t] for t in toks]
            panel_set = set(panel_idx)
            order = sorted(panel_idx, key=lambda i: (-scores[i], pinned[i]))
            rank = {i: r1 for r1, i in enumerate(order, start=1)}
            v = [1.0 if (i in panel_set and rank[i] <= 10) else 0.0
                 for i in test_positions]
            cov = sum(1 for i in test_positions if i in panel_set)
            results["panels"].setdefault(gen, {})[str(k)] = {
                "coverage": cov / n_test,
                "two_stage_recall10": float(np.mean(v)),
                "two_stage_ci": list(bootstrap_ci(np.array(v)))}
    args.out.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    for gen, ks in sorted(results["panels"].items()):
        row = " ".join(f"K{k}: cov {v['coverage']:.3f} 2stage {v['two_stage_recall10']:.4f}"
                       for k, v in sorted(ks.items()))
        print(f"{gen:16s} {row}", flush=True)
    print(f"[{utcnow_iso()}] wrote {args.out}", flush=True)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p_b = sub.add_parser("build")
    p_b.add_argument("--exports-dir", type=Path, default=DEFAULT_EXPORTS)
    p_b.set_defaults(fn=cmd_build)
    p_s = sub.add_parser("stage")
    p_s.add_argument("--scores", required=True,
                     help="reports/m2c/ot_panels/score_vectors/*.scores.json")
    p_s.add_argument("--exports-dir", type=Path, default=DEFAULT_EXPORTS)
    p_s.add_argument("--out", type=Path, default=OUT_DIR / "twostage.json")
    p_s.set_defaults(fn=cmd_stage)
    args = parser.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
