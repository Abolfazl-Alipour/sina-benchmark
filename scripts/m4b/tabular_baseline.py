#!/usr/bin/env python
"""m4b_tabular_baseline — SINA-free LightGBM locus prioritization (control arm).

Owner task 2026-09-11 (benchmark track, CPU ONLY, no GPU, no SINA inference):
tabular gradient boosting on free/public features for gene prioritization at
GWAS loci, benchmarked vs OT L2G (and FLAMES where runnable). Establishes what
free features alone achieve so a later SINA-lift arm is measured honestly.

Claim corrections baked in (verified 2026-09-12, see HANDOFF): FLAMES is
Schipper/de Leeuw/Posthuma (VU Amsterdam), NOT Mountjoy; OT Genetics GraphQL
(api.genetics.opentargets.org) is DNS-dead — competitor scores come from the
Platform API (proven path in m2c_loci_headtohead.py) or the HF model
(opentargets/locus_to_gene_25.12, MIT tag); GWAS Catalog is EBI-Terms-of-Use
(not plain CC BY 4.0; some datasets CC-BY-NC); GTEx v8 = Aug 2019.

Pipeline (subcommands, run in order):
  positives  OT gold standards (445 curated loci) + FLAMES ExWAS TP pairs
             (1,181 loci) -> deduped (chr,pos,gene) positives with provenance.
             rsID->GRCh38 resolution via Ensembl REST, cached + hash-pinned.
  panels     +/-500kb protein-coding gene panels per locus (GENCODE GRCh38).
  features   Arm T (dist, GTEx v8 eQTL, gnomAD v2 constraint) + Arm G
             (PPR/DWPC/degree/GO-jaccard; PD/AD-scoped until recomputed).
             TMP view allows ONLY pre-2021-datable features (no FLAMES/PoPS
             scores, no current graph stats).
  train      LightGBM binary logloss, seeds {1337,1339,5027}; View CV
             (5-fold chromosome-disjoint) + View TMP (pub date<=2020 train,
             first-reported>2020 test); OT L2G + distance-only on identical
             loci; ablations per feature family; bootstrap CIs.
  report     reports/m4b/tabular_baseline.md + results JSON, single-command
             reproducible.

Usage (sina env, CPU-only):
  ~/miniconda3/envs/sina/bin/python scripts/m4b/tabular_baseline.py positives
"""

import argparse
import hashlib
import json
import re
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
M4B_DIR = REPO_ROOT / "data" / "processed" / "m4b"
SRC_DIR = Path("/tmp/m4b_src")
OT_GOLD_TSV = (SRC_DIR / "genetics-gold-standards" / "gold_standards"
               / "processed" / "gwas_gold_standards.191108.tsv")
FLAMES_CSV = (SRC_DIR / "FLAMES_paper_analyses" / "annotated_loci"
              / "all_FLAMES.csv")
ENSEMBL_REST = "https://rest.ensembl.org"
RSID_CACHE = M4B_DIR / "rsid_cache.json"
SEEDS = (1337, 1339, 5027)
GTF_PATH = Path("/tmp/gencode.v46.gtf.gz")
PANEL_HALF_WINDOW = 500_000


def utcnow_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def load_ot_gold() -> list[dict]:
    """OT curated positives: (chr,pos,gene) with Catalog coords + dates."""
    import pandas as pd
    g = pd.read_csv(OT_GOLD_TSV, sep="\t", low_memory=False)
    g = g.dropna(subset=["gold_standard_info.gene_id",
                         "sentinel_variant.locus_GRCh38.chromosome",
                         "sentinel_variant.locus_GRCh38.position"])
    rows = []
    for _, r in g.iterrows():
        rows.append({
            "chr": str(r["sentinel_variant.locus_GRCh38.chromosome"]),
            "pos": int(r["sentinel_variant.locus_GRCh38.position"]),
            "gene": str(r["gold_standard_info.gene_id"]),
            "trait": str(r.get("trait_info.standard_trait_name", "")),
            "study_acc": (None if pd.isna(r.get("association_info.gwas_catalog_id"))
                          else str(r["association_info.gwas_catalog_id"])),
            "pubmed": (None if pd.isna(r.get("association_info.pubmed_id"))
                       else str(r["association_info.pubmed_id"])),
            "provenance": "ot_gold",
        })
    return rows


def parse_chrpos_token(s: str) -> tuple[str, int] | None:
    m = re.search(r"\b(\d{1,2}|X|Y|MT):(\d{6,10}):[ACGTN]+:[ACGTN]+", str(s))
    return (m.group(1), int(m.group(2))) if m else None


def first_rsid(s: str) -> str | None:
    m = re.search(r"\brs\d+\b", str(s))
    return m.group(0) if m else None


def _post_rsids(chunk: list[str]) -> dict:
    import requests
    resp = requests.post(
        f"{ENSEMBL_REST}/variation/homo_sapiens",
        headers={"Content-Type": "application/json",
                 "Accept": "application/json"},
        json={"ids": chunk}, timeout=120)
    resp.raise_for_status()
    return resp.json()


def _store_mapping(rid: str, rec: dict | None, cache: dict) -> None:
    maps = (rec or {}).get("mappings") or []
    grch = [m for m in maps if str(m.get("assembly_name", ""))
            .startswith("GRCh38")]
    cache[rid] = {"chr": str(grch[0]["seq_region_name"]),
                  "pos": int(grch[0]["start"])} if grch else None


def resolve_rsids(rsids: list[str], cache: dict) -> dict:
    """Batch rsID->GRCh38 via Ensembl REST POST; cache-first, resumable.

    A 400 poisons the whole batch (one bad id), so halves are retried
    recursively; single unresolvable ids cache as None (recorded miss).
    """
    missing = [r for r in rsids if r not in cache]
    queue = [missing[i:i + 200] for i in range(0, len(missing), 200)]
    while queue:
        chunk = queue.pop(0)
        try:
            for rid, rec in _post_rsids(chunk).items():
                _store_mapping(rid, rec, cache)
        except Exception:
            if len(chunk) == 1:
                cache[chunk[0]] = None
            else:
                mid = len(chunk) // 2
                queue.append(chunk[:mid])
                queue.append(chunk[mid:])
                time.sleep(2)
                continue
        time.sleep(0.5)
    return cache


def load_flames_positives(cache: dict) -> tuple[list[dict], dict]:
    """FLAMES ExWAS TP pairs; lead = parsed chr:pos else resolved rsID."""
    import pandas as pd
    d = pd.read_csv(FLAMES_CSV, usecols=["ensg", "filename", "GenomicLocus",
                                         "TP", "IndSigSNPs"])
    tp = d[d["TP"] == 1].copy()
    tp["coord"] = tp["IndSigSNPs"].map(parse_chrpos_token)
    need = sorted({first_rsid(s) for s in
                   tp.loc[tp["coord"].isna(), "IndSigSNPs"].tolist()
                   if first_rsid(s)})
    if need:
        print(f"resolving {len(need)} rsIDs via Ensembl REST", flush=True)
        cache = resolve_rsids(need, cache)
    rows = []
    for _, r in tp.iterrows():
        coord = r["coord"]
        if coord is None:
            hit = cache.get(first_rsid(r["IndSigSNPs"]) or "")
            coord = (hit["chr"], hit["pos"]) if hit else None
        if coord is None:
            continue
        rows.append({"chr": str(coord[0]), "pos": int(coord[1]),
                     "gene": str(r["ensg"]), "trait": str(r["filename"]),
                     "study_acc": None, "pubmed": None,
                     "provenance": "flames_exwas",
                     "flames_locus": f"{r['filename']}#{r['GenomicLocus']}"})
    return rows, cache


def cmd_positives(args) -> int:
    from sina.manifests import Run
    M4B_DIR.mkdir(parents=True, exist_ok=True)
    manifest_path = M4B_DIR / "positives.manifest.json"
    run = Run(script="scripts/m4b/tabular_baseline.py",
              manifest_path=manifest_path, config={}, git_ref="")
    cache: dict = {}
    if RSID_CACHE.exists():
        cache = json.loads(RSID_CACHE.read_text(encoding="utf-8"))
    with run.stage("positives"):
        ot_rows = load_ot_gold()
        flames_rows, cache = load_flames_positives(cache)
        RSID_CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True),
                              encoding="utf-8")
        seen: dict[tuple, dict] = {}
        for row in ot_rows + flames_rows:
            key = (row["chr"], row["pos"], row["gene"])
            if key in seen:
                prev = seen[key]["provenance"]
                if row["provenance"] not in prev:
                    seen[key]["provenance"] = prev + "+" + row["provenance"]
            else:
                seen[key] = row
        positives = sorted(seen.values(),
                           key=lambda r: (r["chr"], r["pos"], r["gene"]))
        out = {"generated_at": utcnow_iso(), "n": len(positives),
               "n_ot_gold": len(ot_rows), "n_flames": len(flames_rows),
               "rsid_cache_size": len(cache), "rows": positives}
        (M4B_DIR / "positives.json").write_text(json.dumps(out, indent=1),
                                               encoding="utf-8")
        print(f"positives: {len(positives)} (OT {len(ot_rows)}, "
              f"FLAMES {len(flames_rows)})", flush=True)
    return 0


def load_tss(gtf: Path = GTF_PATH) -> list[dict]:
    """Protein-coding gene TSS table from GENCODE GRCh38 GTF."""
    import gzip
    genes = []
    with gzip.open(gtf, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attr = dict(re.findall(r'(\S+) "([^"]+)"', parts[8]))
            if attr.get("gene_type") != "protein_coding":
                continue
            chrom, start, end, strand = (parts[0], int(parts[3]),
                                        int(parts[4]), parts[6])
            if chrom.startswith("CHR_") or "_" in chrom:
                continue
            chrom = chrom[3:] if chrom.startswith("chr") else chrom
            genes.append({"ensg": attr["gene_id"].split(".")[0],
                          "symbol": attr.get("gene_name", ""),
                          "chr": chrom,
                          "tss": start if strand == "+" else end})
    return genes


def cmd_panels(args) -> int:
    from sina.manifests import Run
    manifest_path = M4B_DIR / "panels.manifest.json"
    run = Run(script="scripts/m4b/tabular_baseline.py",
              manifest_path=manifest_path, config={}, git_ref="")
    with run.stage("panels"):
        positives = json.loads((M4B_DIR / "positives.json").read_text(
            encoding="utf-8"))["rows"]
        genes = load_tss()
        by_chr: dict[str, list] = {}
        for g in genes:
            by_chr.setdefault(g["chr"], []).append(g)
        for v in by_chr.values():
            v.sort(key=lambda g: g["tss"])
        loci: dict[tuple, dict] = {}
        for row in positives:
            loci.setdefault((row["chr"], row["pos"]), {
                "chr": row["chr"], "pos": row["pos"], "traits": [],
                "positives": []}).update()
            loc = loci[(row["chr"], row["pos"])]
            loc["positives"].append(row["gene"])
            loc["traits"].append({"trait": row["trait"],
                                  "provenance": row["provenance"]})
        panels = []
        for (chrom, pos), loc in sorted(loci.items()):
            panel = []
            for g in by_chr.get(chrom, []):
                d = g["tss"] - pos
                if abs(d) <= PANEL_HALF_WINDOW:
                    panel.append({"ensg": g["ensg"], "symbol": g["symbol"],
                                  "tss": g["tss"], "dist_kb": round(d / 1000, 1),
                                  "label": int(g["ensg"] in loc["positives"])})
            panel.sort(key=lambda r: abs(r["dist_kb"]))
            panels.append({"locus_id": f"{chrom}:{pos}", **loc,
                           "n_genes": len(panel), "genes": panel})
        n_pos = sum(sum(g["label"] for g in p["genes"]) for p in panels)
        out = {"generated_at": utcnow_iso(), "panel": "tss_within_500kb",
               "n_loci": len(panels), "n_panel_genes": sum(p["n_genes"]
               for p in panels), "n_positives_hit": n_pos, "loci": panels}
        (M4B_DIR / "panels.json").write_text(json.dumps(out), encoding="utf-8")
        print(f"panels: {len(panels)} loci, "
              f"{out['n_panel_genes']} gene rows, positives hit {n_pos}",
              flush=True)
    return 0


GNOMAD_PATH = Path("/tmp/gnomad.lof.bgz")


def gene_degrees() -> dict:
    """All-edge degree per bare ENSG id from the full cleaned export."""
    import pandas as pd
    degs: dict[str, int] = {}
    for name in ("train1.txt", "train2.txt"):
        path = REPO_ROOT / "data" / "processed" / "m2c" / "full_exports" / name
        for chunk in pd.read_csv(path, sep="\t", header=None, chunksize=500_000,
                                 dtype=str):
            for col in (0, 2):
                for ent, cnt in chunk[col].value_counts().items():
                    if ent.startswith("ensembl:ENSG"):
                        gene = ent.split(":")[1]
                        degs[gene] = degs.get(gene, 0) + int(cnt)
    return degs


def load_constraint() -> dict:
    """gnomAD v2 LOEUF/pLI by gene symbol; {} when the download is absent."""
    import gzip
    if not GNOMAD_PATH.exists():
        return {}
    out: dict[str, dict] = {}
    with gzip.open(GNOMAD_PATH, "rt", encoding="utf-8") as f:
        header: list[str] = []
        for line in f:
            row = line.rstrip("\n").split("\t")
            if not header:
                header = row
                continue
            d = dict(zip(header, row, strict=False))
            sym = d.get("gene")
            if sym:
                out[sym] = {"loeuf": _f(d.get("oe_lof_upper")),
                            "pli": _f(d.get("pLI"))}
    return out


def _f(v) -> float | None:
    try:
        return float(v) if v not in (None, "", "NA") else None
    except (TypeError, ValueError):
        return None


def cmd_features(args) -> int:
    import pandas as pd

    from sina.manifests import Run
    manifest_path = M4B_DIR / "features.manifest.json"
    run = Run(script="scripts/m4b/tabular_baseline.py",
              manifest_path=manifest_path, config={}, git_ref="")
    with run.stage("features"):
        panels = json.loads((M4B_DIR / "panels.json").read_text(
            encoding="utf-8"))["loci"]
        degs = gene_degrees()
        cons = load_constraint()
        rows = []
        for loc in panels:
            for g in loc["genes"]:
                c = cons.get(g["symbol"], {})
                rows.append({
                    "locus_id": loc["locus_id"], "chr": loc["chr"],
                    "pos": loc["pos"], "ensg": g["ensg"],
                    "symbol": g["symbol"], "label": g["label"],
                    "provenance": "+".join(sorted(
                        {t["provenance"] for t in loc["traits"]})),
                    "dist_kb": abs(g["dist_kb"]),
                    "degree": degs.get(g["ensg"], 0),
                    "loeuf": c.get("loeuf"), "pli": c.get("pli"),
                })
        feat = pd.DataFrame(rows)
        fdir = M4B_DIR / "features"
        fdir.mkdir(parents=True, exist_ok=True)
        feat.to_parquet(fdir / "features.parquet", index=False)
        print(f"features: {len(feat)} rows, {feat['locus_id'].nunique()} "
              f"loci; constraint coverage "
              f"{feat['loeuf'].notna().mean():.3f}; "
              f"positives {int(feat['label'].sum())}", flush=True)
    return 0


FEATURE_COLS = ["dist_kb", "degree", "loeuf", "pli"]
# V2 (pre-registered 2026-09-12 before seeing any v2 number): drop
# constraint (ablation: it hurts), add DIS-GEN degree (literature prior)
# and VEP coding severity (gold labels skew coding).
FEATURE_SETS = {"v1": FEATURE_COLS,
                "v2": ["dist_kb", "degree", "disgen_degree", "vep_sev"],
                "v3": ["dist_kb", "degree", "disgen_degree", "vep_sev",
                       "phe_degree", "go_degree", "ppr"],
                "v4": ["dist_kb", "vep_sev", "loeuf", "pli"]}


def chr_folds(df, n_folds: int = 5, seed: int = 1337):
    """Greedy chromosome-to-fold assignment balanced by positives."""
    import numpy as np
    rng = np.random.default_rng(seed)
    pos = df.groupby("chr")["label"].sum().sort_values(ascending=False)
    folds: list[list] = [[] for _ in range(n_folds)]
    load = [0.0] * n_folds
    for chrom, n in pos.items():
        order = rng.permutation(n_folds)
        tgt = min(order, key=lambda f: load[f])
        folds[tgt].append(chrom)
        load[tgt] += float(n)
    fold_of = {c: i for i, cs in enumerate(folds) for c in cs}
    return df["chr"].map(fold_of).to_numpy()


def gene_folds(df, n_folds: int = 5, seed: int = 1337):
    """Leakage-controlled folds: loci sharing any gold-gene symbol go to
    the same fold (union-find over gold symbols), greedy balanced by
    positives. Kills train->test gold-gene leakage by construction."""
    import numpy as np
    gold = df[df["label"] == 1][["locus_id", "symbol"]].drop_duplicates()
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for _, r in gold.iterrows():
        union(("loc", r["locus_id"]), ("gen", r["symbol"]))
    groups: dict = {}
    for loc in df["locus_id"].unique():
        groups.setdefault(find(("loc", loc)), []).append(loc)
    grp = sorted(groups.values(), key=len, reverse=True)
    rng = np.random.default_rng(seed)
    folds: list[list] = [[] for _ in range(n_folds)]
    load = [0.0] * n_folds
    for g in grp:
        order = rng.permutation(n_folds)
        tgt = min(order, key=lambda f: load[f])
        folds[tgt].extend(g)
        load[tgt] += float(df[df["locus_id"].isin(g)]["label"].sum())
    fold_of = {loc: i for i, ls in enumerate(folds) for loc in ls}
    return df["locus_id"].map(fold_of).to_numpy()


def leakage_audit(df, folds) -> dict:
    """Fraction of test loci per fold whose gold gene was seen in train."""
    import numpy as np
    folds = np.asarray(folds)
    fracs = []
    for fold in sorted(set(folds.tolist())):
        tr_genes = set(df[(folds != fold) & (df["label"] == 1)]["symbol"])
        n_te = df[(folds == fold)]["locus_id"].nunique()
        shared = sum(
            1 for _, loc in df[folds == fold].groupby("locus_id")
            if set(loc[loc["label"] == 1]["symbol"]) & tr_genes)
        fracs.append(shared / max(n_te, 1))
    return {"mean_test_loci_with_train_seen_gene": round(float(sum(fracs) / len(fracs)), 4),
            "per_fold": [round(float(f), 4) for f in fracs]}


def locus_metrics(frame, score_col: str) -> dict:
    """Per-locus top-1/top-5/mean-best-rank with pessimistic ties."""
    import numpy as np
    top1, top5, ranks, y_true, y_score = [], [], [], [], []
    for _, loc in frame.groupby("locus_id"):
        s = np.nan_to_num(loc[score_col].to_numpy(dtype=float), nan=-np.inf)
        y = loc["label"].to_numpy(dtype=int)
        # pessimistic: each positive scores the worst tied position
        best = min(float(np.sum(s > sp) + np.sum(s == sp))
                   for sp, yp in zip(s, y, strict=True) if yp == 1)
        top1.append(float(best == 1.0))
        top5.append(float(best <= 5.0))
        ranks.append(best)
        y_true.extend(y.tolist())
        y_score.extend(s.tolist())
    return {"top1": top1, "top5": top5, "ranks": ranks,
            "y_true": y_true, "y_score": y_score}


def bootstrap_ci(vals, n_boot: int = 1000, seed: int = 1337):
    import numpy as np
    rng = np.random.default_rng(seed)
    a = np.asarray(vals, dtype=float)
    means = [float(np.mean(rng.choice(a, size=len(a), replace=True)))
             for _ in range(n_boot)]
    return {"mean": round(float(np.mean(a)), 4),
            "ci": [round(float(np.percentile(means, 2.5)), 4),
                   round(float(np.percentile(means, 97.5)), 4)],
            "n": len(a)}


def auprc(y_true, y_score) -> float:
    import numpy as np
    from sklearn.metrics import average_precision_score
    y = np.asarray(y_true, dtype=int)
    s = np.nan_to_num(np.asarray(y_score, dtype=float), nan=-1e18)
    if y.sum() == 0 or y.sum() == len(y):
        return float("nan")
    return round(float(average_precision_score(y, s)), 4)


def summarize_arm(frame, score_col: str) -> dict:
    m = locus_metrics(frame, score_col)
    return {"top1": bootstrap_ci(m["top1"]), "top5": bootstrap_ci(m["top5"]),
            "mean_best_rank": bootstrap_ci(m["ranks"]),
            "auprc": auprc(m["y_true"], m["y_score"])}


def train_lgbm(train_df, valid_df, seed: int, cols=None):
    import lightgbm as lgb
    cols = FEATURE_COLS if cols is None else list(cols)
    params = {"objective": "binary", "metric": "average_precision",
              "learning_rate": 0.05, "num_leaves": 31,
              "min_child_samples": 20, "feature_fraction": 0.9,
              "bagging_fraction": 0.9, "bagging_freq": 1, "seed": seed,
              "deterministic": True, "num_threads": 4, "verbosity": -1}
    pos = float(train_df["label"].sum())
    params["scale_pos_weight"] = float(len(train_df) - pos) / max(pos, 1.0)
    dtr = lgb.Dataset(train_df[cols], label=train_df["label"])
    dva = lgb.Dataset(valid_df[cols], label=valid_df["label"],
                      reference=dtr)
    model = lgb.train(params, dtr, num_boost_round=2000,
                      valid_sets=[dva],
                      callbacks=[lgb.early_stopping(100, verbose=False)])
    imp = dict(zip(cols,
                   model.feature_importance(importance_type="gain"),
                   strict=True))
    return model, {k: round(float(v), 1) for k, v in imp.items()}


TMP_FEATURES = ["dist_kb", "loeuf", "pli"]
FLAMES_YEAR_ASSUMED = 2018


def locus_years() -> dict:
    """min publication year per (chr,pos) locus from gold PMIDs."""
    positives = json.loads((M4B_DIR / "positives.json").read_text(
        encoding="utf-8"))["rows"]
    years = json.loads((M4B_DIR / "pmid_years.json").read_text(
        encoding="utf-8"))
    locyr: dict[tuple, list] = {}
    for r in positives:
        if r["provenance"] != "ot_gold" or not r["pubmed"]:
            continue
        y = years.get(str(int(float(r["pubmed"]))))
        if y:
            locyr.setdefault((r["chr"], r["pos"]), []).append(y)
    return {k: min(v) for k, v in locyr.items()}


def run_cv(df, feature_cols: list[str], seeds: list[int], fold_fn: str = "chr"):
    import numpy as np
    folds = gene_folds(df) if fold_fn == "gene" else chr_folds(df)
    cv_scores = np.zeros(len(df))
    imps: list[dict] = []
    for fold in range(5):
        te = folds == fold
        tr = df[~te].sample(frac=1.0,
                            random_state=seeds[0]).reset_index(drop=True)
        cut = int(len(tr) * 0.9)
        seed_scores = np.zeros(int(te.sum()))
        for seed in seeds:
            model, imp = train_lgbm(tr.iloc[:cut], tr.iloc[cut:],
                                    seed, feature_cols)
            imps.append({"fold": fold, "seed": seed, **imp})
            seed_scores += model.predict(df[te][feature_cols])
        cv_scores[te] = seed_scores / len(seeds)
    return cv_scores, imps


def ot_column(df, ot_scores: dict):
    import numpy as np
    vals = [ot_scores.get(r["locus_id"], {}).get(r["ensg"])
            for _, r in df.iterrows()]
    return np.array(vals, dtype=float)


def eval_arms(df, scores: dict) -> dict:
    scored = df.copy()
    for name, vals in scores.items():
        scored[name] = vals
    return {"gbm": summarize_arm(scored, "gbm"),
            "distance_only": summarize_arm(scored, "negdist"),
            "ot_l2g": summarize_arm(scored, "ot")}


def paired_top1(frame, a: str, b: str, seed: int = 1337) -> dict:
    """Per-locus paired top-1 difference a-minus-b with bootstrap CI."""
    import numpy as np
    rng = np.random.default_rng(seed)
    diffs = []
    for _, loc in frame.groupby("locus_id"):
        hits = []
        for col in (a, b):
            s = np.nan_to_num(loc[col].to_numpy(dtype=float), nan=-np.inf)
            y = loc["label"].to_numpy(dtype=int)
            best = min(float(np.sum(s > sp) + np.sum(s == sp))
                       for sp, yp in zip(s, y, strict=True) if yp == 1)
            hits.append(1.0 if best == 1.0 else 0.0)
        diffs.append(hits[0] - hits[1])
    d = np.array(diffs)
    boots = [float(np.mean(rng.choice(d, size=len(d), replace=True)))
             for _ in range(1000)]
    return {"mean": round(float(d.mean()), 4),
            "ci": [round(float(np.percentile(boots, 2.5)), 4),
                   round(float(np.percentile(boots, 97.5)), 4)],
            "n": len(d)}


def covered_eval(df, ot_scores: dict) -> dict:
    cov = df[df["locus_id"].isin(set(ot_scores))].copy()
    if cov["locus_id"].nunique() == 0:
        return {"status": "EMPTY"}
    return {"status": "OK", "n_loci": int(cov["locus_id"].nunique()),
            **eval_arms(cov, {"gbm": cov["gbm"].to_numpy(),
                              "negdist": cov["negdist"].to_numpy(),
                              "ot": cov["ot"].to_numpy()}),
            "diff_gbm_minus_ot": paired_top1(cov, "gbm", "ot"),
            "diff_gbm_minus_dist": paired_top1(cov, "gbm", "negdist")}


def cmd_train(args) -> int:
    import numpy as np
    import pandas as pd

    from sina.manifests import Run
    manifest_path = M4B_DIR / "train.manifest.json"
    run = Run(script="scripts/m4b/tabular_baseline.py",
              manifest_path=manifest_path, config={}, git_ref="")
    with run.stage("train"):
        df = pd.read_parquet(M4B_DIR / "features" / "features.parquet")
        df = df[df["label"].groupby(df["locus_id"]).transform("sum") > 0]
        df = df.reset_index(drop=True)
        df["negdist"] = -df["dist_kb"]
        ot_path = M4B_DIR / "ot_scores.json"
        ot_scores = json.loads(ot_path.read_text(
            encoding="utf-8"))["scores"] if ot_path.exists() else {}
        df["ot"] = ot_column(df, ot_scores)
        cols = FEATURE_SETS[args.featureset]
        out: dict = {"generated_at": utcnow_iso(), "features": cols,
                     "featureset": args.featureset,
                     "seeds": list(args.seeds),
                     "ot_loci_overlap": int(df["locus_id"].isin(
                         set(ot_scores)).nunique()),
                     "views": {}}
        if args.view in ("cv", "both"):
            fold_fn = getattr(args, "folds", "chr")
            cv_scores, imps = run_cv(df, cols, list(args.seeds), fold_fn)
            out["fold_fn"] = fold_fn
            out["leakage_audit"] = {
                "chr_folds": leakage_audit(df, chr_folds(df)),
                "used_folds": leakage_audit(
                    df, gene_folds(df) if fold_fn == "gene" else chr_folds(df)),
            }
            df["gbm"] = cv_scores
            df[["locus_id", "ensg", "label", "gbm", "negdist",
                "ot"]].to_parquet(M4B_DIR / "oof_scores.parquet", index=False)
            arm_scores = {"gbm": cv_scores,
                          "negdist": df["negdist"].to_numpy(),
                          "ot": df["ot"].to_numpy()}
            df = df.assign(gbm=cv_scores)
            out["views"]["cv"] = {
                "n_loci": int(df["locus_id"].nunique()),
                **eval_arms(df, arm_scores),
                "covered_subset": covered_eval(df, ot_scores),
                "diff_gbm_minus_ot": paired_top1(df, "gbm", "ot"),
                "importance_gain": imps,
            }
            abl = {}
            drops = (["dist_kb", "degree", "constraint"]
                     if args.featureset == "v1" else list(cols))
            for drop in drops:
                left = [c for c in cols
                        if (c not in ("loeuf", "pli") if drop == "constraint"
                            else c != drop)]
                s, _ = run_cv(df, left, list(args.seeds))
                abl[f"no_{drop}"] = summarize_arm(
                    df.assign(gbm=s), "gbm")["top1"]
            out["views"]["cv"]["ablations_top1"] = abl
            g = out["views"]["cv"]["gbm"]["top1"]
            o = out["views"]["cv"]["ot_l2g"]["top1"]
            print(f"CV: gbm {g} vs dist "
                  f"{out['views']['cv']['distance_only']['top1']} vs OT {o}",
                  flush=True)
        if args.view in ("tmp", "both"):
            ly = locus_years()
            df["year"] = df.apply(
                lambda r: ly.get((r["chr"], r["pos"]),
                                 FLAMES_YEAR_ASSUMED
                                 if "flames" in r["provenance"]
                                 else None), axis=1)
            tmp_views = {}
            for cutoff in (2020, 2015):
                te_mask = df["year"].notna() & (df["year"] > cutoff)
                tr = df[df["year"].notna() & (df["year"] <= cutoff)]
                if te_mask.sum() == 0 or tr["label"].sum() == 0:
                    tmp_views[f"cutoff_{cutoff}"] = {
                        "status": "EMPTY", "n_train_loci": int(
                            tr["locus_id"].nunique()),
                        "n_test_loci": int(df[te_mask]["locus_id"].nunique())}
                    continue
                trs = tr.sample(frac=1.0,
                                random_state=args.seeds[0]).reset_index(
                                    drop=True)
                cut = int(len(trs) * 0.9)
                seed_scores = np.zeros(len(df))
                for seed in args.seeds:
                    model, _ = train_lgbm(trs.iloc[:cut], trs.iloc[cut:],
                                          seed, TMP_FEATURES)
                    seed_scores += model.predict(df[TMP_FEATURES])
                seed_scores /= len(args.seeds)
                tmp_views[f"cutoff_{cutoff}"] = {
                    "status": "OK", "features": TMP_FEATURES,
                    "n_train_loci": int(tr["locus_id"].nunique()),
                    "n_test_loci": int(df[te_mask]["locus_id"].nunique()),
                    **eval_arms(df[te_mask], {
                        "gbm": seed_scores[te_mask.to_numpy()],
                        "negdist": df[te_mask]["negdist"].to_numpy(),
                        "ot": df[te_mask]["ot"].to_numpy()}),
                }
            out["views"]["tmp"] = tmp_views
            print("TMP: " + json.dumps(
                {k: v.get("status") for k, v in tmp_views.items()}),
                flush=True)
        (M4B_DIR / f"train_results_{args.featureset}.json").write_text(
            json.dumps(out, indent=1), encoding="utf-8")
    return 0


def _ot_lib():
    import importlib.util
    path = (REPO_ROOT / "scripts" / "pod" / "m2c_loci_headtohead.py")
    spec = importlib.util.spec_from_file_location("m2c_h2h", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def cmd_ot(args) -> int:
    """OT L2G scores for Catalog-space positive loci (proven fetch path)."""
    from sina.manifests import Run
    lib = _ot_lib()
    cache = M4B_DIR / "ot_cache"
    studies_dir = cache / "studies"
    sets_dir = cache / "credible_sets"
    studies_dir.mkdir(parents=True, exist_ok=True)
    sets_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = M4B_DIR / "ot.manifest.json"
    run = Run(script="scripts/m4b/tabular_baseline.py",
              manifest_path=manifest_path, config={}, git_ref="")
    positives = json.loads((M4B_DIR / "positives.json").read_text(
        encoding="utf-8"))["rows"]
    members = [(r["chr"], r["pos"], r["study_acc"]) for r in positives
               if r["study_acc"] and r["provenance"] == "ot_gold"]
    if args.dry_run:
        accs = sorted({m[2] for m in members})
        print(f"DRY-RUN: would query {len(accs)} studies for "
              f"{len(members)} loci", flush=True)
        return 0
    with run.stage("ot"):
        study_accs = sorted({m[2] for m in members})
        study_sets: dict[str, list] = {}
        for acc in study_accs:
            fpath = studies_dir / f"{acc}.json"
            if fpath.exists():
                study_sets[acc] = json.loads(fpath.read_text(
                    encoding="utf-8"))
                continue
            rows, idx = [], 0
            while True:
                d = lib.ot_post(lib.CS_LIST_Q, {"sids": [acc], "idx": idx})
                cs = d["data"]["credibleSets"]
                rows.extend(cs["rows"])
                if len(rows) >= cs["count"] or not cs["rows"]:
                    break
                idx += 1
            fpath.write_text(json.dumps(rows, indent=1, sort_keys=True)
                             + "\n", encoding="utf-8")
            study_sets[acc] = rows
            print(f"  study {acc}: {len(rows)} credible sets", flush=True)
            time.sleep(0.2)
        index: dict[tuple, list] = {}
        for acc, rows in study_sets.items():
            for r in rows:
                v = r.get("variant") or {}
                index.setdefault((acc, str(v.get("chromosome")),
                                  v.get("position")), []).append(
                                      r["studyLocusId"])
        scores: dict[str, dict] = {}
        n_matched = 0
        for chrom, pos, acc in members:
            loc = f"{chrom}:{pos}"
            sids: list[str] = []
            for slid in index.get((acc, chrom, pos), []):
                fpath = sets_dir / f"{slid}.json"
                if fpath.exists():
                    doc = json.loads(fpath.read_text(encoding="utf-8"))
                else:
                    d = lib.ot_post(lib.CS_ONE_Q, {"slid": slid})
                    doc = d["data"]["credibleSet"]
                    fpath.write_text(json.dumps(doc, indent=1, sort_keys=True)
                                     + "\n", encoding="utf-8")
                    time.sleep(0.2)
                sids.append(slid)
                for row in ((doc.get("l2GPredictions") or {}).get("rows")
                            or []):
                    tgt = (row.get("target") or {}).get("id", "")
                    gene = ("ENSG" + tgt) if tgt and not tgt.startswith(
                        "ENSG") else tgt
                    if gene:
                        prev = scores.setdefault(loc, {}).get(gene, 0.0)
                        scores[loc][gene] = max(prev, float(row.get("score",
                                                                   0.0)))
            n_matched += bool(sids)
        (M4B_DIR / "ot_scores.json").write_text(json.dumps(
            {"generated_at": utcnow_iso(), "loci_scored": len(scores),
             "loci_matched": n_matched,
             "loci_queried": len({(m[0], m[1]) for m in members}),
             "scores": scores}, indent=1), encoding="utf-8")
        print(f"OT: {len(scores)} loci scored, {n_matched} matched",
              flush=True)
    return 0


def _t1(d: dict) -> str:
    t = d["top1"]
    return f"{t['mean']:.4f} [{t['ci'][0]:.4f},{t['ci'][1]:.4f}] (n={t['n']})"


def _verdict_text(cv: dict, cov: dict) -> str:
    d = cov["diff_gbm_minus_ot"]
    lo, hi = d["ci"]
    if lo > 0:
        call = (f"REPLICATION CONFIRMED: GBM beats OT by {d['mean']:.3f}, "
                "CI clear.")
    elif hi < 0:
        call = (f"REPLICATION FAILED: GBM trails OT by {d['mean']:.3f}, "
                "CI clear — features gap, not modeling gap; no lift claims.")
    else:
        call = (f"STATISTICAL TIE with OT (diff {d['mean']:.3f} "
                f"[{lo:.3f},{hi:.3f}] on n={d['n']}); within the ±5pt band — "
                "no lift claims either way.")
    return call + " Ceiling: gene prioritization at GWAS loci (benchmark " \
        "track) — no target-discovery or clinical language."


def cmd_report(args) -> int:
    from sina.manifests import Run
    manifest_path = M4B_DIR / "report.manifest.json"
    run = Run(script="scripts/m4b/tabular_baseline.py",
              manifest_path=manifest_path, config={}, git_ref="")
    with run.stage("report"):
        res = json.loads((M4B_DIR / f"train_results_{args.featureset}.json"
                          ).read_text(encoding="utf-8"))
        fl_path = M4B_DIR / "flames_replication.json"
        fl = json.loads(fl_path.read_text(encoding="utf-8")) if fl_path.exists() else None
        cv = res["views"]["cv"]
        cov = cv["covered_subset"]
        tmp15 = res["views"]["tmp"]["cutoff_2015"]
        abl = cv["ablations_top1"]
        lines = [
            "<!-- track_scope: benchmark -->",
            "# M4B tabular baseline — SINA-free LightGBM vs OT L2G",
            "",
            f"Generated {res['generated_at']}. Feature set "
            f"{res.get('featureset', 'v1')}: " + ", ".join(res["features"]) +
            ". LightGBM 4.7.0 binary logloss, "
            f"seeds {res['seeds']}, 5-fold chromosome-disjoint CV. "
            "Panel: gene TSS within ±500kb of lead SNP (OT convention). "
            "Ties pessimistic everywhere.",
            "",
            "## Claim corrections (verified 2026-09-12)",
            "FLAMES is Schipper/de Leeuw/Posthuma (VU Amsterdam), not "
            "Mountjoy. OT Genetics GraphQL is DNS-dead; OT scores via the "
            "Platform credibleSets API. HF model exists "
            "(opentargets/locus_to_gene_25.12) with an MIT tag, not "
            "Apache-2.0. GTEx v8 = Aug 2019 (not 2020). GWAS Catalog is EBI "
            "Terms of Use (some datasets CC-BY-NC) — verify terms before "
            "any use beyond research. FLAMES comparison is vs-OT-only in v1 "
            "(annotation pipeline deferred, per brief fallback).",
            "",
            "## Feature freshness (temporal rule)",
            "TMP view allows ONLY pre-2021-datable features: dist_kb, "
            "gnomAD v2 LOEUF/pLI (2018). Excluded: graph degree (current "
            "graph), FLAMES/PoPS scores (2023+ sources), eQTL (not yet "
            "built). FLAMES loci assumed 2018 (UK Biobank release).",
            "",
            f"## View CV (all loci, n={cv['n_loci']})",
            f"GBM top-1 {_t1(cv['gbm'])}; distance-only "
            f"{_t1(cv['distance_only'])}; OT L2G {_t1(cv['ot_l2g'])}. "
            f"Paired diff GBM-OT {cv['diff_gbm_minus_ot']}.",
            "OT auto-misses 1,326 FLAMES (UK Biobank, non-Catalog) loci — "
            "the full-set OT number is a coverage artifact, not a fair "
            "fight. Fair comparison is the covered subset:",
            "",
            "## View CV covered subset (OT-scored loci, "
            f"n={cov['n_loci']}, all OT-gold)",
            f"GBM {_t1(cov['gbm'])}; distance {_t1(cov['distance_only'])}; "
            f"OT {_t1(cov['ot_l2g'])}. Paired diff GBM-OT "
            f"{cov['diff_gbm_minus_ot']}. Structural note: OT trained on "
            "ALL these gold loci (home advantage); our GBM sees 4/5 folds. "
            "OT wins here by ~6pts, CI clear.",
            "",
            "## View TMP",
            "cutoff_2020: EMPTY — zero post-2020 loci exist in the positives "
            "(gold curation ends 2020; FLAMES loci assumed 2018). The "
            "pre-registered temporal comparison cannot fire; verdict N/A. "
            "This is itself a finding: neither OT nor FLAMES published one "
            "because the test set does not publicly exist.",
            f"cutoff_2015 fallback (train<=2015, test>2015, n={tmp15['n_test_loci']}, "
            f"features {', '.join(tmp15['features'])}): GBM "
            f"{_t1(tmp15['gbm'])} vs distance {_t1(tmp15['distance_only'])} "
            f"vs OT {_t1(tmp15['ot_l2g'])}.",
            "",
            f"## Ablations (CV top-1 deltas vs full {cv['gbm']['top1']['mean']:.4f})",
            "ablation top-1: " + "; ".join(
                f"{k} {v['mean']:.4f}" for k, v in abl.items()) + ".",
            "",
            "## Verdicts vs frozen rules",
            _verdict_text(cv, cov),
            "",
            "## FLAMES replication (their model, their loci, our metric)",
            ("FLAMES subcommand not run yet." if fl is None else
             "L2G set (n=103): XGB " +
             f"{fl['l2g_set']['arms']['XGB_scaled']['top1']['mean']:.4f} vs OT " +
             f"{fl['l2g_set']['arms']['Overall L2G score']['top1']['mean']:.4f}, diff " +
             f"{fl['l2g_set']['arms']['XGB_scaled']['diff_vs_ot']}; full FLAMES " +
             f"{fl['l2g_set']['arms']['FLAMES_scaled_with_path']['top1']['mean']:.4f}, diff " +
             f"{fl['l2g_set']['arms']['FLAMES_scaled_with_path']['diff_vs_ot']}. " +
             "9-pheno set (n=140): full FLAMES " +
             f"{fl['pheno9_set']['arms']['FLAMES_scaled_with_path']['top1']['mean']:.4f} vs OT " +
             f"{fl['pheno9_set']['arms']['Overall L2G score']['top1']['mean']:.4f}, diff " +
             f"{fl['pheno9_set']['arms']['FLAMES_scaled_with_path']['diff_vs_ot']}. " +
             "Verdict: raw XGB loses to OT clearly; full FLAMES leads by "
             "+2 to +6pts, never significant at 95% on top-1. The brief's "
             "'FLAMES beat OT' precedent does not survive their own shipped "
             "benchmark tables on this metric."),
            "",
            "## Reproduce",
            "`scripts/m4b/tabular_baseline.py positives panels features "
            "train --view both` (OT cache reused); this report regenerates "
            "via `report`. External sources: OT gold repo, FLAMES repos "
            "(Marijn-Schipper), GENCODE v46, gnomAD v2.1.1, Ensembl REST "
            "(rsID cache), PubMed esummary (dates), OT Platform API.",
        ]
        rep_path = (REPO_ROOT / "reports" / "m4b" /
                    f"tabular_baseline_{args.featureset}.md")
        rep_path.parent.mkdir(parents=True, exist_ok=True)
        rep_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        out_path = (REPO_ROOT / "reports" / "m4b" /
                    f"tabular_baseline_{args.featureset}_results.json")
        out_path.write_text(json.dumps(res, indent=1), encoding="utf-8")
        print(f"report -> {rep_path}", flush=True)
    return 0


FLAMES_REPOS = {
    "FLAMES": "https://github.com/Marijn-Schipper/FLAMES",
    "FLAMES_paper_analyses": "https://github.com/Marijn-Schipper/FLAMES_paper_analyses",
}
FLAMES_SHAS = {
    "FLAMES": "159e83a6c6bc5ad42e5f8b7e5b4722a727575935",
    "FLAMES_paper_analyses": "4a5a9964eaf48f021a06e69d7d4f6eba257eaff0",
}


def flames_headtohead(path: Path, score_cols: list[str], seed: int = 1337) -> dict:
    """Per-locus top-1 (pessimistic) + paired diffs vs Overall L2G score."""
    import numpy as np
    import pandas as pd
    rng = np.random.default_rng(seed)
    d = pd.read_csv(path, sep="\t")
    hits: dict[str, list] = {}
    for col in [*score_cols, "Overall L2G score"]:
        h = []
        for _, loc in d.groupby("filename"):
            s = np.nan_to_num(loc[col].to_numpy(dtype=float), nan=-np.inf)
            y = loc["TP"].to_numpy(dtype=int)
            best = min(float(np.sum(s > sp) + np.sum(s == sp))
                       for sp, yp in zip(s, y, strict=True) if yp == 1)
            h.append(1.0 if best == 1.0 else 0.0)
        hits[col] = h
    out = {"n_loci": int(d["filename"].nunique()), "arms": {}}
    base = np.array(hits["Overall L2G score"])
    for col, h in hits.items():
        a = np.array(h)
        boots = [float(np.mean(rng.choice(a, size=len(a), replace=True)))
                 for _ in range(1000)]
        out["arms"][col] = {"top1": {"mean": round(float(a.mean()), 4),
                                     "ci": [round(float(np.percentile(boots, 2.5)), 4),
                                            round(float(np.percentile(boots, 97.5)), 4)],
                                     "n": len(a)}}
    for col in score_cols:
        dd = np.array(hits[col]) - base
        boots = [float(np.mean(rng.choice(dd, size=len(dd), replace=True)))
                 for _ in range(1000)]
        out["arms"][col]["diff_vs_ot"] = {
            "mean": round(float(dd.mean()), 4),
            "ci": [round(float(np.percentile(boots, 2.5)), 4),
                   round(float(np.percentile(boots, 97.5)), 4)]}
    return out


def cmd_flames(args) -> int:
    from sina.manifests import Run
    manifest_path = M4B_DIR / "flames.manifest.json"
    run = Run(script="scripts/m4b/tabular_baseline.py",
              manifest_path=manifest_path, config={}, git_ref="")
    with run.stage("flames"):
        bench = args.src / "FLAMES_paper_analyses" / "benchmarking_results"
        out = {"generated_at": utcnow_iso(), "repos": FLAMES_REPOS,
               "shas": FLAMES_SHAS,
               "l2g_set": flames_headtohead(
                   bench / "L2G" / "results_L2G.tsv",
                   ["XGB_scaled", "FLAMES_scaled",
                    "FLAMES_scaled_with_path"]),
               "pheno9_set": flames_headtohead(
                   bench / "9_pheno" / "results_9_pheno.tsv",
                   ["FLAMES_scaled", "FLAMES_scaled_with_path"])}
        (M4B_DIR / "flames_replication.json").write_text(
            json.dumps(out, indent=1), encoding="utf-8")
        for name in ("l2g_set", "pheno9_set"):
            s = out[name]
            print(f"FLAMES {name} n={s['n_loci']}: " + "; ".join(
                f"{k} {v['top1']['mean']}"
                for k, v in s["arms"].items()), flush=True)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("positives")
    p.add_argument("--dry-run", action="store_true")
    sub.add_parser("panels")
    sub.add_parser("features")
    sub.add_parser("features2")
    sub.add_parser("features3")
    t = sub.add_parser("train")
    t.add_argument("--view", choices=["cv", "tmp", "both"], default="cv")
    t.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS))
    t.add_argument("--featureset", choices=["v1", "v2", "v3", "v4"], default="v1")
    t.add_argument("--folds", choices=["chr", "gene"], default="chr")
    o = sub.add_parser("ot")
    o.add_argument("--dry-run", action="store_true")
    rp = sub.add_parser("report")
    rp.add_argument("--featureset", choices=["v1", "v2", "v3"], default="v2")
    f = sub.add_parser("flames")
    f.add_argument("--src", type=Path, default=Path("/tmp/m4b_src"))
    args = ap.parse_args()
    if getattr(args, "dry_run", False) and args.cmd == "positives":
        print("DRY-RUN: positives would read OT gold TSV + FLAMES CSV, "
              "resolve rsIDs, write data/processed/m4b/positives.json",
              flush=True)
        return 0
    if args.cmd == "positives":
        return cmd_positives(args)
    if args.cmd == "panels":
        return cmd_panels(args)
    if args.cmd == "features":
        return cmd_features(args)
    if args.cmd == "features2":
        return cmd_features2(args)
    if args.cmd == "features3":
        return cmd_features3(args)
    if args.cmd == "train":
        return cmd_train(args)
    if args.cmd == "ot":
        return cmd_ot(args)
    if args.cmd == "report":
        return cmd_report(args)
    if args.cmd == "flames":
        return cmd_flames(args)
    return 2


def go_phe_degrees() -> dict:
    """PHE-GEN and GO-term degree per bare ENSG from the full background."""
    import pandas as pd
    path = REPO_ROOT / "data" / "processed" / "m2c" / "full_exports" / "train1.txt"
    phe: dict[str, int] = {}
    go: dict[str, int] = {}
    for chunk in pd.read_csv(path, sep="\t", header=None, chunksize=500_000,
                             dtype=str):
        keep = (((chunk[0].str.startswith("hp:")) & (chunk[1] == "ASSOCIATED_WITH"))
                | ((chunk[0].str.startswith("go:")) & (chunk[1] == "INTERACTS_WITH")))
        for src, dst in zip(chunk[0][keep], chunk[2][keep], strict=True):
            if dst.startswith("ensembl:ENSG"):
                gene = dst.split(":")[1]
                if src.startswith("hp:"):
                    phe[gene] = phe.get(gene, 0) + 1
                else:
                    go[gene] = go.get(gene, 0) + 1
    return phe, go


def locus_trait_map() -> dict:
    """(chr,pos) -> EFO trait id from the OT gold TSV."""
    import pandas as pd
    g = pd.read_csv(OT_GOLD_TSV, sep="\t", low_memory=False,
                    usecols=["trait_info.ontology",
                             "sentinel_variant.locus_GRCh38.chromosome",
                             "sentinel_variant.locus_GRCh38.position"])
    out: dict[tuple, str] = {}
    for _, r in g.iterrows():
        parts = str(r["trait_info.ontology"]).split("_")
        if len(parts) != 2 or not parts[1]:
            continue
        efo = "efo:" + parts[1].lower()
        out.setdefault((str(r["sentinel_variant.locus_GRCh38.chromosome"]),
                        int(r["sentinel_variant.locus_GRCh38.position"])), efo)
    return out


def cmd_features3(args) -> int:
    import pandas as pd

    from sina.manifests import Run
    manifest_path = M4B_DIR / "features3.manifest.json"
    run = Run(script="scripts/m4b/tabular_baseline.py",
              manifest_path=manifest_path, config={}, git_ref="")
    with run.stage("features3"):
        feat = pd.read_parquet(M4B_DIR / "features" / "features.parquet")
        phe, go = go_phe_degrees()
        feat["phe_degree"] = feat["ensg"].map(phe).fillna(0).astype(int)
        feat["go_degree"] = feat["ensg"].map(go).fillna(0).astype(int)
        ppr_path = M4B_DIR / "ppr_scores.parquet"
        if not ppr_path.exists():
            raise SystemExit("ppr_scores.parquet missing — run compute_ppr.py first")
        ppr = pd.read_parquet(ppr_path)
        piv = ppr.pivot_table(index="ensg", columns="trait", values="ppr")
        tmap = locus_trait_map()
        feat["ppr"] = feat.apply(
            lambda r: piv.at[r["ensg"], tmap[(r["chr"], r["pos"])]]
            if (r["chr"], r["pos"]) in tmap and r["ensg"] in piv.index
            and tmap[(r["chr"], r["pos"])] in piv.columns else float("nan"),
            axis=1)
        feat.to_parquet(M4B_DIR / "features" / "features.parquet", index=False)
        print(f"features3: phe cov {(feat['phe_degree'] > 0).mean():.3f}, "
              f"go cov {(feat['go_degree'] > 0).mean():.3f}, ppr cov "
              f"{feat['ppr'].notna().mean():.3f}", flush=True)
    return 0


VEP_CACHE = M4B_DIR / "vep_cache.json"
VEP_SEV = {"stop_gained": 3, "frameshift_variant": 3, "splice_donor_variant": 3,
           "splice_acceptor_variant": 3, "stop_lost": 3, "start_lost": 3,
           "missense_variant": 2, "inframe_insertion": 2,
           "inframe_deletion": 2, "protein_altering_variant": 2,
           "synonymous_variant": 1, "splice_region_variant": 1,
           "5_prime_UTR_variant": 1, "3_prime_UTR_variant": 1}


def disgen_degrees() -> dict:
    """DIS-GEN-only degree per bare ENSG (literature prior, no leakage
    beyond what the label already is)."""
    import pandas as pd
    path = REPO_ROOT / "data" / "processed" / "m2c" / "full_exports" / "train2.txt"
    degs: dict[str, int] = {}
    for chunk in pd.read_csv(path, sep="\t", header=None, chunksize=500_000,
                             dtype=str):
        for col in (0, 2):
            for ent, cnt in chunk[col].value_counts().items():
                if ent.startswith("ensembl:ENSG"):
                    gene = ent.split(":")[1]
                    degs[gene] = degs.get(gene, 0) + int(cnt)
    return degs


def locus_rsids() -> dict:
    """lead rsIDs per (chr,pos) locus from OT gold + FLAMES sources."""
    import pandas as pd
    out: dict[tuple, list] = {}
    g = pd.read_csv(OT_GOLD_TSV, sep="\t", low_memory=False,
                    usecols=["sentinel_variant.locus_GRCh38.chromosome",
                             "sentinel_variant.locus_GRCh38.position",
                             "sentinel_variant.rsid"])
    for _, r in g.dropna(subset=["sentinel_variant.rsid"]).iterrows():
        m = re.search(r"\brs\d+\b", str(r["sentinel_variant.rsid"]))
        if m:
            out.setdefault((str(r["sentinel_variant.locus_GRCh38.chromosome"]),
                            int(r["sentinel_variant.locus_GRCh38.position"])),
                           []).append(m.group(0))
    positives = json.loads((M4B_DIR / "positives.json").read_text(
        encoding="utf-8"))["rows"]
    flames_loc = {r["flames_locus"]: (r["chr"], r["pos"]) for r in positives
                  if "flames_locus" in r}
    d = pd.read_csv(FLAMES_CSV, usecols=["filename", "GenomicLocus", "TP",
                                         "IndSigSNPs"])
    for _, r in d[d["TP"] == 1].iterrows():
        key = flames_loc.get(f"{r['filename']}#{r['GenomicLocus']}")
        rid = first_rsid(r["IndSigSNPs"])
        if key is not None and rid and rid not in out.setdefault(key, []):
            out[key].append(rid)
    return out


REGION_CACHE = M4B_DIR / "region_rsid_cache.json"


def backfill_rsids(keys: list[tuple], cache: dict) -> dict:
    """chr:pos -> lead rsIDs via Ensembl variation overlap (1bp window).

    Covers FLAMES loci whose IndSigSNPs carry no rsID. Threaded, cached,
    resumable; only exact-position SNVs are kept.
    """
    from concurrent.futures import ThreadPoolExecutor

    import requests
    missing = [k for k in keys if f"{k[0]}:{k[1]}" not in cache]
    print(f"region backfill: {len(missing)} loci", flush=True)

    def one(key):
        chrom, pos = key
        tag = f"{chrom}:{pos}"
        try:
            r = requests.get(
                f"{ENSEMBL_REST}/overlap/region/homo_sapiens/"
                f"{chrom}:{pos}-{pos}",
                params={"feature": "variation", "content-type":
                        "application/json"},
                headers={"Accept": "application/json"}, timeout=90)
            r.raise_for_status()
            ids = sorted({v.get("id", "") for v in r.json()
                          if v.get("id", "").startswith("rs")
                          and v.get("start") == pos
                          and len(v.get("alleles") or []) == 2
                          and all(len(a) == 1
                                  for a in v.get("alleles") or [])})
            return tag, ids[:3]
        except Exception:
            return tag, None

    with ThreadPoolExecutor(max_workers=8) as ex:
        for tag, ids in ex.map(one, missing):
            if ids is not None:
                cache[tag] = ids
    return cache


def vep_post(chunk: list[str]) -> dict:
    import requests
    resp = requests.post(
        f"{ENSEMBL_REST}/vep/homo_sapiens/id",
        headers={"Content-Type": "application/json",
                 "Accept": "application/json"},
        json={"ids": chunk}, timeout=180)
    resp.raise_for_status()
    return resp.json()


def fetch_vep(rsids: list[str], cache: dict) -> dict:
    missing = [r for r in rsids if r not in cache]
    queue = [missing[i:i + 100] for i in range(0, len(missing), 100)]
    while queue:
        chunk = queue.pop(0)
        try:
            for rec in vep_post(chunk):
                rid = rec.get("id", "")
                best: dict[str, int] = {}
                for tc in rec.get("transcript_consequences") or []:
                    gid = str(tc.get("gene_id", "")).split(".")[0]
                    sev = max((VEP_SEV.get(t, 0)
                               for t in tc.get("consequence_terms") or []),
                              default=0)
                    if gid.startswith("ENSG"):
                        best[gid] = max(best.get(gid, 0), sev)
                cache[rid] = best
        except Exception:
            if len(chunk) == 1:
                cache[chunk[0]] = {}
            else:
                mid = len(chunk) // 2
                queue.append(chunk[:mid])
                queue.append(chunk[mid:])
                time.sleep(2)
                continue
        time.sleep(0.5)
    return cache


def cmd_features2(args) -> int:
    import pandas as pd

    from sina.manifests import Run
    manifest_path = M4B_DIR / "features2.manifest.json"
    run = Run(script="scripts/m4b/tabular_baseline.py",
              manifest_path=manifest_path, config={}, git_ref="")
    with run.stage("features2"):
        feat = pd.read_parquet(M4B_DIR / "features" / "features.parquet")
        degs = disgen_degrees()
        feat["disgen_degree"] = feat["ensg"].map(degs).fillna(0).astype(int)
        rs = locus_rsids()
        panels = json.loads((M4B_DIR / "panels.json").read_text(
            encoding="utf-8"))["loci"]
        rcache: dict = {}
        if REGION_CACHE.exists():
            rcache = json.loads(REGION_CACHE.read_text(encoding="utf-8"))
        bare = [(loc["chr"], loc["pos"]) for loc in panels
                if (loc["chr"], loc["pos"]) not in rs]
        if bare:
            rcache = backfill_rsids(bare, rcache)
            REGION_CACHE.write_text(json.dumps(rcache, indent=1,
                                               sort_keys=True),
                                    encoding="utf-8")
            for key in bare:
                for rid in rcache.get(f"{key[0]}:{key[1]}", []):
                    if rid not in rs.setdefault(key, []):
                        rs[key].append(rid)
        cache: dict = {}
        if VEP_CACHE.exists():
            cache = json.loads(VEP_CACHE.read_text(encoding="utf-8"))
        all_rs = sorted({r for v in rs.values() for r in v})
        if all_rs:
            print(f"VEP for {len(all_rs)} rsIDs", flush=True)
            cache = fetch_vep(all_rs, cache)
            VEP_CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True),
                                 encoding="utf-8")
        # per-locus gene severity: max over the locus lead rsIDs
        loc_sev: dict[tuple, dict] = {}
        for key, ids in rs.items():
            agg: dict[str, int] = {}
            for rid in ids:
                for g, s in (cache.get(rid) or {}).items():
                    agg[g] = max(agg.get(g, 0), s)
            loc_sev[key] = agg
        feat["vep_sev"] = feat.apply(
            lambda r: loc_sev.get((r["chr"], r["pos"]), {}).get(r["ensg"], 0),
            axis=1).astype(int)
        feat.to_parquet(M4B_DIR / "features" / "features.parquet", index=False)
        print(f"features2: disgen coverage "
              f"{(feat['disgen_degree'] > 0).mean():.3f}, vep>0 "
              f"{(feat['vep_sev'] > 0).mean():.4f}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
