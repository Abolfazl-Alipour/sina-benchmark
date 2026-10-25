"""Phase 3: pre-registered v5 runs (v4 base + track features).

Arms are fixed here before any v5 number is seen:
  base : dist_kb, vep_sev, loeuf, pli
  +F   : + fm_max_pp, fm_sum_pp, fm_n_vars, fm_lead_in_window
  +C   : + cv_gene_count, cv_gene_maxsub, cv_win_count, cv_win_maxsub
  +H   : + reg_ccre_count, reg_pls_flag, reg_any
  +E   : + eq_best_pp, eq_n_tissues, eq_any (skipped if tracke missing)
  all  : base + every available track
Deciding readout: covered PR AUC diff (arm minus OT) with locus-blocked
bootstrap CI. Top-1 paired diff is secondary. Both fold schemes reported.
"""

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts" / "m4b"))
from tabular_baseline import (  # noqa: E402
    auprc,
    ot_column,
    paired_top1,
    run_cv,
    summarize_arm,
)

M4B_DIR = REPO_ROOT / "data" / "processed" / "m4b"
M5B_DIR = REPO_ROOT / "m5b"
SEEDS = [1337, 1339, 5027]
N_BOOT = 1000

BASE = ["dist_kb", "vep_sev", "loeuf", "pli"]
TRACKS = {
    "F": ("trackf", ["fm_max_pp", "fm_sum_pp", "fm_n_vars",
                     "fm_lead_in_window"]),
    "C": ("trackc", ["cv_gene_count", "cv_gene_maxsub", "cv_win_count",
                     "cv_win_maxsub"]),
    "H": ("trackh", ["reg_ccre_count", "reg_pls_flag", "reg_any"]),
    "E": ("tracke", ["eq_best_pp", "eq_n_tissues", "eq_any"]),
}


def prauc_diff_ci(df, a: str, b: str, n_boot: int = N_BOOT,
                  seed: int = 1337) -> dict:
    """Locus-blocked bootstrap CI on pooled PR AUC difference a-minus-b."""
    import numpy as np

    groups = [(loc, g["label"].to_numpy(dtype=int),
               (g[a].to_numpy(dtype=float), g[b].to_numpy(dtype=float)))
              for loc, g in df.groupby("locus_id")]
    ids = list(range(len(groups)))
    rng = np.random.default_rng(seed)
    diffs = []
    for _ in range(n_boot):
        sample = rng.choice(ids, size=len(ids), replace=True)
        yt, sa, sb = [], [], []
        for i in sample:
            _, y, (xa, xb) = groups[i]
            yt.extend(y.tolist())
            sa.extend(np.nan_to_num(xa, nan=-1e18).tolist())
            sb.extend(np.nan_to_num(xb, nan=-1e18).tolist())
        yarr = np.asarray(yt, dtype=int)
        diffs.append(auprc(yarr, np.asarray(sa)) - auprc(yarr, np.asarray(sb)))
    d = np.array(diffs)
    return {"mean": round(float(np.mean(d)), 4),
            "ci": [round(float(np.percentile(d, 2.5)), 4),
                   round(float(np.percentile(d, 97.5)), 4)],
            "n": len(groups)}


def load_matrix():
    import pandas as pd

    df = pd.read_parquet(M4B_DIR / "features" / "features.parquet")
    df = df[df["label"].groupby(df["locus_id"]).transform("sum") > 0]
    df = df.reset_index(drop=True)
    df["negdist"] = -df["dist_kb"]
    ot_scores = json.loads((M4B_DIR / "ot_scores.json").read_text(
        encoding="utf-8"))["scores"]
    df["ot"] = ot_column(df, ot_scores)
    avail = {}
    for tag, (stem, cols) in TRACKS.items():
        p = M5B_DIR / "cache" / f"{stem}.parquet"
        if p.exists():
            t = pd.read_parquet(p)
            df = df.merge(t, on=["locus_id", "ensg"], how="left")
            for c in cols:
                df[c] = df[c].fillna(0.0)
            avail[tag] = cols
    return df, ot_scores, avail


def run_arm(df, cols, ot_scores, fold_fn):
    scores, imps = run_cv(df, cols, SEEDS, fold_fn)
    df = df.assign(gbm=scores)
    full = {"gbm": summarize_arm(df, "gbm"),
            "distance_only": summarize_arm(df, "negdist"),
            "ot_l2g": summarize_arm(df, "ot"),
            "prauc_diff_gbm_minus_ot": prauc_diff_ci(df, "gbm", "ot"),
            "prauc_diff_gbm_minus_dist": prauc_diff_ci(df, "gbm", "negdist"),
            "top1_diff_gbm_minus_ot": paired_top1(df, "gbm", "ot")}
    cov = df[df["locus_id"].isin(set(ot_scores))].copy()
    covered = {"n_loci": int(cov["locus_id"].nunique()),
               "gbm": summarize_arm(cov, "gbm"),
               "distance_only": summarize_arm(cov, "negdist"),
               "ot_l2g": summarize_arm(cov, "ot"),
               "prauc_diff_gbm_minus_ot": prauc_diff_ci(cov, "gbm", "ot"),
               "prauc_diff_gbm_minus_dist": prauc_diff_ci(
                   cov, "gbm", "negdist"),
               "top1_diff_gbm_minus_ot": paired_top1(cov, "gbm", "ot")}
    return full, covered, imps


def main() -> int:
    df, ot_scores, avail = load_matrix()
    arms = {"base": list(BASE)}
    for tag in ("F", "C", "H", "E"):
        if tag in avail:
            arms[f"base+{tag}"] = BASE + avail[tag]
    arms["all"] = BASE + [c for t in avail.values() for c in t]
    out = {"seeds": SEEDS, "tracks_available": sorted(avail),
           "arms": arms, "runs": {}}
    for fold_fn in ("chr", "gene"):
        for name, cols in arms.items():
            full, covered, imps = run_arm(df, cols, ot_scores, fold_fn)
            out["runs"][f"{fold_fn}:{name}"] = {
                "cols": cols, "full": full, "covered": covered,
                "importance_gain": imps,
                "n_loci_full": int(df["locus_id"].nunique())}
            c = covered["prauc_diff_gbm_minus_ot"]
            print(f"{fold_fn}:{name} covered PR AUC diff "
                  f"{c['mean']} {c['ci']} | top1 diff "
                  f"{covered['top1_diff_gbm_minus_ot']['mean']} "
                  f"{covered['top1_diff_gbm_minus_ot']['ci']}", flush=True)
    (M5B_DIR / "reports" / "v5_results.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
