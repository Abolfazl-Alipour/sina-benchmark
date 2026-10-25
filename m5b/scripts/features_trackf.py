"""Track F: fine-mapping overlap features per (locus, gene).

Pre-registered definitions (fixed before seeing any v5 number):
- Window: credible-set variant position within 100kb of the gene TSS.
- fm_max_pp: max posteriorProbability among in-window variants, pooled
  over all credible sets of the locus (0.0 when none).
- fm_sum_pp: sum of in-window posteriorProbabilities.
- fm_n_vars: count of in-window variants with posteriorProbability > 0.01.
- fm_lead_in_window: 1.0 if any set's lead variant falls in-window.

Inputs are OT fine-mapping variant lists (public evidence), never OT's
gene scores. Output: m5b/cache/trackf.parquet + committed manifest.
"""

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

M4B_DIR = REPO_ROOT / "data" / "processed" / "m4b"
M5B_DIR = REPO_ROOT / "m5b"
WINDOW = 100_000
PP_MIN = 0.01


def overlap_features(variants: list, tss: int) -> dict:
    """Pure overlap math over (position, posterior) pairs. Unit-tested."""
    in_win = [(p, q) for p, q in variants if abs(p - tss) <= WINDOW]
    sig = [q for _, q in in_win if q > PP_MIN]
    return {"fm_max_pp": round(max([q for _, q in in_win] + [0.0]), 6),
            "fm_sum_pp": round(sum(q for _, q in in_win), 6),
            "fm_n_vars": len(sig)}


def main() -> int:
    import pandas as pd

    feat = pd.read_parquet(M4B_DIR / "features" / "features.parquet",
                           columns=["locus_id", "chr", "ensg"])
    loci = json.loads((M4B_DIR / "panels.json").read_text(
        encoding="utf-8"))["loci"]
    panels = {loc["locus_id"]: (loc["chr"], loc["pos"]) for loc in loci}
    tss_of = {}
    for loc in loci:
        for g in loc["genes"]:
            tss_of.setdefault((str(loc["chr"]), g["ensg"]), g["tss"])
    sets_by_locus: dict = {}
    for f in (M5B_DIR / "cache" / "finemap").glob("*.json"):
        if f.name == "manifest.json":
            continue
        sets_by_locus[f.stem] = json.loads(f.read_text(encoding="utf-8"))
    # map set -> locus via lead-variant (chr,pos) == panel (chr,pos)
    lead_of = {}
    for f in (M4B_DIR / "ot_cache" / "credible_sets").glob("*.json"):
        d = json.loads(f.read_text(encoding="utf-8"))
        v = d["variant"]
        lead_of[f.stem] = (str(v["chromosome"]), v["position"], v["id"])
    pos2loc = {(c, p): lid for lid, (c, p) in panels.items()}
    locus_sets: dict = {}
    for stem, (c, p, _vid) in lead_of.items():
        lid = pos2loc.get((c, p))
        if lid:
            locus_sets.setdefault(lid, []).append(stem)

    rows = []
    for _, r in feat.drop_duplicates(["locus_id", "ensg"]).iterrows():
        tss = tss_of.get((str(r["chr"]), r["ensg"]))
        agg = {"fm_max_pp": 0.0, "fm_sum_pp": 0.0, "fm_n_vars": 0,
               "fm_lead_in_window": 0.0}
        if tss is not None:
            variants, leads = [], []
            for stem in locus_sets.get(r["locus_id"], []):
                d = sets_by_locus.get(stem)
                if not d:
                    continue
                for row in d["rows"]:
                    v = row["variant"]
                    if str(v["chromosome"]) == str(r["chr"]):
                        variants.append((v["position"],
                                         row["posteriorProbability"]))
                leads.append(lead_of[stem])
            agg = overlap_features(variants, tss)
            agg["fm_lead_in_window"] = float(any(
                str(lv[0]) == str(r["chr"]) and abs(lv[1] - tss) <= WINDOW
                for lv in leads))
        rows.append({"locus_id": r["locus_id"], "ensg": r["ensg"], **agg})
    out = pd.DataFrame(rows)
    (M5B_DIR / "cache").mkdir(parents=True, exist_ok=True)
    out.to_parquet(M5B_DIR / "cache" / "trackf.parquet", index=False)
    mani = {"n_loci_sets": len(locus_sets),
            "n_rows": len(out),
            "frac_rows_with_signal": round(float((out["fm_max_pp"] > 0).mean()),
                                          4),
            "window": WINDOW, "pp_min": PP_MIN}
    (M5B_DIR / "reports" / "trackf.manifest.json").write_text(
        json.dumps(mani, indent=1), encoding="utf-8")
    print(f"trackf: {len(out)} rows, signal frac {mani['frac_rows_with_signal']}",
          flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
