"""Track H: regulatory-element overlap features per (locus, gene).

Pre-registered definitions (fixed before seeing any v5 number):
- Source: ENCODE Registry V3 cCREs, GRCh38 (bulk download, cached).
- Credible variants: OT fine-mapping rows with posteriorProbability > 0.01
  from the locus's sets (cached).
- Window: variant within 100kb of gene TSS AND overlapping a cCRE.
- reg_ccre_count: # such variants.
- reg_pls_flag: 1.0 if any overlapping cCRE is promoter-like (PLS).
- reg_any: 1.0 if reg_ccre_count > 0.
"""

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
M4B_DIR = REPO_ROOT / "data" / "processed" / "m4b"
M5B_DIR = REPO_ROOT / "m5b"
WINDOW = 100_000
PP_MIN = 0.01


def load_ccres(path: Path) -> dict:
    """chr -> sorted list of (start, end, is_pls)."""
    out: dict = {}
    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            try:
                s, e = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            chrom = parts[0].removeprefix("chr")
            is_pls = "PLS" in parts[5].split(",")
            out.setdefault(chrom, []).append((s, e, is_pls))
    for v in out.values():
        v.sort(key=lambda t: t[0])
    return out


def overlap_flags(positions: list, elements: list) -> tuple:
    """Pure overlap math: positions in any element? in a PLS element?
    Returns (count_any, any_pls). Unit-tested."""
    n_any = 0
    any_pls = False
    for p in positions:
        for s, e, is_pls in elements:
            if s > p:
                break
            if s <= p <= e:
                n_any += 1
                any_pls = any_pls or is_pls
                break
    return n_any, any_pls


def main() -> int:
    import pandas as pd

    ccre = load_ccres(M5B_DIR / "cache" / "GRCh38-cCREs.bed")
    loci = json.loads((M4B_DIR / "panels.json").read_text(
        encoding="utf-8"))["loci"]
    panels = {loc["locus_id"]: (str(loc["chr"]), loc["pos"])
              for loc in loci}
    pos2loc = {(c, p): lid for lid, (c, p) in panels.items()}
    tss_of = {}
    for loc in loci:
        for g in loc["genes"]:
            tss_of.setdefault((str(loc["chr"]), g["ensg"]), g["tss"])
    file_locus = {}
    for f in (M4B_DIR / "ot_cache" / "credible_sets").glob("*.json"):
        d = json.loads(f.read_text(encoding="utf-8"))
        v = d["variant"]
        lid = pos2loc.get((str(v["chromosome"]), v["position"]))
        if lid:
            file_locus[f.stem] = lid
    cred_by_locus: dict = {}
    for f in (M5B_DIR / "cache" / "finemap").glob("*.json"):
        if f.name == "manifest.json":
            continue
        lid = file_locus.get(f.stem)
        if lid is None:
            continue
        bucket = cred_by_locus.setdefault(lid, [])
        seen = {p for p, _ in bucket}
        for row in json.loads(f.read_text(encoding="utf-8"))["rows"]:
            v = row["variant"]
            if row["posteriorProbability"] <= PP_MIN:
                continue
            if str(v["chromosome"]) != panels[lid][0]:
                continue
            if v["position"] not in seen:
                bucket.append((v["position"], row["posteriorProbability"]))
                seen.add(v["position"])
    feat = pd.read_parquet(M4B_DIR / "features" / "features.parquet",
                           columns=["locus_id", "chr", "ensg"])
    rows = []
    for _, r in feat.drop_duplicates(["locus_id", "ensg"]).iterrows():
        tss = tss_of.get((str(r["chr"]), r["ensg"]))
        n_any, any_pls = 0, False
        if tss is not None:
            pos = [p for p, _ in cred_by_locus.get(r["locus_id"], [])
                   if abs(p - tss) <= WINDOW]
            n_any, any_pls = overlap_flags(
                pos, ccre.get(str(r["chr"]), []))
        rows.append({"locus_id": r["locus_id"], "ensg": r["ensg"],
                     "reg_ccre_count": n_any,
                     "reg_pls_flag": float(any_pls),
                     "reg_any": float(n_any > 0)})
    out = pd.DataFrame(rows)
    out.to_parquet(M5B_DIR / "cache" / "trackh.parquet", index=False)
    mani = {"n_rows": len(out),
            "frac_rows_with_signal": round(float((out["reg_any"] > 0).mean()),
                                          4)}
    (M5B_DIR / "reports" / "trackh.manifest.json").write_text(
        json.dumps(mani, indent=1), encoding="utf-8")
    print(f"trackh: {len(out)} rows, signal frac "
          f"{mani['frac_rows_with_signal']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
