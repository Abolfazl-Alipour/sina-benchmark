"""Track C: ClinVar pathogenicity features per (locus, gene).

Pre-registered definitions (fixed before seeing any v5 number):
- Source: ClinVar variant_summary (bulk download, cached), GRCh38 rows only.
- Pathogenic set: ClinicalSignificance in {Pathogenic, Likely pathogenic,
  Pathogenic/Likely pathogenic}.
- cv_gene_count: # P/LP variants attributed to this gene symbol.
- cv_gene_maxsub: max NumberSubmitters among them (evidence weight).
- cv_win_count: # P/LP variants positionally within 100kb of gene TSS.
- cv_win_maxsub: max NumberSubmitters among those.
"""

import gzip
import json
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
M4B_DIR = REPO_ROOT / "data" / "processed" / "m4b"
M5B_DIR = REPO_ROOT / "m5b"
WINDOW = 100_000
PATHO = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
USECOLS = ["GeneSymbol", "ClinicalSignificance", "Assembly", "Chromosome",
           "Start", "Stop", "NumberSubmitters"]


def parse_pathogenic(src: Path) -> pd.DataFrame:
    """Stream ClinVar summary -> GRCh38 pathogenic variant table."""
    chunks = []
    with gzip.open(src, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = [header.index(c) for c in USECOLS]
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx):
                continue
            if parts[idx[2]] != "GRCh38":
                continue
            if parts[idx[1]] not in PATHO:
                continue
            try:
                nsub = int(parts[idx[6]])
                start, stop = int(parts[idx[4]]), int(parts[idx[5]])
            except ValueError:
                continue
            chunks.append((parts[idx[0]], parts[idx[3]], start, stop, nsub))
    return pd.DataFrame(chunks, columns=["symbol", "chr", "start", "stop",
                                         "nsub"])


def gene_features(patho: pd.DataFrame, symbol: str, chrom: str,
                  tss: int) -> dict:
    """Pure per-gene aggregation. Unit-tested."""
    g = patho[patho["symbol"] == symbol]
    w = patho[(patho["chr"] == chrom) & (patho["start"] <= tss + WINDOW)
              & (patho["stop"] >= tss - WINDOW)]
    return {"cv_gene_count": len(g),
            "cv_gene_maxsub": int(g["nsub"].max()) if len(g) else 0,
            "cv_win_count": len(w),
            "cv_win_maxsub": int(w["nsub"].max()) if len(w) else 0}


def main() -> int:
    patho = parse_pathogenic(M5B_DIR / "cache" / "variant_summary.txt.gz")
    patho.to_parquet(M5B_DIR / "cache" / "clinvar_patho.parquet", index=False)
    loci = json.loads((M4B_DIR / "panels.json").read_text(
        encoding="utf-8"))["loci"]
    rows = []
    for loc in loci:
        for g in loc["genes"]:
            rows.append({"locus_id": loc["locus_id"], "ensg": g["ensg"],
                         **gene_features(patho, g["symbol"], str(loc["chr"]),
                                         g["tss"])})
    out = pd.DataFrame(rows)
    out.to_parquet(M5B_DIR / "cache" / "trackc.parquet", index=False)
    mani = {"n_patho_variants": len(patho),
            "frac_rows_with_signal": round(
                float((out["cv_win_count"] > 0).mean()), 4)}
    (M5B_DIR / "reports" / "trackc.manifest.json").write_text(
        json.dumps(mani, indent=1), encoding="utf-8")
    print(f"trackc: {len(patho)} patho variants, {len(out)} rows, "
          f"signal frac {mani['frac_rows_with_signal']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
