"""Track E: eQTL overlap features per (locus, gene).

Pre-registered definitions (fixed before seeing any v5 number):
- Disease side: OT credible-set variants + posteriorProbability (cached).
- Expression side: gene's GTEx v8 independent-eQTL variant (chr,pos) list,
  any tissue (cached per gene).
- Overlap key: exact (chr, pos) match.
- eq_best_pp: max disease posterior among overlapping variants (0.0 none).
- eq_n_tissues: # distinct tissues among overlapping eQTLs.
- eq_any: 1.0 if any overlap else 0.0.
Genes without eQTL cache (outside covered loci) get zeros.
"""

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
M4B_DIR = REPO_ROOT / "data" / "processed" / "m4b"
M5B_DIR = REPO_ROOT / "m5b"


def overlap_stats(dis_variants: list, eqtl_variants: list) -> dict:
    """Pure overlap math. dis: [((chr,pos), pp)], eQTL: [((chr,pos), tissue)].
    Unit-tested."""
    eq_by_key: dict = {}
    for key, tissue in eqtl_variants:
        eq_by_key.setdefault(key, set()).add(tissue)
    hits = [(pp, key) for key, pp in dis_variants if key in eq_by_key]
    tissues = set()
    for _, key in hits:
        tissues |= eq_by_key[key]
    return {"eq_best_pp": round(max([pp for pp, _ in hits] + [0.0]), 6),
            "eq_n_tissues": len(tissues),
            "eq_any": float(bool(hits))}


def main() -> int:
    import pandas as pd

    feat = pd.read_parquet(M4B_DIR / "features" / "features.parquet",
                           columns=["locus_id", "chr", "ensg"])
    panels = {loc["locus_id"]: (str(loc["chr"]), loc["pos"])
              for loc in json.loads((M4B_DIR / "panels.json").read_text(
                  encoding="utf-8"))["loci"]}
    pos2loc = {(c, p): lid for lid, (c, p) in panels.items()}
    # finemap file -> locus via lead-variant match
    file_locus = {}
    for f in (M4B_DIR / "ot_cache" / "credible_sets").glob("*.json"):
        d = json.loads(f.read_text(encoding="utf-8"))
        v = d["variant"]
        lid = pos2loc.get((str(v["chromosome"]), v["position"]))
        if lid:
            file_locus[f.stem] = lid
    dis_by_locus: dict = {}
    for f in (M5B_DIR / "cache" / "finemap").glob("*.json"):
        if f.name == "manifest.json":
            continue
        lid = file_locus.get(f.stem)
        if lid is None:
            continue
        lchr = panels[lid][0]
        bucket = dis_by_locus.setdefault(lid, [])
        seen = {k for k, _ in bucket}
        for row in json.loads(f.read_text(encoding="utf-8"))["rows"]:
            v = row["variant"]
            key = (str(v["chromosome"]), v["position"])
            if key[0] != lchr or key in seen:
                continue
            bucket.append((key, row["posteriorProbability"]))
            seen.add(key)
    eq_by_gene: dict = {}
    for f in (M5B_DIR / "cache" / "eqtl").glob("*.json"):
        if f.name in ("manifest.json", "gene_ids.json"):
            continue
        d = json.loads(f.read_text(encoding="utf-8"))
        eq_by_gene[d["ensg"]] = [((e["chr"], e["pos"]), e["tissue"])
                                 for e in d["eqtls"]]
    rows = []
    for _, r in feat.drop_duplicates(["locus_id", "ensg"]).iterrows():
        rows.append({"locus_id": r["locus_id"], "ensg": r["ensg"],
                     **overlap_stats(
                         dis_by_locus.get(r["locus_id"], []),
                         eq_by_gene.get(r["ensg"], []))})
    out = pd.DataFrame(rows)
    out.to_parquet(M5B_DIR / "cache" / "tracke.parquet", index=False)
    mani = {"n_rows": len(out), "n_genes_with_eqtl": len(eq_by_gene),
            "frac_rows_with_signal": round(float((out["eq_any"] > 0).mean()),
                                          4)}
    (M5B_DIR / "reports" / "tracke.manifest.json").write_text(
        json.dumps(mani, indent=1), encoding="utf-8")
    print(f"tracke: {len(out)} rows, signal frac "
          f"{mani['frac_rows_with_signal']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
