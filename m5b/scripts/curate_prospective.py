"""Track P: curate post-2020 prospective loci from GWAS Catalog split files.

Pre-registered curation rules (fixed before any evaluation number):
- Rows: P-VALUE <= 5e-8, MAPPED_GENE non-empty, CHR_ID/POS valid integers,
  autosomes + X only.
- Locus clustering: within one MAPPED_TRAIT_URI, sort by P-VALUE ascending,
  greedily absorb rows within +/-500kb of the picked lead (mirror of M4B
  panel logic). Lead = lowest p-value row.
- Labels: distinct MAPPED_GENE values of member rows.
- Year: from DATE column; keep > 2020.
Output: m5b/reports/trackp_panels.json + manifest (counts only, no scores).
"""

import json
import zipfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
M5B_DIR = REPO_ROOT / "m5b"
ARCHIVE = M5B_DIR / "cache" / "gwas-associations-split.zip"
P_CUTOFF = 5e-8
HALF_WINDOW = 500_000


def parse_pval(text: str) -> float | None:
    try:
        v = float(text)
    except (TypeError, ValueError):
        return None
    return v if v > 0 else None


def split_genes(text: str) -> list:
    """MAPPED_GENE holds comma lists; intergenic pairs use ' - '.
    Returns clean symbols. Unit-tested."""
    out = []
    for chunk in text.split(","):
        for part in chunk.split(" - "):
            sym = part.strip()
            if sym and sym != "-":
                out.append(sym)
    return out


def cluster_locus(rows: list) -> list:
    """Pure greedy lead-clustering. rows: [(pos, pval, payload)].
    Unit-tested."""
    ordered = sorted(rows, key=lambda r: r[1])
    loci, used = [], set()
    for i, (pos, _p, _pl) in enumerate(ordered):
        if i in used:
            continue
        members = [i]
        for j, (p2, _q, _l) in enumerate(ordered):
            if j not in used and j != i and abs(p2 - pos) <= HALF_WINDOW:
                members.append(j)
                used.add(j)
        used.add(i)
        loci.append([ordered[m] for m in members])
    return loci


def read_hits():
    z = zipfile.ZipFile(ARCHIVE)
    for name in z.namelist():
        with z.open(name) as f:
            head = f.readline().decode(errors="replace").rstrip("\n"
                                                                ).split("\t")
            ix = {c: head.index(c) for c in
                  ("DATE", "CHR_ID", "CHR_POS", "MAPPED_GENE", "P-VALUE",
                   "MAPPED_TRAIT", "MAPPED_TRAIT_URI", "STUDY ACCESSION",
                   "DISEASE/TRAIT")}
            for raw in f:
                parts = raw.decode(errors="replace").rstrip("\n"
                                                            ).split("\t")
                if len(parts) < len(head):
                    continue
                def g(c, _parts=parts, _ix=ix):
                    return _parts[_ix[c]].strip()

                year = (g("DATE") or "")[:4]
                if not year.isdigit() or int(year) <= 2020:
                    continue
                p = parse_pval(g("P-VALUE"))
                if p is None or p > P_CUTOFF:
                    continue
                if not g("MAPPED_GENE") or g("MAPPED_GENE") == "-":
                    continue
                try:
                    chrom, pos = g("CHR_ID"), int(float(g("CHR_POS")))
                except ValueError:
                    continue
                if chrom not in {str(c) for c in range(1, 23)} | {"X"}:
                    continue
                yield {"chr": chrom, "pos": pos, "p": p,
                       "genes": split_genes(g("MAPPED_GENE")),
                       "trait": g("MAPPED_TRAIT") or g("DISEASE/TRAIT"),
                       "trait_uri": g("MAPPED_TRAIT_URI"),
                       "study_acc": g("STUDY ACCESSION"), "year": int(year)}


def main() -> int:
    by_trait: dict = {}
    n_rows = 0
    for h in read_hits():
        n_rows += 1
        by_trait.setdefault(h["trait_uri"] or h["trait"], []).append(h)
    loci = []
    by_chr: dict = {}
    for _trait, hits in by_trait.items():
        for h in hits:
            by_chr.setdefault(h["chr"], []).append((h["pos"], h["p"], h))
    for chr_, rows in by_chr.items():
        for members in cluster_locus(rows):
            lead = min(members, key=lambda m: m[1])
            genes = sorted({g for m in members for g in m[2]["genes"]})
            loci.append({"chr": chr_, "pos": lead[0],
                         "traits": sorted({m[2]["trait"] for m in members}),
                         "lead_p": lead[1],
                         "year_min": min(m[2]["year"] for m in members),
                         "n_rows": len(members), "genes": genes,
                         "studies": sorted({m[2]["study_acc"]
                                            for m in members})})
    loci.sort(key=lambda loc: (loc["chr"], loc["pos"]))
    out = {"n_hits": n_rows, "n_traits": len(by_trait), "n_loci": len(loci),
           "loci": loci}
    (M5B_DIR / "reports" / "trackp_panels.json").write_text(
        json.dumps(out), encoding="utf-8")
    (M5B_DIR / "reports" / "trackp.manifest.json").write_text(json.dumps(
        {k: v for k, v in out.items() if k != "loci"}, indent=1),
        encoding="utf-8")
    print(f"trackp: {n_rows} hits, {len(by_trait)} traits, "
          f"{len(loci)} loci", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
