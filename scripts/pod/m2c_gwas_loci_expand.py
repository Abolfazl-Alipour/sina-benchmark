#!/usr/bin/env python
"""m2c_gwas_loci_expand — GWAS loci panels for expansion diseases (CPU).

Mirrors scripts/pod/m2c_gwas_loci_build.py (same TSV parsing, p-thresh,
+/-500kb locus merge, output schema) with three differences for the
single-model product track: (1) --units takes any MONDO ids (no hardcoded
PD/AD); (2) symbol map comes from the local GENCODE GTF (no OT raw pulls);
(3) gene universe is the full-export vocab (--entity-types, type 1 = gene),
recorded as genes_universe (NOT the 17,448 pinned list).

Usage (sina env, CPU):
    ~/miniconda3/envs/sina/bin/python scripts/pod/m2c_gwas_loci_expand.py \
        --tsv-dir /tmp/gwas_annot --out-dir reports/m2c \
        --units MONDO_0005148 MONDO_0005090 MONDO_0007254 MONDO_0005575 \
        --entity-types data/processed/m2c/full_exports/entity_types.txt
"""

import argparse
import csv
import gzip
import json
import re
import time
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
GTF_PATH = Path("/tmp/gencode.v46.gtf.gz")
P_GWAS = 5e-8
LOCUS_WINDOW = 500_000
PARTS = ["pre2023", "2023", "2024", "2025", "2026"]


def utcnow_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def load_symbol_map() -> dict:
    sym2ens: dict[str, str] = {}
    with gzip.open(GTF_PATH, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attr = dict(re.findall(r'(\S+) "([^"]+)"', parts[8]))
            sym, eid = attr.get("gene_name"), attr.get("gene_id", "").split(".")[0]
            if sym and eid.startswith("ENSG"):
                if attr.get("gene_type") == "protein_coding":
                    sym2ens[sym] = eid
                else:
                    sym2ens.setdefault(sym, eid)
    return sym2ens


def split_genes(cell: str) -> list:
    if not cell:
        return []
    parts = re.split(r"[;,]| - | x | / ", cell)
    return [p.strip() for p in parts if p.strip() and p.strip() not in {"NR", "NA", "-"}]


def parse_pval(cell: str):
    try:
        return float(cell)
    except (TypeError, ValueError):
        return None


def parse_year(cell: str):
    m = re.match(r"(\d{4})-\d{2}-\d{2}", (cell or "").strip())
    return int(m.group(1)) if m else None


def parse_int(cell: str):
    try:
        return int((cell or "").strip())
    except ValueError:
        return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--units", nargs="+", required=True,
                    help="MONDO ids, e.g. MONDO_0005148")
    ap.add_argument("--entity-types", type=Path, required=True)
    ap.add_argument("--p-thresh", type=float, default=P_GWAS)
    args = ap.parse_args()

    sym2ens = load_symbol_map()
    print(f"symbol map: {len(sym2ens)} symbols", flush=True)
    universe = set()
    with args.entity_types.open(encoding="utf-8") as f:
        for line in f:
            eid, typ = line.rstrip("\n").split("\t")
            if typ == "1" and eid.startswith("ensembl:ENSG"):
                universe.add(eid)
    print(f"universe genes: {len(universe)}", flush=True)
    units = {("mondo__" + u.split("_")[1].lower()): u for u in args.units}

    rows_by_unit: dict[str, list] = defaultdict(list)
    for part in PARTS:
        path = args.tsv_dir / f"gwas-catalog-download-associations-alt.{part}.tsv"
        if not path.exists():
            print(f"missing part {part}, skipping", flush=True)
            continue
        with path.open(encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                uri = row.get("MAPPED_TRAIT_URI") or ""
                for unit, mondo in units.items():
                    if mondo in uri:
                        rows_by_unit[unit].append((part, row))

    for unit, rows in rows_by_unit.items():
        print(f"{unit}: {len(rows)} trait rows", flush=True)
        assocs = []
        unresolved: dict[str, int] = defaultdict(int)
        for part, row in rows:
            pval = parse_pval((row.get("P-VALUE") or "").strip())
            year = parse_year(row.get("DATE")) or parse_year(row.get("DATE ADDED TO CATALOG"))
            genes: dict[str, set] = {}
            for cell in split_genes(row.get("SNP_GENE_IDS") or ""):
                if cell.startswith("ENSG"):
                    genes.setdefault("ensembl:" + cell, set()).add("snp_gene")
            for role, col in (("upstream", "UPSTREAM_GENE_ID"),
                              ("downstream", "DOWNSTREAM_GENE_ID")):
                for cell in split_genes(row.get(col) or ""):
                    if cell.startswith("ENSG"):
                        genes.setdefault("ensembl:" + cell, set()).add(role)
            for colname, role in (("MAPPED_GENE", "mapped"), ("REPORTED GENE(S)", "reported")):
                for sym in split_genes(row.get(colname) or ""):
                    base = re.sub(r"\s*[\\(\\[].*$", "", sym).strip()
                    if base in sym2ens:
                        genes.setdefault("ensembl:" + sym2ens[base], set()).add(role)
                    else:
                        unresolved[sym] += 1
            trait_uris = [u.strip() for u in
                          (row.get("MAPPED_TRAIT_URI") or "").split(",")
                          if u.strip()]
            unit_uri = "http://purl.obolibrary.org/obo/" + units[unit]
            assocs.append({
                "part": part, "pvalue": pval, "trait_uris": trait_uris,
                "direct_trait": trait_uris == [unit_uri],
                "genomewide": pval is not None and pval <= args.p_thresh,
                "year": year, "chr": (row.get("CHR_ID") or "").strip(),
                "pos": parse_int(row.get("CHR_POS")),
                "snps": (row.get("SNPS") or "").strip(),
                "pubmed": (row.get("PUBMEDID") or "").strip(),
                "study_acc": (row.get("STUDY ACCESSION") or "").strip(),
                "genes": {g: sorted(r) for g, r in genes.items()},
            })

        gw = [a for a in assocs if a["genomewide"] and a["pos"] is not None]
        gw.sort(key=lambda a: (a["chr"], a["pos"]))
        loci = []
        for a in gw:
            placed = False
            for loc in loci:
                if loc["chr"] == a["chr"] and abs(loc["pos_max"] - a["pos"]) <= LOCUS_WINDOW \
                        and a["pos"] >= loc["pos_min"] - LOCUS_WINDOW:
                    loc["members"].append(a)
                    loc["pos_min"] = min(loc["pos_min"], a["pos"])
                    loc["pos_max"] = max(loc["pos_max"], a["pos"])
                    placed = True
                    break
            if not placed:
                loci.append({"chr": a["chr"], "pos_min": a["pos"],
                             "pos_max": a["pos"], "members": [a]})
        for i, loc in enumerate(loci):
            genes: dict[str, set] = defaultdict(set)
            years = []
            pmin = None
            for a in loc["members"]:
                for g, roles in a["genes"].items():
                    genes[g].update(roles)
                if a["year"]:
                    years.append(a["year"])
                if a["pvalue"] is not None:
                    pmin = a["pvalue"] if pmin is None else min(pmin, a["pvalue"])
            loc["id"] = f"{unit}_locus{i:03d}"
            loc["n_assoc"] = len(loc["members"])
            loc["genes"] = {g: sorted(r) for g, r in genes.items()}
            loc["genes_universe"] = sorted(g for g in genes if g in universe)
            loc["year_min"] = min(years) if years else None
            loc["year_max"] = max(years) if years else None
            loc["p_min"] = pmin

        first_year_direct, first_year_all = {}, {}
        for a in assocs:
            if a["year"] is None or not a["genomewide"]:
                continue
            for g in a["genes"]:
                if g not in universe:
                    continue
                if g not in first_year_all or a["year"] < first_year_all[g]:
                    first_year_all[g] = a["year"]
                if a["direct_trait"] and (g not in first_year_direct
                                          or a["year"] < first_year_direct[g]):
                    first_year_direct[g] = a["year"]

        out = {"unit": unit, "generated_at": utcnow_iso(),
               "source": "EBI GWAS Catalog bulk "
                         "gwas-catalog-associations_ontology-annotated-split.zip",
               "p_thresh": args.p_thresh, "locus_window": LOCUS_WINDOW,
               "gene_universe": "full_exports vocab (NOT pinned 17448)",
               "n_trait_rows": len(rows),
               "n_genomewide": sum(1 for a in assocs if a["genomewide"]),
               "n_loci": len(loci),
               "n_unresolved_symbols": len(unresolved),
               "unresolved_top": sorted(unresolved.items(), key=lambda kv: -kv[1])[:20],
               "gene_first_year_direct": first_year_direct,
               "gene_first_year_all": first_year_all,
               "loci": loci}
        args.out_dir.mkdir(parents=True, exist_ok=True)
        opath = args.out_dir / f"gwas_loci_{unit}.json"
        opath.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        print(f"{unit}: {len(loci)} loci, wrote {opath}", flush=True)
    if not rows_by_unit:
        print("no trait rows matched any unit", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
