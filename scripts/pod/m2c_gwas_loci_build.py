#!/usr/bin/env python
"""m2c_gwas_loci_build — PD/AD GWAS loci from the Catalog bulk TSV (CPU, no GPU).

Owner directive 2026-09-06: temporal-split + per-locus tracks. Input is the
EBI FTP bulk file (gwas-catalog-associations_ontology-annotated-split.zip,
year-split parts; DisGeNET bulk is login-walled, OT parquet is the fallback).
Parses the parts, keeps rows whose MAPPED_TRAIT_URI names the unit MONDO id,
and emits per-unit locus JSON used by both downstream probes:

- temporal: per-gene first-reported year (publication DATE, fallback DATE
  ADDED) -> pre/post-cutoff gene sets for the prospective null test.
- per-locus: chr/pos-grouped loci with candidate gene sets (mapped + reported
  + up/downstream, resolved to ensembl:ENSG where possible) for the
  rank-in-locus probe (our checkpoint score vs nearest-gene baseline).

Gene resolution order per row: SNP_GENE_IDS (already Ensembl) >
UPSTREAM/DOWNSTREAM_GENE_ID (already Ensembl) > MAPPED_GENE / REPORTED GENE(S)
symbols via the OT approvedSymbol map (14,212 symbols from the committed
ot_panels raw pulls). Unresolved symbols are recorded, not dropped silently.

Usage (sina env, CPU):
    ~/miniconda3/envs/sina/bin/python scripts/pod/m2c_gwas_loci_build.py \
        --tsv-dir /tmp/gwas_annot --out-dir reports/m2c
"""

import argparse
import csv
import json
import re
import time
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
OT_RAW_DIR = REPO_ROOT / "data" / "processed" / "m2c" / "ot_panels"
PINNED_PATH = REPO_ROOT / "data" / "processed" / "m2c" / "pinned_genes_17448.txt"

UNITS = {"mondo__0005180": "MONDO_0005180", "mondo__0004975": "MONDO_0004975"}
PARTS = ["pre2023", "2023", "2024", "2025", "2026"]
P_GWAS = 5e-8
LOCUS_WINDOW = 500_000  # same-chr associations within +/-500kb merge


def utcnow_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def load_symbol_map() -> dict:
    sym2ens = {}
    for unit in UNITS:
        path = OT_RAW_DIR / f"raw_{unit}.jsonl"
        with path.open(encoding="utf-8") as f:
            for line in f:
                row = json.loads(line)
                tgt = row.get("target") or {}
                s, i = tgt.get("approvedSymbol"), tgt.get("id")
                if s and i:
                    if not i.startswith("ENSG"):
                        i = "ENSG" + i
                    sym2ens.setdefault(s, i)
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--p-thresh", type=float, default=P_GWAS)
    args = ap.parse_args()

    sym2ens = load_symbol_map()
    print(f"symbol map: {len(sym2ens)} symbols", flush=True)
    pinned = set(PINNED_PATH.read_text(encoding="utf-8").strip().split())
    assert len(pinned) == 17448

    rows_by_unit = defaultdict(list)
    for part in PARTS:
        path = args.tsv_dir / f"gwas-catalog-download-associations-alt.{part}.tsv"
        if not path.exists():
            print(f"missing part {part}, skipping", flush=True)
            continue
        with path.open(encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                uri = row.get("MAPPED_TRAIT_URI") or ""
                for unit, mondo in UNITS.items():
                    if mondo in uri:
                        rows_by_unit[unit].append((part, row))

    for unit, rows in rows_by_unit.items():
        print(f"{unit}: {len(rows)} trait rows", flush=True)
        assocs = []
        unresolved = defaultdict(int)
        for part, row in rows:
            pval = parse_pval((row.get("P-VALUE") or "").strip())
            year_pub = parse_year(row.get("DATE"))
            year_added = parse_year(row.get("DATE ADDED TO CATALOG"))
            year = year_pub or year_added
            try:
                pos = int((row.get("CHR_POS") or "").strip())
            except ValueError:
                pos = None
            genes = {}  # ensembl token -> role set
            for cell in split_genes(row.get("SNP_GENE_IDS") or ""):
                tok = cell if cell.startswith("ENSG") else None
                if tok:
                    genes.setdefault("ensembl:" + tok, set()).add("snp_gene")
            for role, col in (("upstream", "UPSTREAM_GENE_ID"),
                              ("downstream", "DOWNSTREAM_GENE_ID")):
                for cell in split_genes(row.get(col) or ""):
                    if cell.startswith("ENSG"):
                        genes.setdefault("ensembl:" + cell, set()).add(role)
            try:
                up_d = int((row.get("UPSTREAM_GENE_DISTANCE") or "").strip())
            except ValueError:
                up_d = None
            try:
                dn_d = int((row.get("DOWNSTREAM_GENE_DISTANCE") or "").strip())
            except ValueError:
                dn_d = None
            sym_genes = []
            for sym in split_genes(row.get("MAPPED_GENE") or ""):
                if sym in sym2ens:
                    genes.setdefault("ensembl:" + sym2ens[sym], set()).add("mapped")
                    sym_genes.append(sym)
                else:
                    unresolved[sym] += 1
            rep_genes = []
            for sym in split_genes(row.get("REPORTED GENE(S)") or ""):
                base = re.sub(r"\s*[\(\[].*$", "", sym).strip()
                if base in sym2ens:
                    genes.setdefault("ensembl:" + sym2ens[base], set()).add("reported")
                    rep_genes.append(base)
                else:
                    unresolved[sym] += 1
            trait_uris = [u.strip() for u in (row.get("MAPPED_TRAIT_URI") or "").split(",")
                          if u.strip()]
            unit_uri = "http://purl.obolibrary.org/obo/" + UNITS[unit]
            assocs.append({
                "part": part, "pvalue": pval,
                "trait_uris": trait_uris,
                "direct_trait": trait_uris == [unit_uri],
                "genomewide": pval is not None and pval <= args.p_thresh,
                "year_pub": year_pub, "year_added": year_added, "year": year,
                "chr": (row.get("CHR_ID") or "").strip(), "pos": pos,
                "snps": (row.get("SNPS") or "").strip(),
                "pubmed": (row.get("PUBMEDID") or "").strip(),
                "study_acc": (row.get("STUDY ACCESSION") or "").strip(),
                "mapped_symbols": sym_genes, "reported_symbols": rep_genes,
                "upstream_id": (row.get("UPSTREAM_GENE_ID") or "").strip(),
                "upstream_dist": up_d,
                "downstream_id": (row.get("DOWNSTREAM_GENE_ID") or "").strip(),
                "downstream_dist": dn_d,
                "intergenic": (row.get("INTERGENIC") or "").strip(),
                "genes": {g: sorted(r) for g, r in genes.items()},
            })

        # Locus grouping: same chr, positions chained within +/-500kb.
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
            genes = defaultdict(set)
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
            loc["genes_pinned"] = sorted(g for g in genes if g in pinned)
            loc["year_min"] = min(years) if years else None
            loc["year_max"] = max(years) if years else None
            loc["p_min"] = pmin
            for a in loc["members"]:
                a["genes"] = {g: sorted(r) for g, r in a["genes"].items()}

        # Temporal gene sets: first-reported year per mapped/reported/snp gene,
        # twice — direct-trait rows only (sole MONDO id in MAPPED_TRAIT_URI) and
        # all rows (incl. proxy/biomarker traits). Each assoc records its trait set.
        first_year_direct = {}
        first_year_all = {}
        for a in assocs:
            if a["year"] is None or not a["genomewide"]:
                continue
            for g in a["genes"]:
                if g not in pinned:
                    continue
                if g not in first_year_all or a["year"] < first_year_all[g]:
                    first_year_all[g] = a["year"]
                if a["direct_trait"] and (g not in first_year_direct
                                          or a["year"] < first_year_direct[g]):
                    first_year_direct[g] = a["year"]

        out = {"unit": unit, "generated_at": utcnow_iso(),
               "source": "EBI GWAS Catalog bulk "
                         "gwas-catalog-associations_ontology-annotated-split.zip "
                         "(ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/)",
               "p_thresh": args.p_thresh, "locus_window": LOCUS_WINDOW,
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
        opath.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n",
                         encoding="utf-8")
        n_pinned_genes = len(first_year_all)
        print(f"{unit}: {len(loci)} loci, {n_pinned_genes} pinned genes with years, "
              f"unresolved symbols {len(unresolved)}, wrote {opath}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
