#!/usr/bin/env python
"""m2c_vep_culprits — functional gold-standard culprits per GWAS locus (CPU).

Owner directive 2026-09-06 (final per-locus test). Curated gold standards are
too thin for PD/AD (1 PD + 6 drug-class AD genes) and bulk L2G is ~20GB, so:
for each GWAS locus lead SNP (min-p genome-wide row, one rsID per row —
verified), query Ensembl REST VEP and name culprit loci where EXACTLY ONE
pinned gene carries a strict-coding consequence (missense or LoF). This
mirrors L2G's pathogenicity feature class and the gold-standards
"functional" evidence class: lab-grade causality, independent of our train2.

Strict coding = missense_variant, stop_gained, stop_lost, start_lost,
frameshift_variant, splice_donor_variant, splice_acceptor_variant,
protein_altering_variant, inframe_insertion, inframe_deletion,
transcript_ablation. Everything else (incl. splice_region, UTR, intron,
regulatory) is non-culprit evidence.

Novelty (the point of the test): culprit must be in NEITHER the train2 NOR
valid splits of the unit — the checkpoint never saw it as a positive.

Output: reports/m2c/vep_culprits_<unit>.json with culprit loci, VEP cache
notes, novelty counts, and a power verdict (pooled novel n>=15 to run the
top-1 test, per owner rule).

Usage (sina env, CPU; network to rest.ensembl.org):
    ~/miniconda3/envs/sina/bin/python scripts/pod/m2c_vep_culprits.py \
        --loci-json reports/m2c/gwas_loci_mondo__0005180.json \
        --unit mondo__0005180 --cache /tmp/vep_cache.json \
        --out reports/m2c/vep_culprits_mondo__0005180.json
"""

import argparse
import json
import time
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
PINNED_PATH = REPO_ROOT / "data" / "processed" / "m2c" / "pinned_genes_17448.txt"
VEP_URL = "https://rest.ensembl.org/vep/human/id"

STRICT_CODING = {
    "missense_variant", "stop_gained", "stop_lost", "start_lost",
    "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant",
    "protein_altering_variant", "inframe_insertion", "inframe_deletion",
    "transcript_ablation",
}


def utcnow_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def vep_batch(ids: list[str]) -> list[dict]:
    data = json.dumps({"ids": ids}).encode()
    last = None
    for attempt in range(4):
        try:
            req = urllib.request.Request(
                VEP_URL, data=data,
                headers={"Content-Type": "application/json",
                         "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception as e:  # retry transient HTTP/network
            last = e
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"VEP failed after retries: {last}")


def base_ensg(gene_id: str) -> str:
    return "ensembl:" + gene_id.split(".")[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--unit", required=True)
    ap.add_argument("--loci-json", type=Path, required=True)
    ap.add_argument("--exports-dir", type=Path,
                    default=REPO_ROOT / "data" / "processed" / "m2c" / "completion_exports")
    ap.add_argument("--cache", type=Path, default=Path("/tmp/vep_cache.json"))
    ap.add_argument("--batch", type=int, default=150)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    unit, disease = args.unit, args.unit.replace("__", ":")
    pinned = set(PINNED_PATH.read_text(encoding="utf-8").strip().split())
    assert len(pinned) == 17448

    loci = json.loads(args.loci_json.read_text(encoding="utf-8"))
    seen = set()
    for split in ("train2.txt", "valid.txt"):
        with (args.exports_dir / unit / split).open(encoding="utf-8") as f:
            for line in f:
                h, _, t = line.rstrip("\n").split("\t")
                if h == disease:
                    seen.add(t)

    # Lead SNP per locus = min-p member rsID.
    lead_of_locus = {}
    want = []
    for loc in loci["loci"]:
        best = None
        for a in loc["members"]:
            snp = (a.get("snps") or "").strip()
            if not snp.startswith("rs") or a.get("pvalue") is None:
                continue
            if best is None or a["pvalue"] < best[1]:
                best = (snp, a["pvalue"])
        if best:
            lead_of_locus[loc["id"]] = best[0]
            if best[0] not in want:
                want.append(best[0])
    print(f"{unit}: {len(lead_of_locus)} loci with lead rsID, {len(want)} unique",
          flush=True)

    cache = {}
    if args.cache.exists():
        cache = json.loads(args.cache.read_text(encoding="utf-8"))
    missing = [rs for rs in want if rs not in cache]
    print(f"VEP cache hits {len(want) - len(missing)}/{len(want)}", flush=True)
    for i in range(0, len(missing), args.batch):
        chunk = missing[i:i + args.batch]
        for rec in vep_batch(chunk):
            cache[rec.get("input", "")] = rec
        print(f"VEP {min(i + args.batch, len(missing))}/{len(missing)}", flush=True)
        time.sleep(0.5)
    args.cache.write_text(json.dumps(cache) + "\n", encoding="utf-8")

    culprits = []
    for loc in loci["loci"]:
        rs = lead_of_locus.get(loc["id"])
        if not rs or rs not in cache:
            continue
        rec = cache[rs]
        if rec.get("error"):
            continue
        coding = {}
        for tc in rec.get("transcript_consequences") or []:
            gid = tc.get("gene_id") or ""
            if not gid.startswith("ENSG"):
                continue
            tok = base_ensg(gid)
            if tok not in pinned:
                continue
            terms = set(tc.get("consequence_terms") or [])
            strict = sorted(terms & STRICT_CODING)
            if strict:
                coding.setdefault(tok, set()).update(strict)
        genes = loc["genes_pinned"]
        if len(genes) >= 2 and len(coding) == 1:
            culprit = next(iter(coding))
            culprits.append({
                "locus": loc["id"], "chr": loc["chr"], "lead_snp": rs,
                "n_genes": len(genes), "culprit": culprit,
                "consequences": sorted(coding[culprit]),
                "novel": culprit not in seen,
            })

    novel = [c for c in culprits if c["novel"]]
    out = {"unit": unit, "generated_at": utcnow_iso(),
           "vep": VEP_URL, "strict_coding": sorted(STRICT_CODING),
           "n_loci": len(loci["loci"]), "n_lead_snps": len(want),
           "n_single_coding_loci": len(culprits),
           "n_novel_culprit_loci": len(novel),
           "power_ok": len(novel) >= 15,
           "culprits": culprits}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(f"{unit}: {len(culprits)} single-coding loci, {len(novel)} novel, "
          f"power_ok={len(novel) >= 15}; wrote {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
