"""Track E fetch: versioned GTEx v8 IDs + independent eQTL variant lists
for every candidate gene at covered loci. Cache-first + resume-safe."""

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
M4B_DIR = REPO_ROOT / "data" / "processed" / "m4b"
M5B_DIR = REPO_ROOT / "m5b"
DST = M5B_DIR / "cache" / "eqtl"
BASE = "https://gtexportal.org/api/v2"


def get(path: str, params: dict) -> dict:
    url = BASE + path + "?" + urllib.parse.urlencode(params)
    last = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=90) as resp:
                return json.loads(resp.read().decode())
        except Exception as exc:
            last = exc
            time.sleep(2 + attempt * 2)
    raise RuntimeError(f"GTEx fetch failed {path} {params}: {last}")


def gene_version(ensg: str) -> str | None:
    d = get("/reference/gene", {"geneId": ensg, "datasetId": "gtex_v8"})
    rows = d.get("data", [])
    return rows[0]["gencodeId"] if rows else None


def independent_eqtls(vid: str) -> list:
    out, page = [], 0
    while True:
        d = get("/association/independentEqtl",
                {"gencodeId": vid, "datasetId": "gtex_v8",
                 "page": page, "itemsPerPage": 250})
        out.extend(d.get("data", []))
        info = d.get("paging_info", {})
        if page + 1 >= info.get("numberOfPages", 1):
            break
        page += 1
    return out


def main() -> int:
    import pandas as pd

    DST.mkdir(parents=True, exist_ok=True)
    feat = pd.read_parquet(M4B_DIR / "features" / "features.parquet",
                           columns=["locus_id", "ensg"])
    ot = json.loads((M4B_DIR / "ot_scores.json").read_text(
        encoding="utf-8"))["scores"]
    genes = sorted(feat[feat["locus_id"].isin(ot)]["ensg"].unique().tolist())
    ids_path = DST / "gene_ids.json"
    ids = json.loads(ids_path.read_text(encoding="utf-8")) if \
        ids_path.exists() else {}
    n_hit = n_new = 0
    for i, ensg in enumerate(genes):
        if ensg in ids:
            n_hit += 1
            continue
        ids[ensg] = gene_version(ensg)
        n_new += 1
        time.sleep(0.15)
        if (i + 1) % 200 == 0:
            ids_path.write_text(json.dumps(ids), encoding="utf-8")
            print(f"eqtl ids: {i + 1}/{len(genes)}", flush=True)
    ids_path.write_text(json.dumps(ids), encoding="utf-8")
    n_q = n_e = n_err = 0
    for i, ensg in enumerate(genes):
        vid = ids.get(ensg)
        if not vid:
            continue
        fpath = DST / f"{ensg}.json"
        if fpath.exists():
            continue
        try:
            rows = independent_eqtls(vid)
        except RuntimeError as exc:
            fpath.write_text(json.dumps({"ensg": ensg, "gencodeId": vid,
                                         "eqtls": [], "error": str(exc)[:200]}),
                             encoding="utf-8")
            n_err += 1
            print(f"eqtl SKIP {ensg} {vid}: {str(exc)[:120]}", flush=True)
            time.sleep(1.0)
            continue
        keep = [{"chr": str(r["chromosome"]).removeprefix("chr"),
                 "pos": r["pos"], "tissue": r["tissueSiteDetailId"]}
                for r in rows]
        fpath.write_text(json.dumps({"ensg": ensg, "gencodeId": vid,
                                     "eqtls": keep}), encoding="utf-8")
        n_q += 1
        n_e += len(keep)
        time.sleep(0.15)
        if (i + 1) % 200 == 0:
            print(f"eqtl lists: {i + 1}/{len(genes)}", flush=True)
    (DST / "manifest.json").write_text(json.dumps(
        {"n_genes": len(genes), "n_mapped_ids": sum(1 for v in ids.values()
                                                   if v),
         "n_eqtls": n_e + sum(len(json.loads(
             (DST / f"{g}.json").read_text(
                 encoding="utf-8"))["eqtls"]) for g in genes
             if (DST / f"{g}.json").exists()),
         "note": "cumulative counts include cache hits"}, indent=1),
        encoding="utf-8")
    print(f"eqtl done: {len(genes)} genes, {n_q} newly pulled, "
          f"{n_err} skipped with errors", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
