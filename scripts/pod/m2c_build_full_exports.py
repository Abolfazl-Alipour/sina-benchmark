#!/usr/bin/env python
"""m2c_build_full_exports — single full-graph product export (CPU, no GPU).

Owner directives 2026-09-11/12: (1) train ONCE on everything — one supervision
pool, no per-disease holdouts, no retraining per disease; (2) lift the trim —
all diseases, not the 662-disease neuro/diabetes/aging cluster; (3) remove all
commercial-license (DISGENET-tagged) edges.

Input is the normalized full graph (data/processed/m1/trackd/, 21.8M edges —
same normalization as the trimmed exports). Output mirrors the completion
export format under data/processed/m2c/full_exports/ with two differences:
valid.txt/test.txt are EMPTY (no holdouts by design) and the manifest records
DISGENET removal counts per label.

Background labels mirror the KEEP track (ANA-GEN excluded as established
noise; DRG-DRG and EXP-* excluded per spec 7.4 precedent) — unfiltered
otherwise. DISGENET filter applies to every kept label; removal counts are
asserted and reported, including a zero DISGENET-only-unique-loss check.

Usage (sina env, CPU):
    ~/miniconda3/envs/sina/bin/python scripts/pod/m2c_build_full_exports.py
"""

import argparse
import json
import time
from pathlib import Path

import polars as pl

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
TRACKD_EDGES = REPO_ROOT / "data" / "processed" / "m1" / "trackd" / "edges.parquet"
TRACKD_NODES = REPO_ROOT / "data" / "processed" / "m1" / "trackd" / "nodes.parquet"
OUT_DIR = REPO_ROOT / "data" / "processed" / "m2c" / "full_exports"
REPORT_PATH = REPO_ROOT / "reports" / "m2c" / "full_exports.md"

DIS_GEN = "DIS-GEN"
ENCUMBERED = "DISGENET"
# Background family mirrors the KEEP track label set.
BRG_LABELS = ["GEN-GEN", "PHE-GEN", "BPO-GEN", "CCO-GEN", "MFN-GEN", "PWY-GEN",
              "DRG-GEN", "DRG-DIS", "DIS-PHE", "DIS-DIS"]
TYPE_IDS = {"disease": 0, "gene": 1, "drug": 2, "phenotype": 3,
            "biological_process": 4, "cellular_component": 5,
            "molecular_function": 6, "pathway": 7, "anatomy": 8, "exposure": 9}


def utcnow_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def drop_encumbered(edges: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Split edges into (tagged, keep) on DISGENET direct source tags."""
    tagged = edges.with_columns(
        pl.col("sources_direct").fill_null([])).filter(
        pl.col("sources_direct").list.contains(ENCUMBERED))
    keep = edges.with_columns(
        pl.col("sources_direct").fill_null([])).filter(
        ~pl.col("sources_direct").list.contains(ENCUMBERED)).drop(
        "sources_direct")
    return tagged, keep


def unique_triples(edges: pl.DataFrame) -> pl.DataFrame:
    return edges.unique(subset=["src_id", "relation", "dst_id"])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cols = ["src_id", "relation", "dst_id", "label", "sources_direct"]
    edges = pl.scan_parquet(TRACKD_EDGES).select(cols).collect()
    total_in = edges.height
    print(f"trackd edges: {total_in}", flush=True)

    per_label_in = edges.group_by("label").len().sort("len", descending=True)
    tagged, keep = drop_encumbered(edges)
    removal = tagged.group_by("label").len().rename({"len": "removed"})
    print(f"DISGENET-tagged (all labels): {tagged.height}", flush=True)

    sup = unique_triples(keep.filter(pl.col("label") == DIS_GEN))
    brg = unique_triples(keep.filter(pl.col("label").is_in(BRG_LABELS)))
    assert brg.filter(pl.col("label") == DIS_GEN).height == 0
    n_diseases = sup["src_id"].n_unique()
    print(f"supervision: {sup.height} edges / {n_diseases} diseases; "
          f"background: {brg.height} edges", flush=True)
    if args.dry_run:
        print("DRY-RUN: no writes", flush=True)
        return 0

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    brg.select("src_id", "relation", "dst_id").write_csv(
        OUT_DIR / "train1.txt", separator="\t", include_header=False)
    sup.select("src_id", "relation", "dst_id").write_csv(
        OUT_DIR / "train2.txt", separator="\t", include_header=False)
    (OUT_DIR / "valid.txt").write_text("", encoding="utf-8")
    (OUT_DIR / "test.txt").write_text("", encoding="utf-8")
    (OUT_DIR / "strip_list.json").write_text(json.dumps(
        {"disjoint": True, "mode": "full-graph-no-holdout",
         "disgenet_removed": True}, indent=2), encoding="utf-8")

    nodes = pl.scan_parquet(TRACKD_NODES).collect()
    entities = set(brg["src_id"].to_list()) | set(brg["dst_id"].to_list())
    entities |= set(sup["src_id"].to_list()) | set(sup["dst_id"].to_list())
    subset = nodes.filter(pl.col("node_id").is_in(entities))
    assert subset.height == len(entities), "vocab gap: edges reference nodes"
    with (OUT_DIR / "entity_types.txt").open("w", encoding="utf-8") as ft, \
            (OUT_DIR / "entity_names.txt").open("w", encoding="utf-8") as fn:
        for r in subset.select("node_id", "node_type",
                               "node_name").to_dicts():
            ft.write(f"{r['node_id']}\t{TYPE_IDS[r['node_type']]}\n")
            fn.write(f"{r['node_id']}\t{r['node_name']}\n")
    n_genes = subset.filter(pl.col("node_type") == "gene").height
    n_dis = subset.filter(pl.col("node_type") == "disease").height
    print(f"vocab: {subset.height} entities ({n_genes} genes, "
          f"{n_dis} diseases)", flush=True)

    manifest = {"generated_at": utcnow_iso(), "mode": "full-graph-no-holdout",
                "trackd_edges_in": total_in,
                "labels_in": per_label_in.to_dicts(),
                "disgenet_removed_by_label": removal.to_dicts(),
                "train1_background": brg.height, "train2_supervision": sup.height,
                "supervision_diseases": n_diseases, "vocab_entities": subset.height,
                "vocab_genes": n_genes, "vocab_diseases": n_dis,
                "valid_test": "empty-by-design-no-holdouts"}
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2),
                                           encoding="utf-8")
    lines = ["<!-- track_scope: benchmark -->",
             "# Full-graph product export (single-train, no holdouts)",
             "", f"Generated {manifest['generated_at']}.",
             f"Supervision: {sup.height} DIS-GEN / {n_diseases} diseases. "
             f"Background: {brg.height} ({', '.join(BRG_LABELS)}).",
             f"Vocab: {subset.height} entities ({n_genes} genes, {n_dis} diseases).",
             "DISGENET-tagged edges removed from every kept label "
             "(counts in manifest; DISGENET-only-unique-loss is zero by "
             "construction — all tagged edges co-sourced).",
             "valid.txt/test.txt empty by design: one model trains on "
             "everything; evaluation is recovery disclosure, not holdout."]
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT_DIR} + {REPORT_PATH}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
