#!/usr/bin/env python
"""m2_build_completion_exports — within-disease 95/5 completion exports (Phase 1).

Adapts m2_build_lodo_exports.py but does NOT edit it in place.
Completion semantics: for each main disease (MONDO:0004975 AD, MONDO:0005180 PD),
split its DIS-GEN edges 95/5 by evidence_group_id (seed 1337, twins co-assigned),
then split the 95% train portion 15/85 into valid/train_keep. train2 = other
diseases (full 191,953 pool minus target's test+valid) + target's train_keep.
Valid/test are the target's own held-out genes only (so metrics are for the
target disease). Disjointness asserts mirror LODO: train2∩train1=0, etc.

Output mirrors LODO format under data/processed/m2c/completion_exports/<unit>/
with train1.txt (1,431,945 KEEP), train2.txt, valid.txt, test.txt,
entity_names/types, strip_list.json (empty for completion), manifest.parquet.

Usage: python scripts/m2_build_completion_exports.py --config configs/m2c_completion_exports_keep.yaml [--dry-run] [--resume]
"""

import argparse
import json
import sys
from pathlib import Path

import polars as pl
from pydantic import BaseModel

from sina.config import load_config
from sina.logging import configure_logging, get_logger
from sina.manifests import Run
from sina.splits import GROUP_COL, assert_split_disjoint, group_assignment

log = get_logger(__name__)

DIS_GEN = "DIS-GEN"

class UnitEntry(BaseModel):
    curie: str
    name: str
    tier: str

class CompletionExportConfig(BaseModel):
    supervision: Path
    brg: Path
    trackd_nodes: Path
    output_dir: Path
    report_path: Path
    manifest_path: Path
    split_seed: int = 1337
    valid_fraction: float = 0.15
    completion_test_fraction: float = 0.05
    drop_brg_labels: list[str] = []
    units: list[UnitEntry]

def curie_safe(curie: str) -> str:
    return curie.replace(":", "__")

def write_triples(edges: pl.DataFrame, path: Path) -> None:
    edges.select("src_id", "relation", "dst_id").write_csv(path, separator="\t", include_header=False)

def write_entity_files(nodes: pl.DataFrame, frames: list[pl.DataFrame], out: Path) -> None:
    entities: set[str] = set()
    for f in frames:
        entities |= set(f["src_id"].to_list()) | set(f["dst_id"].to_list())
    type_ids = {"disease":0,"gene":1,"drug":2,"phenotype":3,"biological_process":4,"cellular_component":5,"molecular_function":6,"pathway":7,"anatomy":8,"exposure":9}
    subset = nodes.filter(pl.col("node_id").is_in(entities))
    with (out / "entity_types.txt").open("w", encoding="utf-8") as ft, (out / "entity_names.txt").open("w", encoding="utf-8") as fn:
        for r in subset.select("node_id", "node_type", "node_name").to_dicts():
            ft.write(f"{r['node_id']}\t{type_ids[r['node_type']]}\n")
            fn.write(f"{r['node_id']}\t{r['node_name']}\n")

def completion_split(unit_edges: pl.DataFrame, test_fraction: float, valid_fraction: float, seed: int):
    # First split: 95% train vs 5% test on evidence_group_id
    assign_test = group_assignment(unit_edges[GROUP_COL], {"train": 1-test_fraction, "test": test_fraction}, seed)
    joined = unit_edges.join(assign_test, on=GROUP_COL, how="left")
    train_groups = joined.filter(pl.col("partition") == "train").drop("partition")
    test = joined.filter(pl.col("partition") == "test").drop("partition")
    # Second split: within train, 15% valid vs 85% keep
    assign_valid = group_assignment(train_groups[GROUP_COL], {"valid": valid_fraction, "train_keep": 1-valid_fraction}, seed)
    joined2 = train_groups.join(assign_valid, on=GROUP_COL, how="left")
    valid = joined2.filter(pl.col("partition") == "valid").drop("partition")
    train_keep = joined2.filter(pl.col("partition") == "train_keep").drop("partition")
    # Disjointness on groups
    assert_split_disjoint({"valid": valid, "test": test})
    assert_split_disjoint({"valid": valid, "train_keep": train_keep})
    assert_split_disjoint({"test": test, "train_keep": train_keep})
    return train_keep, valid, test

def build_exports(cfg: CompletionExportConfig, state: dict) -> None:
    sup, brg_full = state["sup"], state["brg"]
    brg = brg_full.filter(~pl.col("label").is_in(cfg.drop_brg_labels))
    assert brg.filter(pl.col("label") == DIS_GEN).height == 0
    pge = brg.filter(pl.col("label") == "PHE-GEN").height
    bpo = brg.filter(pl.col("label") == "BPO-GEN").height
    assert pge > 500_000, f"BRG PHE-GEN {pge} <=500k"
    assert bpo > 100_000, f"BRG BPO-GEN {bpo} <=100k"
    # BRG is identical across units
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    # Prepare manifest rows
    manifest_rows = []
    for u in cfg.units:
        log.info("building_unit", curie=u.curie)
        u_edges = sup.filter(pl.col("src_id") == u.curie)
        assert u_edges.height >= 200, f"{u.curie} has {u_edges.height} <200, kill trigger"
        train_keep, valid, test = completion_split(u_edges, cfg.completion_test_fraction, cfg.valid_fraction, cfg.split_seed)
        # train2 = other diseases + train_keep
        other = sup.filter(pl.col("src_id") != u.curie)
        train2 = pl.concat([other, train_keep])
        # Disjointness asserts (tuples)
        def _tuples(df):
            return set(zip(df["src_id"].to_list(), df["relation"].to_list(), df["dst_id"].to_list(), strict=False))
        assert _tuples(train2).isdisjoint(_tuples(valid)), f"{u.curie}: train2∩valid"
        assert _tuples(train2).isdisjoint(_tuples(test)), f"{u.curie}: train2∩test"
        assert _tuples(valid).isdisjoint(_tuples(test)), f"{u.curie}: valid∩test"
        # train2 ∩ train1 should be 0 (train1 has no DIS-GEN, so automatic) but assert
        assert train2.filter(pl.col("label") == DIS_GEN).height == train2.height  # all are DIS-GEN
        assert brg.filter(pl.col("label") == DIS_GEN).height == 0
        out = cfg.output_dir / curie_safe(u.curie)
        out.mkdir(parents=True, exist_ok=True)
        write_triples(brg, out / "train1.txt")
        write_triples(train2, out / "train2.txt")
        write_triples(valid, out / "valid.txt")
        write_triples(test, out / "test.txt")
        write_entity_files(state["nodes"], [brg, train2, valid, test], out)
        # strip_list empty for completion
        (out / "strip_list.json").write_text(json.dumps({"disjoint": True, "mode": "completion_95_5"}, indent=2))
        manifest_rows.append({"unit": u.curie, "name": u.name, "train2": train2.height, "valid": valid.height, "test": test.height, "brg": brg.height})
        log.info("unit_done", curie=u.curie, train2=train2.height, valid=valid.height, test=test.height)
    # Write manifest parquet + report
    import pyarrow as pa, pyarrow.parquet as pq
    tbl = pa.Table.from_pylist(manifest_rows)
    pq.write_table(tbl, cfg.output_dir / "manifest.parquet")
    # Report
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)
    with cfg.report_path.open("w", encoding="utf-8") as f:
        f.write("<!-- track_scope: benchmark -->\n# Completion exports (Phase 1 — 95/5 within-disease)\n\n")
        f.write(f"BRG KEEP {brg.height} edges; supervision 191953 total; seed {cfg.split_seed}; valid {cfg.valid_fraction} of train portion; test {cfg.completion_test_fraction} of unit.\n\n")
        f.write("| unit | name | train2 | valid | test | brg |\n|---|---|---:|---:|---:|---:|\n")
        for r in manifest_rows:
            f.write(f"| {r['unit']} | {r['name']} | {r['train2']} | {r['valid']} | {r['test']} | {r['brg']} |\n")
        f.write("\nDisjointness asserts PASS (train2∩valid, train2∩test, valid∩test, train2∩train1 via label). Twins co-assigned via evidence_group_id.\n")

def load_inputs(cfg: CompletionExportConfig) -> dict:
    sup = pl.read_parquet(cfg.supervision)
    brg = pl.read_parquet(cfg.brg)
    nodes = pl.read_parquet(cfg.trackd_nodes)
    return {"sup": sup, "brg": brg, "nodes": nodes}

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    configure_logging()
    cfg = load_config(args.config, CompletionExportConfig)
    if args.dry_run:
        print(f"would build {len(cfg.units)} completion exports → {cfg.output_dir}")
        return 0
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)
    with Run(script="m2_build_completion_exports", manifest_path=cfg.manifest_path, config={"config": str(args.config)}) as run:
        with run.stage("build", outputs=[cfg.report_path, cfg.output_dir / "manifest.parquet"]):
            state = load_inputs(cfg)
            build_exports(cfg, state)
    log.info("completion_exports_done")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
