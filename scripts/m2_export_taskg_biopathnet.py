#!/usr/bin/env python
"""m2_export_taskg_biopathnet — export Task-G data in BioPathNet format (spec §7.4).

BioPathNet's `biomedical` dataset expects, per dataset dir:
  train1.txt  — BRG message-passing triples (never supervised)
  train2.txt  — supervision train triples (DIS-GEN)
  valid.txt / test.txt — eval triples
  entity_types.txt — entity<TAB>type_id for EVERY entity in the graph
  entity_names.txt — entity<TAB>display name

BRG edges referencing held-out family diseases are excluded at training time
(split-aware BRG builder, spec §7.4(2b)).

Usage: python scripts/m2_export_taskg_biopathnet.py --config configs/m2b_taskg_export.yaml
           [--dry-run] [--resume]
"""

import argparse
import sys
from pathlib import Path

import polars as pl
from pydantic import BaseModel

from sina.config import load_config
from sina.logging import configure_logging, get_logger
from sina.manifests import Run

log = get_logger(__name__)


class SplitSpec(BaseModel):
    relation_label: str


class ExportConfig(BaseModel):
    taskg_dir: Path
    splits_dir: Path
    trackd_nodes: Path
    families_config: Path
    splits: dict[str, SplitSpec]
    type_ids: dict[str, int]
    output_dir: Path
    manifest_path: Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, ExportConfig)

    if args.dry_run:
        print(f"would export Task-G BioPathNet data for {list(cfg.splits)} → {cfg.output_dir}")
        return 0

    cfg.output_dir.mkdir(parents=True, exist_ok=True)

    with Run(
        script="m2_export_taskg_biopathnet",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        for split_name in cfg.splits:
            if args.resume and split_name in run.done_stages:
                log.info("split_skipped", split=split_name)
                continue
            with run.stage(split_name):
                export_split(cfg, split_name)

    log.info("export_done")
    return 0


def export_split(cfg: ExportConfig, split_name: str) -> None:
    spec = cfg.splits[split_name]
    out = cfg.output_dir / split_name
    out.mkdir(parents=True, exist_ok=True)

    # Held-out family nodes: the disease side of the split's valid+test eval
    # edges (DIS-GEN is disease→gene; spec §7.4(2b) excludes only BRG edges
    # referencing held-out family diseases, not the candidate genes).
    valid = pl.read_parquet(cfg.splits_dir / split_name / "valid.parquet")
    test = pl.read_parquet(cfg.splits_dir / split_name / "test.parquet")
    eval_edges = pl.concat([valid, test])
    family_nodes = set(eval_edges["src_id"].to_list())

    # train1: BRG message passing. CHANGES/0013 (owner-approved 2026-08-01):
    # retain low-leakage context edges (DIS-PHE, DIS-DIS, hierarchies, …)
    # for held-out diseases so query nodes stay connected to the graph;
    # only DRG-DIS / EXP-DIS edges touching held-out diseases are excluded
    # (close 2-hop proxies of the DIS-GEN supervision target).
    brg = pl.read_parquet(cfg.taskg_dir / "taskg_message_passing.parquet")
    touches_family = pl.col("src_id").is_in(family_nodes) | pl.col("dst_id").is_in(family_nodes)
    brg = brg.filter(~touches_family | ~pl.col("label").is_in(["DRG-DIS", "EXP-DIS"]))
    write_triples(brg, out / "train1.txt")

    # train2: supervision train DIS-GEN edges from the split.
    train = pl.read_parquet(cfg.splits_dir / split_name / "train.parquet")
    sup = train.filter(pl.col("label") == spec.relation_label)
    write_triples(sup, out / "train2.txt")

    # valid/test.
    write_triples(valid.filter(pl.col("label") == spec.relation_label), out / "valid.txt")
    write_triples(test.filter(pl.col("label") == spec.relation_label), out / "test.txt")

    # entity_types + entity_names for every entity in all four files.
    all_edges = pl.concat(
        [
            brg.select("src_id", "dst_id"),
            sup.select("src_id", "dst_id"),
            valid.select("src_id", "dst_id"),
            test.select("src_id", "dst_id"),
        ]
    )
    entities = set(all_edges["src_id"].to_list()) | set(all_edges["dst_id"].to_list())
    nodes = pl.read_parquet(cfg.trackd_nodes)
    write_entity_files(cfg, nodes, entities, out)

    log.info(
        "split_exported",
        split=split_name,
        train1=brg.height,
        train2=sup.height,
        valid=valid.height,
        test=test.height,
        entities=len(entities),
        family_excluded=len(family_nodes),
    )


def write_triples(edges: pl.DataFrame, path: Path) -> None:
    edges.select("src_id", "relation", "dst_id").write_csv(
        path, separator="\t", include_header=False
    )


def write_entity_files(
    cfg: ExportConfig, nodes: pl.DataFrame, entities: set[str], out: Path
) -> None:
    type_map = {t: i for t, i in cfg.type_ids.items()}
    subset = nodes.filter(pl.col("node_id").is_in(entities))
    missing = entities - set(subset["node_id"].to_list())
    if missing:
        raise ValueError(f"{len(missing)} entities missing from nodes table")
    with (
        (out / "entity_types.txt").open("w", encoding="utf-8") as ft,
        (out / "entity_names.txt").open("w", encoding="utf-8") as fn,
    ):
        for r in subset.select("node_id", "node_type", "node_name").to_dicts():
            tid = type_map.get(r["node_type"])
            if tid is None:
                raise ValueError(f"unmapped node_type {r['node_type']} for {r['node_id']}")
            ft.write(f"{r['node_id']}\t{tid}\n")
            fn.write(f"{r['node_id']}\t{r['node_name']}\n")


if __name__ == "__main__":
    sys.exit(main())
