#!/usr/bin/env python
"""m2_build_taskg — build the curated Task-G edge subset on OptimusKG (spec §7.4).

Stages (resumable; --dry-run prints the plan):
  1. load — canonical Track-D edges + relation curation config
  2. build — supervised DIS-GEN task edges + message-passing context subset
  3. report — counts per role → reports/m2a/taskg_report.md

The output is the graph BioPathNet Task-G will train/evaluate on. ANA-GEN is
excluded by default (expression-ubiquity noise); EXP-*/DRG-DRG are excluded
per spec §7.4. Relation curation itself is owner sign-off gated — this script
consumes `configs/m2_trackd_relations.yaml` as-is.

Usage: python scripts/m2_build_taskg.py --config configs/m2_trackd_relations.yaml
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


class RoleEntry(BaseModel):
    label: str
    n_edges: int
    note: str = ""
    relations: list[str] | None = None
    role: str | None = None
    default_weight: float | None = None


class TaskGConfig(BaseModel):
    status: str
    supervised: list[RoleEntry]
    message_passing: list[RoleEntry]
    downweighted: list[RoleEntry]
    excluded: list[RoleEntry]


DEFAULT_EDGES = Path("data/processed/m1/trackd/edges.parquet")
DEFAULT_OUT_DIR = Path("data/processed/m2/taskg")
DEFAULT_REPORT = Path("reports/m2a/taskg_report.md")
DEFAULT_MANIFEST = Path("data/raw/manifests/m2_build_taskg.json")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--edges", type=Path, default=DEFAULT_EDGES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report-path", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--manifest-path", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, TaskGConfig)

    if args.dry_run:
        n_sup = sum(e.n_edges for e in cfg.supervised)
        n_mp = sum(e.n_edges for e in cfg.message_passing)
        print(
            f"would build Task-G subset: {n_sup:,} supervised + {n_mp:,} "
            f"message-passing edges → {args.output_dir}"
        )
        return 0

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report_path.parent.mkdir(parents=True, exist_ok=True)

    with Run(
        script="m2_build_taskg",
        manifest_path=args.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        if not (args.resume and "build" in run.done_stages):
            with run.stage("build"):
                counts = build_subset(cfg, args.edges, args.output_dir)
        else:
            counts = None

        with run.stage("report", outputs=[args.report_path]):
            write_report(cfg, args, counts)

    log.info("taskg_done")
    return 0


def build_subset(cfg: TaskGConfig, edges_path: Path, out_dir: Path) -> dict:
    edges = pl.read_parquet(edges_path)
    sup_labels = {e.label for e in cfg.supervised}
    mp_labels = {e.label for e in cfg.message_passing}
    down_labels = {e.label for e in cfg.downweighted}

    supervised = edges.filter(pl.col("label").is_in(sup_labels))
    message = edges.filter(pl.col("label").is_in(mp_labels))
    downweighted = edges.filter(pl.col("label").is_in(down_labels))

    supervised.write_parquet(out_dir / "taskg_supervised.parquet")
    message.write_parquet(out_dir / "taskg_message_passing.parquet")
    downweighted.write_parquet(out_dir / "taskg_downweighted.parquet")

    counts = {
        "supervised": supervised.height,
        "message_passing": message.height,
        "downweighted": downweighted.height,
        "total": edges.height,
    }
    log.info("taskg_subset_built", **counts)
    return counts


def write_report(cfg: TaskGConfig, args: argparse.Namespace, counts: dict | None) -> None:
    if counts is None:
        counts = {}
        for role in ("supervised", "message_passing", "downweighted"):
            p = args.output_dir / f"taskg_{role}.parquet"
            counts[role] = pl.read_parquet(p).height if p.exists() else 0

    def rows(entries: list[RoleEntry]) -> str:
        return "\n".join(f"| {e.label} | {e.n_edges:,} | {e.note} |" for e in entries)

    text = f"""<!-- track_scope: benchmark -->
# M2B Task-G relation curation on OptimusKG (spec §7.4)

**Status: {cfg.status}** — the curation list below is a DRAFT. BioPathNet
Task-G runs must not start until the owner signs off (§0.3.9).

## Supervised task edges (prediction target)

| label | edges | note |
|---|---|---|
{rows(cfg.supervised)}

## Message-passing context (never supervised)

| label | edges | note |
|---|---|---|
{rows(cfg.message_passing)}

## Downweighted by default (ablation only)

| label | edges | note |
|---|---|---|
{rows(cfg.downweighted)}

## Excluded

| label | edges | note |
|---|---|---|
{rows(cfg.excluded)}

## Built subset

- Supervised: {counts.get("supervised", 0):,} edges → `taskg_supervised.parquet`
- Message passing: {counts.get("message_passing", 0):,} edges → `taskg_message_passing.parquet`
- Downweighted (kept aside): {counts.get("downweighted", 0):,} edges → `taskg_downweighted.parquet`

## Reproduction

- Config: `{args.config}`
- Run: `python scripts/m2_build_taskg.py --config {args.config} --resume`
"""
    args.report_path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
