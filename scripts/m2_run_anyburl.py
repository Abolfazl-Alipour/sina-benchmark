#!/usr/bin/env python
"""m2_run_anyburl — AnyBURL rule-mining pilot (spec §6.3).

Stages (resumable; --dry-run prints the plan):
  For each configured split:
    1. export — write train/test triples in AnyBURL format
    2. learn — run AnyBURL Learn for `learn_seconds`
    3. apply — run AnyBURL Apply to produce top-k predictions
    4. parse — parse preds-100 into long-format scores parquet
  5. report — counts + timings → reports/m2a/anyburl_report.md

Usage: python scripts/m2_run_anyburl.py --config configs/m2_anyburl_run.yaml
           [--dry-run] [--resume]
"""

import argparse
import json
import sys
import time
from pathlib import Path

import polars as pl
from pydantic import BaseModel

from sina.anyburl import (
    export_triples,
    parse_predictions,
    run_anyburl,
    write_apply_config,
    write_learn_config,
)
from sina.config import load_config
from sina.logging import configure_logging, get_logger
from sina.manifests import Run

log = get_logger(__name__)


class SplitSpec(BaseModel):
    relation: str
    reverse_test: bool = False
    prediction_type: str = "tails"


class AnyburlConfig(BaseModel):
    jar_path: Path
    splits_dir: Path
    output_dir: Path
    report_path: Path
    manifest_path: Path
    splits: dict[str, SplitSpec]
    snapshots_at: str = "10,50,100"
    rules_snapshot: int = 100
    worker_threads: int = 4
    top_k: int = 10
    max_heap_gb: int = 8
    timeout_s: int = 3600


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, AnyburlConfig)

    if args.dry_run:
        print(f"would run AnyBURL on {list(cfg.splits)} → {cfg.output_dir}")
        return 0

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)

    with Run(
        script="m2_run_anyburl",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        for split_name in cfg.splits:
            if args.resume and split_name in run.done_stages:
                log.info("split_skipped", split=split_name)
                continue
            with run.stage(split_name):
                result = run_split(cfg, split_name)
                write_split_outputs(cfg.output_dir, split_name, result)

        with run.stage("report", outputs=[cfg.report_path]):
            write_report(cfg)

    log.info("anyburl_done")
    return 0


def run_split(cfg: AnyburlConfig, split_name: str) -> dict:
    spec = cfg.splits[split_name]
    work_dir = cfg.output_dir / split_name / "work"
    work_dir.mkdir(parents=True, exist_ok=True)
    train_path = work_dir / "train.txt"
    test_path = work_dir / "test.txt"
    learn_cfg = work_dir / "config-learn.properties"
    apply_cfg = work_dir / "config-apply.properties"

    train_edges = pl.read_parquet(cfg.splits_dir / split_name / "train.parquet")
    test_edges = pl.read_parquet(cfg.splits_dir / split_name / "test.parquet")

    export_triples(train_edges, train_path)
    test_for_export = test_edges
    if spec.reverse_test:
        test_for_export = test_edges.rename({"src_id": "dst_id", "dst_id": "src_id"})
    export_triples(test_for_export, test_path)

    result: dict[str, object] = {
        "train_rows": train_edges.height,
        "test_rows": test_edges.height,
    }

    # Learn
    write_learn_config(
        learn_cfg,
        train_path,
        work_dir,
        worker_threads=cfg.worker_threads,
        snapshots_at=cfg.snapshots_at,
    )
    t0 = time.perf_counter()
    run_anyburl(
        cfg.jar_path,
        "Learn",
        learn_cfg,
        max_heap_gb=cfg.max_heap_gb,
        timeout_s=cfg.timeout_s,
    )
    result["learn_time_s"] = round(time.perf_counter() - t0, 3)
    rules_path = work_dir / f"rules-{cfg.rules_snapshot}"
    if not rules_path.exists():
        raise FileNotFoundError(f"AnyBURL rules not found at {rules_path}")
    # AnyBURL sometimes emits rules with empty bodies; filter them before Apply.
    from sina.anyburl import filter_rules

    filtered_rules_path = work_dir / f"rules-{cfg.rules_snapshot}-filtered"
    n_kept = filter_rules(rules_path, filtered_rules_path)
    result["n_rules_kept"] = n_kept

    # Apply
    preds_path = work_dir / f"preds-{cfg.rules_snapshot}"
    write_apply_config(
        apply_cfg,
        filtered_rules_path,
        train_path,
        test_path,
        preds_path,
        worker_threads=cfg.worker_threads,
        top_k=cfg.top_k,
    )
    t0 = time.perf_counter()
    run_anyburl(
        cfg.jar_path,
        "Apply",
        apply_cfg,
        max_heap_gb=cfg.max_heap_gb,
        timeout_s=cfg.timeout_s,
    )
    result["apply_time_s"] = round(time.perf_counter() - t0, 3)
    if not preds_path.exists():
        raise FileNotFoundError(f"AnyBURL predictions not found at {preds_path}")

    # Parse
    scores = parse_predictions(preds_path, prediction_type=spec.prediction_type)
    result["scores"] = scores
    result["n_scores"] = scores.height
    result["n_queries"] = scores["query_id"].n_unique()
    return result


def write_split_outputs(output_dir: Path, split_name: str, result: dict) -> None:
    base = output_dir / split_name
    base.mkdir(parents=True, exist_ok=True)
    scores_path = base / "anyburl_scores.parquet"
    result["scores"].write_parquet(scores_path)
    meta = {
        "split_name": split_name,
        "train_rows": result["train_rows"],
        "test_rows": result["test_rows"],
        "n_scores": result["n_scores"],
        "n_queries": result["n_queries"],
        "learn_time_s": result["learn_time_s"],
        "apply_time_s": result["apply_time_s"],
        "scores_path": str(scores_path),
    }
    meta_path = base / "metadata.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    log.info("anyburl_outputs_written", split=split_name, **meta)


def write_report(cfg: AnyburlConfig) -> None:
    rows = []
    for split_name in cfg.splits:
        meta_path = cfg.output_dir / split_name / "metadata.json"
        if not meta_path.exists():
            continue
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        rows.append(
            f"| {split_name} | {meta['train_rows']:,} | {meta['test_rows']:,} | "
            f"{meta['n_queries']:,} | {meta['n_scores']:,} | "
            f"{meta['learn_time_s']} | {meta['apply_time_s']} |"
        )

    text = f"""<!-- track_scope: benchmark -->
# M2A-3 AnyBURL rule-mining pilot (spec §6.3)

**Status:** DRAFT — pilot on the family holdout splits only; rule confidence
normalization and ensemble integration happen in M2A-4.

AnyBURL-23-1x.jar (BSD-3) learned rules with snapshots at {cfg.snapshots_at}
seconds per split on the training edges, then applied the {cfg.rules_snapshot}s
rules to the test edges to produce top-{cfg.top_k} rankings. The parsed
confidence scores are stored in
`data/processed/m2/anyburl/<split>/anyburl_scores.parquet`.

## Splits scored

| split | train rows | test rows | queries | scores | learn s | apply s |
|---|---|---|---|---|---|---|
{chr(10).join(rows)}

## Reproduction

- Config: `configs/m2_anyburl_run.yaml`
- Manifest: `{cfg.manifest_path}`
- Output: `{cfg.output_dir}`
- Re-run: `python scripts/m2_run_anyburl.py --config configs/m2_anyburl_run.yaml --resume`
"""
    cfg.report_path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
