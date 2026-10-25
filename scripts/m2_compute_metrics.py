#!/usr/bin/env python
"""m2_compute_metrics — ranking metrics + statistics (spec §6.4/§6.5).

Stages (resumable; --dry-run prints the plan):
  For each configured split:
    1. load split edges + baseline score tables
    2. compute raw + filtered MRR/Hits@k for each baseline
    3. compute degree-stratified metrics
    4. compute null-control metrics
    5. write metrics JSON
  6. report — summary table → reports/m2a/metrics_report.md

Usage: python scripts/m2_compute_metrics.py --config configs/m2_metrics.yaml
           [--dry-run] [--resume]
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
from sina.metrics import degree_stratified_metrics, link_prediction_metrics

log = get_logger(__name__)


class SplitSpec(BaseModel):
    track: str
    query_col: str
    answer_col: str
    query_type: str
    filter_train: bool = True
    filter_valid: bool = True


class MetricsConfig(BaseModel):
    splits_dir: Path
    baselines_dir: Path
    trackr_nodes: Path
    trackd_nodes: Path
    output_dir: Path
    report_path: Path
    manifest_path: Path
    splits: dict[str, SplitSpec]
    baselines: list[str]
    null_controls: list[str]
    ks: list[int]
    n_degree_bins: int = 4


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, MetricsConfig)

    if args.dry_run:
        print(f"would compute metrics for {list(cfg.splits)} → {cfg.output_dir}")
        return 0

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)

    with Run(
        script="m2_compute_metrics",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        for split_name in cfg.splits:
            if args.resume and split_name in run.done_stages:
                log.info("split_skipped", split=split_name)
                continue
            with run.stage(split_name):
                result = compute_split_metrics(cfg, split_name)
                write_metrics(cfg.output_dir, split_name, result)

        with run.stage("report", outputs=[cfg.report_path]):
            write_report(cfg)

    log.info("metrics_done")
    return 0


def compute_split_metrics(cfg: MetricsConfig, split_name: str) -> dict:
    spec = cfg.splits[split_name]
    split_dir = cfg.splits_dir / split_name
    train = pl.read_parquet(split_dir / "train.parquet")
    valid = pl.read_parquet(split_dir / "valid.parquet")
    test = pl.read_parquet(split_dir / "test.parquet")

    # Keep only edges whose query node is a disease (disease→gene direction).
    nodes_path = cfg.trackr_nodes if spec.track == "R" else cfg.trackd_nodes
    disease_nodes = set(
        pl.read_parquet(nodes_path)
        .filter(pl.col("node_type") == spec.query_type)["node_id"]
        .to_list()
    )
    test = test.filter(pl.col(spec.query_col).is_in(disease_nodes))

    filter_frames = []
    if spec.filter_train:
        filter_frames.append(train)
    if spec.filter_valid:
        filter_frames.append(valid)
    filter_frames.append(test)
    filter_edges = pl.concat(filter_frames) if filter_frames else None
    # Restrict filter set to test queries for speed.
    if filter_edges is not None:
        test_queries = test[spec.query_col].unique().implode()
        filter_edges = filter_edges.filter(pl.col(spec.query_col).is_in(test_queries))

    result: dict[str, object] = {
        "split_name": split_name,
        "n_test_edges": test.height,
        "baselines": {},
        "nulls": {},
        "degree_stratified": {},
    }

    # Baselines
    for baseline in cfg.baselines:
        path = cfg.baselines_dir / split_name / f"{baseline}.parquet"
        if not path.exists():
            log.warning("baseline_missing", split=split_name, baseline=baseline)
            continue
        scores = pl.read_parquet(path)
        m = link_prediction_metrics(
            scores,
            test,
            filter_edges,
            query_col=spec.query_col,
            answer_col=spec.answer_col,
            ks=tuple(cfg.ks),
        )
        result["baselines"][baseline] = m
        # degree stratified
        strat = degree_stratified_metrics(
            scores,
            test,
            train,
            filter_edges,
            query_col=spec.query_col,
            answer_col=spec.answer_col,
            ks=tuple(cfg.ks),
            n_bins=cfg.n_degree_bins,
        )
        result["degree_stratified"][baseline] = strat

    # Null controls
    for null in cfg.null_controls:
        path = cfg.baselines_dir / split_name / f"{null}.parquet"
        if not path.exists():
            log.warning("null_missing", split=split_name, null=null)
            continue
        null_edges = pl.read_parquet(path)
        if null == "degree_matched_nulls":
            # columns: src_id, true_dst_id, null_dst_id, label
            if spec.track == "R":
                # disease→gene: query = true_dst_id (disease), answer = null_dst_id (gene)
                null_edges = null_edges.select(
                    pl.col("true_dst_id").alias(spec.query_col),
                    pl.col("null_dst_id").alias(spec.answer_col),
                )
            else:
                # Track-D already stores (disease, gene); nulls replace the gene.
                null_edges = null_edges.select(
                    pl.col("src_id").alias(spec.query_col),
                    pl.col("null_dst_id").alias(spec.answer_col),
                )
        deg_path = cfg.baselines_dir / split_name / "degree_scores.parquet"
        if not deg_path.exists():
            continue
        deg_scores = pl.read_parquet(deg_path)
        m = link_prediction_metrics(
            deg_scores,
            null_edges,
            filter_edges,
            query_col=spec.query_col,
            answer_col=spec.answer_col,
            ks=tuple(cfg.ks),
        )
        result["nulls"][null] = m

    return result


def write_metrics(output_dir: Path, split_name: str, result: dict) -> None:
    base = output_dir / split_name
    base.mkdir(parents=True, exist_ok=True)
    path = base / "metrics.json"
    path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    log.info("metrics_written", split=split_name)


def write_report(cfg: MetricsConfig) -> None:
    rows = []
    for split_name in cfg.splits:
        path = cfg.output_dir / split_name / "metrics.json"
        if not path.exists():
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        for name, m in data["baselines"].items():
            rows.append(
                f"| {split_name} | {name} | {m['n_queries']:,} | "
                f"{m['raw']['mrr']:.4f} | {m['filtered']['mrr']:.4f} | "
                f"{m['raw']['hits@1']:.4f} | {m['filtered']['hits@1']:.4f} | "
                f"{m['raw']['hits@10']:.4f} | {m['filtered']['hits@10']:.4f} |"
            )
        for name, m in data["nulls"].items():
            rows.append(
                f"| {split_name} | {name} | {m['n_queries']:,} | "
                f"{m['raw']['mrr']:.4f} | {m['filtered']['mrr']:.4f} | "
                f"{m['raw']['hits@1']:.4f} | {m['filtered']['hits@1']:.4f} | "
                f"{m['raw']['hits@10']:.4f} | {m['filtered']['hits@10']:.4f} |"
            )

    text = f"""<!-- track_scope: benchmark -->
# M2A-4 Metrics + statistics (spec §6.4/§6.5)

**Status:** DRAFT — covers degree/PPR/nulls only; KGE and AnyBURL (unusable on
zero-shot family) are excluded. Multi-seed statistics and paired Wilcoxon tests
require KGE/AnyBURL scores across seeds and are pending.

Metrics are computed with `sina.metrics.link_prediction_metrics` (filtered and
raw MRR/Hits@k; filter set = train+valid+test truths) and
`degree_stratified_metrics` (degree-quantile bins).

**Interpretation note (mandatory, §6.5):** zero-shot family holdout runs are
expected to score 50-80% below transductive settings. The degree baseline
captures popularity bias; the degree-matched null shows what a popularity-matched
random ranking achieves. A model must beat both to be considered useful.

## Headline metrics

| split | model | queries | raw MRR | filt MRR | raw h@1 | filt h@1 | raw h@10 | filt h@10 |
|---|---|---|---|---|---|---|---|---|
{chr(10).join(rows)}

## Reproduction

- Config: `configs/m2_metrics.yaml`
- Manifest: `{cfg.manifest_path}`
- Output: `{cfg.output_dir}`
- Re-run: `python scripts/m2_compute_metrics.py --config configs/m2_metrics.yaml --resume`
"""
    cfg.report_path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
