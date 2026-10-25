#!/usr/bin/env python
"""m2_run_baselines — degree/PPR/null-control pilots (spec §6.3).

Stages (resumable; --dry-run prints the plan):
  For each configured split:
    1. load split train/test edges + candidate node pool
    2. degree baseline — score every candidate gene for every query disease
    3. PPR baseline — personalized PageRank from each query disease
    4. degree-matched nulls — sampled negatives per true edge
    5. relation-shuffled nulls — permuted relation labels
  6. report — manifest + counts → reports/m2a/baselines_report.md

Usage: python scripts/m2_run_baselines.py --config configs/m2_baselines.yaml
           [--dry-run] [--resume]
"""

import argparse
import json
import sys
import time
from pathlib import Path

import polars as pl
from pydantic import BaseModel

from sina.baselines import (
    degree_baseline,
    degree_matched_nulls,
    ppr_baseline,
    relation_shuffled_nulls,
)
from sina.config import load_config
from sina.logging import configure_logging, get_logger
from sina.manifests import Run

log = get_logger(__name__)


class SplitSpec(BaseModel):
    track: str
    query_col: str
    query_type: str
    candidate_type: str
    eval_relation: str | None = None
    relation_col: str | None = None
    eval_label: str | None = None
    label_col: str | None = None


class PprConfig(BaseModel):
    alpha: float = 0.15
    max_iter: int = 100
    tol: float = 1.0e-6


class DegreeMatchedConfig(BaseModel):
    n_negatives: int = 10
    bins: int = 20
    seed: int = 1337


class RelationShuffledConfig(BaseModel):
    seed: int = 1337


class BaselinesConfig(BaseModel):
    splits_dir: Path
    trackr_nodes: Path
    trackd_nodes: Path
    splits: dict[str, SplitSpec]
    baselines: list[str]
    null_controls: list[str]
    ppr: PprConfig
    degree_matched: DegreeMatchedConfig
    relation_shuffled: RelationShuffledConfig
    output_dir: Path
    report_path: Path
    manifest_path: Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, BaselinesConfig)

    if args.dry_run:
        print(f"would score splits {list(cfg.splits)} with {cfg.baselines} → {cfg.output_dir}")
        return 0

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)

    track_nodes = {
        "R": pl.read_parquet(cfg.trackr_nodes),
        "D": pl.read_parquet(cfg.trackd_nodes),
    }

    with Run(
        script="m2_run_baselines",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        for split_name in cfg.splits:
            if args.resume and split_name in run.done_stages:
                log.info("split_skipped", split=split_name)
                continue
            with run.stage(split_name):
                spec = cfg.splits[split_name]
                result = score_split(cfg, spec, split_name, track_nodes[spec.track])
                write_split_outputs(cfg.output_dir, split_name, result)

        with run.stage("report", outputs=[cfg.report_path]):
            write_report(cfg)

    log.info("baselines_done")
    return 0


def score_split(
    cfg: BaselinesConfig,
    spec: SplitSpec,
    split_name: str,
    nodes: pl.DataFrame,
) -> dict:
    split_dir = cfg.splits_dir / split_name
    train = pl.read_parquet(split_dir / "train.parquet")
    test = pl.read_parquet(split_dir / "test.parquet")

    candidates = set(nodes.filter(pl.col("node_type") == spec.candidate_type)["node_id"].to_list())

    # Reduce test to eval edges only (family/txgnn already produce eval-only test,
    # but be explicit).
    if spec.eval_relation and spec.relation_col:
        test = test.filter(pl.col(spec.relation_col) == spec.eval_relation)
    if spec.eval_label and spec.label_col:
        test = test.filter(pl.col(spec.label_col) == spec.eval_label)

    result: dict[str, object] = {
        "train_rows": train.height,
        "test_rows": test.height,
        "n_candidates": len(candidates),
        "n_queries_in_test": test[spec.query_col].n_unique(),
    }

    # Degree baseline
    if "degree" in cfg.baselines:
        t0 = time.perf_counter()
        deg = degree_baseline(train, candidates)
        # Cross-join query diseases in test with candidate degrees
        query_ids = test[spec.query_col].unique().to_frame("query_id")
        deg_renamed = deg.rename({"node_id": "candidate_id", "degree": "score"})
        deg_scores = query_ids.join(deg_renamed, how="cross")
        result["degree_time_s"] = round(time.perf_counter() - t0, 3)
        result["degree_scores"] = deg_scores

    # PPR baseline
    if "ppr" in cfg.baselines:
        t0 = time.perf_counter()
        query_nodes = test[spec.query_col].unique().to_list()
        ppr_scores = ppr_baseline(
            train,
            query_nodes=query_nodes,
            candidate_nodes=candidates,
            alpha=cfg.ppr.alpha,
            max_iter=cfg.ppr.max_iter,
            tol=cfg.ppr.tol,
        )
        result["ppr_time_s"] = round(time.perf_counter() - t0, 3)
        result["ppr_scores"] = ppr_scores.rename({"ppr_score": "score"})

    # Null controls
    if "degree_matched" in cfg.null_controls:
        t0 = time.perf_counter()
        nulls = degree_matched_nulls(
            train,
            test,
            candidate_pool=candidates,
            n_negatives=cfg.degree_matched.n_negatives,
            bins=cfg.degree_matched.bins,
            seed=cfg.degree_matched.seed,
        )
        result["degree_matched_time_s"] = round(time.perf_counter() - t0, 3)
        result["degree_matched_nulls"] = nulls

    if "relation_shuffled" in cfg.null_controls:
        t0 = time.perf_counter()
        shuffled = relation_shuffled_nulls(test, seed=cfg.relation_shuffled.seed)
        result["relation_shuffled_time_s"] = round(time.perf_counter() - t0, 3)
        result["relation_shuffled_nulls"] = shuffled

    return result


def write_split_outputs(output_dir: Path, split_name: str, result: dict) -> None:
    base = output_dir / split_name
    base.mkdir(parents=True, exist_ok=True)
    meta: dict[str, object] = {"split_name": split_name}

    for key in ("degree_scores", "ppr_scores"):
        if key in result:
            path = base / f"{key}.parquet"
            result[key].write_parquet(path)
            meta[key] = {"rows": result[key].height, "path": str(path)}

    for key in ("degree_matched_nulls", "relation_shuffled_nulls"):
        if key in result:
            path = base / f"{key}.parquet"
            result[key].write_parquet(path)
            meta[key] = {"rows": result[key].height, "path": str(path)}

    for timing_key in (
        "degree_time_s",
        "ppr_time_s",
        "degree_matched_time_s",
        "relation_shuffled_time_s",
    ):
        if timing_key in result:
            meta[timing_key] = result[timing_key]

    meta["train_rows"] = result["train_rows"]
    meta["test_rows"] = result["test_rows"]
    meta["n_candidates"] = result["n_candidates"]
    meta["n_queries_in_test"] = result["n_queries_in_test"]

    meta_path = base / "metadata.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    log.info("baseline_outputs_written", split=split_name, **meta)


def write_report(cfg: BaselinesConfig) -> None:
    rows = []
    timings = []
    for split_name in cfg.splits:
        meta_path = cfg.output_dir / split_name / "metadata.json"
        if not meta_path.exists():
            continue
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        rows.append(
            f"| {split_name} | {meta['train_rows']:,} | {meta['test_rows']:,} | "
            f"{meta['n_candidates']:,} | {meta['n_queries_in_test']:,} | "
            f"{'degree_scores' in meta} | {'ppr_scores' in meta} | "
            f"{'degree_matched_nulls' in meta} | {'relation_shuffled_nulls' in meta} |"
        )
        timings.append(
            f"| {split_name} | "
            f"{meta.get('degree_time_s', '—')} | "
            f"{meta.get('ppr_time_s', '—')} | "
            f"{meta.get('degree_matched_time_s', '—')} | "
            f"{meta.get('relation_shuffled_time_s', '—')} |"
        )

    text = f"""<!-- track_scope: benchmark -->
# M2A-3 Baseline + null-control pilots (spec §6.3)

**Status:** DRAFT — KGE (TransE/ComplEx/RotatE/InductiveNodePiece) and AnyBURL
pilots are not yet run. This report covers degree, PPR, and the two null
controls on the small-query splits (family holdouts + TxGNN neurodegenerative
area) to measure popularity bias and structural baselines.

## Splits scored

| split | train rows | test rows | candidates | queries | degree | PPR | deg-null | rel-null |
|---|---|---|---|---|---|---|---|---|
{chr(10).join(rows)}

## Timing (seconds)

| split | degree | PPR | deg-matched null | rel-shuffled null |
|---|---|---|---|---|
{chr(10).join(timings)}

## Design notes

- **Degree baseline:** ranks candidate genes by total training degree. Expected
  to be a strong popularity-biased baseline on dense disease-gene graphs.
- **PPR baseline:** personalized PageRank from each query disease on the
  symmetric-normalized training graph (restart alpha = {cfg.ppr.alpha}).
- **Degree-matched null:** for each true edge, sample {cfg.degree_matched.n_negatives}
  negative gene(s) from the same degree bin as the true answer.
- **Relation-shuffled null:** permutes relation labels on the test edge set,
  preserving the structural degree distribution but breaking semantics.

## Reproduction

- Config: `{cfg.output_dir.parent.parent / "configs/m2_baselines.yaml"}`
- Manifest: `{cfg.manifest_path}`
- Output: `{cfg.output_dir}`
- Re-run: `python scripts/m2_run_baselines.py --config configs/m2_baselines.yaml --resume`
"""
    cfg.report_path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
