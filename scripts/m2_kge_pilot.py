#!/usr/bin/env python
"""m2_kge_pilot — PyKEEN KGE pilot on Track-R family holdout (spec §6.3).

Stages (resumable; --dry-run prints the plan):
  1. load — split train/test edges + candidate node pool
  2. train — train the configured KGE model for a small number of epochs
  3. score — score candidate genes for each query disease
  4. report — timing + score counts → reports/m2a/kge_pilot_report.md

Usage: python scripts/m2_kge_pilot.py --config configs/m2_kge_pilot.yaml
           [--dry-run] [--resume]
"""

import argparse
import json
import sys
import time
from pathlib import Path

import polars as pl
from pydantic import BaseModel

from sina.config import load_config
from sina.kge import train_and_score_kge
from sina.logging import configure_logging, get_logger
from sina.manifests import Run

log = get_logger(__name__)


class SplitSpec(BaseModel):
    query_col: str
    query_type: str
    candidate_type: str
    eval_relation: str
    relation_col: str


class ModelSpec(BaseModel):
    name: str
    embedding_dim: int
    epochs: int
    seed: int


class KgePilotConfig(BaseModel):
    splits_dir: Path
    output_dir: Path
    report_path: Path
    manifest_path: Path
    splits: dict[str, SplitSpec]
    models: list[ModelSpec]
    device: str = "cpu"
    rss_limit_gb: float = 40.0
    subsample_prefix: str | None = None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, KgePilotConfig)

    if args.dry_run:
        print(f"would run KGE pilot on {list(cfg.splits)} → {cfg.output_dir}")
        return 0

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)

    with Run(
        script="m2_kge_pilot",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
        rss_limit_gb=cfg.rss_limit_gb,
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

    log.info("kge_pilot_done")
    return 0


def run_split(cfg: KgePilotConfig, split_name: str) -> dict:
    spec = cfg.splits[split_name]
    split_dir = cfg.splits_dir / split_name
    train = pl.read_parquet(split_dir / "train.parquet")
    if cfg.subsample_prefix:
        train = train.filter(pl.col("evidence_group_id").str.starts_with(cfg.subsample_prefix))
    test = pl.read_parquet(split_dir / "test.parquet")
    nodes_path = Path("data/processed/m1/trackr/nodes.parquet")
    nodes = pl.read_parquet(nodes_path)

    # Keep only eval edges with disease query.
    test = test.filter(pl.col(spec.relation_col) == spec.eval_relation)
    query_nodes = test[spec.query_col].unique().to_list()
    candidates = set(nodes.filter(pl.col("node_type") == spec.candidate_type)["node_id"].to_list())

    result: dict[str, object] = {
        "split_name": split_name,
        "train_rows": train.height,
        "test_rows": test.height,
        "n_queries": len(query_nodes),
        "n_candidates": len(candidates),
        "models": {},
    }

    for model_spec in cfg.models:
        t0 = time.perf_counter()
        scores = train_and_score_kge(
            train,
            query_nodes=query_nodes,
            candidate_nodes=candidates,
            relation=spec.eval_relation,
            model=model_spec.name,
            embedding_dim=model_spec.embedding_dim,
            epochs=model_spec.epochs,
            seed=model_spec.seed,
            device=cfg.device,
        )
        elapsed = time.perf_counter() - t0
        model_key = f"{model_spec.name}_d{model_spec.embedding_dim}_e{model_spec.epochs}"
        result["models"][model_key] = {
            "scores": scores,
            "time_s": round(elapsed, 3),
            "n_scores": scores.height,
        }
        log.info(
            "kge_model_done",
            split=split_name,
            model=model_key,
            time_s=elapsed,
            n_scores=scores.height,
        )
    return result


def write_split_outputs(output_dir: Path, split_name: str, result: dict) -> None:
    base = output_dir / split_name
    base.mkdir(parents=True, exist_ok=True)
    meta = {
        "split_name": split_name,
        "train_rows": result["train_rows"],
        "test_rows": result["test_rows"],
        "n_queries": result["n_queries"],
        "n_candidates": result["n_candidates"],
        "models": {},
    }
    for model_key, m in result["models"].items():
        path = base / f"{model_key}_scores.parquet"
        m["scores"].write_parquet(path)
        meta["models"][model_key] = {
            "time_s": m["time_s"],
            "n_scores": m["n_scores"],
            "path": str(path),
        }
    meta_path = base / "metadata.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    log.info("kge_outputs_written", split=split_name)


def write_report(cfg: KgePilotConfig) -> None:
    rows = []
    for split_name in cfg.splits:
        meta_path = cfg.output_dir / split_name / "metadata.json"
        if not meta_path.exists():
            continue
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        for model_key, m in meta["models"].items():
            rows.append(
                f"| {split_name} | {model_key} | {meta['n_queries']:,} | "
                f"{meta['n_candidates']:,} | {m['n_scores']:,} | {m['time_s']} |"
            )

    text = f"""<!-- track_scope: benchmark -->
# M2A-3 KGE pilot report (spec §6.3)

**Status:** DRAFT — single-epoch pilot on Track-R family holdout to measure
PyKEEN throughput; full-scale multi-epoch runs are estimated from this pilot.

## Pilot results

| split | model | queries | candidates | scores | time s |
|---|---|---|---|---|---|
{chr(10).join(rows)}

## Full-scale estimate

- Pilot: {cfg.models[0].epochs} epoch(s), embedding_dim={cfg.models[0].embedding_dim}.
- Full-scale target: TransE-256/ComplEx-256/RotatE-512, ≥5 seeds, ≥50 epochs.
- Estimated time per full run: pilot_time x (256/64) x 50 x 5 seeds.
- At pilot throughput, full-scale KGE is a multi-hour job and should be run
  on a scheduled long-running job or cloud CPU.

## Reproduction

- Config: `configs/m2_kge_pilot.yaml`
- Manifest: `{cfg.manifest_path}`
- Output: `{cfg.output_dir}`
- Re-run: `python scripts/m2_kge_pilot.py --config configs/m2_kge_pilot.yaml --resume`
"""
    cfg.report_path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
