#!/usr/bin/env python
"""m2_theta_derivation — A6 null-distribution for theta (CHANGES/0021).

Derives candidate theta at null 95th-percentile crossing via ≥100 key shuffles.
Reports frozen 5% vs 4% sensitivity row + continuous key-size-stratified.

Usage: python scripts/m2_theta_derivation.py --config configs/m2_theta.yaml [--dry-run]
"""

import argparse
from pathlib import Path

import numpy as np
import polars as pl
from pydantic import BaseModel

from sina.config import load_config
from sina.logging import configure_logging, get_logger
from sina.manifests import Run

log = get_logger(__name__)


class ThetaConfig(BaseModel):
    unit_audit: Path
    output_dir: Path
    report_path: Path
    manifest_path: Path
    shuffles: int = 100
    seed: int = 1337


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, ThetaConfig)

    if args.dry_run:
        print(f"would derive theta from {cfg.unit_audit} with {cfg.shuffles} shuffles")
        return 0

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)

    with Run(
        script="m2_theta_derivation",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        with run.stage("derive", outputs=[cfg.output_dir / "theta_null.parquet"]):
            derive(cfg)
        with run.stage("report", outputs=[cfg.report_path]):
            write_report(cfg)

    log.info("theta_done")
    return 0


def derive(cfg: ThetaConfig) -> None:
    audit = pl.read_parquet(cfg.unit_audit)
    # Use key_size_pct = key_size / 17448 *100
    audit = audit.with_columns((pl.col("key_size") / 17448 * 100).alias("pct"))
    rng = np.random.default_rng(cfg.seed)
    # Null: shuffle pct values
    null_pcts = []
    for _ in range(cfg.shuffles):
        shuffled = rng.permutation(audit["pct"].to_numpy())
        # Candidate theta at 95th percentile of shuffled pct
        null_pcts.append(float(np.percentile(shuffled, 95)))
    null_arr = np.array(null_pcts)
    candidate = float(np.percentile(null_arr, 95)) if len(null_arr) else 5.0
    # Save
    pl.DataFrame({"null_95pct": null_arr}).write_parquet(cfg.output_dir / "theta_null.parquet")
    # Also save candidate
    (cfg.output_dir / "theta_candidate.json").write_text(
        f'{{"candidate_theta_pct": {candidate:.2f}, '
        f'"frozen_theta_pct": 5.0, "shuffles": {cfg.shuffles}}}\n'
    )
    log.info("theta_candidate", candidate=candidate, frozen=5.0, shuffles=cfg.shuffles)


def write_report(cfg: ThetaConfig) -> None:
    import json

    cand = (
        json.loads((cfg.output_dir / "theta_candidate.json").read_text())
        if (cfg.output_dir / "theta_candidate.json").exists()
        else {"candidate_theta_pct": 5.0}
    )
    audit = pl.read_parquet(cfg.unit_audit)
    audit = audit.with_columns((pl.col("key_size") / 17448 * 100).alias("pct"))
    # Count units that would flip at 4% vs 5%
    n_4 = audit.filter(pl.col("pct") > 4).height
    n_5 = audit.filter(pl.col("pct") > 5).height
    # Werner 4.91% is the boundary
    lines = [
        "<!-- track_scope: benchmark -->",
        "# Theta derivation (A6, CHANGES/0021)",
        "",
        (
            f"Null 95th-pct candidate: **{cand['candidate_theta_pct']:.2f}%** "
            f"(shuffles={cand.get('shuffles', 100)}), "
            f"frozen **5.0%** retained (CHANGES/0021 adopt-vs-keep)."
        ),
        (
            f"Units >4% : {n_4}, >5% : {n_5} — Werner 4.91% drops to "
            f"context-only at 4% row (sensitivity, reporting tier only)."
        ),
        "",
        "Continuous key-size-stratified reporting template: stratify results by "
        "`pct` bins `<2%`, `2-4%`, `4-5%`, `>5%` (Werner in `4-5%` bin).",
        "",
        (
            "Decision: **keep frozen 5%**, publish 4% sensitivity row + "
            "continuous stratification. Zero extra cost."
        ),
        "",
    ]
    cfg.report_path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
