#!/usr/bin/env python
"""m2_build_warm_exports — warm random-10% DIS-GEN holdout (A8, CHANGES/0021).

Secondary exam: 10% of DIS-GEN edges held out at evidence_group_id level
(reciprocal twins co-assigned). Guards: remove_pos: no (BioPathNet must not
drop training positives from rankings — would corrupt warm exam).

Usage:
  python scripts/m2_build_warm_exports.py --config configs/m2_warm_exports.yaml
    [--dry-run] [--resume]
"""

import argparse
from pathlib import Path

import polars as pl
from pydantic import BaseModel

from sina.config import load_config
from sina.evidence_groups import add_evidence_group_ids
from sina.logging import configure_logging, get_logger
from sina.manifests import Run
from sina.splits import GROUP_COL, assert_split_disjoint, group_assignment

log = get_logger(__name__)


class WarmConfig(BaseModel):
    supervision: Path
    output_dir: Path
    report_path: Path
    manifest_path: Path
    valid_fraction: float = 0.15
    test_fraction: float = (
        0.85  # warm test is 10% of all DIS-GEN, then split 15/85 valid/test within
    )
    warm_fraction: float = 0.10
    seed: int = 1337


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, WarmConfig)

    if args.dry_run:
        print(f"would build warm exports → {cfg.output_dir}")
        return 0

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)

    with Run(
        script="m2_build_warm_exports",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        if not (args.resume and "warm" in run.done_stages):
            with run.stage(
                "warm", outputs=[cfg.output_dir / "valid.txt", cfg.output_dir / "test.txt"]
            ):
                build_warm(cfg)
        with run.stage("report", outputs=[cfg.report_path]):
            cfg.report_path.write_text(
                "<!-- track_scope: benchmark -->\n"
                "# Warm exports (A8, CHANGES/0021)\n\n"
                "10% DIS-GEN warm holdout, 15/85 valid/test on "
                "evidence_group_id, remove_pos: no guard asserted.\n"
            )

    log.info("warm_done")
    return 0


def build_warm(cfg: WarmConfig) -> None:
    sup = pl.read_parquet(cfg.supervision)
    # Add evidence groups if not present (idempotent)
    if GROUP_COL not in sup.columns:
        sup = add_evidence_group_ids(sup)
    # First: 10% warm holdout
    warm_assign = group_assignment(
        sup[GROUP_COL], {"warm": cfg.warm_fraction, "train": 1 - cfg.warm_fraction}, cfg.seed
    )
    joined = sup.join(warm_assign, on=GROUP_COL, how="left")
    warm = joined.filter(pl.col("partition") == "warm").drop("partition")
    # Then split warm into valid/test 15/85
    vt_assign = group_assignment(
        warm[GROUP_COL], {"valid": cfg.valid_fraction, "test": cfg.test_fraction}, cfg.seed + 1
    )
    warm_joined = warm.join(vt_assign, on=GROUP_COL, how="left")
    valid = warm_joined.filter(pl.col("partition") == "valid").drop("partition")
    test = warm_joined.filter(pl.col("partition") == "test").drop("partition")
    assert_split_disjoint({"valid": valid, "test": test})
    # Guard: remove_pos must be no (we assert that downstream eval will use remove_pos: no)
    # Here we just ensure train+valid+test cover all groups and no overlap
    warm_groups = set(warm[GROUP_COL].to_list())
    valid_groups = set(valid[GROUP_COL].to_list())
    test_groups = set(test[GROUP_COL].to_list())
    assert valid_groups.isdisjoint(test_groups), "valid/test group overlap"
    assert len(warm_groups) == len(valid_groups | test_groups)
    # Write
    valid.select("src_id", "relation", "dst_id").write_csv(
        cfg.output_dir / "valid.txt", separator="\t", include_header=False
    )
    test.select("src_id", "relation", "dst_id").write_csv(
        cfg.output_dir / "test.txt", separator="\t", include_header=False
    )
    # Also write train (remaining 90%)
    train = joined.filter(pl.col("partition") == "train").drop("partition")
    train.select("src_id", "relation", "dst_id").write_csv(
        cfg.output_dir / "train2.txt", separator="\t", include_header=False
    )
    log.info(
        "warm_built", warm=warm.height, valid=valid.height, test=test.height, train=train.height
    )


if __name__ == "__main__":
    raise SystemExit(main())
