#!/usr/bin/env python
"""m2_build_external_eval — build the pre-registered external eval set (spec §6.2).

**Current status:** skeleton. The external data source is not yet chosen
(PENDING_DATA_SOURCE). This script validates the disjointness contract and
produces an empty placeholder manifest until a source is registered.

Usage: python scripts/m2_build_external_eval.py --config configs/m2_external_eval.yaml
           [--dry-run] [--resume]
"""

import argparse
import json
import sys
from pathlib import Path

from pydantic import BaseModel

from sina.config import load_config
from sina.logging import configure_logging, get_logger
from sina.manifests import Run

log = get_logger(__name__)


class ExternalEvalConfig(BaseModel):
    name: str
    status: str
    description: str
    source: str | None
    license_verbatim: str | None
    cutoff_date: str | None
    relations: dict[str, str]
    output_dir: Path
    manifest_path: Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, ExternalEvalConfig)

    if args.dry_run:
        print(f"would build external eval set '{cfg.name}' (status: {cfg.status})")
        return 0

    cfg.output_dir.mkdir(parents=True, exist_ok=True)

    with Run(
        script="m2_build_external_eval",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        with run.stage("validate"):
            if cfg.status != "READY":
                log.warning(
                    "external_eval_not_ready",
                    status=cfg.status,
                    message="data source not registered; producing placeholder manifest",
                )

        with run.stage("build"):
            if cfg.status == "READY":
                build_external_eval(cfg)
            else:
                write_placeholder(cfg)

    log.info("external_eval_done")
    return 0


def build_external_eval(cfg: ExternalEvalConfig) -> None:
    """Build the external eval set once a source is registered.

    TODO(owner): register the external source, then:
      1. download + hash into REGISTRY
      2. parse into canonical (src_id, relation, dst_id)
      3. add evidence_group_ids
      4. assert disjointness vs. M1 train edges at evidence_group level
      5. write parquet + metadata
    """
    raise NotImplementedError("external data source not yet registered")


def write_placeholder(cfg: ExternalEvalConfig) -> None:
    """Write an empty placeholder manifest documenting the pending source."""
    placeholder = {
        "name": cfg.name,
        "status": cfg.status,
        "description": cfg.description,
        "source": cfg.source,
        "license_verbatim": cfg.license_verbatim,
        "cutoff_date": cfg.cutoff_date,
        "relations": cfg.relations,
        "note": "External eval set is pending owner data-source registration.",
    }
    path = cfg.output_dir / "placeholder.json"
    path.write_text(json.dumps(placeholder, indent=2), encoding="utf-8")
    log.info("placeholder_written", path=str(path))


if __name__ == "__main__":
    sys.exit(main())
