#!/usr/bin/env python
"""m2_download_anyburl — register + download the AnyBURL jar (spec §6.3).

Usage:
    python scripts/m2_download_anyburl.py --config configs/m2_anyburl.yaml [--dry-run] [--resume]

Same download contract as M1: registry-driven, resumable, SHA-256-verified.
"""

import argparse
import sys
from pathlib import Path

from pydantic import BaseModel

from sina.config import load_config
from sina.download import download_asset
from sina.logging import configure_logging, get_logger
from sina.manifests import Run
from sina.registry import RegistryAsset, load_registry, save_registry

log = get_logger(__name__)


class DownloadConfig(BaseModel):
    track: str
    registry_path: Path
    dest_dir: Path
    manifest_path: Path
    assets: list[RegistryAsset]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, DownloadConfig)

    registry = load_registry(cfg.registry_path)
    with Run(
        script="m2_download_anyburl",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config), "track": cfg.track},
    ) as run:
        for asset in cfg.assets:
            asset.track = cfg.track
            stage = f"download_{asset.name}"
            if args.resume and stage in run.done_stages:
                log.info("stage_skipped", stage=stage)
                continue
            with run.stage(stage, outputs=[cfg.dest_dir / asset.filename, cfg.registry_path]):
                registry.upsert(asset)
                save_registry(registry, cfg.registry_path)
                stored = registry.get(asset.name)
                download_asset(stored, cfg.dest_dir, dry_run=args.dry_run)
                registry.upsert(stored)
                save_registry(registry, cfg.registry_path)
    log.info("download_run_done", track=cfg.track)
    return 0


if __name__ == "__main__":
    sys.exit(main())
