#!/usr/bin/env python
"""m1_download — register + download assets for one track (spec §0.3.5/§0.3.6).

Usage:
    python scripts/m1_download.py --config configs/m1_trackr.yaml [--dry-run] [--resume]

Asset metadata comes from the YAML config; entries are upserted into
data/raw/REGISTRY.json, then downloaded via sina.download.download_asset
(resumable, SHA-256-verified, registry-updating). Never download science data
any other way (HANDOFF §3).
"""

import argparse
import sys
from pathlib import Path

from pydantic import BaseModel

from sina.config import load_config
from sina.download import download_asset
from sina.logging import configure_logging, get_logger
from sina.manifests import Run
from sina.registry import RegistryAsset, apply_pin, load_registry, save_registry

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
        script="m1_download",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config), "track": cfg.track},
    ) as run:
        for asset in cfg.assets:
            asset.track = cfg.track
            asset = apply_pin(asset, registry)
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
