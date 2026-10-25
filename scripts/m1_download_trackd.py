#!/usr/bin/env python
"""m1_download_trackd — OptimusKG via the official client (spec §5.2.1).

Downloads the stratified per-type Parquet files with the pinned `optimuskg`
PyPI client (1.0.0), then checksums every file into data/raw/REGISTRY.json
(spec §0.3.6: SHA-256-verified assets, license recorded verbatim).

Usage: python scripts/m1_download_trackd.py --config configs/m1_trackd.yaml
           [--dry-run] [--resume]
"""

import argparse
import sys
from pathlib import Path

import optimuskg
from optimuskg import _dataverse
from pydantic import BaseModel

from sina.config import load_config
from sina.logging import configure_logging, get_logger
from sina.manifests import Run
from sina.registry import (
    RegistryAsset,
    apply_pin,
    load_registry,
    mark_verified,
    save_registry,
    sha256_file,
)

log = get_logger(__name__)


class TrackDDownloadConfig(BaseModel):
    track: str
    registry_path: Path
    cache_dir: Path
    manifest_path: Path
    dataset_version: str
    license_verbatim: str
    files: list[str]  # Dataverse-relative paths, e.g. "nodes/gene.parquet"


def asset_name(rel_path: str) -> str:
    return "optimuskg_" + rel_path.removesuffix(".parquet").replace("/", "_")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, TrackDDownloadConfig)
    optimuskg.set_cache_dir(cfg.cache_dir)

    registry = load_registry(cfg.registry_path)
    with Run(
        script="m1_download_trackd",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config), "track": cfg.track},
    ) as run:
        for rel in cfg.files:
            name = asset_name(rel)
            stage = f"download_{name}"
            if args.resume and stage in run.done_stages:
                log.info("stage_skipped", stage=stage)
                continue
            with run.stage(stage, outputs=[cfg.registry_path]):
                entry = _dataverse.resolve(rel)
                asset = RegistryAsset(
                    name=name,
                    version=cfg.dataset_version,
                    url=f"{optimuskg.get_server()}/api/access/datafile/{entry.id}",
                    filename=rel,
                    license_verbatim=cfg.license_verbatim,
                    track=cfg.track,
                )
                if args.dry_run:
                    log.info("dry_run_would_download", name=name, url=asset.url)
                    continue
                local = optimuskg.get_file(rel)
                asset = apply_pin(asset, registry)
                if asset.sha256 is not None and sha256_file(local) != asset.sha256:
                    raise ValueError(
                        f"SHA-256 mismatch for {name!r}: cached {local} does not match "
                        "the registry pin — refusing to mark verified"
                    )
                mark_verified(asset, local)
                registry.upsert(asset)
                save_registry(registry, cfg.registry_path)
                log.info("asset_verified", name=name, sha256=asset.sha256)
    log.info("download_run_done", track=cfg.track, n_files=len(cfg.files))
    return 0


if __name__ == "__main__":
    sys.exit(main())
