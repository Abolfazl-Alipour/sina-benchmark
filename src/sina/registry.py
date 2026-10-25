"""Asset registry: data/raw/REGISTRY.json schema, hashing, verification.

Plan §0.3.6: SHA-256-verified downloads vs data/raw/REGISTRY.json; licenses verbatim.
"""

import hashlib
import json
import os
import tempfile
from datetime import UTC, datetime
from enum import StrEnum
from pathlib import Path

from pydantic import BaseModel, Field

CHUNK_SIZE = 8 * 1024 * 1024  # 8 MiB


class AssetStatus(StrEnum):
    PENDING = "pending_download"
    VERIFIED = "verified"


class RegistryAsset(BaseModel):
    """One downloadable asset tracked by the registry."""

    name: str
    version: str
    url: str
    filename: str
    sha256: str | None = None
    size_bytes: int | None = None
    license_verbatim: str = ""
    track: str = "shared"  # "R" | "D" | "C" | "shared"
    status: AssetStatus = AssetStatus.PENDING
    downloaded_at: datetime | None = None
    notes: str = ""


class Registry(BaseModel):
    schema_version: int = 1
    assets: list[RegistryAsset] = Field(default_factory=list)

    def get(self, name: str) -> RegistryAsset:
        for asset in self.assets:
            if asset.name == name:
                return asset
        raise KeyError(f"asset {name!r} not in registry")

    def upsert(self, asset: RegistryAsset) -> None:
        self.assets = [a for a in self.assets if a.name != asset.name]
        self.assets.append(asset)


def apply_pin(asset: RegistryAsset, registry: "Registry") -> RegistryAsset:
    """Carry a previously verified registry pin (sha256/size) onto a config asset.

    Config files list assets without hashes (hashes are recorded by the download
    scripts on first fetch). A bare upsert would clobber the stored pin with
    None, silently disabling verification on re-runs — this let an HTTP 202
    bot-check page get marked "verified" as a 0-byte PrimeKG on 2026-08-15.
    The pin is only carried when name AND version match and the incoming asset
    does not declare its own hash.
    """
    if asset.sha256 is not None:
        return asset
    try:
        existing = registry.get(asset.name)
    except KeyError:
        return asset
    if existing.version == asset.version and existing.sha256 is not None:
        asset.sha256 = existing.sha256
        asset.size_bytes = existing.size_bytes
    return asset


def sha256_file(path: Path | str, chunk_size: int = CHUNK_SIZE) -> str:
    """Streaming SHA-256 of a file."""
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_registry(path: Path | str) -> Registry:
    path = Path(path)
    if not path.exists():
        return Registry()
    with path.open("r", encoding="utf-8") as fh:
        return Registry.model_validate(json.load(fh))


def save_registry(registry: Registry, path: Path | str) -> None:
    """Atomic write: tmp file + os.replace, so a crash never truncates the registry."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = registry.model_dump_json(indent=2)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=path.name, suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(payload)
            fh.write("\n")
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def verify_asset(asset: RegistryAsset, file_path: Path | str) -> bool:
    """Check a downloaded file against the registry entry's SHA-256.

    Returns True when the asset has no recorded hash yet (first download —
    the hash gets recorded from this file) or when the hash matches.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        return False
    if asset.sha256 is None:
        return True
    return sha256_file(file_path) == asset.sha256


def mark_verified(asset: RegistryAsset, file_path: Path | str) -> RegistryAsset:
    """Fill in hash/size/timestamp after a successful download + verification."""
    file_path = Path(file_path)
    asset.sha256 = sha256_file(file_path)
    asset.size_bytes = file_path.stat().st_size
    asset.downloaded_at = datetime.now(UTC)
    asset.status = AssetStatus.VERIFIED
    return asset
