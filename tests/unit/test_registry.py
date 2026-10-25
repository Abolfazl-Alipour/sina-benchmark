"""Unit tests for sina.registry — hashing, roundtrip, verification."""

import hashlib

import pytest

from sina.registry import (
    AssetStatus,
    Registry,
    RegistryAsset,
    load_registry,
    mark_verified,
    save_registry,
    sha256_file,
    verify_asset,
)


def _asset(name: str = "primekg", sha256: str | None = None) -> RegistryAsset:
    return RegistryAsset(
        name=name,
        version="v1",
        url="https://example.org/data.csv",
        filename="data.csv",
        sha256=sha256,
        license_verbatim="CC BY-NC-SA 4.0",
        track="R",
    )


def test_sha256_file_known_content(tmp_path):
    p = tmp_path / "x.bin"
    p.write_bytes(b"hello sina")
    assert sha256_file(p) == hashlib.sha256(b"hello sina").hexdigest()


def test_registry_roundtrip_is_atomic(tmp_path):
    reg_path = tmp_path / "REGISTRY.json"
    reg = Registry()
    reg.upsert(_asset())
    save_registry(reg, reg_path)
    loaded = load_registry(reg_path)
    assert loaded.schema_version == 1
    assert loaded.get("primekg").license_verbatim == "CC BY-NC-SA 4.0"
    # no temp files left behind
    assert list(tmp_path.glob("*.tmp")) == []


def test_load_missing_registry_returns_empty(tmp_path):
    assert load_registry(tmp_path / "nope.json").assets == []


def test_upsert_replaces_same_name(tmp_path):
    reg = Registry()
    reg.upsert(_asset(sha256=None))
    reg.upsert(_asset(sha256="abc"))
    assert len(reg.assets) == 1
    assert reg.get("primekg").sha256 == "abc"


def test_get_missing_raises(tmp_path):
    with pytest.raises(KeyError):
        Registry().get("nope")


def test_verify_asset_cases(tmp_path):
    p = tmp_path / "data.csv"
    p.write_bytes(b"payload")
    digest = hashlib.sha256(b"payload").hexdigest()

    assert verify_asset(_asset(sha256=None), p) is True  # first download: no hash yet
    assert verify_asset(_asset(sha256=digest), p) is True
    assert verify_asset(_asset(sha256="0" * 64), p) is False
    assert verify_asset(_asset(sha256=digest), tmp_path / "missing.csv") is False


def test_mark_verified_fills_metadata(tmp_path):
    p = tmp_path / "data.csv"
    p.write_bytes(b"payload")
    asset = mark_verified(_asset(), p)
    assert asset.status is AssetStatus.VERIFIED
    assert asset.sha256 == hashlib.sha256(b"payload").hexdigest()
    assert asset.size_bytes == len(b"payload")
    assert asset.downloaded_at is not None
