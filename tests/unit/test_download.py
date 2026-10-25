"""Unit tests for sina.download — resume, skip, hash enforcement.

Network is never touched: a fake fetcher simulates Range-resumable transfers.
"""

import hashlib
from pathlib import Path

import pytest

from sina.download import download_asset
from sina.registry import AssetStatus, RegistryAsset

PAYLOAD = b"sina-payload-" * 4096  # ~53 KB
DIGEST = hashlib.sha256(PAYLOAD).hexdigest()


def _asset(sha256: str | None = None) -> RegistryAsset:
    return RegistryAsset(
        name="t",
        version="v1",
        url="https://example.org/blob",
        filename="blob.bin",
        sha256=sha256,
    )


def make_fetcher(payload: bytes):
    """Fake Range-aware fetcher: appends payload[offset:] to the part file."""

    def fetch(url: str, offset: int, part_file: Path) -> None:
        assert url == "https://example.org/blob"
        with part_file.open("ab") as fh:
            fh.write(payload[offset:])

    return fetch


def test_fresh_download_records_hash(tmp_path):
    asset = _asset()
    dest = download_asset(asset, tmp_path, fetcher=make_fetcher(PAYLOAD))
    assert dest.read_bytes() == PAYLOAD
    assert asset.status is AssetStatus.VERIFIED
    assert asset.sha256 == DIGEST
    assert not (tmp_path / "blob.bin.part").exists()


def test_resume_partial_download(tmp_path):
    part = tmp_path / "blob.bin.part"
    cut = len(PAYLOAD) // 3
    part.write_bytes(PAYLOAD[:cut])  # interrupted earlier run
    dest = download_asset(_asset(sha256=DIGEST), tmp_path, fetcher=make_fetcher(PAYLOAD))
    assert dest.read_bytes() == PAYLOAD


def test_hash_mismatch_aborts_and_cleans_part(tmp_path):
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        download_asset(_asset(sha256="0" * 64), tmp_path, fetcher=make_fetcher(PAYLOAD))
    assert not (tmp_path / "blob.bin").exists()
    assert not (tmp_path / "blob.bin.part").exists()


def test_existing_verified_file_skips_fetch(tmp_path):
    (tmp_path / "blob.bin").write_bytes(PAYLOAD)

    def fail_fetch(url, offset, part_file):  # pragma: no cover - must not be called
        raise AssertionError("fetcher called despite verified file")

    dest = download_asset(_asset(sha256=DIGEST), tmp_path, fetcher=fail_fetch)
    assert dest.read_bytes() == PAYLOAD


def test_existing_mismatched_file_refuses_overwrite(tmp_path):
    (tmp_path / "blob.bin").write_bytes(b"corrupted")
    with pytest.raises(ValueError, match="refusing to overwrite"):
        download_asset(_asset(sha256=DIGEST), tmp_path, fetcher=make_fetcher(PAYLOAD))
    assert (tmp_path / "blob.bin").read_bytes() == b"corrupted"  # untouched


def test_dry_run_downloads_nothing(tmp_path):
    dest = download_asset(_asset(), tmp_path, fetcher=make_fetcher(PAYLOAD), dry_run=True)
    assert not dest.exists()
    assert not (tmp_path / "blob.bin.part").exists()


def test_zero_byte_download_aborts(tmp_path):
    """A 0-byte 'success' (e.g. bot-check empty body) must never verify."""

    def empty_fetch(url, offset, part_file):
        part_file.touch()

    with pytest.raises(ValueError, match="zero-byte download"):
        download_asset(_asset(), tmp_path, fetcher=empty_fetch)
    assert not (tmp_path / "blob.bin").exists()
    assert not (tmp_path / "blob.bin.part").exists()


def test_non_payload_status_aborts(tmp_path, monkeypatch):
    """HTTP 202 (bot-check challenge page) is not a payload — must raise."""
    import contextlib

    import httpx

    from sina import download as download_mod

    class FakeResponse:
        status_code = 202
        content = b""

        def raise_for_status(self):  # 202 does not raise — that is the trap
            return None

        def iter_bytes(self, chunk_size):
            return iter([b"<html>challenge</html>"])

    @contextlib.contextmanager
    def fake_stream(method, url, headers=None, follow_redirects=False, timeout=None):
        yield FakeResponse()

    monkeypatch.setattr(httpx, "stream", fake_stream)
    from tenacity import wait_none

    part = tmp_path / "blob.bin.part"
    with pytest.raises(ValueError, match="unexpected HTTP 202"):
        download_mod._httpx_fetch.retry_with(wait=wait_none(), reraise=True)(
            "https://example.org/blob", 0, part
        )


def test_apply_pin_preserves_registry_hash(tmp_path):
    from sina.registry import Registry, apply_pin

    registry = Registry(assets=[_asset(sha256=DIGEST)])
    pinned = apply_pin(_asset(), registry)  # config asset without hash
    assert pinned.sha256 == DIGEST


def test_apply_pin_not_carried_across_versions(tmp_path):
    from sina.registry import Registry, apply_pin

    registry = Registry(assets=[_asset(sha256=DIGEST)])
    newer = _asset()
    newer.version = "v2"
    assert apply_pin(newer, registry).sha256 is None


def test_apply_pin_keeps_config_declared_hash(tmp_path):
    from sina.registry import Registry, apply_pin

    registry = Registry(assets=[_asset(sha256=DIGEST)])
    own = _asset(sha256="1" * 64)
    assert apply_pin(own, registry).sha256 == "1" * 64
