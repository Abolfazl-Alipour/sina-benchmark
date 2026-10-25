"""Resumable, hash-verified downloads (plan §0.3.5/§0.3.6).

Every science-path download goes through download_asset():
- skips when the destination already exists and matches the registry hash;
- resumes a partial .part file via HTTP Range when the server supports it;
- verifies SHA-256 before promoting .part to the final name;
- updates REGISTRY.json (hash, size, timestamp) after verification.

The HTTP fetch is injectable so unit tests never touch the network
(no dummy data in the science path — fakes live in tests only).
"""

from collections.abc import Callable
from pathlib import Path

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential

from sina.logging import get_logger
from sina.registry import RegistryAsset, mark_verified, sha256_file, verify_asset

log = get_logger(__name__)

# fetcher(url, offset_bytes, part_file) appends the response body to part_file
Fetcher = Callable[[str, int, Path], None]


@retry(stop=stop_after_attempt(4), wait=wait_exponential(multiplier=2, min=2, max=30))
def _httpx_fetch(url: str, offset: int, part_file: Path) -> None:
    headers = {"Range": f"bytes={offset}-"} if offset else {}
    with (
        httpx.stream("GET", url, headers=headers, follow_redirects=True, timeout=120.0) as resp,
        part_file.open("ab") as fh,
    ):
        resp.raise_for_status()
        # raise_for_status lets 1xx/2xx through: a 202 bot-check page is NOT a
        # payload (seen 2026-08-15: Dataverse WAF challenge from an EC2 IP).
        if resp.status_code not in (200, 206):
            raise ValueError(
                f"unexpected HTTP {resp.status_code} for {url} — refusing non-payload response"
            )
        # 200 on a resume request means the server ignored Range — restart cleanly.
        if offset and resp.status_code == 200:
            log.warning("range_not_supported_restarting", url=url)
            fh.truncate(0)
        for chunk in resp.iter_bytes(chunk_size=8 * 1024 * 1024):
            fh.write(chunk)


def download_asset(
    asset: RegistryAsset,
    dest_dir: Path | str,
    fetcher: Fetcher = _httpx_fetch,
    dry_run: bool = False,
) -> Path:
    """Download + verify one registry asset. Idempotent and resumable."""
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / asset.filename
    part = dest.with_name(dest.name + ".part")

    if dest.exists():
        if verify_asset(asset, dest):
            log.info("asset_present_verified", name=asset.name, path=str(dest))
            return dest
        raise ValueError(
            f"{dest} exists but does not match registry SHA-256 for {asset.name!r}; "
            "refusing to overwrite — investigate or delete manually"
        )

    if dry_run:
        log.info("dry_run_would_download", name=asset.name, url=asset.url, dest=str(dest))
        return dest

    offset = part.stat().st_size if part.exists() else 0
    if offset:
        log.info("resuming_partial_download", name=asset.name, offset=offset)
    log.info("download_start", name=asset.name, url=asset.url)
    fetcher(asset.url, offset, part)

    if not part.exists() or part.stat().st_size == 0:
        part.unlink(missing_ok=True)
        raise ValueError(
            f"zero-byte download for {asset.name!r} — refusing to verify an empty payload"
        )
    digest = sha256_file(part)
    if asset.sha256 is not None and digest != asset.sha256:
        part.unlink()
        raise ValueError(
            f"SHA-256 mismatch for {asset.name!r}: got {digest}, expected {asset.sha256} "
            "(partial file removed; safe to retry)"
        )
    part.rename(dest)
    mark_verified(asset, dest)
    log.info("download_verified", name=asset.name, sha256=asset.sha256, size=asset.size_bytes)
    return dest
