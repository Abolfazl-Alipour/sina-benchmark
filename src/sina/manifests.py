"""Run/stage manifests with peak-RSS tracking.

Plan §0.3.3: peak RSS < 40 GB asserted. Plan §0.3.5: checkpoints ≤30 min,
--resume everywhere. Every long script writes a manifest JSON recording its
stages, output hashes, and peak memory so runs are auditable and resumable.
"""

import threading
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path

import psutil
from pydantic import BaseModel, Field

from sina.logging import get_logger
from sina.registry import sha256_file

log = get_logger(__name__)

MAX_RSS_GB = 40.0  # plan §0.3.3


def _utcnow() -> datetime:
    return datetime.now(UTC)


class RssGuardError(RuntimeError):
    """Raised when a run exceeds the plan's peak-RSS budget."""


class RssSampler:
    """Background-thread peak-RSS sampler for the current process (incl. children)."""

    def __init__(self, interval_s: float = 5.0, limit_gb: float = MAX_RSS_GB) -> None:
        self.interval_s = interval_s
        self.limit_gb = limit_gb
        self.peak_gb = 0.0
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def _sample_once(self) -> float:
        proc = psutil.Process()
        total = proc.memory_info().rss
        for child in proc.children(recursive=True):
            try:
                total += child.memory_info().rss
            except psutil.NoSuchProcess:
                continue
        return total / 1e9

    def _loop(self) -> None:
        while not self._stop.is_set():
            self.peak_gb = max(self.peak_gb, self._sample_once())
            self._stop.wait(self.interval_s)

    def __enter__(self) -> "RssSampler":
        self._thread = threading.Thread(target=self._loop, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *exc: object) -> None:
        self._stop.set()
        if self._thread:
            self._thread.join(timeout=2)
        self.peak_gb = max(self.peak_gb, self._sample_once())

    def check(self) -> None:
        """Assert the plan's RSS budget (call at end of run)."""
        if self.peak_gb > self.limit_gb:
            msg = (
                f"peak RSS {self.peak_gb:.1f} GB exceeded "
                f"plan budget {self.limit_gb:.1f} GB (§0.3.3)"
            )
            raise RssGuardError(msg)


class StageRecord(BaseModel):
    name: str
    started_at: datetime = Field(default_factory=_utcnow)
    finished_at: datetime | None = None
    status: str = "running"  # running | done | failed
    outputs: dict[str, str] = Field(default_factory=dict)  # path -> sha256
    detail: str = ""


class RunManifest(BaseModel):
    """One script execution: identity, config echo, stages, peak RSS."""

    run_id: str
    script: str
    config: dict = Field(default_factory=dict)
    git_ref: str = ""
    started_at: datetime = Field(default_factory=_utcnow)
    finished_at: datetime | None = None
    stages: list[StageRecord] = Field(default_factory=list)
    peak_rss_gb: float = 0.0
    status: str = "running"  # running | done | failed

    def stage_names_done(self) -> set[str]:
        return {s.name for s in self.stages if s.status == "done"}


class Run:
    """Context manager tying an RssSampler to a RunManifest with stage checkpoints."""

    def __init__(
        self,
        script: str,
        manifest_path: Path | str,
        config: dict | None = None,
        git_ref: str = "",
        rss_limit_gb: float = MAX_RSS_GB,
    ) -> None:
        self.manifest_path = Path(manifest_path)
        self.sampler = RssSampler(limit_gb=rss_limit_gb)
        if self.manifest_path.exists():
            # --resume: reload prior manifest, keep completed stages.
            self.manifest = RunManifest.model_validate_json(
                self.manifest_path.read_text(encoding="utf-8")
            )
            log.info("run_resumed", script=script, run_id=self.manifest.run_id)
        else:
            self.manifest = RunManifest(
                run_id=_utcnow().strftime("%Y%m%dT%H%M%SZ"),
                script=script,
                config=config or {},
                git_ref=git_ref,
            )

    @property
    def done_stages(self) -> set[str]:
        return self.manifest.stage_names_done()

    @contextmanager
    def stage(self, name: str, outputs: list[Path] | None = None):
        """Checkpoint one stage; skipped by callers when name in run.done_stages."""
        record = StageRecord(name=name)
        self.manifest.stages.append(record)
        self._save()
        try:
            yield record
        except Exception:
            record.status = "failed"
            self.manifest.status = "failed"
            self._save()
            raise
        record.status = "done"
        record.finished_at = _utcnow()
        for out in outputs or []:
            out = Path(out)
            if out.exists() and out.is_file():
                record.outputs[str(out)] = sha256_file(out)
        self._save()  # checkpoint after every stage
        log.info("stage_done", stage=name)

    def _save(self) -> None:
        self.manifest_path.parent.mkdir(parents=True, exist_ok=True)
        self.manifest_path.write_text(self.manifest.model_dump_json(indent=2), encoding="utf-8")

    def __enter__(self) -> "Run":
        self.sampler.__enter__()
        self._save()
        return self

    def __exit__(self, exc_type, *exc) -> None:
        self.sampler.__exit__()
        self.manifest.peak_rss_gb = round(self.sampler.peak_gb, 3)
        self.manifest.finished_at = _utcnow()
        self.manifest.status = "failed" if exc_type else "done"
        self._save()
        if exc_type is None:
            self.sampler.check()  # RSS budget assertion (§0.3.3)
