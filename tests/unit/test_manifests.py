"""Unit tests for sina.manifests — run/stage checkpoints, resume, RSS guard."""

import json

import pytest

from sina.manifests import RssGuardError, RssSampler, Run, RunManifest


def test_run_records_stages_and_outputs(tmp_path):
    out = tmp_path / "result.txt"
    mp = tmp_path / "manifest.json"
    with Run("script_x", mp, config={"a": 1}) as run, run.stage("build", outputs=[out]):
        out.write_text("done")
    data = json.loads(mp.read_text())
    assert data["status"] == "done"
    assert data["stages"][0]["name"] == "build"
    assert data["stages"][0]["status"] == "done"
    assert str(out) in data["stages"][0]["outputs"]
    assert data["peak_rss_gb"] > 0


def test_run_resume_keeps_completed_stages(tmp_path):
    mp = tmp_path / "manifest.json"
    with Run("script_x", mp) as run, run.stage("s1"):
        pass
    # simulate crash: new Run over same manifest
    with Run("script_x", mp) as run:
        assert "s1" in run.done_stages
        with run.stage("s2"):
            pass
    manifest = RunManifest.model_validate_json(mp.read_text())
    assert manifest.stage_names_done() == {"s1", "s2"}


def test_failed_stage_marked_and_reraised(tmp_path):
    mp = tmp_path / "manifest.json"
    ctx = pytest.raises(RuntimeError, match="boom")
    with ctx, Run("script_x", mp) as run, run.stage("fragile"):
        raise RuntimeError("boom")
    data = json.loads(mp.read_text())
    assert data["status"] == "failed"
    assert data["stages"][0]["status"] == "failed"


def test_rss_sampler_tracks_peak():
    with RssSampler(interval_s=0.01) as sampler:
        _ = [0] * 1000
    assert sampler.peak_gb > 0
    sampler.check()  # under the 40 GB plan budget — must not raise


def test_rss_guard_raises_over_budget():
    sampler = RssSampler(limit_gb=1e-9)
    sampler.peak_gb = 1.0
    with pytest.raises(RssGuardError):
        sampler.check()
