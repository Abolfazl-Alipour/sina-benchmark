"""Integration tests over real M2A-3 baseline outputs (spec §6.3).

These run locally where `data/processed/m2/baselines/` exists; they SKIP in CI.
"""

from pathlib import Path

import polars as pl
import pytest

BASELINES_DIR = Path("data/processed/m2/baselines")

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(not BASELINES_DIR.exists(), reason="M2 baselines not built"),
]


def test_baseline_outputs_exist():
    splits = [d.name for d in BASELINES_DIR.iterdir() if d.is_dir()]
    assert splits
    for split in splits:
        base = BASELINES_DIR / split
        assert (base / "metadata.json").exists()
        for name in (
            "degree_scores",
            "degree_matched_nulls",
            "relation_shuffled_nulls",
        ):
            assert (base / f"{name}.parquet").exists()


def test_degree_scores_have_reasonable_shape():
    for split_dir in BASELINES_DIR.iterdir():
        if not split_dir.is_dir():
            continue
        meta = __import__("json").loads((split_dir / "metadata.json").read_text(encoding="utf-8"))
        deg = pl.read_parquet(split_dir / "degree_scores.parquet")
        # one row per (query, candidate)
        assert deg.height == meta["n_queries_in_test"] * meta["n_candidates"]
        assert deg["score"].is_not_null().all()


def test_nulls_are_non_trivial():
    for split_dir in BASELINES_DIR.iterdir():
        if not split_dir.is_dir():
            continue
        nulls = pl.read_parquet(split_dir / "degree_matched_nulls.parquet")
        assert nulls.height > 0
        assert (nulls["null_dst_id"] != nulls["true_dst_id"]).all()
