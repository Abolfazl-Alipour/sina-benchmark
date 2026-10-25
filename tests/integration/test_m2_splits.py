"""Integration tests over real M2A-2 split artifacts (spec §6.2).

These run locally where `data/processed/m2/splits/` exists; they SKIP in CI.
Every split is checked for the §5.4 leakage contract and basic shape.
"""

from pathlib import Path

import polars as pl
import pytest

from sina.splits import SplitLeakageError, assert_split_disjoint

SPLITS_DIR = Path("data/processed/m2/splits")

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(not SPLITS_DIR.exists(), reason="M2 splits not built"),
]

SPLITS = [
    "trackr_random",
    "trackd_random",
    "trackr_family",
    "trackd_family",
    "trackr_family_alzheimer",
    "trackr_family_parkinson",
    "trackd_family_alzheimer",
    "trackd_family_parkinson",
    "trackr_txgnn_neurodegenerative",
    "trackd_temporal",
]


def _load_split(name: str) -> dict[str, pl.DataFrame]:
    base = SPLITS_DIR / name
    return {part: pl.read_parquet(base / f"{part}.parquet") for part in ("train", "valid", "test")}


@pytest.mark.parametrize("name", SPLITS)
def test_split_exists(name):
    for part in ("train", "valid", "test"):
        assert (SPLITS_DIR / name / f"{part}.parquet").exists()
    assert (SPLITS_DIR / name / "metadata.json").exists()


@pytest.mark.parametrize("name", SPLITS)
def test_split_leakage_guard(name):
    split = _load_split(name)
    assert_split_disjoint(split)


@pytest.mark.parametrize("name", SPLITS)
def test_split_partitions_cover_groups(name):
    split = _load_split(name)
    groups = set()
    for part in ("train", "valid", "test"):
        groups.update(split[part]["evidence_group_id"].unique().to_list())
    meta = __import__("json").loads(
        (SPLITS_DIR / name / "metadata.json").read_text(encoding="utf-8")
    )
    meta_groups = sum(meta["partitions"][p]["groups"] for p in ("train", "valid", "test"))
    assert len(groups) == meta_groups


def test_random_splits_are_non_trivial():
    for name in ("trackr_random", "trackd_random"):
        split = _load_split(name)
        assert split["test"].height > 0
        assert split["valid"].height > 0
        assert split["train"].height > split["test"].height


def test_family_splits_hide_family_neighborhood():
    # Track-D family eval edges are DIS-GEN; Track-R are disease_protein.
    # We just assert eval sets are non-empty and train is the bulk.
    for name in ("trackr_family", "trackd_family"):
        split = _load_split(name)
        assert split["test"].height > 0
        assert split["train"].height > split["test"].height


def test_txgnn_split_uses_indications():
    split = _load_split("trackr_txgnn_neurodegenerative")
    assert split["test"].height > 0
    rels = set(split["test"]["relation"].unique().to_list())
    assert rels <= {"indication", "contraindication"}


def test_temporal_split_is_secondary_but_nonempty():
    split = _load_split("trackd_temporal")
    assert split["test"].height > 0
    # test edges all have year_initial > cutoff (2018)
    test_years = split["test"]["year_initial"].cast(pl.Int32, strict=False)
    assert test_years.min() > 2018


def test_disjointness_raises_on_corrupted_split():
    split = _load_split("trackr_random")
    bad = {
        "train": split["train"],
        "test": pl.concat([split["test"], split["train"].head(1)]),
    }
    with pytest.raises(SplitLeakageError):
        assert_split_disjoint(bad)
