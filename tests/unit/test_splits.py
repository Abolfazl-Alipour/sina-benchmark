"""Unit tests for sina.splits — leakage-critical split primitives (spec §6.2).

Synthetic fixtures only — no science data in tests (spec §0.2.2).
"""

import polars as pl
import pytest

from sina.splits import (
    SplitLeakageError,
    assert_split_disjoint,
    family_holdout_split,
    random_edge_split,
    temporal_split,
    txgnn_area_split,
)


@pytest.fixture()
def toy_edges() -> pl.DataFrame:
    """Four groups (g1..g4); g1 has reciprocal twins, g4 has duplicates."""
    return pl.DataFrame(
        {
            "src_id": ["a", "b", "c", "c", "d", "d"],
            "dst_id": ["b", "a", "e", "f", "g", "g"],
            "relation": ["ppi", "ppi", "touches", "touches", "inhibits", "inhibits"],
            "evidence_group_id": ["g1", "g1", "g2", "g3", "g4", "g4"],
        }
    )


@pytest.fixture()
def family_edges() -> pl.DataFrame:
    """Disease family split fixture.

    family = {dis1, dis2}; eval relation = disease_protein.
    Reciprocal twin of a family eval edge shares g_eval so it must land in the
    same partition (leakage guard).
    """
    return pl.DataFrame(
        {
            "src_id": ["dis1", "prot:1", "dis2", "prot:2", "other", "dis1"],
            "dst_id": ["prot:1", "dis1", "prot:2", "dis2", "prot:3", "prot:3"],
            "relation": [
                "disease_protein",
                "disease_protein",
                "disease_protein",
                "disease_protein",
                "ppi",
                "ppi",
            ],
            "label": ["", "", "", "", "", ""],
            "evidence_group_id": ["g_eval", "g_eval", "g2", "g2", "g3", "g4"],
        }
    )


@pytest.fixture()
def temporal_edges() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "src_id": ["d1", "d2", "d3", "d4", "d5", "x1"],
            "dst_id": ["g1", "g2", "g3", "g4", "g5", "y1"],
            "relation": ["DIS-GEN"] * 5 + ["ppi"],
            "label": ["DIS-GEN"] * 5 + ["ppi"],
            "year_initial": [2015, 2017, 2018, 2019, 2017, None],
            "year_final": [2016, 2019, 2018, 2020, 2022, None],
            "evidence_group_id": ["t1", "t2", "t3", "t4", "t5", "nt"],
        }
    )


def test_assert_split_disjoint_ok(toy_edges):
    split = {
        "train": toy_edges.filter(pl.col("evidence_group_id").is_in(["g1", "g2"])),
        "valid": toy_edges.filter(pl.col("evidence_group_id") == "g3"),
        "test": toy_edges.filter(pl.col("evidence_group_id") == "g4"),
    }
    assert_split_disjoint(split)


def test_assert_split_disjoint_raises_on_leakage(toy_edges):
    split = {
        "train": toy_edges.filter(pl.col("evidence_group_id") == "g1"),
        "test": toy_edges.filter(pl.col("evidence_group_id").is_in(["g1", "g2"])),
    }
    with pytest.raises(SplitLeakageError, match="g1"):
        assert_split_disjoint(split)


def test_random_edge_split_disjoint_and_deterministic(toy_edges):
    s1 = random_edge_split(toy_edges, {"train": 0.7, "valid": 0.2, "test": 0.1}, seed=7)
    s2 = random_edge_split(toy_edges, {"train": 0.7, "valid": 0.2, "test": 0.1}, seed=7)
    assert_split_disjoint(s1)
    assert_split_disjoint(s2)
    for part in ("train", "valid", "test"):
        g1 = set(s1[part]["evidence_group_id"].unique().to_list())
        g2 = set(s2[part]["evidence_group_id"].unique().to_list())
        assert g1 == g2


def test_random_edge_split_preserves_all_groups(toy_edges):
    split = random_edge_split(toy_edges, {"train": 0.8, "valid": 0.1, "test": 0.1}, seed=9)
    seen = (
        set(split["train"]["evidence_group_id"].unique().to_list())
        | set(split["valid"]["evidence_group_id"].unique().to_list())
        | set(split["test"]["evidence_group_id"].unique().to_list())
    )
    assert seen == {"g1", "g2", "g3", "g4"}


def test_family_holdout_zero_shot_and_neighborhood_removed(family_edges):
    family = {"dis1", "dis2"}
    split = family_holdout_split(
        family_edges,
        family,
        eval_relation_col="relation",
        eval_relation_values={"disease_protein"},
        seed=11,
    )
    assert_split_disjoint(split)
    train_groups = set(split["train"]["evidence_group_id"].unique().to_list())
    assert "g_eval" not in train_groups
    assert "g2" not in train_groups  # family eval edge removed from train
    assert "g4" not in train_groups  # non-eval family edge removed (neighborhood)
    assert "g3" in train_groups
    # reciprocal twin of eval edge co-assigned
    eval_groups = set(split["valid"]["evidence_group_id"].unique().to_list()) | set(
        split["test"]["evidence_group_id"].unique().to_list()
    )
    assert "g_eval" in eval_groups


def test_family_holdout_reciprocal_twins_share_partition(family_edges):
    """CI guard: reciprocal twins sharing evidence_group_id must not straddle."""
    split = family_holdout_split(
        family_edges,
        {"dis1", "dis2"},
        eval_relation_col="relation",
        eval_relation_values={"disease_protein"},
        seed=13,
    )
    for part in ("train", "valid", "test"):
        gids = split[part]["evidence_group_id"].unique().to_list()
        assert gids.count("g_eval") <= 1


def test_txgnn_area_split_protocol():
    edges = pl.DataFrame(
        {
            "src_id": ["d1", "p1", "d1", "d1", "x1"],
            "dst_id": ["p1", "d1", "d2", "t1", "y1"],
            "relation": ["indication", "indication", "ppi", "treatment", "ppi"],
            "evidence_group_id": ["i1", "i1", "n1", "n2", "o1"],
        }
    )
    split = txgnn_area_split(
        edges,
        area_disease_ids={"d1"},
        indication_relations={"indication", "contraindication"},
        seed=17,
    )
    assert_split_disjoint(split)
    assert set(split["test"]["evidence_group_id"].unique().to_list()) == {"i1"}
    assert "o1" in set(split["train"]["evidence_group_id"].unique().to_list())


def test_temporal_split_cutoff_and_straddlers_dropped(temporal_edges):
    split = temporal_split(temporal_edges, cutoff_year=2018)
    assert_split_disjoint(split)
    # t1 (final 2016) and t3 (final 2018) train; t4 (init 2019) test
    train_groups = set(split["train"]["evidence_group_id"].unique().to_list())
    test_groups = set(split["test"]["evidence_group_id"].unique().to_list())
    assert {"t1", "t3", "nt"}.issubset(train_groups)
    assert {"t4"}.issubset(test_groups)
    # t2 straddles 2018 (init 2017, final 2019) → dropped entirely
    assert "t2" not in train_groups and "t2" not in test_groups
    # t5 has init 2017/final 2022; straddles, dropped
    assert "t5" not in train_groups and "t5" not in test_groups
