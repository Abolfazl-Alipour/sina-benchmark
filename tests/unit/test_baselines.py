"""Unit tests for sina.baselines — deterministic, leakage-free scoring."""

import polars as pl
import pytest

from sina.baselines import (
    degree_baseline,
    degree_matched_nulls,
    ppr_baseline,
    relation_shuffled_nulls,
)


@pytest.fixture()
def toy_graph() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "src_id": ["a", "b", "c", "d", "e"],
            "dst_id": ["b", "c", "d", "e", "a"],
            "relation": ["r1", "r1", "r2", "r2", "r1"],
        }
    )


def test_degree_baseline_counts_undirected(toy_graph):
    deg = degree_baseline(toy_graph, {"a", "b", "c", "d", "e"})
    by_node = {r["node_id"]: r["degree"] for r in deg.to_dicts()}
    assert by_node["a"] == 2
    assert by_node["b"] == 2
    assert by_node["c"] == 2
    assert by_node["d"] == 2
    assert by_node["e"] == 2


def test_degree_baseline_candidates_only(toy_graph):
    deg = degree_baseline(toy_graph, {"a", "c"})
    assert set(deg["node_id"].to_list()) == {"a", "c"}


def test_ppr_baseline_returns_query_candidate_pairs(toy_graph):
    scores = ppr_baseline(
        toy_graph,
        query_nodes=["a"],
        candidate_nodes=["b", "c", "d", "e"],
        alpha=0.15,
        max_iter=50,
    )
    assert scores.height == 4
    assert set(scores["query_id"].unique().to_list()) == {"a"}
    # PPR scores are non-negative and sum roughly to restart weight
    assert scores["ppr_score"].sum() <= 1.0 + 1e-6
    assert scores["ppr_score"].min() >= 0.0


def test_ppr_baseline_query_is_excluded(toy_graph):
    scores = ppr_baseline(
        toy_graph,
        query_nodes=["a"],
        candidate_nodes=["a", "b"],
        alpha=0.15,
        max_iter=50,
    )
    assert "a" not in set(scores["candidate_id"].to_list())


def test_degree_matched_nulls_preserves_query(toy_graph):
    true_edges = pl.DataFrame({"src_id": ["a"], "dst_id": ["b"]})
    nulls = degree_matched_nulls(
        toy_graph,
        true_edges,
        candidate_pool=["a", "b", "c", "d", "e"],
        n_negatives=2,
        bins=2,
        seed=7,
    )
    assert nulls.height == 2
    assert set(nulls["src_id"].to_list()) == {"a"}
    assert "b" not in set(nulls["null_dst_id"].to_list())


def test_degree_matched_nulls_deterministic(toy_graph):
    true_edges = pl.DataFrame({"src_id": ["a"], "dst_id": ["b"]})
    n1 = degree_matched_nulls(toy_graph, true_edges, ["a", "b", "c", "d", "e"], 2, 2, seed=7)
    n2 = degree_matched_nulls(toy_graph, true_edges, ["a", "b", "c", "d", "e"], 2, 2, seed=7)
    assert n1.equals(n2)


def test_relation_shuffled_nulls_is_permutation(toy_graph):
    shuffled = relation_shuffled_nulls(toy_graph, seed=11)
    assert shuffled.height == toy_graph.height
    assert (
        shuffled["relation"]
        .value_counts()
        .sort("relation")
        .equals(toy_graph["relation"].value_counts().sort("relation"))
    )
