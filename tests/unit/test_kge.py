"""Unit tests for sina.kge — tiny deterministic training + scoring."""

import polars as pl
import pytest

from sina.kge import train_and_score_kge


@pytest.fixture()
def toy_edges():
    return pl.DataFrame(
        {
            "src_id": ["a", "b", "c", "d", "a", "b"],
            "dst_id": ["b", "c", "d", "e", "c", "d"],
            "relation": ["r", "r", "r", "r", "s", "s"],
        }
    )


def test_kge_train_and_score(toy_edges):
    scores = train_and_score_kge(
        toy_edges,
        query_nodes=["a"],
        candidate_nodes=["b", "c", "d", "e"],
        relation="r",
        model="TransE",
        embedding_dim=8,
        epochs=5,
        seed=42,
        device="cpu",
    )
    assert scores.height == 4
    assert set(scores["query_id"].unique().to_list()) == {"a"}
    assert scores["score"].is_not_null().all()


def test_kge_unseen_query_returns_empty(toy_edges):
    scores = train_and_score_kge(
        toy_edges,
        query_nodes=["z"],
        candidate_nodes=["b", "c"],
        relation="r",
        model="TransE",
        embedding_dim=8,
        epochs=5,
        seed=42,
        device="cpu",
    )
    assert scores.height == 0
