"""Unit tests for sina.stats — multi-seed aggregation and comparisons."""

import polars as pl

from sina.stats import aggregate_seeds, compare_models, per_disease_mrr


def test_per_disease_mrr():
    scores = pl.DataFrame(
        {
            "query_id": ["q1", "q1", "q2", "q2"],
            "candidate_id": ["a", "b", "a", "b"],
            "score": [0.9, 0.1, 0.2, 0.8],
        }
    )
    test = pl.DataFrame({"src_id": ["q1", "q2"], "dst_id": ["a", "b"]})
    mrr = per_disease_mrr(scores, test)
    assert mrr["q1"] == 1.0
    assert mrr["q2"] == 1.0


def test_aggregate_seeds():
    seed_metrics = [
        {"mrr": 0.1, "hits@1": 0.05, "n_queries": 10},
        {"mrr": 0.2, "hits@1": 0.15, "n_queries": 10},
        {"mrr": 0.3, "hits@1": 0.25, "n_queries": 10},
    ]
    agg = aggregate_seeds(seed_metrics)
    assert agg["n_seeds"] == 3
    assert abs(agg["mrr"]["mean"] - 0.2) < 1e-9
    assert agg["mrr"]["sd"] > 0


def test_compare_models():
    scores_a = [
        pl.DataFrame(
            {
                "query_id": ["q1", "q1", "q2", "q2"],
                "candidate_id": ["a", "b", "a", "b"],
                "score": [0.9, 0.1, 0.2, 0.8],
            }
        )
    ]
    scores_b = [
        pl.DataFrame(
            {
                "query_id": ["q1", "q1", "q2", "q2"],
                "candidate_id": ["a", "b", "a", "b"],
                "score": [0.1, 0.9, 0.8, 0.2],
            }
        )
    ]
    test = pl.DataFrame({"src_id": ["q1", "q2"], "dst_id": ["a", "b"]})
    result = compare_models(scores_a, scores_b, test)
    assert result["n_seeds"] == 1
    assert result["median_p"] is not None
