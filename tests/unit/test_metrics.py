"""Unit tests for sina.metrics — deterministic ranking metrics and statistics."""

import polars as pl
import pytest

from sina.metrics import (
    degree_stratified_metrics,
    link_prediction_metrics,
    mean_sd_ci,
    paired_wilcoxon,
)


@pytest.fixture()
def toy_scores():
    return pl.DataFrame(
        {
            "query_id": ["q1", "q1", "q1", "q2", "q2", "q2"],
            "candidate_id": ["a", "b", "c", "a", "b", "c"],
            "score": [0.9, 0.8, 0.7, 0.1, 0.3, 0.2],
        }
    )


@pytest.fixture()
def toy_test():
    return pl.DataFrame({"src_id": ["q1", "q2"], "dst_id": ["a", "b"]})


def test_link_prediction_metrics_raw(toy_scores, toy_test):
    m = link_prediction_metrics(toy_scores, toy_test, ks=(1, 3))
    assert m["n_queries"] == 2
    # q1: a is rank 1; q2: b is rank 1
    assert m["raw"]["hits@1"] == 1.0
    assert m["raw"]["mrr"] == 1.0


def test_link_prediction_metrics_filtered(toy_scores, toy_test):
    # add a filter edge that would outrank the true answer for q1
    filter_edges = pl.DataFrame({"src_id": ["q1"], "dst_id": ["c"]})
    scores = pl.concat(
        [
            toy_scores,
            pl.DataFrame({"query_id": ["q1"], "candidate_id": ["c"], "score": [0.95]}),
        ]
    )
    m = link_prediction_metrics(scores, toy_test, filter_edges, ks=(1, 3))
    # raw: a is rank 2 (c outranks); filtered: c removed, a rank 1
    assert m["raw"]["hits@1"] == 0.5
    assert m["filtered"]["hits@1"] == 1.0


def test_degree_stratified_metrics(toy_scores, toy_test):
    train = pl.DataFrame({"src_id": ["a", "b"], "dst_id": ["b", "c"]})
    m = degree_stratified_metrics(toy_scores, toy_test, train, ks=(1,))
    assert m  # at least one bin


def test_mean_sd_ci():
    mean, sd, (lo, hi) = mean_sd_ci([1.0, 2.0, 3.0, 4.0, 5.0])
    assert mean == 3.0
    assert sd > 0
    assert lo < mean < hi


def test_paired_wilcoxon():
    x = [1, 2, 3, 4, 5]
    y = [2, 3, 4, 5, 6]
    _, p = paired_wilcoxon(x, y)
    assert 0 <= p <= 1
