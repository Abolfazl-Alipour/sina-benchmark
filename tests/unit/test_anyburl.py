"""Unit tests for sina.anyburl — triple export, config generation, parsing."""

import polars as pl

from sina.anyburl import (
    export_triples,
    filter_rules,
    parse_predictions,
    write_apply_config,
    write_learn_config,
)


def test_export_triples(tmp_path):
    edges = pl.DataFrame(
        {
            "src_id": ["a", "b"],
            "relation": ["r1", "off-label use"],
            "dst_id": ["c", "d"],
        }
    )
    path = tmp_path / "triples.txt"
    export_triples(edges, path)
    assert path.read_text() == "a\tr1\tc\nb\toff-label_use\td\n"


def test_write_learn_config(tmp_path):
    path = tmp_path / "learn.properties"
    write_learn_config(
        path,
        tmp_path / "train.txt",
        tmp_path / "out",
        worker_threads=2,
        snapshots_at="10,50",
    )
    text = path.read_text()
    assert "PATH_TRAINING" in text
    assert "WORKER_THREADS = 2" in text
    assert "SNAPSHOTS_AT = 10,50" in text


def test_write_apply_config(tmp_path):
    path = tmp_path / "apply.properties"
    write_apply_config(
        path,
        tmp_path / "rules-100",
        tmp_path / "train.txt",
        tmp_path / "test.txt",
        tmp_path / "preds-100",
        valid_path=tmp_path / "valid.txt",
        worker_threads=2,
        top_k=5,
    )
    text = path.read_text()
    assert "PATH_RULES" in text
    assert "PATH_TEST" in text
    assert "PATH_VALID" in text
    assert "TOP_K_OUTPUT = 5" in text


def test_filter_rules_drops_empty_body(tmp_path):
    src = tmp_path / "rules"
    dst = tmp_path / "rules-filtered"
    src.write_text(
        "13\t13\t1.0\tr(X,Y) <= r(Y,X)\n"
        "267\t840\t0.0003\tdrug_drug(X,Y) <= \n"
        "5\t5\t1.0\tr2(X,Y) <= r3(X,Y)\n"
    )
    kept = filter_rules(src, dst)
    assert kept == 2
    assert "drug_drug" not in dst.read_text()


def test_parse_predictions_tails(tmp_path):
    pred = tmp_path / "preds-100"
    pred.write_text(
        "a r1 c\nHeads: \nTails: d\t0.9\te\t0.7\nb r2 d\nHeads: \nTails: c\t0.5\n\nbadline\n"
    )
    scores = parse_predictions(pred, prediction_type="tails")
    assert scores.height == 3
    assert set(scores["query_id"].unique().to_list()) == {"a", "b"}
    assert set(scores["candidate_id"].unique().to_list()) == {"d", "e", "c"}
    assert scores["score"].max() == 0.9


def test_parse_predictions_heads(tmp_path):
    pred = tmp_path / "preds-100"
    pred.write_text("a r1 c\nHeads: d\t0.9\te\t0.7\nTails: \nb r2 d\nHeads: c\t0.5\nTails: \n")
    scores = parse_predictions(pred, prediction_type="heads")
    assert scores.height == 3
    assert set(scores["query_id"].unique().to_list()) == {"c", "d"}
    assert set(scores["candidate_id"].unique().to_list()) == {"d", "e", "c"}
