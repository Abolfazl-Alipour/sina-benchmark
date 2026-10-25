"""Unit tests for Track-D evidence-group semantics (spec §5.2.3).

Synthetic fixtures only.
"""

import polars as pl

from sina.evidence_groups import add_evidence_group_ids


def _edges() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "src_id": ["mondo:1", "ensembl:A", "mondo:2", "mondo:2"],
            "dst_id": ["ensembl:A", "mondo:1", "ensembl:B", "ensembl:B"],
            "relation": ["ASSOCIATED_WITH"] * 4,
            "undirected": [True, True, True, True],
            "sources_direct": [["d1"], ["d1"], ["d2"], ["d1"]],
            "sources_indirect": [[], [], [], []],
        }
    )


def _with_srcset(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        _s=pl.struct(direct=pl.col("sources_direct"), indirect=pl.col("sources_indirect"))
    )


def _group_ids(df: pl.DataFrame) -> pl.Series:
    out = add_evidence_group_ids(_with_srcset(df), undirected_col="undirected", sources_col="_s")
    return out["evidence_group_id"]


def test_undirected_twins_share_id_when_source_sets_match():
    twins = _group_ids(_edges()).head(2)
    assert twins.n_unique() == 1


def test_same_triple_different_source_set_gets_distinct_id():
    dup = _group_ids(_edges()).tail(2)
    # same (from, relation, to) but different source-set → different group (§5.4 granularity)
    assert dup.n_unique() == 2


def test_source_order_does_not_change_key():
    df = pl.DataFrame(
        {
            "src_id": ["a", "a"],
            "dst_id": ["b", "b"],
            "relation": ["r", "r"],
            "undirected": [False, False],
            "sources_direct": [["x", "y"], ["y", "x"]],
            "sources_indirect": [["z"], ["z"]],
        }
    )
    out = add_evidence_group_ids(_with_srcset(df), undirected_col="undirected", sources_col="_s")
    assert out["evidence_group_id"].n_unique() == 1


def test_directed_row_reverse_not_grouped():
    df = pl.DataFrame(
        {
            "src_id": ["a", "b"],
            "dst_id": ["b", "a"],
            "relation": ["r", "r"],
            "undirected": [False, False],
        }
    )
    out = add_evidence_group_ids(df, undirected_col="undirected")
    assert out["evidence_group_id"].n_unique() == 2
