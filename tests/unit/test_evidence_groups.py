"""Unit tests for sina.evidence_groups (spec §5.4 leakage-critical contract).

Synthetic fixtures only — no science data in tests (spec §0.2.2).
"""

import polars as pl
import pytest

from sina.evidence_groups import (
    add_evidence_group_ids,
    canonical_edge_key,
    evidence_group_id,
    reciprocal_audit,
    symmetric_relations,
)


@pytest.fixture()
def toy_edges() -> pl.DataFrame:
    # g1<->g2 are reciprocal twins (symmetric relation); g2->g3 is directed;
    # two duplicate rows of the same directed assertion (multi-source duplicate).
    return pl.DataFrame(
        {
            "src_id": ["ncbigene:1", "ncbigene:2", "ncbigene:2", "ncbigene:2"],
            "dst_id": ["ncbigene:2", "ncbigene:1", "ncbigene:3", "ncbigene:3"],
            "relation": ["ppi", "ppi", "inhibits", "inhibits"],
        }
    )


def test_canonical_key_order_insensitive_only_when_symmetric():
    assert canonical_edge_key("a", "ppi", "b", symmetric=True) == canonical_edge_key(
        "b", "ppi", "a", symmetric=True
    )
    assert canonical_edge_key("a", "inhibits", "b", symmetric=False) != canonical_edge_key(
        "b", "inhibits", "a", symmetric=False
    )


def test_group_id_deterministic():
    a = evidence_group_id("ncbigene:1", "ppi", "ncbigene:2", symmetric=True)
    b = evidence_group_id("ncbigene:1", "ppi", "ncbigene:2", symmetric=True)
    assert a == b and len(a) == 16


def test_reciprocal_audit_finds_symmetric_relation(toy_edges):
    audit = reciprocal_audit(toy_edges)
    frac = {r["relation"]: r["reverse_fraction"] for r in audit.to_dicts()}
    assert frac["ppi"] == 1.0
    assert frac["inhibits"] == 0.0
    assert symmetric_relations(audit, threshold=0.99) == {"ppi"}


def test_reciprocal_twins_share_group_id(toy_edges):
    """CI guard (§13.1): reciprocal twins MUST share one evidence_group_id."""
    out = add_evidence_group_ids(toy_edges, symmetric={"ppi"})
    twins = out.filter(pl.col("relation") == "ppi")["evidence_group_id"].unique()
    assert twins.len() == 1


def test_multi_source_duplicates_share_group_id(toy_edges):
    """CI guard (§13.1): duplicate assertions MUST share one evidence_group_id."""
    out = add_evidence_group_ids(toy_edges, symmetric={"ppi"})
    dups = out.filter(pl.col("relation") == "inhibits")["evidence_group_id"].unique()
    assert dups.len() == 1


def test_directed_reverse_gets_distinct_group():
    edges = pl.DataFrame(
        {
            "src_id": ["a", "b"],
            "dst_id": ["b", "a"],
            "relation": ["inhibits", "inhibits"],
        }
    )
    out = add_evidence_group_ids(edges, symmetric=set())
    assert out["evidence_group_id"].n_unique() == 2
