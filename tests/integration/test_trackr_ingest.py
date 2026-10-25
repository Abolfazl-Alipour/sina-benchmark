"""Integration tests over the real Track-R canonical tables (spec §5.5).

These run locally where `data/processed/m1/trackr/` exists; they SKIP in CI
(no data artifacts in git). Subsampling is deterministic: rows whose
evidence_group_id's first byte is 0x00..0x02 (≈1% of groups) — group-level, so
the reciprocal co-assignment guard holds inside the subsample too.
"""

from pathlib import Path

import polars as pl
import pytest

TRACKR_DIR = Path("data/processed/m1/trackr")
NODES = TRACKR_DIR / "nodes.parquet"
EDGES = TRACKR_DIR / "edges.parquet"

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(not EDGES.exists(), reason="Track-R canonical tables not built"),
]

PUBLISHED_NODES = 129_375
PUBLISHED_EDGES = 8_100_498
TOLERANCE = 0.01


@pytest.fixture(scope="module")
def nodes() -> pl.DataFrame:
    return pl.read_parquet(NODES)


@pytest.fixture(scope="module")
def edges() -> pl.DataFrame:
    return pl.read_parquet(EDGES)


def test_counts_within_tolerance(nodes, edges):
    assert abs(nodes.height - PUBLISHED_NODES) / PUBLISHED_NODES <= TOLERANCE
    assert abs(edges.height - PUBLISHED_EDGES) / PUBLISHED_EDGES <= TOLERANCE


def test_fk_integrity(nodes, edges):
    node_ids = set(nodes["node_id"].to_list())
    dangling = edges.filter(~pl.col("src_id").is_in(node_ids) | ~pl.col("dst_id").is_in(node_ids))
    assert dangling.height == 0


def test_reciprocal_coassignment_on_subsample(edges):
    """1% group-level subsample: reciprocal twins of symmetric relations must
    either both be in or both be out (co-assignment guard, §13.1)."""
    sub = edges.filter(pl.col("evidence_group_id").str.slice(0, 2).is_in(["00", "01", "02"]))
    assert sub.height > 0
    audit_relations = {r["relation"] for r in sub.group_by("relation").agg(n=pl.len()).to_dicts()}
    for rel in audit_relations:
        pairs = sub.filter(pl.col("relation") == rel)
        # group id uniqueness per unordered pair: twins share an id by construction;
        # check no group id maps to two different unordered endpoint sets
        per_group = pairs.group_by("evidence_group_id").agg(
            lo=pl.min_horizontal("src_id", "dst_id").n_unique(),
            hi=pl.max_horizontal("src_id", "dst_id").n_unique(),
        )
        assert per_group.filter((pl.col("lo") > 1) | (pl.col("hi") > 1)).height == 0


def test_node_id_curie_prefixes(nodes):
    prefixes = {n.split(":", 1)[0] for n in nodes["node_id"].head(1000).to_list()}
    assert all(":" in n for n in nodes["node_id"].head(1000).to_list())
    assert prefixes  # CURIE form type:id
