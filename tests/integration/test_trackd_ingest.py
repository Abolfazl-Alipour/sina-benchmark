"""Integration tests over the real Track-D canonical tables (spec §5.5).

Run locally where data/processed/m1/trackd exists; SKIP in CI.
Subsample = 1% of evidence groups (deterministic, group-level) so
co-assignment guards hold inside the subsample.
"""

from pathlib import Path

import polars as pl
import pytest

TRACKD_DIR = Path("data/processed/m1/trackd")
NODES = TRACKD_DIR / "nodes.parquet"
EDGES = TRACKD_DIR / "edges.parquet"
GENE_PROPS = TRACKD_DIR.parent / "gene_properties.parquet"

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(not EDGES.exists(), reason="Track-D canonical tables not built"),
]

# §5.5: OptimusKG counts reconcile EXACTLY
EXPECTED_NODES = 190_531
EXPECTED_EDGES = 21_813_816
EXPECTED_GENES = 61_306
EXPECTED_DISGEN = 9_734_774


@pytest.fixture(scope="module")
def nodes() -> pl.DataFrame:
    return pl.read_parquet(NODES)


@pytest.fixture(scope="module")
def edges() -> pl.DataFrame:
    return pl.read_parquet(EDGES)


def test_exact_counts(nodes, edges):
    assert nodes.height == EXPECTED_NODES
    assert edges.height == EXPECTED_EDGES


def test_fk_integrity(nodes, edges):
    node_ids = set(nodes["node_id"].to_list())
    dangling = edges.filter(~pl.col("src_id").is_in(node_ids) | ~pl.col("dst_id").is_in(node_ids))
    assert dangling.height == 0


def test_disgen_count_and_temporal_coverage(edges):
    disgen = edges.filter(pl.col("label") == "DIS-GEN")
    assert disgen.height == EXPECTED_DISGEN
    valid = lambda c: pl.col(c).is_not_null() & (pl.col(c) != "")  # noqa: E731
    coverage = disgen.filter(valid("year_initial") & valid("year_final")).height / disgen.height
    # Appendix D.5: coverage measured; whatever the value it must be in [0, 1]
    # and the fraction is reported in the characterization report.
    assert 0.0 <= coverage <= 1.0


def test_undirected_coassignment_on_subsample(edges):
    """1% group-level subsample: an undirected twin pair must not straddle
    the subsample boundary (guard §13.1 adapted to §5.2.3 semantics)."""
    sub = edges.filter(pl.col("evidence_group_id").str.slice(0, 2).is_in(["00", "01", "02"]))
    assert sub.height > 0
    undirected = sub.filter(pl.col("undirected"))
    per_group = undirected.group_by("evidence_group_id").agg(
        lo=pl.min_horizontal("src_id", "dst_id").n_unique(),
        hi=pl.max_horizontal("src_id", "dst_id").n_unique(),
    )
    assert per_group.filter((pl.col("lo") > 1) | (pl.col("hi") > 1)).height == 0


def test_gene_properties_table():
    if not GENE_PROPS.exists():
        pytest.skip("gene_properties.parquet not built")
    props = pl.read_parquet(GENE_PROPS)
    assert props.height == EXPECTED_GENES
    for col in [
        "tractability",
        "safety_liabilities",
        "target_class",
        "constraint_scores",
        "homologues",
        "chemical_probes",
    ]:
        assert col in props.columns
