"""Evidence group IDs — the leakage-critical contract (spec §5.4).

M2 splits may ONLY split on `evidence_group_id`. Therefore:

- **Reciprocal twins** of a symmetric (undirected) relation MUST share one ID.
- **Multi-source duplicates** (same assertion, several provenance rows) MUST
  share one ID.

Group ID = SHA-1 hex (first 16 chars) of the canonical edge key
``relation|a|b`` where ``(a, b)`` is sorted for symmetric relations and kept
ordered for directed ones. Which relations count as symmetric is decided
empirically per graph by :func:`reciprocal_audit` (a relation is symmetric when
essentially every edge's reverse twin is also present) — never hard-coded.
"""

import hashlib
from collections.abc import Iterable

import polars as pl

ID_LEN = 16


def canonical_edge_key(src: str, relation: str, dst: str, *, symmetric: bool) -> str:
    """Order-insensitive key for symmetric relations, ordered key otherwise."""
    a, b = sorted((src, dst)) if symmetric else (src, dst)
    return f"{relation}|{a}|{b}"


def evidence_group_id(src: str, relation: str, dst: str, *, symmetric: bool) -> str:
    """Stable group ID for one assertion."""
    key = canonical_edge_key(src, relation, dst, symmetric=symmetric)
    return hashlib.sha1(key.encode("utf-8")).hexdigest()[:ID_LEN]


def reciprocal_audit(
    edges: pl.DataFrame,
    src_col: str = "src_id",
    dst_col: str = "dst_id",
    rel_col: str = "relation",
) -> pl.DataFrame:
    """Per relation: how many edges have their reverse twin also present.

    Returns a frame with columns:
    ``relation, n_edges, n_with_reverse, reverse_fraction``.
    """
    present = edges.select(rel_col, src_col, dst_col).unique()
    reversed_view = present.rename({src_col: dst_col, dst_col: src_col})
    with_flag = present.join(
        reversed_view.with_columns(_has_reverse=pl.lit(True)),
        on=[rel_col, src_col, dst_col],
        how="left",
    ).with_columns(pl.col("_has_reverse").fill_null(False))
    return (
        with_flag.group_by(rel_col)
        .agg(
            n_edges=pl.len(),
            n_with_reverse=pl.col("_has_reverse").sum(),
        )
        .with_columns(reverse_fraction=pl.col("n_with_reverse") / pl.col("n_edges"))
        .sort(rel_col)
    )


def symmetric_relations(audit: pl.DataFrame, threshold: float = 0.99) -> set[str]:
    """Relations whose reverse-twin fraction meets the threshold."""
    return set(audit.filter(pl.col("reverse_fraction") >= threshold)["relation"].to_list())


def _source_set_key(sources_col: str) -> pl.Expr:
    """Canonical source-set key from a struct{direct, indirect} column.

    The whole-column cast comes first: parquet/test frames may type an
    all-empty list column as List(Null), on which ``list.join`` fails.
    """
    empty = pl.lit([], dtype=pl.List(pl.String))

    def canon(field: str) -> pl.Expr:
        return (
            pl.col(sources_col)
            .struct.field(field)
            .cast(pl.List(pl.String))  # List(Null) columns break list.join
            .fill_null(empty)
            .list.sort()
            .list.join(";")
        )

    return pl.concat_str([canon("direct"), canon("indirect")], separator="#")


def add_evidence_group_ids(
    edges: pl.DataFrame,
    symmetric: Iterable[str] = (),
    src_col: str = "src_id",
    dst_col: str = "dst_id",
    rel_col: str = "relation",
    undirected_col: str | None = None,
    sources_col: str | None = None,
) -> pl.DataFrame:
    """Append an ``evidence_group_id`` column (duplicates + twins share an ID).

    Symmetry is either per-relation (``symmetric`` set — Track-R, where the file
    materializes both directions of every edge) or per-row (``undirected_col``
    boolean — Track-D, spec §5.2.3). When ``sources_col`` names a
    struct{direct, indirect} column, the canonical source-set is part of the key
    (assertion+provenance granularity, spec §5.4).
    """
    sym = set(symmetric)
    is_sym = pl.col(undirected_col) if undirected_col else pl.col(rel_col).is_in(sorted(sym))
    lo = pl.min_horizontal(src_col, dst_col)
    hi = pl.max_horizontal(src_col, dst_col)
    key = (
        pl.when(is_sym)
        .then(pl.concat_str([rel_col, lo, hi], separator="|"))
        .otherwise(pl.concat_str([rel_col, src_col, dst_col], separator="|"))
    )
    if sources_col:
        key = pl.concat_str([key, _source_set_key(sources_col)], separator="|")
    keyed = edges.with_columns(_key=key)
    key_to_id = (
        keyed.select("_key")
        .unique()
        .with_columns(
            evidence_group_id=pl.col("_key").map_elements(
                lambda k: hashlib.sha1(k.encode("utf-8")).hexdigest()[:ID_LEN],
                return_dtype=pl.String,
            )
        )
    )
    return keyed.join(key_to_id, on="_key", how="left").drop("_key")
