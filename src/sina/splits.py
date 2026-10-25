"""Split primitives (spec §6.2) — the leakage-critical contract in code.

HARD RULE (§5.4): splits partition **evidence_group_id** sets, never edge rows.
Every function here returns {"train", "valid", "test"} edge-row frames whose
group-id sets are pairwise disjoint — enforced by :func:`assert_split_disjoint`,
which is also wired into CI guard tests (§13.1).

All randomness is seeded and deterministic (same inputs + seed → same split).
"""

import hashlib

import polars as pl

GROUP_COL = "evidence_group_id"


class SplitLeakageError(RuntimeError):
    """A reciprocal/source-duplicate twin straddles a split boundary (§13.1)."""


def assert_split_disjoint(split: dict[str, pl.DataFrame]) -> None:
    """Fail if any evidence_group_id appears in more than one partition."""
    seen: dict[str, str] = {}
    for name, frame in split.items():
        for gid in frame[GROUP_COL].unique().to_list():
            if gid in seen:
                raise SplitLeakageError(f"evidence group {gid} in both {seen[gid]!r} and {name!r}")
            seen[gid] = name


def group_assignment(group_ids: pl.Series, fractions: dict[str, float], seed: int) -> pl.DataFrame:
    """Deterministically assign each group id to a partition by hashing.

    Hashing (sha1 of seed+gid) instead of shuffling keeps the assignment
    stable across runs and independent of row order.
    """
    parts = list(fractions.items())
    total = sum(f for _, f in parts)

    def assign(gid: str) -> str:
        h = int.from_bytes(hashlib.sha1(f"{seed}|{gid}".encode()).digest()[:8], "big") / 2**64
        acc = 0.0
        for name, frac in parts:
            acc += frac / total
            if h < acc:
                return name
        return parts[-1][0]

    # One unique() call only: Series.unique() row order is nondeterministic
    # across calls, so building the two columns from two separate unique()
    # calls can misalign group IDs with their partitions.
    uniq = group_ids.unique().sort()
    return pl.DataFrame(
        {
            GROUP_COL: uniq,
            "partition": uniq.map_elements(assign, return_dtype=pl.String),
        }
    )


def _by_partition(edges: pl.DataFrame, assignment: pl.DataFrame) -> dict[str, pl.DataFrame]:
    joined = edges.join(assignment, on=GROUP_COL, how="left")
    split = {
        name: joined.filter(pl.col("partition") == name).drop("partition")
        for name in ("train", "valid", "test")
    }
    assert_split_disjoint(split)
    return split


def random_edge_split(
    edges: pl.DataFrame,
    fractions: dict[str, float] | None = None,
    seed: int = 1337,
) -> dict[str, pl.DataFrame]:
    """Random group-level split — comparability only, 'inflated by design' (§6.2.4)."""
    fractions = fractions or {"train": 0.9, "valid": 0.05, "test": 0.05}
    return _by_partition(edges, group_assignment(edges[GROUP_COL], fractions, seed))


def family_holdout_split(
    edges: pl.DataFrame,
    family_node_ids: set[str],
    eval_relation_col: str,
    eval_relation_values: set[str],
    seed: int = 1337,
    valid_fraction: float = 0.1,
) -> dict[str, pl.DataFrame]:
    """SINA disease-family holdout (§6.2.2) — strict zero-shot on one family.

    - eval edges (disease→gene relations) touching a family disease go to
      test/valid (by group, seeded);
    - ALL other edges touching a family disease are removed from train
      (their molecular neighborhood is hidden, TxGNN-style);
    - everything else is train.
    """
    touches_family = pl.col("src_id").is_in(family_node_ids) | pl.col("dst_id").is_in(
        family_node_ids
    )
    is_eval = touches_family & pl.col(eval_relation_col).is_in(eval_relation_values)
    eval_edges = edges.filter(is_eval)
    neighborhood_removed = edges.filter(touches_family & ~is_eval)
    train = edges.filter(~touches_family)
    assignment = group_assignment(
        eval_edges[GROUP_COL], {"valid": valid_fraction, "test": 1 - valid_fraction}, seed
    )
    joined = eval_edges.join(assignment, on=GROUP_COL, how="left")
    split = {
        "train": train,
        "valid": joined.filter(pl.col("partition") == "valid").drop("partition"),
        "test": joined.filter(pl.col("partition") == "test").drop("partition"),
    }
    assert_split_disjoint(split)
    split["_removed_neighborhood"] = neighborhood_removed
    return split


def txgnn_area_split(
    edges: pl.DataFrame,
    area_disease_ids: set[str],
    indication_relations: set[str],
    relation_col: str = "relation",
    neighborhood_remove_fraction: float = 0.95,
    seed: int = 1337,
) -> dict[str, pl.DataFrame]:
    """TxGNN-protocol split (§6.2.1), one disease area.

    Held-out diseases: their indication/contraindication edges → test; ~95% of
    their remaining (neighborhood) edges are removed from train; the rest of
    the graph is train. All at evidence-group granularity.
    """
    touches = pl.col("src_id").is_in(area_disease_ids) | pl.col("dst_id").is_in(area_disease_ids)
    is_indication = touches & pl.col(relation_col).is_in(indication_relations)
    test = edges.filter(is_indication)
    neighborhood = edges.filter(touches & ~is_indication)
    keep = group_assignment(
        neighborhood[GROUP_COL],
        {"keep": 1 - neighborhood_remove_fraction, "drop": neighborhood_remove_fraction},
        seed,
    ).filter(pl.col("partition") == "keep")
    train = edges.filter(~touches)
    split = {"train": train, "valid": edges.head(0), "test": test}
    assert_split_disjoint(split)
    split["_kept_neighborhood"] = neighborhood.join(
        keep.select(GROUP_COL), on=GROUP_COL, how="semi"
    )
    return split


def temporal_split(
    edges: pl.DataFrame,
    cutoff_year: int,
    label_value: str = "DIS-GEN",
) -> dict[str, pl.DataFrame]:
    """Temporal split (§6.2.3 — SECONDARY on Track-D v2.0, see Appendix D.5).

    train: temporal edges with year_final ≤ cutoff; test: year_initial > cutoff.
    Non-temporal edges (no year fields) stay in train. DIS-GEN edges without
    years are excluded from both partitions' eval role but remain in train.
    """
    # Compute group-id sets on the small DIS-GEN subset first, then filter the
    # full frame once — avoids holding several full-size copies at peak RSS.
    temporal = edges.filter(pl.col("label") == label_value)
    year_init = pl.col("year_initial").cast(pl.Int32, strict=False)
    year_final = pl.col("year_final").cast(pl.Int32, strict=False)
    has_years = year_init.is_not_null()
    test_gids = temporal.filter(has_years & (year_init > cutoff_year))[GROUP_COL].unique()
    # train: year_final ≤ cutoff, or no temporal annotation at all; boundary
    # straddlers (year_initial ≤ cutoff < year_final) land in neither set and
    # are dropped entirely.
    train_gids = temporal.filter(~has_years | (year_final <= cutoff_year))[GROUP_COL].unique()
    del temporal
    test_gid_frame = test_gids.to_frame(GROUP_COL)
    test = edges.join(test_gid_frame, on=GROUP_COL, how="semi")
    train = edges.filter(
        (pl.col("label") != label_value) | pl.col(GROUP_COL).is_in(train_gids.implode())
    ).join(test_gid_frame, on=GROUP_COL, how="anti")
    split = {
        "train": train,
        "valid": edges.head(0),
        "test": test,
    }
    assert_split_disjoint(split)
    return split
