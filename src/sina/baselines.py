"""Baseline scorers and null controls for SINA link-prediction evaluation.

All functions operate on the training graph only and return a long-format
DataFrame of (query_id, candidate_id, score) for downstream ranking and
metrics (spec §6.3). They are deterministic: same inputs + seed → same scores.
"""

from collections import Counter
from collections.abc import Iterable

import numpy as np
import polars as pl
import scipy.sparse as sp


def degree_baseline(
    train_edges: pl.DataFrame,
    candidate_nodes: Iterable[str],
    node_cols: tuple[str, str] = ("src_id", "dst_id"),
) -> pl.DataFrame:
    """Score every candidate by its total degree in the training graph."""
    src_col, dst_col = node_cols
    candidates = set(candidate_nodes)
    counts = Counter()
    for col in (src_col, dst_col):
        for node in train_edges[col].to_list():
            counts[node] += 1
    rows = [(n, float(counts[n])) for n in candidates]
    return pl.DataFrame({"node_id": [r[0] for r in rows], "degree": [r[1] for r in rows]})


def _build_indexed_graph(
    edges: pl.DataFrame, node_cols: tuple[str, str]
) -> tuple[sp.csr_matrix, dict[str, int], list[str]]:
    """Build an undirected, unweighted, symmetric-normalized CSR graph."""
    src_col, dst_col = node_cols
    nodes = sorted(set(edges[src_col].unique().to_list()) | set(edges[dst_col].unique().to_list()))
    idx = {n: i for i, n in enumerate(nodes)}
    src = edges[src_col].replace_strict(idx).to_numpy()
    dst = edges[dst_col].replace_strict(idx).to_numpy()
    n = len(nodes)
    data = np.ones(len(src) * 2, dtype=np.float64)
    rows = np.concatenate([src, dst])
    cols = np.concatenate([dst, src])
    adj = sp.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    # Symmetric normalization D^{-1/2} A D^{-1/2}
    deg = np.asarray(adj.sum(axis=1)).ravel()
    with np.errstate(divide="ignore"):
        d_inv_sqrt = np.power(deg, -0.5)
    d_inv_sqrt[np.isinf(d_inv_sqrt)] = 0.0
    d_mat = sp.diags(d_inv_sqrt)
    return d_mat @ adj @ d_mat, idx, nodes


def ppr_baseline(
    train_edges: pl.DataFrame,
    query_nodes: Iterable[str],
    candidate_nodes: Iterable[str],
    alpha: float = 0.15,
    max_iter: int = 100,
    tol: float = 1e-6,
    node_cols: tuple[str, str] = ("src_id", "dst_id"),
) -> pl.DataFrame:
    """Personalized PageRank from each query node on the training graph.

    Returns one row per (query_node, candidate_node) pair. Uses power iteration
    on a symmetric-normalized adjacency matrix (random-walk restart).
    """
    query_nodes = list(query_nodes)
    candidate_nodes = list(candidate_nodes)
    if not query_nodes or not candidate_nodes:
        return pl.DataFrame(
            {"query_id": [], "candidate_id": [], "ppr_score": []},
            schema={"query_id": pl.String, "candidate_id": pl.String, "ppr_score": pl.Float64},
        )

    norm_adj, idx, all_nodes = _build_indexed_graph(train_edges, node_cols)
    n = norm_adj.shape[0]
    candidate_indices = [idx[node] for node in candidate_nodes if node in idx]

    rows: list[tuple[str, str, float]] = []
    for q in query_nodes:
        if q not in idx:
            continue
        q_vec = np.zeros(n)
        q_vec[idx[q]] = 1.0
        p = q_vec.copy()
        for _ in range(max_iter):
            p_new = (1 - alpha) * (norm_adj @ p) + alpha * q_vec
            if np.linalg.norm(p_new - p, ord=1) < tol:
                p = p_new
                break
            p = p_new
        for ci in candidate_indices:
            cand = all_nodes[ci]
            if cand == q:
                continue
            rows.append((q, cand, float(p[ci])))

    return pl.DataFrame(
        {
            "query_id": [r[0] for r in rows],
            "candidate_id": [r[1] for r in rows],
            "ppr_score": [r[2] for r in rows],
        },
        schema={
            "query_id": pl.String,
            "candidate_id": pl.String,
            "ppr_score": pl.Float64,
        },
    )


def degree_matched_nulls(
    train_edges: pl.DataFrame,
    true_edges: pl.DataFrame,
    candidate_pool: Iterable[str],
    n_negatives: int = 1,
    bins: int = 10,
    seed: int = 1337,
    node_cols: tuple[str, str] = ("src_id", "dst_id"),
) -> pl.DataFrame:
    """For each true (query, answer) edge, sample negatives with degree-matched genes.

    The null preserves the query node and replaces the answer node with a gene
    drawn from the same degree bin (populariy-bias control, spec §0.2.4).
    """
    rng = np.random.default_rng(seed)
    src_col, dst_col = node_cols
    candidates = list(candidate_pool)
    deg = degree_baseline(train_edges, candidates, node_cols=node_cols)
    deg = deg.with_columns(
        bin=pl.col("degree").qcut(bins, allow_duplicates=True, include_breaks=False).cast(pl.String)
    )
    bin_to_candidates = {b: rows["node_id"].to_list() for b, rows in deg.group_by("bin")}
    answer_deg = degree_baseline(
        train_edges, set(true_edges[dst_col].unique().to_list()), node_cols=node_cols
    )
    answer_deg = answer_deg.with_columns(
        bin=pl.col("degree").qcut(bins, allow_duplicates=True, include_breaks=False).cast(pl.String)
    )
    answer_bin = {r["node_id"]: r["bin"] for r in answer_deg.to_dicts()}

    rows: list[tuple[str, str, str, int]] = []
    for r in true_edges.select(src_col, dst_col).to_dicts():
        q, a = r[src_col], r[dst_col]
        b = answer_bin.get(a)
        pool = [c for c in bin_to_candidates.get(b, candidates) if c != a]
        if not pool:
            pool = [c for c in candidates if c != a]
        chosen = rng.choice(pool, size=min(n_negatives, len(pool)), replace=False)
        for neg in chosen:
            rows.append((q, a, neg, 1))
    return pl.DataFrame(
        {
            src_col: [r[0] for r in rows],
            "true_" + dst_col: [r[1] for r in rows],
            "null_" + dst_col: [r[2] for r in rows],
            "label": [r[3] for r in rows],
        }
    )


def relation_shuffled_nulls(
    edges: pl.DataFrame,
    relation_col: str = "relation",
    seed: int = 1337,
) -> pl.DataFrame:
    """Return a copy of `edges` with relation labels shuffled (structural null)."""
    rng = np.random.default_rng(seed)
    relations = edges[relation_col].to_list()
    rng.shuffle(relations)
    return edges.with_columns(pl.Series(relation_col, relations))
