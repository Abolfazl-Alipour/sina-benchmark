"""Ranking metrics and statistics for SINA link-prediction evaluation (spec §6.4/§6.5).

All functions are deterministic and operate on long-format score tables
(query_id, candidate_id, score). Filtered metrics remove known true edges
(train+valid+test) from the candidate ranking before computing ranks.
"""

import bisect
from collections import defaultdict

import numpy as np
import polars as pl
from scipy import stats


def link_prediction_metrics(
    scores: pl.DataFrame,
    test_edges: pl.DataFrame,
    filter_edges: pl.DataFrame | None = None,
    query_col: str = "src_id",
    answer_col: str = "dst_id",
    ks: tuple[int, ...] = (1, 3, 10),
) -> dict:
    """Compute raw and filtered MRR + Hits@k for a score table.

    `scores` must have columns query_id, candidate_id, score. Higher score =
    better. For each test edge, the true answer is ranked among all candidates
    that appear in the score table for that query.
    """
    # Build known-true set per query for filtering
    known: dict[str, set[str]] = defaultdict(set)
    if filter_edges is not None:
        known_pairs = set(
            zip(
                filter_edges[query_col].to_list(),
                filter_edges[answer_col].to_list(),
                strict=False,
            )
        )
        for q, a in known_pairs:
            known[q].add(a)

    # Group scores by query
    per_query: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for r in scores.to_dicts():
        per_query[r["query_id"]].append((r["candidate_id"], r["score"]))

    results = {
        "raw": {f"hits@{k}": 0.0 for k in ks} | {"mrr": 0.0},
        "filtered": {f"hits@{k}": 0.0 for k in ks} | {"mrr": 0.0},
        "n_queries": 0,
    }

    # Rank once per unique query, then compute filtered rank per test edge by
    # binary-searching the answer's score in the filtered candidate list.
    raw_cache: dict[str, dict[str, int]] = {}
    filt_scores: dict[str, list[float]] = {}
    score_map: dict[str, dict[str, float]] = {}
    for q, cand_scores in per_query.items():
        sorted_raw = sorted(cand_scores, key=lambda x: x[1], reverse=True)
        raw_cache[q] = {c: i for i, (c, _) in enumerate(sorted_raw, start=1)}
        score_map[q] = {c: s for c, s in cand_scores}
        filtered = [s for c, s in sorted_raw if c not in known.get(q, set())]
        filt_scores[q] = sorted(filtered, reverse=True)

    n = 0
    for r in test_edges.select(query_col, answer_col).to_dicts():
        q, a = r[query_col], r[answer_col]
        if q not in raw_cache:
            continue
        n += 1
        rank_raw = raw_cache[q].get(a)
        answer_score = score_map[q].get(a)
        rank_filt = None
        if answer_score is not None:
            # number of filtered candidates with score > answer_score
            rank_filt = bisect.bisect_left([-s for s in filt_scores[q]], -answer_score) + 1

        results["raw"]["mrr"] += 1.0 / rank_raw if rank_raw else 0.0
        results["filtered"]["mrr"] += 1.0 / rank_filt if rank_filt else 0.0
        for k in ks:
            if rank_raw and rank_raw <= k:
                results["raw"][f"hits@{k}"] += 1.0
            if rank_filt and rank_filt <= k:
                results["filtered"][f"hits@{k}"] += 1.0

    results["n_queries"] = n
    if n:
        results["raw"]["mrr"] /= n
        results["filtered"]["mrr"] /= n
        for k in ks:
            results["raw"][f"hits@{k}"] /= n
            results["filtered"][f"hits@{k}"] /= n
    return results


def _rank_of(ranked: list[tuple[str, float]], target: str) -> int | None:
    """1-based rank of target in a sorted list, or None if absent."""
    for i, (c, _) in enumerate(ranked, start=1):
        if c == target:
            return i
    return None


def degree_stratified_metrics(
    scores: pl.DataFrame,
    test_edges: pl.DataFrame,
    train_edges: pl.DataFrame,
    filter_edges: pl.DataFrame | None = None,
    query_col: str = "src_id",
    answer_col: str = "dst_id",
    ks: tuple[int, ...] = (1, 3, 10),
    n_bins: int = 4,
) -> dict:
    """Compute metrics per degree bin of the true answer (popularity stratification)."""
    # Degree of each answer node in the training graph
    degree: dict[str, int] = defaultdict(int)
    for col in ("src_id", "dst_id"):
        for node in train_edges[col].to_list():
            degree[node] += 1

    test_with_deg = test_edges.with_columns(
        answer_degree=pl.col(answer_col).map_elements(
            lambda x: degree.get(x, 0), return_dtype=pl.Int64
        )
    )
    # Bin by quantile
    deg_values = [degree.get(r[answer_col], 0) for r in test_edges.select(answer_col).to_dicts()]
    if not deg_values:
        return {}
    quantiles = np.quantile(deg_values, np.linspace(0, 1, n_bins + 1)[1:-1])

    def bin_of(d: int) -> int:
        return int(np.searchsorted(quantiles, d, side="right"))

    results = {}
    for b in range(n_bins):
        sub = test_with_deg.filter(
            pl.col("answer_degree").map_elements(
                lambda x, b=b: bin_of(x) == b, return_dtype=pl.Boolean
            )
        )
        if sub.height == 0:
            continue
        m = link_prediction_metrics(scores, sub, filter_edges, query_col, answer_col, ks)
        m["n_edges"] = sub.height
        results[f"degree_bin_{b}"] = m
    return results


def mean_sd_ci(
    values: list[float], confidence: float = 0.95
) -> tuple[float, float, tuple[float, float]]:
    """Mean, sample SD, and (1-alpha) CI for a list of values."""
    arr = np.asarray(values, dtype=float)
    n = len(arr)
    if n == 0:
        return 0.0, 0.0, (0.0, 0.0)
    mean = float(np.mean(arr))
    sd = float(np.std(arr, ddof=1)) if n > 1 else 0.0
    if n < 2:
        return mean, sd, (mean, mean)
    se = sd / np.sqrt(n)
    t = stats.t.ppf((1 + confidence) / 2, df=n - 1)
    return mean, sd, (mean - t * se, mean + t * se)


def paired_wilcoxon(x: list[float], y: list[float]) -> tuple[float, float]:
    """Paired Wilcoxon signed-rank test (two-sided). Returns (statistic, p-value)."""
    stat, p = stats.wilcoxon(np.asarray(x), np.asarray(y), alternative="two-sided")
    return float(stat), float(p)
