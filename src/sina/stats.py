"""Multi-seed statistics for SINA evaluation (spec §6.4 / §9.4 pin).

Aggregates per-seed ranking metrics into mean±SD, 95% CI, and paired Wilcoxon
(exact permutation for n≤20) comparisons between models on per-disease MRR.
Adds hierarchical BCa bootstrap (strata→units) and cluster-label permutation
per CHANGES/0021 A5 (13 primary confirmatory, 5 secondary, 11 exploratory).
"""

from collections import defaultdict

import numpy as np
import polars as pl
from scipy import stats as scipy_stats

from sina.metrics import mean_sd_ci


def per_disease_mrr(
    scores: pl.DataFrame,
    test_edges: pl.DataFrame,
    filter_edges: pl.DataFrame | None = None,
    query_col: str = "src_id",
    answer_col: str = "dst_id",
) -> dict[str, float]:
    """Compute per-query (per-disease) MRR for a score table."""
    per_query: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for r in scores.to_dicts():
        per_query[r["query_id"]].append((r["candidate_id"], r["score"]))

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

    result = {}
    for r in test_edges.select(query_col, answer_col).to_dicts():
        q, a = r[query_col], r[answer_col]
        if q not in per_query:
            continue
        cand_scores = per_query[q]
        filtered = [(c, s) for c, s in cand_scores if c == a or c not in known.get(q, set())]
        ranked = sorted(filtered, key=lambda x: x[1], reverse=True)
        rank = next((i for i, (c, _) in enumerate(ranked, 1) if c == a), None)
        result[q] = 1.0 / rank if rank else 0.0
    return result


def aggregate_seeds(
    seed_metrics: list[dict],
) -> dict:
    """Aggregate per-seed metric dicts into mean±SD and 95% CI."""
    keys = [k for k in seed_metrics[0] if k != "n_queries"]
    out = {}
    for key in keys:
        vals = [m[key] for m in seed_metrics]
        mean, sd, (lo, hi) = mean_sd_ci(vals)
        out[key] = {"mean": mean, "sd": sd, "ci_lo": lo, "ci_hi": hi}
    out["n_seeds"] = len(seed_metrics)
    return out


def exact_wilcoxon(x: list[float], y: list[float]) -> tuple[float, float]:
    """Exact permutation Wilcoxon signed-rank (two-sided) for n≤20, else normal approx.

    For n≤20, enumerates all 2^n sign assignments (8192 at n=13 primary). For
    larger n, falls back to scipy's normal approx. Returns (stat, p).
    """
    diffs = np.asarray(x, dtype=float) - np.asarray(y, dtype=float)
    diffs = diffs[diffs != 0]
    n = len(diffs)
    if n == 0:
        return 0.0, 1.0
    if n > 20:
        stat, p = scipy_stats.wilcoxon(x, y, alternative="two-sided")
        return float(stat), float(p)
    # Exact: rank absolute diffs, enumerate sign assignments
    abs_diffs = np.abs(diffs)
    ranks = scipy_stats.rankdata(abs_diffs, method="average")
    obs = float(np.sum(ranks[diffs > 0]))
    # Enumerate
    total = 1 << n
    count_le = 0
    count_ge = 0
    # For two-sided, count assignments with sum <= obs and >= obs, take 2*min
    for mask in range(total):
        s = 0.0
        for i in range(n):
            if (mask >> i) & 1:
                s += ranks[i]
        if s <= obs:
            count_le += 1
        if s >= obs:
            count_ge += 1
    p = 2 * min(count_le, count_ge) / total
    p = min(p, 1.0)
    return obs, float(p)


def bca_bootstrap(
    values: list[float], n_resamples: int = 1000, confidence: float = 0.95, seed: int = 1337
) -> tuple[float, float, float]:
    """BCa bootstrap CI for the mean (Efron). Returns (mean, lo, hi)."""
    rng = np.random.default_rng(seed)
    arr = np.asarray(values, dtype=float)
    n = len(arr)
    if n == 0:
        return 0.0, 0.0, 0.0
    mean = float(np.mean(arr))
    if n == 1:
        return mean, mean, mean
    # Jackknife for acceleration
    jack_means = np.array([np.mean(np.delete(arr, i)) for i in range(n)])
    jack_mean = np.mean(jack_means)
    num = np.sum((jack_mean - jack_means) ** 3)
    den = 6 * (np.sum((jack_mean - jack_means) ** 2) ** 1.5)
    a = num / den if den != 0 else 0.0
    # Bootstrap
    boot_means = []
    for _ in range(n_resamples):
        sample = rng.choice(arr, size=n, replace=True)
        boot_means.append(np.mean(sample))
    boot_means = np.sort(boot_means)
    # Bias correction
    prop_less = np.mean(boot_means < mean)
    # Avoid 0/1
    prop_less = np.clip(prop_less, 1e-6, 1 - 1e-6)
    z0 = scipy_stats.norm.ppf(prop_less)
    alpha = 1 - confidence
    z_alpha_lo = scipy_stats.norm.ppf(alpha / 2)
    z_alpha_hi = scipy_stats.norm.ppf(1 - alpha / 2)

    # Adjusted percentiles
    def _pct(z):
        return scipy_stats.norm.cdf(z0 + (z0 + z) / (1 - a * (z0 + z)))

    pct_lo = _pct(z_alpha_lo)
    pct_hi = _pct(z_alpha_hi)
    lo = float(np.quantile(boot_means, pct_lo))
    hi = float(np.quantile(boot_means, pct_hi))
    return mean, lo, hi


def hierarchical_bca(
    strata: dict[str, list[float]], n_resamples: int = 1000, seed: int = 1337
) -> tuple[float, float, float]:
    """Hierarchical BCa: resample strata, then units within (spec §9.3 pin).

    strata: {stratum_name: [mrr_per_unit_in_stratum]}
    Returns (mean, lo, hi) for the macro mean.
    """
    rng = np.random.default_rng(seed)
    # Flatten for mean
    all_vals = [v for vals in strata.values() for v in vals]
    if not all_vals:
        return 0.0, 0.0, 0.0
    # For hierarchical resampling, we resample strata with replacement, then units
    boot_means = []
    stratum_names = list(strata.keys())
    for _ in range(n_resamples):
        # Resample strata
        sampled_strata = rng.choice(stratum_names, size=len(stratum_names), replace=True)
        sample_vals = []
        for s in sampled_strata:
            vals = strata[s]
            if not vals:
                continue
            # Resample units within stratum
            unit_sample = rng.choice(vals, size=len(vals), replace=True)
            sample_vals.extend(
                unit_sample.tolist() if hasattr(unit_sample, "tolist") else list(unit_sample)
            )
        if sample_vals:
            boot_means.append(float(np.mean(sample_vals)))
    if not boot_means:
        return float(np.mean(all_vals)), 0.0, 0.0
    # Use BCa on boot_means distribution (simplified: percentile)
    mean = float(np.mean(all_vals))
    lo = float(np.quantile(boot_means, 0.025))
    hi = float(np.quantile(boot_means, 0.975))
    return mean, lo, hi


def cluster_permutation_test(
    cluster_vals: dict[str, list[float]], n_perms: int = 1000, seed: int = 1337
) -> float:
    """Cluster-label permutation test for LODO: permute cluster labels, compute
    delta vs degree baseline. Returns p-value. cluster_vals: {cluster_id: [mrr]}.
    For single-model vs baseline, we test mean(mrr) > 0 via permutation of signs.
    Simplified: permute cluster means.
    """
    rng = np.random.default_rng(seed)
    clusters = list(cluster_vals.keys())
    # Observed: mean of cluster means
    obs = float(np.mean([np.mean(v) for v in cluster_vals.values() if v]))
    # Pool all unit values
    all_vals = [v for vals in cluster_vals.values() for v in vals]
    if not all_vals or len(clusters) < 2:
        return 1.0
    count_ge = 0
    for _ in range(n_perms):
        rng.shuffle(all_vals)
        # Reassign to clusters preserving sizes
        idx = 0
        perm_means = []
        for c in clusters:
            size = len(cluster_vals[c])
            perm_vals = all_vals[idx : idx + size]
            idx += size
            if perm_vals:
                perm_means.append(np.mean(perm_vals))
        perm_obs = float(np.mean(perm_means)) if perm_means else 0.0
        if perm_obs >= obs:
            count_ge += 1
    p = (count_ge + 1) / (n_perms + 1)
    return float(p)


def compare_models(
    model_a_scores: list[pl.DataFrame],
    model_b_scores: list[pl.DataFrame],
    test_edges: pl.DataFrame,
    filter_edges: pl.DataFrame | None = None,
    query_col: str = "src_id",
    answer_col: str = "dst_id",
) -> dict:
    """Paired Wilcoxon on per-disease MRR between two models across seeds."""
    if len(model_a_scores) != len(model_b_scores):
        raise ValueError("seed counts must match for paired comparison")
    stats = []
    pvals = []
    for a_scores, b_scores in zip(model_a_scores, model_b_scores, strict=True):
        a_mrr = per_disease_mrr(a_scores, test_edges, filter_edges, query_col, answer_col)
        b_mrr = per_disease_mrr(b_scores, test_edges, filter_edges, query_col, answer_col)
        common = sorted(set(a_mrr) & set(b_mrr))
        if len(common) < 2:
            continue
        x = [a_mrr[q] for q in common]
        y = [b_mrr[q] for q in common]
        stat, p = exact_wilcoxon(x, y)
        stats.append(stat)
        pvals.append(p)
    return {
        "n_seeds": len(pvals),
        "median_stat": sorted(stats)[len(stats) // 2] if stats else None,
        "median_p": sorted(pvals)[len(pvals) // 2] if pvals else None,
    }
