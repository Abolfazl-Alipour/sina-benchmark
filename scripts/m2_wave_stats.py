"""m2_wave_stats — pre-registered multi-seed stats for Phase D waves (CHANGES/0024).

Joins the per-seed Wave results JSONs (run_wave1_benchmark.py) with the
KEEP-graph baselines parquet and runs the frozen stats pins (CHANGES/0020
§13 / spec §6.4, §9.4) over the 13 primary units:

- exact Wilcoxon signed-rank (2^13 enumeration) of model MRR vs the best
  same-graph null (max of degree/PPR MRR per unit)
- hierarchical BCa bootstrap (>=1000 resamples; strata = frozen co-holdout
  clusters from configs/m2_lodo_exports.yaml)
- cluster-label permutation test

The frozen degenerate rule (A3) would exclude flagged units from the primary
gate. Because the KEEP triviality re-run flags all 13 primaries under the
≤2-hop prong while no baseline reaches MRR>=0.99 (conflict on record,
CHANGES/0024 §(b) — owner decision pending), this script reports BOTH the
all-13-primaries stats and the degenerate-excluded stats, clearly labeled;
it does not resolve the conflict.

Usage:
  python scripts/m2_wave_stats.py --config configs/m2_wave_stats.yaml
"""

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import polars as pl
import yaml
from pydantic import BaseModel

from sina.config import load_config
from sina.logging import configure_logging, get_logger
from sina.stats import (
    aggregate_seeds,
    cluster_permutation_test,
    exact_wilcoxon,
    hierarchical_bca,
)

log = get_logger(__name__)


class WaveStatsConfig(BaseModel):
    results_jsons: list[Path]  # one per seed (wave1 seed 1337, then 1339/5027)
    lodo_baselines: Path
    triviality: Path
    exports_config: Path  # configs/m2_lodo_exports.yaml — frozen clusters source
    report_path: Path
    results_json: Path
    manifest_path: Path
    n_resamples: int = 1000
    seed: int = 1337


def utcnow_iso() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, WaveStatsConfig)

    # Per-seed per-unit model MRR (primaries only — results files may carry
    # other tiers in later waves; the stats pin covers the 13 primaries).
    per_seed: list[dict[str, float]] = []
    seeds = []
    for path in cfg.results_jsons:
        payload = json.loads(path.read_text(encoding="utf-8"))
        seeds.append(payload.get("seed"))
        per_seed.append({r["curie"]: float(r["test_mrr"]) for r in payload["results"]})

    baselines = {r["unit"]: r for r in pl.read_parquet(cfg.lodo_baselines).to_dicts()}
    triv = {r["unit"]: bool(r["degenerate"]) for r in pl.read_parquet(cfg.triviality).to_dicts()}
    frozen = yaml.safe_load(cfg.exports_config.read_text(encoding="utf-8"))
    unit_cluster = {u["curie"]: u["cluster"] for u in frozen["units"]}

    primaries = sorted(u["curie"] for u in frozen["units"] if u["tier"] == "primary")
    missing = [u for u in primaries if any(u not in s for s in per_seed)]
    if missing:
        raise SystemExit(f"results missing primary units: {missing}")

    # Mean MRR per unit across seeds + aggregate_seeds over per-seed macro MRR
    model_mrr = {u: sum(s[u] for s in per_seed) / len(per_seed) for u in primaries}
    null_mrr = {
        u: max(
            float(baselines[u].get("degree_filtered_mrr", 0.0)),
            float(baselines[u].get("ppr_filtered_mrr", 0.0)),
        )
        for u in primaries
    }
    seed_macro = [
        {"mrr": sum(s[u] for u in primaries) / len(primaries), "n_queries": len(primaries)}
        for s in per_seed
    ]
    macro = aggregate_seeds(seed_macro)["mrr"]

    def gate_stats(units: list[str]) -> dict:
        x = [model_mrr[u] for u in units]
        y = [null_mrr[u] for u in units]
        stat, p = exact_wilcoxon(x, y)
        strata: dict[str, list[float]] = {}
        deltas: dict[str, list[float]] = {}
        for u in units:
            strata.setdefault(unit_cluster[u], []).append(model_mrr[u])
            deltas.setdefault(unit_cluster[u], []).append(model_mrr[u] - null_mrr[u])
        mean, lo, hi = hierarchical_bca(strata, n_resamples=cfg.n_resamples, seed=cfg.seed)
        perm_p = cluster_permutation_test(deltas, n_perms=cfg.n_resamples, seed=cfg.seed)
        return {
            "n_units": len(units),
            "units": units,
            "exact_wilcoxon_stat": stat,
            "exact_wilcoxon_p": p,
            "hierarchical_bca_mean": mean,
            "hierarchical_bca_ci": [lo, hi],
            "cluster_permutation_p": perm_p,
            "n_resamples": cfg.n_resamples,
        }

    all13 = gate_stats(primaries)
    non_degen = [u for u in primaries if not triv.get(u, False)]
    excluded = gate_stats(non_degen) if len(non_degen) >= 2 else {"note": "<2 non-degenerate units"}

    results = {
        "created_at": utcnow_iso(),
        "seeds": seeds,
        "n_seeds": len(seeds),
        "macro_mrr_over_seeds": macro,
        "per_unit": {
            u: {
                "model_mrr_mean": model_mrr[u],
                "best_null_mrr": null_mrr[u],
                "cluster": unit_cluster[u],
                "degenerate_keep_a3": triv.get(u, False),
            }
            for u in primaries
        },
        "primary_gate_all_13": all13,
        "primary_gate_degenerate_excluded": excluded,
    }
    cfg.results_json.parent.mkdir(parents=True, exist_ok=True)
    results_text = json.dumps(results, indent=2) + "\n"
    cfg.results_json.write_text(results_text, encoding="utf-8")

    def fmt_gate(g: dict) -> list[str]:
        if "note" in g:
            return [f"- {g['note']}"]
        return [
            f"- units: {g['n_units']}",
            f"- exact Wilcoxon (2^n enumeration): stat {g['exact_wilcoxon_stat']:.1f}, "
            f"p = {g['exact_wilcoxon_p']:.4f}",
            f"- hierarchical BCa (n={g['n_resamples']}): mean {g['hierarchical_bca_mean']:.5f}, "
            f"95% CI [{g['hierarchical_bca_ci'][0]:.5f}, {g['hierarchical_bca_ci'][1]:.5f}]",
            f"- cluster-label permutation: p = {g['cluster_permutation_p']:.4f}",
        ]

    lines = [
        "<!-- track_scope: benchmark -->",
        "# Phase D wave stats — pre-registered pins (CHANGES/0020 §13, CHANGES/0024)",
        "",
        (
            f"Model vs best same-graph null (degree/PPR) on per-unit test MRR, "
            f"{len(seeds)} seed(s) {seeds}, 13 primary LODO units, KEEP graph "
            f"(train1 = 1,431,945). Macro test MRR over seeds: "
            f"{macro['mean']:.5f} ± {macro['sd']:.5f} "
            f"(95% CI [{macro['ci_lo']:.5f}, {macro['ci_hi']:.5f}]). "
            "Disease-associated gene recovery benchmark; MRR exact."
        ),
        "",
        (
            "**Degenerate-flag conflict (CHANGES/0024 §(b), owner decision pending):** "
            "the frozen A3 ≤2-hop prong flags all 13 primaries on the KEEP graph. "
            "Both variants are reported; neither preempts the owner call."
        ),
        "",
        "## Primary gate — all 13 primaries",
        *fmt_gate(all13),
        "",
        "## Primary gate — degenerate-excluded (frozen A3, literal)",
        *fmt_gate(excluded),
        "",
        "## Per-unit",
        "",
        "| unit | cluster | model MRR (seed mean) | best null MRR | A3 degenerate (KEEP) |",
        "|---|---|---|---|---|",
        *(
            f"| {u} | {unit_cluster[u]} | {model_mrr[u]:.5f} | {null_mrr[u]:.4f} | "
            f"{triv.get(u, False)} |"
            for u in sorted(primaries, key=lambda k: -model_mrr[k])
        ),
        "",
        (
            "Reproduction: `python scripts/m2_wave_stats.py --config "
            "configs/m2_wave_stats.yaml`. With 1 seed (Wave 1 only) all stats are "
            "directional; the pre-registered read requires seeds {1337,1339,5027}."
        ),
    ]
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)
    text = "\n".join(lines) + "\n"
    cfg.report_path.write_text(text, encoding="utf-8")

    cfg.manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = {
        "script": "m2_wave_stats",
        "created_at": utcnow_iso(),
        "status": "done",
        "outputs": {
            str(cfg.report_path): hashlib.sha256(text.encode("utf-8")).hexdigest(),
            str(cfg.results_json): hashlib.sha256(results_text.encode("utf-8")).hexdigest(),
        },
    }
    cfg.manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    log.info("wave_stats_done", n_seeds=len(seeds), wilcoxon_p=all13["exact_wilcoxon_p"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
