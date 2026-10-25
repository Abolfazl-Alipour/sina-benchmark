"""m2_wave1_report — Phase D Wave 1 leaderboard report (CHANGES/0024).

Joins the Wave 1 benchmark results JSON (run_wave1_benchmark.py,
13 primary LODO units x seed 1337 x 1 pass on the KEEP 1,431,945-edge
exports) with the KEEP-graph baselines parquet and the KEEP triviality
audit, and emits reports/m2b/wave1_benchmark.md.

Framing (frozen): Wave 1 is directional only (1 seed, no significance);
island-band expectation MRR 0.01-0.03 — beating the degree/PPR nulls IS
the result. AUROC/AP are batch-averaged approximations (streaming-eval
patch) and are reported as such; MRR/Hits@k are exact. Claim ceiling:
disease-associated gene recovery only.

Usage:
  python scripts/m2_wave1_report.py --config configs/m2_wave1_report.yaml
"""

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import polars as pl
from pydantic import BaseModel

from sina.config import load_config
from sina.logging import configure_logging, get_logger

log = get_logger(__name__)


class Wave1ReportConfig(BaseModel):
    results_json: Path
    lodo_baselines: Path
    triviality: Path
    report_path: Path
    manifest_path: Path


def utcnow_iso() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, Wave1ReportConfig)

    payload = json.loads(cfg.results_json.read_text(encoding="utf-8"))
    results = payload["results"]
    baselines = {r["unit"]: r for r in pl.read_parquet(cfg.lodo_baselines).to_dicts()}
    triv = {r["unit"]: r for r in pl.read_parquet(cfg.triviality).to_dicts()}

    lines = [
        "<!-- track_scope: benchmark -->",
        "# Phase D Wave 1 — 13 primary LODO units (seed 1337, 1 pass, KEEP graph)",
        "",
        (
            f"BioPathNet (NBFNet full-graph message passing, CHANGES/0020 pinned config) "
            f"trained per unit on the keep_disphe exports (train1 = 1,431,945 edges incl. "
            f"DIS-PHE), batch {payload.get('batch_size')}, seed {payload.get('seed')}. "
            f"Run window: {payload.get('started_at')} -> "
            f"{payload.get('finished_at', 'in progress')}."
        ),
        "",
        (
            "Framing (frozen, CHANGES/0024): Wave 1 is **directional only** — 1 seed, no "
            "significance claims. Island-band expectation: MRR 0.01-0.03; beating the "
            "same-graph degree/PPR nulls IS the result. This benchmark measures "
            "disease-associated gene recovery under leave-one-disease-out."
        ),
        "",
        (
            "Caveats: AUROC/AP below are batch-averaged approximations (streaming-eval "
            "patch) — MRR/Hits@k are exact and are the only gate metrics. The frozen A3 "
            "≤2-hop triviality prong flags all 13 units on the KEEP graph "
            "(triviality_audit_keep.md); stats treatment of that conflict is pending "
            "owner decision (CHANGES/0024 §(b))."
        ),
        "",
        "| unit | name | test | model MRR | H@10 | H@100 | degree MRR | PPR MRR | lift | ≤2-hop% |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]

    n_beat = 0
    mrrs = []
    for r in sorted(results, key=lambda x: -x["test_mrr"]):
        b = baselines.get(r["curie"], {})
        t = triv.get(r["curie"], {})
        best_null = max(
            float(b.get("degree_filtered_mrr", 0.0)), float(b.get("ppr_filtered_mrr", 0.0))
        )
        lift = (r["test_mrr"] / best_null) if best_null > 0 else float("nan")
        if r["test_mrr"] > best_null:
            n_beat += 1
        mrrs.append(r["test_mrr"])
        lines.append(
            f"| {r['curie']} | {r['name']} | {b.get('n_test', '?')} | "
            f"{r['test_mrr']:.5f} | {r['test_hits10'] * 100:.2f}% | "
            f"{r['test_hits100'] * 100:.2f}% | {b.get('degree_filtered_mrr', float('nan')):.4f} | "
            f"{b.get('ppr_filtered_mrr', float('nan')):.4f} | {lift:.2f}x | "
            f"{100.0 * float(t.get('pct_le2', 0.0)):.1f}% |"
        )

    n = len(results)
    mean_mrr = sum(mrrs) / n if n else 0.0
    lines += [
        "",
        (
            f"Units evaluated: {n}/13 — {n_beat}/{n} beat the best same-graph null "
            f"(degree/PPR); mean test MRR {mean_mrr:.5f}."
        ),
        "",
        (
            "Reproduction: `python scripts/pod/run_wave1_benchmark.py --batch-size 32 "
            "--seed 1337` (results + manifest: "
            "`reports/m2b/wave1_benchmark_results.json` / `.manifest.json`; checkpoints "
            "gitignored, EBS-persisted). Nulls: `reports/m2b/lodo_baselines_keep.md`; "
            "negative controls: `reports/m2b/lodo_shuffles_keep.md` + GPU collapse runs "
            "(`--unit-set shuffle`)."
        ),
    ]

    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)
    text = "\n".join(lines) + "\n"
    cfg.report_path.write_text(text, encoding="utf-8")

    cfg.manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = {
        "script": "m2_wave1_report",
        "created_at": utcnow_iso(),
        "status": "done",
        "outputs": {str(cfg.report_path): hashlib.sha256(text.encode("utf-8")).hexdigest()},
    }
    cfg.manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    log.info("wave1_report_done", units=n, beat_nulls=n_beat, mean_mrr=round(mean_mrr, 5))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
