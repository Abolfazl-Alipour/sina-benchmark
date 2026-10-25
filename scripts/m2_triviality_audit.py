#!/usr/bin/env python
"""m2_triviality_audit — A3 triviality audit + degenerate handling (CHANGES/0021).

Per LODO unit, compute shortest-path distance in train1 (BRG) from the held-out
disease to each test answer gene. A unit is degenerate (excluded from primary
gate, reported separately) if EITHER:
  - any baseline attains MRR >=0.99 from query isolation (PPR 1.0 observed), OR
  - >10% of its test answers are reachable in <=2 hops in train1.

Outputs:
  - data/processed/m2/lodo_audits/triviality.parquet (per-unit distribution)
  - updates data/processed/m2/lodo_baselines/lodo_baselines.parquet with
    valid_baseline flag + degenerate flag
  - reports/m2b/triviality_audit.md

Usage:
  python scripts/m2_triviality_audit.py --config configs/m2_triviality.yaml
    [--dry-run] [--resume]
"""

import argparse
import json
from collections import deque
from pathlib import Path

import polars as pl
from pydantic import BaseModel

from sina.config import load_config
from sina.logging import configure_logging, get_logger
from sina.manifests import Run

log = get_logger(__name__)


class TrivialityConfig(BaseModel):
    lodo_exports_dir: Path
    lodo_baselines: Path
    output_dir: Path
    report_path: Path
    manifest_path: Path
    train1_edges_label: str = "1,420,095"


def bfs_distances(train1: pl.DataFrame, source: str, targets: set[str]) -> dict[str, int | None]:
    """BFS from source on undirected train1 graph (src_id/dst_id). Returns dist or None."""
    # Build adjacency
    adj: dict[str, set[str]] = {}
    for s, d in zip(train1["src_id"].to_list(), train1["dst_id"].to_list(), strict=False):
        adj.setdefault(s, set()).add(d)
        adj.setdefault(d, set()).add(s)
    dist: dict[str, int] = {source: 0}
    q = deque([source])
    visited = {source}
    while q:
        u = q.popleft()
        for v in adj.get(u, []):
            if v not in visited:
                visited.add(v)
                dist[v] = dist[u] + 1
                q.append(v)
    return {t: dist.get(t) for t in targets}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, TrivialityConfig)

    if args.dry_run:
        print(f"would audit {cfg.lodo_exports_dir} → {cfg.report_path}")
        return 0

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)

    with Run(
        script="m2_triviality_audit",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        if not (args.resume and "audit" in run.done_stages):
            with run.stage("audit", outputs=[cfg.output_dir / "triviality.parquet"]):
                audit(cfg)
        with run.stage("report", outputs=[cfg.report_path]):
            write_report(cfg)

    log.info("triviality_done")
    return 0


def audit(cfg: TrivialityConfig) -> None:
    baselines = pl.read_parquet(cfg.lodo_baselines) if cfg.lodo_baselines.exists() else None
    # Map unit -> degenerate due to baseline MRR >=0.99
    baseline_degen = set()
    if baselines is not None:
        for r in baselines.to_dicts():
            if (
                max(float(r.get("ppr_filtered_mrr", 0)), float(r.get("degree_filtered_mrr", 0)))
                >= 0.99
            ):
                baseline_degen.add(r["unit"])

    unit_dirs = sorted(
        [
            d
            for d in cfg.lodo_exports_dir.iterdir()
            if d.is_dir() and (d / "strip_list.json").exists()
        ]
    )
    rows = []
    for ud in unit_dirs:
        strip = json.loads((ud / "strip_list.json").read_text())
        unit = strip["unit"]
        train1 = pl.read_csv(
            ud / "train1.txt",
            separator="\t",
            has_header=False,
            new_columns=["src_id", "relation", "dst_id"],
        )
        test = pl.read_csv(
            ud / "test.txt",
            separator="\t",
            has_header=False,
            new_columns=["src_id", "relation", "dst_id"],
        )
        answers = set(test["dst_id"].to_list())
        if not answers:
            rows.append(
                {
                    "unit": unit,
                    "tier": strip.get("tier", ""),
                    "n_test": 0,
                    "pct_le2": 0.0,
                    "dist_1": 0,
                    "dist_2": 0,
                    "unreachable": 0,
                    "degenerate_baseline": unit in baseline_degen,
                    "degenerate_trivial": False,
                    "degenerate": unit in baseline_degen,
                }
            )
            continue
        dists = bfs_distances(train1, unit, answers)
        # Count
        c1 = sum(1 for v in dists.values() if v == 1)
        c2 = sum(1 for v in dists.values() if v == 2)
        unreachable = sum(1 for v in dists.values() if v is None)
        pct_le2 = (c1 + c2) / len(answers) if answers else 0
        degen_trivial = pct_le2 > 0.10
        degen = (unit in baseline_degen) or degen_trivial
        rows.append(
            {
                "unit": unit,
                "tier": strip.get("tier", ""),
                "n_test": len(answers),
                "dist_1": c1,
                "dist_2": c2,
                "dist_3": sum(1 for v in dists.values() if v == 3),
                "dist_4": sum(1 for v in dists.values() if v == 4),
                "dist_5": sum(1 for v in dists.values() if v == 5),
                "unreachable": unreachable,
                "pct_le2": float(pct_le2),
                "degenerate_baseline": unit in baseline_degen,
                "degenerate_trivial": bool(degen_trivial),
                "degenerate": bool(degen),
            }
        )
        log.info("unit_triviality", unit=unit, pct_le2=pct_le2, degen=degen)

    df = pl.DataFrame(rows)
    df.write_parquet(cfg.output_dir / "triviality.parquet")
    # Update baselines parquet with valid_baseline flag
    if baselines is not None:
        # valid_baseline = not degenerate
        degen_map = {r["unit"]: r["degenerate"] for r in rows}
        baselines = baselines.with_columns(
            pl.col("unit")
            .map_elements(lambda u: not degen_map.get(u, False), return_dtype=pl.Boolean)
            .alias("valid_baseline")
        )
        # Also add degenerate column
        baselines = baselines.with_columns(
            pl.col("unit")
            .map_elements(lambda u: degen_map.get(u, False), return_dtype=pl.Boolean)
            .alias("degenerate")
        )
        baselines.write_parquet(cfg.lodo_baselines)
        log.info("baselines_updated_with_degenerate", units=df.height)


def write_report(cfg: TrivialityConfig) -> None:
    df = pl.read_parquet(cfg.output_dir / "triviality.parquet")
    lines = [
        "<!-- track_scope: benchmark -->",
        "# Triviality audit (A3, CHANGES/0021)",
        "",
        (
            f"Per-unit shortest-path distance in `train1` (BRG, {cfg.train1_edges_label} edges) "
            "from held-out disease to test answers. Degenerate if baseline "
            "MRR>=0.99 OR >10% <=2 hops."
        ),
        "",
        "| unit | tier | test | ≤1 | ≤2 | ≤2% | 3 | 4 | 5 | unreachable | degenerate |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in df.sort("pct_le2", descending=True).to_dicts():
        lines.append(
            f"| {r['unit']} | {r['tier']} | {r['n_test']} | {r['dist_1']} | "
            f"{r['dist_2']} | {r['pct_le2']:.3f} | {r['dist_3']} | "
            f"{r['dist_4']} | {r['dist_5']} | {r['unreachable']} | "
            f"{r['degenerate']} |"
        )
    degen = df.filter(pl.col("degenerate")).height
    lines += [
        "",
        (
            f"Degenerate units: {degen} / {df.height} "
            f"(excluded from primary gate, reported separately)."
        ),
        "",
    ]
    # Tiny smoke test note
    lines += ["Tiny smoke test: ~10-node graph baseline MRR <1 asserted in unit tests.", ""]
    cfg.report_path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
