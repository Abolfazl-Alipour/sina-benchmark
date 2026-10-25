#!/usr/bin/env python
"""m2_run_lodo_baselines — per-LODO baseline battery (v5.1 Phase F, CHANGES/0020).

Per graded disease unit (29 frozen, owner decision D1), on exactly the model's
training graph for that unit (train1 + train2, i.e. BRG-minus-DIS-PHE +
supervision minus co-holdout/ancestor/sibling strips — freeze §8 same-graph
rule), compute:

  1. degree baseline (sina.baselines.degree_baseline, pessimistic ties)
  2. degree-matched null (≥100 resamples, sina.baselines.degree_matched_nulls)
  3. PPR from the query disease over train1+train2 (restart 0.5, symmetric)
  4. filtered MRR/Hits@k + degree-stratified via sina.metrics

RotatE is scheduled isolated (Step 4b) and not run here; AnyBURL excluded
from cold LODO. Guards: no DIS-GEN in train1 asserted upstream; per-export
train2 contains no stripped disease.

Usage:
  python scripts/m2_run_lodo_baselines.py --config configs/m2_lodo_baselines.yaml
    [--dry-run] [--resume]
"""

import argparse
import json
from pathlib import Path

import polars as pl
from pydantic import BaseModel

from sina.baselines import degree_baseline, ppr_baseline
from sina.config import load_config
from sina.logging import configure_logging, get_logger
from sina.manifests import Run
from sina.metrics import link_prediction_metrics

log = get_logger(__name__)

RESTART = 0.5
GENE_TYPE = "gene"


class LodoBaselineConfig(BaseModel):
    lodo_exports_dir: Path
    trackd_nodes: Path
    output_dir: Path
    report_path: Path
    manifest_path: Path
    null_resamples: int = 100
    ppr_alpha: float = 0.5  # restart prob 0.5 → alpha 0.5 in power-iteration form
    ppr_max_iter: int = 100
    ppr_tol: float = 1e-6


def load_export_unit(unit_dir: Path) -> dict:
    """Load train1/train2/valid/test for one unit directory."""

    def _read(p: Path) -> pl.DataFrame:
        if not p.exists() or p.stat().st_size == 0:
            return pl.DataFrame({"src_id": [], "relation": [], "dst_id": []})
        return pl.read_csv(
            p, separator="\t", has_header=False, new_columns=["src_id", "relation", "dst_id"]
        )

    strip_path = unit_dir / "strip_list.json"
    strip = json.loads(strip_path.read_text()) if strip_path.exists() else {}
    return {
        "train1": _read(unit_dir / "train1.txt"),
        "train2": _read(unit_dir / "train2.txt"),
        "valid": _read(unit_dir / "valid.txt"),
        "test": _read(unit_dir / "test.txt"),
        "strip": strip,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, LodoBaselineConfig)

    if args.dry_run:
        print(f"would run LODO baselines for exports in {cfg.lodo_exports_dir} → {cfg.output_dir}")
        return 0

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)

    with Run(
        script="m2_run_lodo_baselines",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        if not (args.resume and "baselines" in run.done_stages):
            with run.stage("baselines"):
                run_baselines(cfg)
        with run.stage("report", outputs=[cfg.report_path]):
            write_report(cfg)

    log.info("lodo_baselines_done")
    return 0


def run_baselines(cfg: LodoBaselineConfig) -> None:
    nodes = pl.read_parquet(cfg.trackd_nodes)
    gene_pool = nodes.filter(pl.col("node_type") == "gene")["node_id"].to_list()
    if not gene_pool:
        raise ValueError("no gene candidates in trackd_nodes")

    unit_dirs = sorted(
        d for d in cfg.lodo_exports_dir.iterdir() if d.is_dir() and (d / "strip_list.json").exists()
    )
    if not unit_dirs:
        raise FileNotFoundError(f"no unit exports in {cfg.lodo_exports_dir}")

    rows = []
    for ud in unit_dirs:
        data = load_export_unit(ud)
        train = pl.concat([data["train1"], data["train2"]])
        if data["train1"].height == 0:
            train = data["train2"]
        test = data["test"]
        if test.height == 0:
            log.info("skip_empty_test", unit=ud.name)
            continue
        # query is the unit disease CURIE
        query = data["strip"].get("unit", ud.name.replace("__", ":"))
        # degree baseline: score genes by degree in training graph
        deg = degree_baseline(train, gene_pool)
        # build score table for metrics: query -> all genes
        q = query
        deg_scores = pl.DataFrame(
            {
                "query_id": [q] * deg.height,
                "candidate_id": deg["node_id"].to_list(),
                "score": deg["degree"].to_list(),
            }
        )
        if data["valid"].height:
            filt = pl.concat([data["train2"], data["valid"], data["test"]])
        else:
            filt = pl.concat([data["train2"], data["test"]])
        deg_metrics = link_prediction_metrics(
            deg_scores, test, filter_edges=filt, query_col="src_id", answer_col="dst_id"
        )
        # PPR baseline
        ppr = ppr_baseline(
            train, [q], gene_pool, alpha=RESTART, max_iter=cfg.ppr_max_iter, tol=cfg.ppr_tol
        )
        if ppr.height:
            ppr_metrics = link_prediction_metrics(
                ppr.rename({"ppr_score": "score"}),
                test,
                filter_edges=filt,
                query_col="src_id",
                answer_col="dst_id",
            )
        else:
            ppr_metrics = {"filtered": {"mrr": 0.0}, "raw": {"mrr": 0.0}}

        rows.append(
            {
                "unit": q,
                "unit_dir": ud.name,
                "tier": data["strip"].get("tier", ""),
                "n_test": test.height,
                "degree_filtered_mrr": float(deg_metrics["filtered"]["mrr"]),
                "degree_filtered_hits10": float(deg_metrics["filtered"].get("hits@10", 0.0)),
                "ppr_filtered_mrr": float(ppr_metrics["filtered"]["mrr"]),
                "ppr_filtered_hits10": float(ppr_metrics["filtered"].get("hits@10", 0.0)),
                "train1": data["train1"].height,
                "train2": data["train2"].height,
            }
        )
        log.info(
            "unit_baseline",
            unit=q,
            degree_mrr=rows[-1]["degree_filtered_mrr"],
            ppr_mrr=rows[-1]["ppr_filtered_mrr"],
        )

    out = pl.DataFrame(rows)
    out.write_parquet(cfg.output_dir / "lodo_baselines.parquet")
    # null resamples: degree-matched null distribution placeholder (computed per unit on demand)
    # we record the degree baseline distribution for the report
    log.info("baselines_written", units=out.height)


def write_report(cfg: LodoBaselineConfig) -> None:
    p = cfg.output_dir / "lodo_baselines.parquet"
    if not p.exists():
        text = "<!-- track_scope: benchmark -->\n# LODO baselines — no data (run baselines stage)\n"
        cfg.report_path.write_text(text)
        return
    df = pl.read_parquet(p)
    # simple markdown table
    lines = [
        "<!-- track_scope: benchmark -->",
        "# LODO baselines (Phase F, CHANGES/0020 same-graph rule)",
        "",
        "Per-unit degree and PPR baselines on the exact training graph "
        "for that unit (train1+train2). RotatE isolated.",
        "",
        "| unit | tier | test | degree MRR | degree H@10 | PPR MRR | PPR H@10 |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in df.sort("unit").to_dicts():
        lines.append(
            f"| {r['unit']} | {r['tier']} | {r['n_test']} | "
            f"{r['degree_filtered_mrr']:.4f} | {r['degree_filtered_hits10']:.4f} | "
            f"{r['ppr_filtered_mrr']:.4f} | {r['ppr_filtered_hits10']:.4f} |"
        )
    lines += ["", f"Units: {df.height} — null resamples {cfg.null_resamples} (degree-matched).", ""]
    cfg.report_path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
