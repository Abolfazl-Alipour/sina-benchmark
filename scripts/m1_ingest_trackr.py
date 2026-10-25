#!/usr/bin/env python
"""m1_ingest_trackr — PrimeKG → canonical SINA tables (spec §5.1).

Stages (resumable via the manifest; --dry-run prints the plan only):
  1. characterize  — delimiter sniff, row counts, node-type/relation enums,
                     reciprocal audit → characterization report
  2. canonicalize  — nodes.parquet + edges.parquet with evidence_group_id
                     (spec §5.4: reciprocal twins + duplicates share one ID)
  3. validate      — counts vs published (±1%), FK integrity

Usage: python scripts/m1_ingest_trackr.py --config configs/m1_trackr_ingest.yaml
           [--dry-run] [--resume]
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import polars as pl
from pydantic import BaseModel

from sina.config import load_config
from sina.evidence_groups import add_evidence_group_ids, reciprocal_audit, symmetric_relations
from sina.logging import configure_logging, get_logger
from sina.manifests import Run

log = get_logger(__name__)


class IngestConfig(BaseModel):
    track: str
    input_csv: Path
    out_dir: Path
    interim_dir: Path
    report_path: Path
    manifest_path: Path
    published_counts: dict[str, int]
    tolerance: float = 0.01
    reciprocal_threshold: float = 0.99
    source_prefix_map: dict[str, str] = {}


def sniff_delimiter(path: Path) -> str:
    with path.open("r", encoding="utf-8", newline="") as fh:
        sample = fh.read(1 << 16)
    return csv.Sniffer().sniff(sample, delimiters=",\t;|").delimiter


def curie(source: str, local_id: str, prefix_map: dict[str, str]) -> str:
    prefix = prefix_map.get(source, source.lower())
    return f"{prefix}:{local_id}"


def load_kg(path: Path, delimiter: str, prefix_map: dict[str, str]) -> pl.DataFrame:
    df = pl.read_csv(
        path,
        separator=delimiter,
        schema_overrides={"x_id": pl.String, "y_id": pl.String},
    )
    return df.with_columns(
        src_id=pl.struct(["x_source", "x_id"]).map_elements(
            lambda s: curie(s["x_source"], s["x_id"], prefix_map), return_dtype=pl.String
        ),
        dst_id=pl.struct(["y_source", "y_id"]).map_elements(
            lambda s: curie(s["y_source"], s["y_id"], prefix_map), return_dtype=pl.String
        ),
    )


def characterize(df: pl.DataFrame, delimiter: str, cfg: IngestConfig) -> dict:
    x_types = sorted(df["x_type"].unique().to_list())
    y_types = sorted(df["y_type"].unique().to_list())
    relations = sorted(df["relation"].unique().to_list())
    audit = reciprocal_audit(df)
    symmetric = sorted(symmetric_relations(audit, cfg.reciprocal_threshold))
    return {
        "delimiter": repr(delimiter),
        "n_rows": df.height,
        "node_types": sorted(set(x_types) | set(y_types)),
        "n_relations": len(relations),
        "relations": relations,
        "reciprocal_audit": audit.to_dicts(),
        "symmetric_relations": symmetric,
        "reciprocal_threshold": cfg.reciprocal_threshold,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, IngestConfig)

    if args.dry_run:
        print(f"would ingest {cfg.input_csv} → {cfg.out_dir} (report: {cfg.report_path})")
        return 0

    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    cfg.interim_dir.mkdir(parents=True, exist_ok=True)
    nodes_path = cfg.out_dir / "nodes.parquet"
    edges_path = cfg.out_dir / "edges.parquet"
    char_path = cfg.interim_dir / "characterization.json"

    with Run(
        script="m1_ingest_trackr",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config), "track": cfg.track},
    ) as run:
        if args.resume and "characterize" in run.done_stages:
            char = json.loads(char_path.read_text(encoding="utf-8"))
        else:
            with run.stage("characterize", outputs=[char_path]):
                delimiter = sniff_delimiter(cfg.input_csv)
                df = load_kg(cfg.input_csv, delimiter, cfg.source_prefix_map)
                char = characterize(df, delimiter, cfg)
                char_path.write_text(json.dumps(char, indent=2), encoding="utf-8")

        if not (args.resume and "canonicalize" in run.done_stages):
            with run.stage("canonicalize", outputs=[nodes_path, edges_path]):
                df = load_kg(cfg.input_csv, sniff_delimiter(cfg.input_csv), cfg.source_prefix_map)
                edges = df.select(
                    "src_id",
                    "dst_id",
                    pl.col("relation"),
                    pl.col("display_relation"),
                )
                edges = add_evidence_group_ids(edges, char["symmetric_relations"])
                x = df.select(
                    node_id=pl.col("src_id"),
                    node_type=pl.col("x_type"),
                    node_name=pl.col("x_name"),
                    node_source=pl.col("x_source"),
                )
                y = df.select(
                    node_id=pl.col("dst_id"),
                    node_type=pl.col("y_type"),
                    node_name=pl.col("y_name"),
                    node_source=pl.col("y_source"),
                )
                nodes = pl.concat([x, y]).unique(subset=["node_id"], keep="first")
                edges.write_parquet(edges_path)
                nodes.write_parquet(nodes_path)

        with run.stage("validate", outputs=[cfg.report_path]):
            nodes = pl.read_parquet(nodes_path)
            edges = pl.read_parquet(edges_path)
            pub = cfg.published_counts
            node_dev = abs(nodes.height - pub["nodes"]) / pub["nodes"]
            edge_dev = abs(edges.height - pub["edges"]) / pub["edges"]
            if node_dev > cfg.tolerance or edge_dev > cfg.tolerance:
                raise ValueError(
                    f"count reconciliation failed: nodes {nodes.height} vs {pub['nodes']} "
                    f"({node_dev:.4%}), edges {edges.height} vs {pub['edges']} ({edge_dev:.4%})"
                )
            node_ids = set(nodes["node_id"].to_list())
            dangling = edges.filter(
                ~pl.col("src_id").is_in(node_ids) | ~pl.col("dst_id").is_in(node_ids)
            )
            if dangling.height:
                raise ValueError(f"FK integrity failed: {dangling.height} dangling edge endpoints")
            write_report(cfg, char, nodes.height, edges.height)

    log.info("ingest_done", track=cfg.track)
    return 0


def write_report(cfg: IngestConfig, char: dict, n_nodes: int, n_edges: int) -> None:
    audit_rows = "\n".join(
        f"| {r['relation']} | {r['n_edges']} | {r['n_with_reverse']} "
        f"| {r['reverse_fraction']:.4f} |"
        for r in char["reciprocal_audit"]
    )
    pub_nodes, pub_edges = cfg.published_counts["nodes"], cfg.published_counts["edges"]
    node_dev = abs(n_nodes - pub_nodes) / pub_nodes
    edge_dev = abs(n_edges - pub_edges) / pub_edges
    text = f"""<!-- track_scope: benchmark -->
# Track-R (PrimeKG) characterization — M1 (spec §5.1/§5.5)

**Status: DRAFT — reviewer-auditable evidence (Gate 1)**

Source: `data/raw/primekg/kg.csv` (Dataverse doi:10.7910/DVN/IXA7BM v2.1,
sha256 in `data/raw/REGISTRY.json`).

## Empirical facts (from the file, not from the paper)

- Delimiter sniffed: `{char["delimiter"]}`; rows: **{char["n_rows"]:,}**
- Node types ({len(char["node_types"])}): {", ".join(f"`{t}`" for t in char["node_types"])}
- Relations ({len(char["relations"])}): {", ".join(f"`{r}`" for r in char["relations"])}
- Symmetric relations (reverse-fraction ≥ {char["reciprocal_threshold"]}):
  {", ".join(f"`{r}`" for r in char["symmetric_relations"]) or "none"}

## Reciprocal audit

| relation | edges | with reverse twin | fraction |
|---|---|---|---|
{audit_rows}

## Canonical tables

- `data/processed/m1/trackr/nodes.parquet` — **{n_nodes:,}** nodes
  (`node_id` CURIEs: `ncbigene:` / `mondo:` / `drugbank:` / …)
- `data/processed/m1/trackr/edges.parquet` — **{n_edges:,}** edges with
  `evidence_group_id` (reciprocal twins + duplicates share one ID, spec §5.4)

## Count reconciliation (published: {pub_nodes:,} nodes / {pub_edges:,}
directed edges; tolerance ±{cfg.tolerance:.0%})

- nodes: {n_nodes:,} (deviation {node_dev:.4%})
- edges: {n_edges:,} (deviation {edge_dev:.4%})
- FK integrity: all edge endpoints present in `nodes.parquet` — **PASS**
"""
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)
    cfg.report_path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
