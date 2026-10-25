#!/usr/bin/env python
"""m1_ingest_trackd — OptimusKG → canonical SINA tables (spec §5.2).

Stages (resumable; --dry-run prints the plan):
  1. canonicalize_nodes — 10 per-type files → nodes.parquet (CURIE node_id)
  2. canonicalize_edges — 27 per-type files → edges.parquet with
     evidence_group_id (§5.2.3: order-insensitive only when undirected=True,
     key includes the canonical source-set) and properties struct preserved
     (superset struct, Polars-native — no JSON-string parsing)
  3. extract_properties — m1 gene property table (§5.2.4)
  4. validate — EXACT count reconciliation (§5.5), FK integrity, DIS-GEN
     temporal coverage (Appendix D.5) → characterization report

Usage: python scripts/m1_ingest_trackd.py --config configs/m1_trackd_ingest.yaml
           [--dry-run] [--resume]
"""

import argparse
import shutil
import sys
from pathlib import Path

import polars as pl
from pydantic import BaseModel

from sina.config import load_config
from sina.evidence_groups import add_evidence_group_ids
from sina.logging import configure_logging, get_logger
from sina.manifests import Run

log = get_logger(__name__)

# DIS-GEN/PHE-GEN typed property fields lifted to top-level columns (§5.2.4)
DISGEN_FIELDS = [
    "evidence_score",
    "evidence_count",
    "disease_specificity_index",
    "disease_pleiotropy_index",
    "evidence_index",
    "disgenet_score",
    "year_initial",
    "year_final",
    "number_of_pmids",
    "number_of_snps",
]


class IngestConfig(BaseModel):
    track: str
    raw_dir: Path
    out_dir: Path
    interim_dir: Path
    report_path: Path
    manifest_path: Path
    published_counts: dict[str, int]
    tolerance: float = 0.0
    source_prefix_map: dict[str, str] = {}
    temporal_coverage_threshold: float = 0.5


def curie_expr(col: str, prefix_map: dict[str, str]) -> pl.Expr:
    """Vectorized CURIE: 'MONDO_0004975' → 'mondo:0004975', 'ENSG…' → 'ensembl:…'."""
    raw = pl.col(col).str.extract(r"^([A-Za-z]+)")
    prefix = raw.replace_strict(prefix_map, default=raw.str.to_lowercase(), return_dtype=pl.String)
    local = pl.col(col).str.replace(r"^[A-Za-z]+[_:]", "")
    return pl.concat_str([prefix, pl.lit(":"), local])


def extend_struct(df: pl.DataFrame, col: str, all_fields: dict[str, pl.DataType]) -> pl.Expr:
    """Rebuild a struct column with the superset of fields (nulls for missing).

    Present fields are cast to the canonical (first-seen) dtype so per-file
    parts share one schema and can be streamed into a single parquet.
    """
    present = {f.name: f.dtype for f in df.schema[col].fields}
    return pl.struct(
        [
            pl.col(col).struct.field(name).cast(dtype)
            if name in present
            else pl.lit(None, dtype=dtype).alias(name)
            for name, dtype in all_fields.items()
        ]
    ).alias(col)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, IngestConfig)
    if args.dry_run:
        print(f"would ingest {cfg.raw_dir} → {cfg.out_dir} (report: {cfg.report_path})")
        return 0

    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    cfg.interim_dir.mkdir(parents=True, exist_ok=True)
    nodes_path = cfg.out_dir / "nodes.parquet"
    edges_path = cfg.out_dir / "edges.parquet"
    props_path = cfg.out_dir.parent / "gene_properties.parquet"

    with Run(
        script="m1_ingest_trackd",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config), "track": cfg.track},
    ) as run:
        if not (args.resume and "canonicalize_nodes" in run.done_stages):
            with run.stage("canonicalize_nodes", outputs=[nodes_path]):
                frames = []
                for f in sorted(cfg.raw_dir.glob("nodes/*.parquet")):
                    node_type = f.stem
                    df = pl.read_parquet(f)
                    frames.append(
                        df.select(
                            node_id=curie_expr("id", cfg.source_prefix_map),
                            node_type=pl.lit(node_type),
                            # OptimusKG `label` is the node-TYPE label ("DIS");
                            # the human-readable name is properties.name
                            node_name=pl.coalesce(
                                pl.col("properties").struct.field("name"), pl.col("label")
                            ),
                            native_id=pl.col("id"),
                            properties=pl.col("properties"),
                        )
                    )
                nodes = pl.concat(frames, how="diagonal_relaxed")
                nodes = nodes.with_columns(
                    mondo_id=pl.when(pl.col("native_id").str.starts_with("MONDO_"))
                    .then(pl.col("native_id").str.replace("MONDO_", "mondo:"))
                    .otherwise(None)
                )
                nodes.write_parquet(nodes_path)

        if not (args.resume and "canonicalize_edges" in run.done_stages):
            with run.stage("canonicalize_edges", outputs=[edges_path]):
                files = sorted(cfg.raw_dir.glob("edges/*.parquet"))
                parts_dir = cfg.interim_dir / "edges_parts"
                parts_dir.mkdir(parents=True, exist_ok=True)
                # Schema survey (no data read): superset struct fields + canonical
                # DIS-GEN column dtypes so every part shares one schema.
                all_fields: dict[str, pl.DataType] = {}
                schemas: dict[str, dict] = {}
                for f in files:
                    schema = pl.read_parquet_schema(f)
                    schemas[f.stem] = schema
                    for fld in schema["properties"].fields:
                        all_fields.setdefault(fld.name, fld.dtype)
                disgen_dtypes = {
                    fld.name: fld.dtype for fld in schemas["disease_gene"]["properties"].fields
                }
                # Per-file processing keeps peak RSS far below the §0.3.3 budget:
                # transform one file, assign group ids, write a part, release.
                for f in files:
                    df = pl.read_parquet(f)
                    present = {x.name for x in df.schema["properties"].fields}
                    out = df.select(
                        src_id=curie_expr("from", cfg.source_prefix_map),
                        dst_id=curie_expr("to", cfg.source_prefix_map),
                        label=pl.col("label"),
                        relation=pl.col("relation"),
                        undirected=pl.col("undirected"),
                        sources_direct=pl.col("properties")
                        .struct.field("sources")
                        .struct.field("direct"),
                        sources_indirect=pl.col("properties")
                        .struct.field("sources")
                        .struct.field("indirect"),
                        properties=extend_struct(df, "properties", all_fields),
                        **{
                            name: pl.col("properties")
                            .struct.field(name)
                            .cast(disgen_dtypes[name])
                            .alias(name)
                            if name in present
                            else pl.lit(None, dtype=disgen_dtypes[name]).alias(name)
                            for name in DISGEN_FIELDS
                        },
                    )
                    # §5.2.3: group key includes the canonical source-set
                    out = out.with_columns(
                        _srcset=pl.struct(
                            direct=pl.col("sources_direct"),
                            indirect=pl.col("sources_indirect"),
                        )
                    )
                    out = add_evidence_group_ids(
                        out, undirected_col="undirected", sources_col="_srcset"
                    ).drop("_srcset")
                    out.write_parquet(parts_dir / f.name)
                    del df, out
                pl.scan_parquet(parts_dir / "*.parquet").sink_parquet(edges_path)
                shutil.rmtree(parts_dir)

        if not (args.resume and "extract_properties" in run.done_stages):
            with run.stage("extract_properties", outputs=[props_path]):
                genes = pl.read_parquet(cfg.raw_dir / "nodes/gene.parquet")
                props = genes.select(
                    gene_id=curie_expr("id", cfg.source_prefix_map),
                    symbol=pl.col("properties").struct.field("symbol"),
                    name=pl.col("properties").struct.field("name"),
                    biotype=pl.col("properties").struct.field("biotype"),
                    tractability=pl.col("properties").struct.field("tractability"),
                    safety_liabilities=pl.col("properties").struct.field("safety_liabilities"),
                    target_class=pl.col("properties").struct.field("target_class"),
                    constraint_scores=pl.col("properties").struct.field("constraint_scores"),
                    homologues=pl.col("properties").struct.field("homologues"),
                    chemical_probes=pl.col("properties").struct.field("chemical_probes"),
                )
                props.write_parquet(props_path)

        with run.stage("validate", outputs=[cfg.report_path]):
            nodes = pl.read_parquet(nodes_path)
            edges = pl.read_parquet(edges_path)
            pub = cfg.published_counts
            if nodes.height != pub["nodes"] or edges.height != pub["edges"]:
                raise ValueError(
                    f"EXACT reconciliation failed (§5.5): nodes {nodes.height} vs "
                    f"{pub['nodes']}, edges {edges.height} vs {pub['edges']}"
                )
            node_ids = set(nodes["node_id"].to_list())
            dangling = edges.filter(
                ~pl.col("src_id").is_in(node_ids) | ~pl.col("dst_id").is_in(node_ids)
            )
            if dangling.height:
                raise ValueError(f"FK integrity failed: {dangling.height} dangling endpoints")
            disgen = edges.filter(pl.col("label") == "DIS-GEN")
            coverage = temporal_coverage(disgen)
            write_report(cfg, nodes, edges, coverage)
            if coverage["both_years_fraction"] < cfg.temporal_coverage_threshold:
                log.warning(
                    "temporal_coverage_below_threshold",
                    fraction=coverage["both_years_fraction"],
                    note="temporal split downgrades to secondary (Appendix D.5)",
                )

    log.info("ingest_done", track=cfg.track)
    return 0


def temporal_coverage(disgen: pl.DataFrame) -> dict:
    """DIS-GEN year-field coverage (Appendix D assumption 5)."""
    n = disgen.height
    valid = lambda c: pl.col(c).is_not_null() & (pl.col(c) != "")  # noqa: E731
    both = disgen.filter(valid("year_initial") & valid("year_final")).height
    either = disgen.filter(valid("year_initial") | valid("year_final")).height
    return {
        "disgen_edges": n,
        "both_years": both,
        "either_year": either,
        "both_years_fraction": both / n if n else 0.0,
        "either_year_fraction": either / n if n else 0.0,
    }


def write_report(cfg: IngestConfig, nodes: pl.DataFrame, edges: pl.DataFrame, cov: dict) -> None:
    type_rows = "\n".join(
        f"| {r['node_type']} | {r['n']:,} |"
        for r in nodes.group_by("node_type").agg(n=pl.len()).sort("node_type").to_dicts()
    )
    label_rows = "\n".join(
        f"| {r['label']} | {r['n']:,} |"
        for r in edges.group_by("label").agg(n=pl.len()).sort("label").to_dicts()
    )
    n_mondo = nodes.filter(pl.col("mondo_id").is_not_null()).height
    n_groups = edges["evidence_group_id"].n_unique()
    pub = cfg.published_counts
    both_pct = f"{cov['both_years_fraction']:.1%}"
    either_pct = f"{cov['either_year_fraction']:.1%}"
    temporal_verdict = (
        "**PRIMARY**"
        if cov["both_years_fraction"] >= cfg.temporal_coverage_threshold
        else "**DOWNGRADED to secondary**"
    )
    text = f"""<!-- track_scope: benchmark -->
# Track-D (OptimusKG) characterization — M1 (spec §5.2/§5.5)

**Status: DRAFT — reviewer-auditable evidence (Gate 1)**

Source: stratified per-type Parquet, Dataverse doi:10.7910/DVN/IYNGEV v2.0
(37 files, sha256 each in `data/raw/REGISTRY.json`; downloaded via pinned
`optimuskg` 1.0.0 client).

## Canonical tables

- `data/processed/m1/trackd/nodes.parquet` — **{nodes.height:,}** nodes,
  10 types, CURIE `node_id`; MONDO-native diseases: {n_mondo:,}
  (crosswalk gap documented in `CHANGES/0008`)
- `data/processed/m1/trackd/edges.parquet` — **{edges.height:,}** edges,
  `label`/`relation`/`undirected`, provenance `sources_direct`/`sources_indirect`,
  properties superset struct preserved; **{n_groups:,}** evidence groups
  (§5.2.3: order-insensitive only when `undirected=True`)
- `data/processed/m1/gene_properties.parquet` — per-gene property table
  (§5.2.4): tractability, safety_liabilities, target_class, constraint_scores,
  homologues, chemical_probes. NOTE: structure/pocket/druggable flags are NOT
  present in OptimusKG v2.0 — see `CHANGES/0008`.

## Node types

| node_type | nodes |
|---|---|
{type_rows}

## Edge labels

| label | edges |
|---|---|
{label_rows}

## Exact count reconciliation (§5.5 — must match exactly)

- nodes: {nodes.height:,} vs published {pub["nodes"]:,} — **EXACT**
- edges: {edges.height:,} vs published {pub["edges"]:,} — **EXACT**
- FK integrity: all edge endpoints present in nodes.parquet — **PASS**

## DIS-GEN temporal coverage (Appendix D assumption 5)

- DIS-GEN edges: {cov["disgen_edges"]:,}
- with BOTH `year_initial` and `year_final`: {cov["both_years"]:,} ({both_pct})
- with either year field: {cov["either_year"]:,} ({either_pct})
- threshold for a primary temporal split: {cfg.temporal_coverage_threshold:.0%}
  → temporal split is {temporal_verdict}
"""
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)
    cfg.report_path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
