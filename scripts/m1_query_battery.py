#!/usr/bin/env python
"""m1_query_battery — DuckDB views + identical query battery per track (§5.2.5).

Registers one view per canonical table (``<track>_nodes`` / ``<track>_edges``)
and runs the same battery against each track: counts, relation distribution,
degree stats, FK integrity (anti-join), evidence-group multiplicity. Results go
to ``reports/m1/query_battery.md`` (claim-ceiling-linted).

Usage: python scripts/m1_query_battery.py --config configs/m1_query_battery.yaml
           [--dry-run] [--resume]
"""

import argparse
import sys
from pathlib import Path

import duckdb
from pydantic import BaseModel

from sina.config import load_config
from sina.logging import configure_logging, get_logger
from sina.manifests import Run

log = get_logger(__name__)


class TrackEntry(BaseModel):
    label: str
    nodes: Path
    edges: Path


class BatteryConfig(BaseModel):
    report_path: Path
    manifest_path: Path
    top_k: int = 10
    tracks: dict[str, TrackEntry]


def run_battery(con: duckdb.DuckDBPyConnection, track: str, top_k: int) -> dict:
    n = f"{track}_nodes"
    e = f"{track}_edges"
    battery = {
        "node_count": f"SELECT COUNT(*) AS v FROM {n}",
        "edge_count": f"SELECT COUNT(*) AS v FROM {e}",
        "relation_count": f"SELECT COUNT(DISTINCT relation) AS v FROM {e}",
        "node_type_count": f"SELECT COUNT(DISTINCT node_type) AS v FROM {n}",
        "fk_dangling": (
            f"SELECT COUNT(*) AS v FROM {e} "
            f"WHERE src_id NOT IN (SELECT node_id FROM {n}) "
            f"   OR dst_id NOT IN (SELECT node_id FROM {n})"
        ),
        "group_count": f"SELECT COUNT(DISTINCT evidence_group_id) AS v FROM {e}",
        "groups_with_twin": (
            f"SELECT COUNT(*) AS v FROM (SELECT evidence_group_id FROM {e} "
            f"GROUP BY evidence_group_id HAVING COUNT(*) = 2)"
        ),
        "top_relations": (
            f"SELECT relation, COUNT(*) AS n FROM {e} "
            f"GROUP BY relation ORDER BY n DESC LIMIT {top_k}"
        ),
        "top_degree_nodes": (
            f"SELECT deg.node_id, no.node_name, deg.degree FROM ("
            f"  SELECT node_id, SUM(d) AS degree FROM ("
            f"    SELECT src_id AS node_id, COUNT(*) AS d FROM {e} GROUP BY src_id"
            f"    UNION ALL"
            f"    SELECT dst_id AS node_id, COUNT(*) AS d FROM {e} GROUP BY dst_id"
            f"  ) GROUP BY node_id ORDER BY degree DESC LIMIT {top_k}"
            f") deg LEFT JOIN {n} no USING (node_id)"
        ),
    }
    results = {}
    for name, sql in battery.items():
        rows = con.execute(sql).fetchall()
        results[name] = rows
        log.info("battery_query", track=track, query=name)
    return results


def fmt_scalar(rows: list) -> str:
    return f"{rows[0][0]:,}"


def fmt_table(rows: list, headers: tuple[str, ...]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, BatteryConfig)
    if args.dry_run:
        print(f"would run battery over tracks: {', '.join(cfg.tracks)} → {cfg.report_path}")
        return 0

    con = duckdb.connect()
    sections = []
    with Run(
        script="m1_query_battery",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        for track, entry in cfg.tracks.items():
            stage = f"battery_{track}"
            if args.resume and stage in run.done_stages:
                log.info("stage_skipped", stage=stage)
                continue
            with run.stage(stage, outputs=[cfg.report_path]):
                con.execute(
                    f"CREATE OR REPLACE VIEW {track}_nodes AS "
                    f"SELECT * FROM read_parquet('{entry.nodes}')"
                )
                con.execute(
                    f"CREATE OR REPLACE VIEW {track}_edges AS "
                    f"SELECT * FROM read_parquet('{entry.edges}')"
                )
                r = run_battery(con, track, cfg.top_k)
                sections.append(f"""## {entry.label}

- nodes: **{fmt_scalar(r["node_count"])}** ({fmt_scalar(r["node_type_count"])} types)
- edges: **{fmt_scalar(r["edge_count"])}** ({fmt_scalar(r["relation_count"])} relations)
- evidence groups: **{fmt_scalar(r["group_count"])}**
  (groups with exactly 2 rows: {fmt_scalar(r["groups_with_twin"])} —
  on Track-R these are reciprocal twins, because kg.csv stores both directions
  of every edge; on Track-D undirected edges are stored once (§5.2.3), so
  0 is the expected value)
- FK integrity (dangling endpoints): **{fmt_scalar(r["fk_dangling"])}**

Top relations:

{fmt_table(r["top_relations"], ("relation", "edges"))}

Top-degree nodes:

{fmt_table(r["top_degree_nodes"], ("node_id", "name", "degree"))}
""")

    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)
    cfg.report_path.write_text(
        "<!-- track_scope: benchmark -->\n"
        "# M1 DuckDB query battery (spec §5.2.5)\n\n"
        "**Status: DRAFT — reviewer-auditable evidence (Gate 1)**\n\n"
        "Identical battery over DuckDB views of each track's canonical tables.\n\n"
        + "\n".join(sections),
        encoding="utf-8",
    )
    log.info("battery_done", report=str(cfg.report_path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
