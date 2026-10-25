#!/usr/bin/env python
"""m2_build_splits — build leakage-safe splits for Phase 2 (spec §6.2).

Stages (resumable; --dry-run prints the plan):
  1. load — canonical edges, disease families, TxGNN area files
  2. random_r / random_d — random edge split on evidence_group_id
  3. family_r / family_d — disease-family zero-shot holdout (AD/PD families)
  4. txgnn_r — TxGNN-protocol neurodegenerative area split on Track-R
  5. temporal_d — temporal split on Track-D DIS-GEN (secondary; Appendix D.5)
  6. write — parquet partitions + metadata JSON per split
  7. report — counts + leakage-guard status → reports/m2a/splits_report.md

Usage: python scripts/m2_build_splits.py --config configs/m2_splits.yaml
           [--dry-run] [--resume]
"""

import argparse
import gc
import json
import re
import sys
from pathlib import Path

import polars as pl
import yaml
from pydantic import BaseModel

from sina.config import load_config
from sina.logging import configure_logging, get_logger
from sina.manifests import Run
from sina.splits import (
    GROUP_COL,
    SplitLeakageError,
    assert_split_disjoint,
    family_holdout_split,
    group_assignment,
    random_edge_split,
    txgnn_area_split,
)

log = get_logger(__name__)


class TrackConfig(BaseModel):
    edges: Path


class TrackRConfig(TrackConfig):
    relation_col: str = "relation"
    family_eval_relations: list[str]


class TrackDConfig(TrackConfig):
    label_col: str = "label"
    family_eval_label: str = "DIS-GEN"
    temporal_label: str = "DIS-GEN"


class AreaConfig(BaseModel):
    track: str
    file: str
    indication_relations: list[str]


class RandomConfig(BaseModel):
    fractions: dict[str, float]
    seed: int


class FamilyConfig(BaseModel):
    valid_fraction: float
    seed: int


class TxGnnConfig(BaseModel):
    neighborhood_keep_fraction: float
    seed: int


class TemporalConfig(BaseModel):
    cutoff_year: int


class SplitsConfig(BaseModel):
    trackr: TrackRConfig
    trackd: TrackDConfig
    families_config: Path
    txgnn_dir: Path
    areas: dict[str, AreaConfig]
    random: RandomConfig
    family: FamilyConfig
    txgnn: TxGnnConfig
    temporal: TemporalConfig
    output_dir: Path
    report_path: Path
    manifest_path: Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, SplitsConfig)

    if args.dry_run:
        print(
            f"would build splits from {cfg.trackr.edges} and {cfg.trackd.edges} "
            f"→ {cfg.output_dir} (report: {cfg.report_path})"
        )
        return 0

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)

    with Run(
        script="m2_build_splits",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        # Stage: load inputs (small metadata only — edge frames are loaded
        # per stage and freed after, to keep peak RSS under the §0.3.3 budget)
        if not (args.resume and "load" in run.done_stages):
            with run.stage("load"):
                families = load_families(cfg.families_config)
                areas = load_areas(cfg)
                log.info("inputs_loaded", families=list(families.keys()), areas=list(areas.keys()))

        # --- Track-R stages (small graph) ---
        if not (args.resume and "random_r" in run.done_stages):
            with run.stage("random_r"):
                trackr_edges = pl.read_parquet(cfg.trackr.edges)
                random_r = random_edge_split(
                    trackr_edges,
                    fractions=cfg.random.fractions,
                    seed=cfg.random.seed,
                )
                write_split(cfg.output_dir, "trackr_random", random_r)
                del trackr_edges, random_r
                gc.collect()

        if not (args.resume and "family_r" in run.done_stages):
            with run.stage("family_r"):
                trackr_edges = pl.read_parquet(cfg.trackr.edges)
                trackr_nodes = pl.read_parquet(cfg.trackr.edges.parent / "nodes.parquet")
                family_sets_r = {
                    name: expand_family_trackr(members, trackr_nodes)
                    for name, members in families.items()
                }
                # Per-family splits (separate AD/PD programs, spec §7.4(3))
                for fam_name, fam_nodes in family_sets_r.items():
                    fam_split = family_holdout_split(
                        trackr_edges,
                        fam_nodes,
                        eval_relation_col=cfg.trackr.relation_col,
                        eval_relation_values=set(cfg.trackr.family_eval_relations),
                        seed=cfg.family.seed,
                        valid_fraction=cfg.family.valid_fraction,
                    )
                    write_split(cfg.output_dir, f"trackr_family_{fam_name}", fam_split)
                    del fam_split
                    gc.collect()
                # Pooled split retained for continuity with earlier baselines
                family_r = combine_family_splits(
                    trackr_edges,
                    family_sets_r,
                    cfg.trackr.relation_col,
                    set(cfg.trackr.family_eval_relations),
                    cfg.family.seed,
                    cfg.family.valid_fraction,
                )
                write_split(cfg.output_dir, "trackr_family", family_r)
                del trackr_edges, trackr_nodes, family_r
                gc.collect()

        # Stage: TxGNN area split (Track-R only, neurodegenerative in config)
        if not (args.resume and "txgnn_r" in run.done_stages):
            with run.stage("txgnn_r"):
                trackr_edges = pl.read_parquet(cfg.trackr.edges)
                area_cfg = cfg.areas["neurodigenerative"]
                txgnn_r = txgnn_area_split(
                    trackr_edges,
                    area_disease_ids=areas["neurodigenerative"],
                    indication_relations=set(area_cfg.indication_relations),
                    relation_col=cfg.trackr.relation_col,
                    neighborhood_remove_fraction=1.0 - cfg.txgnn.neighborhood_keep_fraction,
                    seed=cfg.txgnn.seed,
                )
                write_split(cfg.output_dir, "trackr_txgnn_neurodegenerative", txgnn_r)
                del trackr_edges, txgnn_r
                gc.collect()

        # --- Track-D stages (21.8M-edge graph; one frame in memory at a time) ---
        if not (args.resume and "random_d" in run.done_stages):
            with run.stage("random_d"):
                # Streaming path: the joined frame is never materialized.
                write_random_split_streaming(
                    cfg.trackd.edges,
                    cfg.output_dir / "trackd_random",
                    fractions=cfg.random.fractions,
                    seed=cfg.random.seed,
                )
                gc.collect()

        if not (args.resume and "family_d" in run.done_stages):
            with run.stage("family_d"):
                family_sets_d = {name: set(members) for name, members in families.items()}
                # Per-family splits (separate AD/PD programs, spec §7.4(3))
                for fam_name, fam_nodes in family_sets_d.items():
                    write_family_split_streaming(
                        cfg.trackd.edges,
                        cfg.output_dir / f"trackd_family_{fam_name}",
                        family_nodes=fam_nodes,
                        label_col=cfg.trackd.label_col,
                        eval_label=cfg.trackd.family_eval_label,
                        seed=cfg.family.seed,
                        valid_fraction=cfg.family.valid_fraction,
                    )
                    gc.collect()
                # Pooled split retained for continuity with earlier baselines
                write_family_split_streaming(
                    cfg.trackd.edges,
                    cfg.output_dir / "trackd_family",
                    family_nodes=set().union(*family_sets_d.values()),
                    label_col=cfg.trackd.label_col,
                    eval_label=cfg.trackd.family_eval_label,
                    seed=cfg.family.seed,
                    valid_fraction=cfg.family.valid_fraction,
                )
                gc.collect()

        # Stage: temporal split (Track-D only; secondary because of low coverage)
        if not (args.resume and "temporal_d" in run.done_stages):
            with run.stage("temporal_d"):
                write_temporal_split_streaming(
                    cfg.trackd.edges,
                    cfg.output_dir / "trackd_temporal",
                    cutoff_year=cfg.temporal.cutoff_year,
                    label_value=cfg.trackd.temporal_label,
                )
                gc.collect()

        # Stage: report
        with run.stage("report", outputs=[cfg.report_path]):
            write_report(cfg, args.config, locals())

    log.info("splits_done")
    return 0


def load_families(path: Path) -> dict[str, list[str]]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    return {k: v["members"] for k, v in data["families"].items()}


def load_areas(cfg: SplitsConfig) -> dict[str, set[str]]:
    areas: dict[str, set[str]] = {}
    for name, acfg in cfg.areas.items():
        df = pl.read_csv(cfg.txgnn_dir / acfg.file)
        # TxGNN CSV node_id is the numeric MONDO local id. Track-R canonical
        # node_id form is "mondo:<id>" with no leading zeros.
        areas[name] = {f"mondo:{int(x)}" for x in df["node_id"].cast(pl.String).to_list()}
    return areas


def expand_family_trackr(members: list[str], trackr_nodes: pl.DataFrame) -> set[str]:
    """Map MONDO family ids to Track-R disease node ids.

    Track-R collapses many MONDO diseases into grouped nodes
    (e.g., mondo_grouped:11913_..._4975_... for Alzheimer disease). We therefore
    match both the exact canonical id and any grouped node that contains a
    family member's numeric id.
    """
    numeric = {re.sub(r"^mondo:0*", "", m) for m in members}
    exact_re = re.compile(
        r"^mondo:(?:" + "|".join(f"(?:0*)?{re.escape(n)}" for n in numeric) + r")$"
    )
    grouped_re = re.compile(
        r"^mondo_grouped:.*(?:^|_)(?:"
        + "|".join(f"(?:0*)?{re.escape(n)}" for n in numeric)
        + r")(?:_|$)"
    )

    def match(node_id: str) -> bool:
        return bool(exact_re.match(node_id) or grouped_re.match(node_id))

    disease_nodes = trackr_nodes.filter(pl.col("node_type") == "disease")
    return set(
        disease_nodes.filter(pl.col("node_id").map_elements(match, return_dtype=pl.Boolean))[
            "node_id"
        ].to_list()
    )


def combine_family_splits(
    edges: pl.DataFrame,
    family_sets: dict[str, set[str]],
    eval_col: str,
    eval_values: set[str],
    seed: int,
    valid_fraction: float,
) -> dict[str, pl.DataFrame]:
    """Hold out multiple families by pooling their nodes, then partition eval edges."""
    pooled = set().union(*family_sets.values())
    split = family_holdout_split(
        edges,
        family_node_ids=pooled,
        eval_relation_col=eval_col,
        eval_relation_values=eval_values,
        seed=seed,
        valid_fraction=valid_fraction,
    )
    return split


def write_split(output_dir: Path, split_name: str, split: dict[str, pl.DataFrame]) -> dict:
    """Persist train/valid/test parquet and a metadata sidecar."""
    base = output_dir / split_name
    base.mkdir(parents=True, exist_ok=True)
    meta: dict[str, object] = {"split_name": split_name, "partitions": {}}
    for part in ("train", "valid", "test"):
        path = base / f"{part}.parquet"
        split[part].write_parquet(path)
        meta["partitions"][part] = {
            "rows": split[part].height,
            "groups": split[part]["evidence_group_id"].n_unique(),
            "path": str(path),
        }
    # leakage-critical sidecars
    for sidecar in ("_removed_neighborhood", "_kept_neighborhood"):
        if sidecar in split:
            path = base / f"{sidecar.lstrip('_')}.parquet"
            split[sidecar].write_parquet(path)
            meta[sidecar] = {
                "rows": split[sidecar].height,
                "groups": split[sidecar]["evidence_group_id"].n_unique(),
                "path": str(path),
            }
    try:
        assert_split_disjoint({k: split[k] for k in ("train", "valid", "test")})
        meta["leakage_guard"] = "PASS"
    except SplitLeakageError as exc:
        meta["leakage_guard"] = f"FAIL: {exc}"
        raise
    meta_path = base / "metadata.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    log.info("split_written", split=split_name, **meta["partitions"])
    return meta


def write_random_split_streaming(
    edges_path: Path,
    out_dir: Path,
    fractions: dict[str, float],
    seed: int,
) -> None:
    """Streaming random group-level split for very large edge tables.

    The gid→partition assignment is built on unique gids only (~small), then
    each partition is written via a lazy join + streaming sink — the joined
    frame is never materialized, keeping peak RSS under the §0.3.3 budget.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    gids = (
        pl.scan_parquet(edges_path)
        .select(GROUP_COL)
        .unique()
        .collect(engine="streaming")[GROUP_COL]
    )
    assignment = group_assignment(gids, fractions, seed).lazy()
    meta: dict[str, object] = {
        "split_name": out_dir.name,
        "partitions": {},
        "streaming": True,
    }
    for part in ("train", "valid", "test"):
        path = out_dir / f"{part}.parquet"
        (
            pl.scan_parquet(edges_path)
            .join(assignment, on=GROUP_COL, how="left")
            .filter(pl.col("partition") == part)
            .drop("partition")
            .sink_parquet(path, engine="streaming")
        )
        stats = (
            pl.scan_parquet(path)
            .select(rows=pl.len(), groups=pl.col(GROUP_COL).n_unique())
            .collect(engine="streaming")
            .to_dicts()[0]
        )
        meta["partitions"][part] = {**stats, "path": str(path)}
    # Leakage guard: pairwise shared-gid count must be zero (checked lazily).
    gid_frames = {
        p: pl.scan_parquet(out_dir / f"{p}.parquet").select(GROUP_COL).unique()
        for p in ("train", "valid", "test")
    }
    for a, b in (("train", "valid"), ("train", "test"), ("valid", "test")):
        shared = gid_frames[a].join(gid_frames[b], on=GROUP_COL).collect(engine="streaming").height
        if shared:
            meta["leakage_guard"] = f"FAIL: {shared} shared gids {a}/{b}"
            raise SplitLeakageError(f"{shared} evidence groups in both {a} and {b}")
    meta["leakage_guard"] = "PASS"
    (out_dir / "metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    log.info("split_written", split=out_dir.name, **meta["partitions"])


def _partition_stats(out_dir: Path, partitions: tuple[str, ...]) -> dict:
    meta: dict[str, object] = {"split_name": out_dir.name, "partitions": {}, "streaming": True}
    for part in partitions:
        path = out_dir / f"{part}.parquet"
        stats = (
            pl.scan_parquet(path)
            .select(rows=pl.len(), groups=pl.col(GROUP_COL).n_unique())
            .collect(engine="streaming")
            .to_dicts()[0]
        )
        meta["partitions"][part] = {**stats, "path": str(path)}
    gid_frames = {
        p: pl.scan_parquet(out_dir / f"{p}.parquet").select(GROUP_COL).unique()
        for p in ("train", "valid", "test")
    }
    for a, b in (("train", "valid"), ("train", "test"), ("valid", "test")):
        shared = gid_frames[a].join(gid_frames[b], on=GROUP_COL).collect(engine="streaming").height
        if shared:
            meta["leakage_guard"] = f"FAIL: {shared} shared gids {a}/{b}"
            raise SplitLeakageError(f"{shared} evidence groups in both {a} and {b}")
    meta["leakage_guard"] = "PASS"
    (out_dir / "metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    log.info("split_written", split=out_dir.name, **meta["partitions"])
    return meta


def write_family_split_streaming(
    edges_path: Path,
    out_dir: Path,
    family_nodes: set[str],
    label_col: str,
    eval_label: str,
    seed: int,
    valid_fraction: float,
) -> None:
    """Streaming family holdout for very large edge tables (spec §6.2.2)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    touches = pl.col("src_id").is_in(family_nodes) | pl.col("dst_id").is_in(family_nodes)
    is_eval = touches & (pl.col(label_col) == eval_label)
    lf = pl.scan_parquet(edges_path)
    eval_edges = lf.filter(is_eval).collect(engine="streaming")
    lf.filter(~touches).sink_parquet(out_dir / "train.parquet", engine="streaming")
    lf.filter(touches & ~is_eval).sink_parquet(
        out_dir / "removed_neighborhood.parquet", engine="streaming"
    )
    assignment = group_assignment(
        eval_edges[GROUP_COL], {"valid": valid_fraction, "test": 1 - valid_fraction}, seed
    )
    joined = eval_edges.join(assignment, on=GROUP_COL, how="left")
    for part in ("valid", "test"):
        joined.filter(pl.col("partition") == part).drop("partition").write_parquet(
            out_dir / f"{part}.parquet"
        )
    _partition_stats(out_dir, ("train", "valid", "test", "removed_neighborhood"))


def write_temporal_split_streaming(
    edges_path: Path,
    out_dir: Path,
    cutoff_year: int,
    label_value: str,
) -> None:
    """Streaming temporal split for very large edge tables (spec §6.2.3)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    year_init = pl.col("year_initial").cast(pl.Int32, strict=False)
    year_final = pl.col("year_final").cast(pl.Int32, strict=False)
    has_years = year_init.is_not_null()
    temporal = (
        pl.scan_parquet(edges_path)
        .filter(pl.col("label") == label_value)
        .collect(engine="streaming")
    )
    test_gids = temporal.filter(has_years & (year_init > cutoff_year))[GROUP_COL].unique()
    train_gids = temporal.filter(~has_years | (year_final <= cutoff_year))[GROUP_COL].unique()
    del temporal
    gc.collect()
    test_gid_frame = test_gids.to_frame(GROUP_COL).lazy()
    lf = pl.scan_parquet(edges_path)
    lf.join(test_gid_frame, on=GROUP_COL, how="semi").sink_parquet(
        out_dir / "test.parquet", engine="streaming"
    )
    lf.filter(
        (pl.col("label") != label_value) | pl.col(GROUP_COL).is_in(train_gids.implode())
    ).join(test_gid_frame, on=GROUP_COL, how="anti").sink_parquet(
        out_dir / "train.parquet", engine="streaming"
    )
    pl.DataFrame(schema=pl.scan_parquet(edges_path).head(0).collect().schema).write_parquet(
        out_dir / "valid.parquet"
    )
    _partition_stats(out_dir, ("train", "valid", "test"))


def write_report(cfg: SplitsConfig, config_path: Path, local_vars: dict) -> None:
    """Generate reports/m2a/splits_report.md."""
    meta_files = sorted(d / "metadata.json" for d in cfg.output_dir.iterdir() if d.is_dir())
    # Per-split summary rows
    summary_lines = []
    for data in sorted(
        [json.loads(m.read_text()) for m in meta_files],
        key=lambda x: x["split_name"],
    ):
        parts = data["partitions"]
        total_groups = sum(parts[p]["groups"] for p in ("train", "valid", "test"))
        total_rows = sum(parts[p]["rows"] for p in ("train", "valid", "test"))
        summary_lines.append(
            f"| {data['split_name']} | {total_rows:,} | {total_groups:,} | "
            f"{parts['train']['groups']:,} | {parts['valid']['groups']:,} | "
            f"{parts['test']['groups']:,} | {data.get('leakage_guard', 'N/A')} |"
        )

    # Track-R family expansion note
    family_sets_r = local_vars.get("family_sets_r", {})
    family_expansion_lines = [
        f"- **{name}**: {len(members)} Track-R disease node(s)"
        for name, members in family_sets_r.items()
    ]

    text = f"""<!-- track_scope: benchmark -->
# M2A-2 Split workstream — leakage-guarded partitions (spec §6.2)

**Status:** DRAFT — split definitions are engineering outputs; biological
sign-off of disease-family membership is in `configs/disease_families.yaml`.

All splits operate on `evidence_group_id` (§5.4). The leakage guard
(`assert_split_disjoint`) was run for every split and is also enforced by CI in
`tests/unit/test_splits.py`.

## Splits produced

| split | total rows | total groups | train groups | valid groups | test groups | leakage guard |
|---|---|---|---|---|---|---|
{chr(10).join(summary_lines)}

## Design notes

- **random_edge_split** (Track-R + Track-D): random group-level partition
  (90/5/5). This split is intentionally inflated by shared evidence groups
  (§6.2.4) and is provided for comparability only.
- **family_holdout_split** (Track-R + Track-D): pooled AD + PD family nodes.
  Eval edges touching a family disease go to valid/test; all other edges
  touching a family disease are removed from train (neighborhood suppression).
- **txgnn_area_split** (Track-R, neurodegenerative area): indication/
  contraindication edges for the 68 TxGNN-area MONDO diseases → test;
  95% of remaining family-neighborhood edges removed from train per the TxGNN
  protocol.
- **temporal_split** (Track-D): train = DIS-GEN edges with `year_final ≤ 2018`
  or no years; test = DIS-GEN edges with `year_initial > 2018`; boundary
  straddlers dropped. **Secondary split:** only ~0.25% of DIS-GEN edges carry
  year annotations (Appendix D.5), so this is a coverage-limited probe.

### Track-R family node expansion

Track-R represents several family-level diseases as `mondo_grouped:*` nodes.
The family config lists individual MONDO ids; for Track-R we additionally
include any grouped disease node whose id contains a family member numeric id:

{chr(10).join(family_expansion_lines)}

## Config / reproduction

- Config: `{config_path}`
- Manifest: `{cfg.manifest_path}`
- Output: `{cfg.output_dir}`
- Re-run: `python scripts/m2_build_splits.py --config {config_path} --resume`
"""
    cfg.report_path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
