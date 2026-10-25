#!/usr/bin/env python
"""m2_build_crosswalk — MONDO crosswalk + disease family drafts (spec §6.2.2).

Stages (resumable; --dry-run prints the plan):
  1. crosswalk — parse mondo.obo xrefs → mondo_crosswalk.parquet; apply to
     Track-D disease nodes (native id → MONDO) and measure coverage
  2. families — is_a descendant closure of the AD/PD roots →
     configs/disease_families.yaml (DRAFT, owner sign-off gate §0.3.9)
  3. report — coverage + family sizes → reports/m2a/mondo_crosswalk.md

Usage: python scripts/m2_build_crosswalk.py --config configs/m2_crosswalk.yaml
           [--dry-run] [--resume]
"""

import argparse
import sys
from pathlib import Path

import polars as pl
import yaml
from pydantic import BaseModel

from sina.config import load_config
from sina.logging import configure_logging, get_logger
from sina.manifests import Run
from sina.obo import build_xref_map, descendant_closure, parse_obo

log = get_logger(__name__)


class FamilySpec(BaseModel):
    key: str
    root: str
    root_name: str


class CrosswalkConfig(BaseModel):
    mondo_obo: Path
    trackd_nodes: Path
    out_dir: Path
    report_path: Path
    manifest_path: Path
    families_config: Path
    xref_prefixes: list[str]
    families: list[FamilySpec]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, CrosswalkConfig)
    if args.dry_run:
        print(f"would parse {cfg.mondo_obo} → {cfg.out_dir} (report: {cfg.report_path})")
        return 0

    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)
    crosswalk_path = cfg.out_dir / "mondo_crosswalk.parquet"
    trackd_map_path = cfg.out_dir / "trackd_disease_mondo.parquet"

    with Run(
        script="m2_build_crosswalk",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        if not (args.resume and "crosswalk" in run.done_stages):
            with run.stage("crosswalk", outputs=[crosswalk_path, trackd_map_path]):
                terms = list(parse_obo(cfg.mondo_obo))
                xref_map = build_xref_map(terms, set(cfg.xref_prefixes))
                # canonical node_id form ("mondo:0004975"), matching nodes.parquet
                xref_map = xref_map.with_columns(mondo_id=pl.col("mondo_id").str.to_lowercase())
                xref_map.write_parquet(crosswalk_path)
                coverage = apply_to_trackd(cfg, xref_map, trackd_map_path)
                (cfg.out_dir / "coverage.yaml").write_text(
                    yaml.safe_dump(coverage), encoding="utf-8"
                )

        if not (args.resume and "families" in run.done_stages):
            # owner sign-off gate (§0.3.9): never regenerate an APPROVED
            # families file — curation edits (e.g., CHANGES/0011) would be lost
            if families_config_approved(cfg.families_config):
                log.warning(
                    "families_config_approved_skip_regeneration",
                    path=str(cfg.families_config),
                )
            else:
                with run.stage("families", outputs=[cfg.families_config]):
                    terms = list(parse_obo(cfg.mondo_obo))
                    release = read_data_version(cfg.mondo_obo)
                    families = {}
                    for fam in cfg.families:
                        # canonical node_id form ("mondo:…") so splits match node_ids
                        members = sorted(m.lower() for m in descendant_closure(terms, fam.root))
                        if len(members) <= 1:
                            raise ValueError(f"family {fam.key}: {fam.root} has no descendants")
                        families[fam.key] = {
                            "root": fam.root.lower(),
                            "root_name": fam.root_name,
                            "n_members": len(members),
                            "members": members,
                        }
                        log.info("family_built", family=fam.key, n_members=len(members))
                    payload = {
                        "status": "DRAFT — owner sign-off required (§0.3.9 semantics gate)",
                        "generated_by": "scripts/m2_build_crosswalk.py",
                        "mondo_release": release,
                        "families": families,
                    }
                    cfg.families_config.write_text(
                        yaml.safe_dump(payload, sort_keys=False), encoding="utf-8"
                    )

        with run.stage("report", outputs=[cfg.report_path]):
            write_report(cfg, crosswalk_path, trackd_map_path)

    log.info("crosswalk_done")
    return 0


def families_config_approved(path: Path) -> bool:
    """True iff the families config exists and carries an APPROVED status."""
    if not path.exists():
        return False
    try:
        payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError:
        return False
    return str(payload.get("status", "")).startswith("APPROVED")


def read_data_version(path: Path) -> str:
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("data-version:"):
                return line.split(":", 1)[1].strip()
            if line == "[Term]":
                break
    return "unknown"


def apply_to_trackd(cfg: CrosswalkConfig, xref_map: pl.DataFrame, out_path: Path) -> dict:
    """Map Track-D disease nodes' native ids to MONDO via the crosswalk."""
    nodes = pl.read_parquet(cfg.trackd_nodes)
    diseases = nodes.filter(pl.col("node_type") == "disease").select(
        "node_id", "native_id", "mondo_id"
    )

    def external_id(native: str) -> str:
        # MONDO_0004975 → MONDO:0004975 ; DOID_0050890 → DOID:0050890 ;
        # Orphanet_1637 → Orphanet:1637 ; OTAR_0000001 → OTAR:0000001
        prefix, _, local = native.partition("_")
        return f"{prefix}:{local}"

    diseases = diseases.with_columns(
        external_id=pl.col("native_id").map_elements(external_id, return_dtype=pl.String)
    )
    mapped = diseases.join(
        xref_map.select("external_id", pl.col("mondo_id").alias("mondo_id_xref")),
        on="external_id",
        how="left",
    )
    mapped = mapped.with_columns(mondo_final=pl.coalesce(["mondo_id", "mondo_id_xref"])).select(
        "node_id", "native_id", "external_id", mondo_id=pl.col("mondo_final")
    )
    mapped.write_parquet(out_path)
    n = diseases.height
    n_native = diseases.filter(pl.col("mondo_id").is_not_null()).height
    n_mapped = mapped.filter(pl.col("mondo_id").is_not_null()).height
    return {
        "trackd_diseases": n,
        "mondo_native": n_native,
        "mondo_after_crosswalk": n_mapped,
        "coverage": n_mapped / n,
        "unmapped": n - n_mapped,
    }


def write_report(cfg: CrosswalkConfig, crosswalk_path: Path, trackd_map_path: Path) -> None:
    coverage = yaml.safe_load((cfg.out_dir / "coverage.yaml").read_text())
    mapped = pl.read_parquet(trackd_map_path).with_columns(
        prefix=pl.col("native_id").str.extract(r"^([A-Za-z]+)")
    )
    fam = yaml.safe_load(cfg.families_config.read_text())
    prefix_rows = "\n".join(
        f"| {r['prefix']} | {r['total']:,} | {r['n_mapped']:,} | {r['n_mapped'] / r['total']:.1%} |"
        for r in mapped.group_by("prefix")
        .agg(total=pl.len(), n_mapped=pl.col("mondo_id").is_not_null().sum())
        .sort("total", descending=True)
        .to_dicts()
    )
    fam_rows = "\n".join(
        f"| {k} | {v['root']} ({v['root_name']}) | {v['n_members']:,} |"
        for k, v in fam["families"].items()
    )
    fam_trackd_rows = []
    for k, v in fam["families"].items():
        n_in = mapped.filter(pl.col("mondo_id").is_in(v["members"])).height
        fam_trackd_rows.append(f"| {k} | {v['n_members']:,} | {n_in:,} |")
    text = f"""<!-- track_scope: benchmark -->
# MONDO crosswalk + disease families — M2A-1 (spec §6.2.2)

**Status: DRAFT — disease families await owner sign-off (§0.3.9)**

MONDO release: `{fam["mondo_release"]}` (`data/raw/mondo/mondo.obo`, sha256 in
REGISTRY). Parser: `src/sina/obo.py` (unit-tested on fixture OBO). All MONDO
ids below use the canonical node_id form (`mondo:…`).

## Crosswalk coverage (closes the CHANGES/0008 gap)

- Track-D disease nodes: **{coverage["trackd_diseases"]:,}**
- MONDO-native (no mapping needed): {coverage["mondo_native"]:,}
- MONDO after crosswalk: **{coverage["mondo_after_crosswalk"]:,}
  ({coverage["coverage"]:.1%})**
- Unmapped: {coverage["unmapped"]:,} — dominated by **OBA** (Ontology of
  Biological Attributes; these are attribute terms, not diseases — they can
  never map) and **EFO** terms without a MONDO xref in this release. The
  mappable disease universe is therefore ~25k nodes, not 36k.

Per-prefix coverage:

| prefix | Track-D nodes | MONDO-mapped | coverage |
|---|---|---|---|
{prefix_rows}

## Disease families (DRAFT — owner sign-off required)

| family | root | members (root + is_a descendants) |
|---|---|---|
{fam_rows}

Track-D presence of each family (MONDO-mapped nodes inside the family):

| family | family size | present in Track-D |
|---|---|---|
{chr(10).join(fam_trackd_rows)}

Config: `configs/disease_families.yaml`. Splits (M2A-2) hold out evidence
groups touching these families; AD and PD are separate holdouts (§6.2.2).
**Limitation for sign-off:** diseases without a MONDO mapping (EFO/OBA-only)
cannot be family-attributed; the owner may expand family membership at
sign-off (e.g., name-matched EFO terms) — any edit re-triggers M2A-2.
"""
    cfg.report_path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
