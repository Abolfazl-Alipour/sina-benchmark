#!/usr/bin/env python
"""m2_lodo_shuffles — B3 dual-shuffles negative-control generator.

Directives v2 Phase B3 + CHANGES/0022 Light Gate 0 permutation control. For
each configured LODO unit and each shuffle mode, rewrite ONLY train2.txt (the
DIS-GEN supervision) of ``data/processed/m2/lodo_exports/<unit>/`` into a
sibling export dir ``<unit>_shuffle_<mode>/``; train1/valid/test/entity files
and strip_list.json are copied unchanged:

  across_disease  popularity null — each gene keeps its exact supervision
                  degree; diseases resampled uniformly from the pooled disease
                  roster (all configured units' train2, config-selectable)
  within_disease  specificity null — each disease keeps its exact supervision
                  degree; genes resampled uniformly from the unit's own train2

Leakage contract (sina.lodo_shuffles): no duplicate tuples, shuffled train2
disjoint from that unit's valid/test tuples (tuple-level, mirroring the A2
guard of m2_build_lodo_exports.py), preserved marginal asserted exactly, fail
loudly after max_attempts resamples. Seeded determinism (config seed, default
1337); sorted output lines give byte-identical train2.txt on re-run. Expected
result when the harness trains on the shuffled supervision: filtered MRR
collapses to the degree baseline (~0.01) — the disease-associated gene
recovery signal is destroyed.

Usage:
  python scripts/m2_lodo_shuffles.py --config configs/m2_lodo_shuffles.yaml
    [--dry-run] [--resume]
"""

import argparse
import json
import shutil
import sys
from pathlib import Path

from pydantic import BaseModel, field_validator

from sina.config import load_config
from sina.lodo_shuffles import (
    MODES,
    across_disease_shuffle,
    curie_safe,
    read_triples,
    verify_shuffle,
    within_disease_shuffle,
    write_triples,
)
from sina.logging import configure_logging, get_logger
from sina.manifests import Run
from sina.registry import sha256_file

log = get_logger(__name__)

COPY_FILES = [
    "train1.txt",
    "valid.txt",
    "test.txt",
    "entity_types.txt",
    "entity_names.txt",
    "strip_list.json",
]


class ShuffleConfig(BaseModel):
    exports_dir: Path
    units: list[str]
    modes: list[str] = list(MODES)
    seed: int = 1337
    disease_pool: str = "pooled"  # pooled | per_unit (across_disease disease roster)
    source_suffix: str = ""  # e.g. "_keep_disphe" to source from KEEP-graph exports
    output_suffix: str = ""  # e.g. "_keep" -> output dirs "<unit>_keep_shuffle_<mode>"
    max_attempts: int = 1000
    manifest_path: Path
    report_path: Path

    @field_validator("modes")
    @classmethod
    def _modes_known(cls, v: list[str]) -> list[str]:
        unknown = set(v) - set(MODES)
        if unknown:
            raise ValueError(f"unknown shuffle modes: {sorted(unknown)}")
        return v

    @field_validator("disease_pool")
    @classmethod
    def _pool_known(cls, v: str) -> str:
        if v not in ("pooled", "per_unit"):
            raise ValueError("disease_pool must be 'pooled' or 'per_unit'")
        return v


def shuffle_out_dir(cfg: ShuffleConfig, unit: str, mode: str) -> Path:
    return cfg.exports_dir / f"{curie_safe(unit)}{cfg.output_suffix}_shuffle_{mode}"


def shuffle_src_dir(cfg: ShuffleConfig, unit: str) -> Path:
    return cfg.exports_dir / f"{curie_safe(unit)}{cfg.source_suffix}"


def pooled_disease_roster(cfg: ShuffleConfig) -> list[str]:
    """Diseases present in train2 across ALL configured units (pooled roster)."""
    diseases: set[str] = set()
    for unit in cfg.units:
        train2 = read_triples(shuffle_src_dir(cfg, unit) / "train2.txt")
        diseases |= {t[0] for t in train2}
    log.info("disease_roster_pooled", units=len(cfg.units), diseases=len(diseases))
    return sorted(diseases)


def shuffle_unit(cfg: ShuffleConfig, unit: str, pooled_roster: list[str]) -> None:
    """Write both shuffle-mode export dirs for one unit (+ verification rows)."""
    src_dir = shuffle_src_dir(cfg, unit)
    train2 = read_triples(src_dir / "train2.txt")
    forbidden = set(read_triples(src_dir / "valid.txt")) | set(read_triples(src_dir / "test.txt"))
    if not train2:
        raise ValueError(f"{unit}: empty train2.txt in {src_dir}")
    roster = pooled_roster if cfg.disease_pool == "pooled" else sorted({t[0] for t in train2})

    for mode in cfg.modes:
        seed = f"{cfg.seed}:{unit}:{mode}"  # deterministic per (seed, unit, mode)
        if mode == "across_disease":
            shuffled = across_disease_shuffle(train2, roster, forbidden, seed, cfg.max_attempts)
        else:
            shuffled = within_disease_shuffle(train2, forbidden, seed, cfg.max_attempts)
        row = verify_shuffle(train2, shuffled, mode, forbidden)

        out_dir = shuffle_out_dir(cfg, unit, mode)
        out_dir.mkdir(parents=True, exist_ok=True)
        for name in COPY_FILES:
            shutil.copy2(src_dir / name, out_dir / name)
        write_triples(shuffled, out_dir / "train2.txt")

        provenance = {
            "unit": unit,
            "source_dir": str(src_dir),
            "output_dir": str(out_dir),
            "mode": mode,
            "seed": cfg.seed,
            "seed_key": seed,
            "disease_pool": cfg.disease_pool if mode == "across_disease" else "n/a",
            "max_attempts": cfg.max_attempts,
            "sha256_source_train2": sha256_file(src_dir / "train2.txt"),
            "sha256_shuffled_train2": sha256_file(out_dir / "train2.txt"),
            "verification": row,
        }
        (out_dir / "shuffle_provenance.json").write_text(
            json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
        )
        log.info(
            "unit_shuffled",
            unit=unit,
            mode=mode,
            edges=row["edges_after"],
            marginal=row["marginal_checked"],
        )


def write_report(cfg: ShuffleConfig) -> None:
    """Assemble reports/m2b/lodo_shuffles.md from on-disk provenance rows."""
    rows = []
    for unit in cfg.units:
        for mode in cfg.modes:
            prov_path = shuffle_out_dir(cfg, unit, mode) / "shuffle_provenance.json"
            prov = json.loads(prov_path.read_text(encoding="utf-8"))
            rows.append(prov)

    table_rows = "\n".join(
        f"| {p['unit']} | {p['mode']} | {p['verification']['edges_before']:,} | "
        f"{p['verification']['edges_after']:,} | {p['verification']['genes']:,} | "
        f"{p['verification']['diseases']:,} | "
        f"{p['verification']['marginal_checked']} | "
        f"{'PASS' if p['verification']['degree_marginal_preserved'] else 'FAIL'} | "
        f"{p['verification']['n_duplicate_tuples']} | "
        f"{'PASS' if p['verification']['valid_test_disjoint'] else 'FAIL'} | "
        f"`{p['sha256_shuffled_train2'][:12]}` |"
        for p in rows
    )
    units_s = ", ".join(f"`{u}`" for u in cfg.units)
    modes_s = ", ".join(f"`{m}`" for m in cfg.modes)

    text = f"""<!-- track_scope: benchmark -->
# LODO supervision shuffles — dual negative controls (Directives v2 B3; CHANGES/0022)

Permutation negative control for the Light Gate 0 gate: the DIS-GEN supervision
(train2.txt) of {len(cfg.units)} LODO export units is rewritten under two
degree-preserving nulls, one sibling export dir per unit and mode
(`<unit>_shuffle_<mode>/` under `{cfg.exports_dir}`). train1/valid/test,
entity files and strip_list.json are copied unchanged; each new dir carries a
`shuffle_provenance.json` (mode, seed, source unit, sha256 of source and
shuffled train2). Units: {units_s}. Modes: {modes_s}. Seed: {cfg.seed}.

- **across_disease (popularity null):** every gene keeps its exact supervision
  degree; its disease edges are reassigned to uniformly sampled diseases
  (roster: {cfg.disease_pool}). Gene popularity survives; disease pairing is
  destroyed.
- **within_disease (specificity null):** every disease keeps its exact
  supervision degree; its edges are reassigned to uniformly sampled genes of
  the same unit's train2. Disease degree survives; gene identity is destroyed.

Expected result when the harness trains on a shuffled export: filtered MRR on
the held-out key collapses to the degree baseline (~0.01). If it does not, the
pipeline leaks supervision signal and disease-associated gene recovery claims
are blocked until the leak is found (CHANGES/0022 permutation gate).

## Verification (all rows asserted at shuffle time)

| unit | mode | before | after | genes | diseases | marginal | check | dups | disjoint | sha256 |
|---|---|---|---|---|---|---|---|---|---|---|
{table_rows}

Leakage contract: no duplicate (disease, relation, gene) tuples in any shuffled
train2; shuffled tuples are disjoint from that unit's valid.txt/test.txt
tuples (tuple-level resampling, max {cfg.max_attempts} attempts per edge, hard
failure beyond); the preserved marginal matches exactly; re-running with seed
{cfg.seed} reproduces byte-identical train2.txt files (sorted lines).

## Reproduction

`python scripts/m2_lodo_shuffles.py --config configs/m2_lodo_shuffles.yaml --resume`
"""
    cfg.report_path.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, ShuffleConfig)

    if args.dry_run:
        print(
            f"would shuffle {len(cfg.units)} units x {len(cfg.modes)} modes "
            f"(seed {cfg.seed}, disease_pool {cfg.disease_pool}) → {cfg.exports_dir} "
            f"(report: {cfg.report_path})"
        )
        for unit in cfg.units:
            for mode in cfg.modes:
                print(f"  {unit} [{mode}] → {shuffle_out_dir(cfg, unit, mode)}")
        return 0

    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)
    pooled_roster = pooled_disease_roster(cfg) if cfg.disease_pool == "pooled" else []

    with Run(
        script="m2_lodo_shuffles",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        for unit in cfg.units:
            stage = f"shuffle_{curie_safe(unit)}"
            if args.resume and stage in run.done_stages:
                continue
            outputs = [shuffle_out_dir(cfg, unit, mode) / "train2.txt" for mode in cfg.modes] + [
                shuffle_out_dir(cfg, unit, mode) / "shuffle_provenance.json" for mode in cfg.modes
            ]
            with run.stage(stage, outputs=outputs):
                shuffle_unit(cfg, unit, pooled_roster)

        with run.stage("report", outputs=[cfg.report_path]):
            write_report(cfg)

    log.info("lodo_shuffles_done", units=len(cfg.units), modes=cfg.modes)
    return 0


if __name__ == "__main__":
    sys.exit(main())
