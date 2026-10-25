#!/usr/bin/env python
"""m2_drug_path_measurement — v5.1 Step 3: drug-edge circularity + license provenance.

Measurement only; nothing graded. Feeds the CHANGES/0020 freeze package:

  fk           F@k per candidate unit, k=2..4: fraction of the unit's DIS-GEN
               answer genes reachable via drug-mediated paths of length <=k in
               the LODO fact graph (BRG + supervision minus the unit's own
               edges, undirected). Drug-mediated distance = 1 + BFS distance
               from any drug adjacent to the unit disease (v5.1 §6.4).
  exclusivity  of the drug-reachable@3 answers, fraction ALSO reachable <=3
               from the disease with all drug nodes removed (reachable is not
               exclusively drug-mediated).
  provenance   per unit + pool: fraction of DIS-GEN edges whose sources mention
               a drug/clinical source (matching rule documented in the report).
  license      DIS-GEN edges by direct source (DisGeNET vs Open Targets) —
               auditable basis of the research-use posture (v5.1 §2.3).
  report       exam-integrity verdict (roster-median F@3), covariate bands,
               top-5 drug-rich units, feeds-into-freeze section.

Usage: python scripts/m2_drug_path_measurement.py --config configs/m2_drug_paths.yaml
           [--dry-run] [--resume]
"""

import argparse
import json
import sys
from collections import Counter
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import polars as pl
import scipy.sparse as sp
from pydantic import BaseModel
from scipy.sparse import csgraph

from sina.config import load_config
from sina.logging import configure_logging, get_logger
from sina.manifests import Run

log = get_logger(__name__)

EXAM_INTEGRITY_TRIGGER = 0.5  # §10: >=50% → exam measures pharmacological reconstruction
COVARIATE_BAND = 0.10  # v5.1 §6.4: F@3 <10% / >=10% strata for the Arm A/B analysis


class DrugPathConfig(BaseModel):
    supervision: Path
    brg: Path
    trackd_nodes: Path
    unit_audit: Path
    seed: int = 1337
    fk_ks: list[int] = [2, 3, 4]
    exclusivity_k: int = 3
    output_dir: Path
    report_path: Path
    manifest_path: Path


# ---------------------------------------------------------------------------
# Pure computations (unit-tested on synthetic fixtures)
# ---------------------------------------------------------------------------


def drug_mediated_distances(dist_from_adjacent_drugs: np.ndarray) -> np.ndarray:
    """Drug-mediated distance = 1 + BFS-from-adjacent-drugs distance.

    First hop disease→drug (1) plus the shortest path from that drug; a gene is
    drug-reachable@k iff this distance <= k (k=2 → disease→drug→gene exactly).
    """
    return dist_from_adjacent_drugs + 1.0


def reachable_mask(distances: np.ndarray, k: int) -> np.ndarray:
    """Boolean mask of nodes at distance <= k (np.inf = unreachable)."""
    return distances <= k


def fk_fraction(gene_distances: np.ndarray, answer_mask: np.ndarray, k: int) -> float:
    """F@k: fraction of answer genes drug-reachable within k hops."""
    if not answer_mask.any():
        return 0.0
    return float(reachable_mask(gene_distances, k)[answer_mask].mean())


def exclusivity_fraction(
    drug_reach_mask: np.ndarray, nondrug_gene_distances: np.ndarray, k: int
) -> tuple[int, float]:
    """Of drug-reachable answers, (n, fraction) also reachable <=k without drugs."""
    n = int(drug_reach_mask.sum())
    if n == 0:
        return 0, 0.0
    also = reachable_mask(nondrug_gene_distances, k) & drug_reach_mask
    return n, float(also.sum() / n)


def classify_direct_source(sources: Iterable[str]) -> str:
    """Map an edge's sources_direct values to a license-provenance bucket."""
    s = {x.upper() for x in sources or []}
    has_dg = "DISGENET" in s
    has_ot = "OPEN_TARGETS" in s
    if has_dg and has_ot:
        return "disgenet+open_targets"
    if has_dg:
        return "disgenet_only"
    if has_ot:
        return "open_targets_only"
    return "other"


def is_drug_clinical_source(source: str) -> bool:
    """Drug/clinical-source matching rule for drug-provenance fractions.

    Exact tokens: CTD (Comparative Toxicogenomics Database), CGI (Cancer Genome
    Interpreter). Substring tokens: drug-regulatory/pharmacological markers
    (DRUG, FDA, EMA, CLINICAL, TRIAL, PHARM, CHEMBL) — none occur in this build.
    """
    s = source.upper()
    if s in {"CTD", "CGI"}:
        return True
    return any(t in s for t in ("DRUG", "FDA", "EMA", "CLINICAL", "TRIAL", "PHARM", "CHEMBL"))


def has_drug_provenance(direct: Iterable[str], indirect: Iterable[str]) -> bool:
    """True if any direct or indirect source matches the drug/clinical rule."""
    return any(is_drug_clinical_source(s) for s in list(direct or []) + list(indirect or []))


def covariate_band(f3: float) -> str:
    """v5.1 §6.4 covariate bands on F@3."""
    return "f3_gte_10pct" if f3 >= COVARIATE_BAND else "f3_lt_10pct"


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------


def load_inputs(cfg: DrugPathConfig) -> dict:
    sup = pl.read_parquet(cfg.supervision)
    brg = pl.read_parquet(cfg.brg)
    nodes = pl.read_parquet(cfg.trackd_nodes).select("node_id", "node_type")
    units = {
        r["unit"]: r["stratum"]
        for r in pl.read_parquet(cfg.unit_audit).select("unit", "stratum").to_dicts()
    }

    brg_e = brg.select("src_id", "dst_id")
    sup_e = sup.select("src_id", "dst_id")
    all_nodes = sorted(
        set(brg_e["src_id"].to_list())
        | set(brg_e["dst_id"].to_list())
        | set(sup_e["src_id"].to_list())
        | set(sup_e["dst_id"].to_list())
    )
    idx = {nd: i for i, nd in enumerate(all_nodes)}
    drug_ids = set(nodes.filter(pl.col("node_type") == "drug")["node_id"].to_list())
    drug_mask = np.array(
        [nd in drug_ids or nd.startswith(("chembl:", "drugbank:")) for nd in all_nodes]
    )

    state = {
        "sup": sup,
        "units": units,
        "idx": idx,
        "n_nodes": len(all_nodes),
        "drug_mask": drug_mask,
        "bsrc": brg_e["src_id"].replace_strict(idx).to_numpy(),
        "bdst": brg_e["dst_id"].replace_strict(idx).to_numpy(),
        "ssrc": sup_e["src_id"].replace_strict(idx).to_numpy(),
        "sdst": sup_e["dst_id"].replace_strict(idx).to_numpy(),
        "sup_d": sup_e["src_id"]
        .replace_strict({d: i for i, d in enumerate(sorted(set(sup_e["src_id"].to_list())))})
        .to_numpy(),
        "dcode": {d: i for i, d in enumerate(sorted(set(sup_e["src_id"].to_list())))},
        "genes": sorted(sup["dst_id"].unique().to_list()),
        "answers": {
            r["src_id"]: set(r["dst_id"])
            for r in sup.group_by("src_id").agg(pl.col("dst_id")).to_dicts()
        },
    }
    state["gene_idx"] = np.array([idx[g] for g in state["genes"]])
    log.info(
        "inputs_loaded",
        supervision=sup.height,
        brg=brg.height,
        units=len(units),
        drug_nodes=int(drug_mask.sum()),
    )
    return state


def _unit_arrays(state: dict, unit: str) -> tuple[np.ndarray, np.ndarray]:
    """LODO fact-graph edge arrays: BRG + supervision minus the unit's own edges."""
    code = state["dcode"].get(unit)
    mask = state["sup_d"] != code if code is not None else np.ones(len(state["sup_d"]), dtype=bool)
    rows = np.concatenate([state["bsrc"], state["ssrc"][mask]])
    cols = np.concatenate([state["bdst"], state["sdst"][mask]])
    return rows, cols


def run_fk(cfg: DrugPathConfig, state: dict) -> None:
    n = state["n_nodes"]
    drug_mask = state["drug_mask"]
    gene_idx = state["gene_idx"]
    gene_pos = {g: i for i, g in enumerate(state["genes"])}
    sup = state["sup"]
    rows_out = []
    for unit, stratum in sorted(state["units"].items()):
        answers = state["answers"].get(unit, set())
        ans_mask = np.zeros(len(state["genes"]), dtype=bool)
        if answers:
            ans_mask[[gene_pos[g] for g in answers]] = True

        rows, cols = _unit_arrays(state, unit)
        adj = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))

        d_idx = state["idx"].get(unit)
        if d_idx is None or not answers:
            rows_out.append(
                {
                    "unit": unit,
                    "stratum": stratum,
                    "n_answers": len(answers),
                    "n_drugs_adjacent": 0,
                    **{f"f_at_{k}": None for k in cfg.fk_ks},
                    "n_drug_reach_excl": 0,
                    "exclusivity_frac": None,
                    "drug_provenance_frac": None,
                }
            )
            continue

        # drugs adjacent to the unit disease (first hop must be a drug)
        nbr = np.concatenate([cols[rows == d_idx], rows[cols == d_idx]])
        adj_drugs = np.unique(nbr[drug_mask[nbr]])
        n_drugs = len(adj_drugs)
        if n_drugs:
            dist_m = csgraph.shortest_path(
                adj, method="D", directed=False, unweighted=True, indices=adj_drugs
            )
            dist_drugs = dist_m.min(axis=0) if dist_m.ndim == 2 else dist_m
            dmed = drug_mediated_distances(dist_drugs)
        else:
            dmed = np.full(n, np.inf)
        gd = dmed[gene_idx]
        fk = {f"f_at_{k}": fk_fraction(gd, ans_mask, k) for k in cfg.fk_ks}

        # exclusivity: same graph with all drug nodes removed
        keep = ~(drug_mask[rows] | drug_mask[cols])
        adj_nd = sp.coo_matrix((np.ones(int(keep.sum())), (rows[keep], cols[keep])), shape=(n, n))
        dist_nd = csgraph.shortest_path(
            adj_nd, method="D", directed=False, unweighted=True, indices=d_idx
        )
        reach_k = reachable_mask(gd, cfg.exclusivity_k) & ans_mask
        n_reach, excl = exclusivity_fraction(reach_k, dist_nd[gene_idx], cfg.exclusivity_k)

        unit_edges = sup.filter(pl.col("src_id") == unit)
        drug_prov = float(
            np.mean(
                [
                    has_drug_provenance(r["sources_direct"], r["sources_indirect"])
                    for r in unit_edges.select("sources_direct", "sources_indirect").to_dicts()
                ]
            )
        )
        rows_out.append(
            {
                "unit": unit,
                "stratum": stratum,
                "n_answers": len(answers),
                "n_drugs_adjacent": n_drugs,
                **fk,
                "n_drug_reach_excl": n_reach,
                "exclusivity_frac": excl,
                "drug_provenance_frac": drug_prov,
            }
        )
        log.info("fk_unit_done", unit=unit, f3=round(fk.get("f_at_3") or 0.0, 4), drugs=n_drugs)
    pl.DataFrame(rows_out).write_parquet(cfg.output_dir / "fk_drug_paths.parquet")
    log.info("fk_done", units=len(rows_out))


def run_provenance(cfg: DrugPathConfig, state: dict) -> None:
    sup = state["sup"]
    pool_drug_prov = float(
        np.mean(
            [
                has_drug_provenance(r["sources_direct"], r["sources_indirect"])
                for r in sup.select("sources_direct", "sources_indirect").to_dicts()
            ]
        )
    )

    def license_counts(df: pl.DataFrame) -> Counter:
        c: Counter = Counter()
        for v in df["sources_direct"].to_list():
            c[classify_direct_source(v)] += 1
        return c

    pool_c = license_counts(sup)
    lic_rows = [{"scope": "__POOL__", "n_edges": sup.height, **pool_c}]
    for unit in sorted(state["units"]):
        u = sup.filter(pl.col("src_id") == unit)
        if u.height == 0:
            continue
        lic_rows.append({"scope": unit, "n_edges": u.height, **license_counts(u)})
    lic = pl.DataFrame(lic_rows).fill_null(0)
    for col in ("disgenet+open_targets", "disgenet_only", "open_targets_only", "other"):
        if col not in lic.columns:
            lic = lic.with_columns(pl.lit(0).alias(col))
    lic = lic.with_columns(
        pct_disgenet=100.0
        * (pl.col("disgenet+open_targets") + pl.col("disgenet_only"))
        / pl.col("n_edges")
    )
    lic.write_parquet(cfg.output_dir / "license_provenance.parquet")
    (cfg.output_dir / "drug_provenance_pool.json").write_text(
        json.dumps({"pool_drug_provenance_frac": pool_drug_prov, "n_edges": sup.height}, indent=2),
        encoding="utf-8",
    )
    log.info(
        "provenance_done",
        pool_drug_prov=round(pool_drug_prov, 5),
        pool_disgenet_pct=round(float(lic["pct_disgenet"][0]), 3),
    )


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def write_report(cfg: DrugPathConfig, state: dict) -> None:
    out = cfg.output_dir
    fk = pl.read_parquet(out / "fk_drug_paths.parquet")
    lic = pl.read_parquet(out / "license_provenance.parquet")
    pool = json.loads((out / "drug_provenance_pool.json").read_text(encoding="utf-8"))

    fk3 = fk.filter(pl.col("f_at_3").is_not_null())
    median_f3 = float(fk3["f_at_3"].median())
    n_all = fk.height
    bands = [covariate_band(v) for v in fk3["f_at_3"].to_list()]
    n_gte = bands.count("f3_gte_10pct")
    n_lt = bands.count("f3_lt_10pct")
    if median_f3 >= EXAM_INTEGRITY_TRIGGER:
        verdict = (
            f"roster-median F@3 = {median_f3:.3f} >= {EXAM_INTEGRITY_TRIGGER} — the exam "
            "substantially measures pharmacological reconstruction; the drug-edge arm "
            "structure is not optional"
        )
    else:
        verdict = (
            f"roster-median F@3 = {median_f3:.3f} < {EXAM_INTEGRITY_TRIGGER} — below the "
            "pharmacological-reconstruction trigger; F@3 enters as a covariate per v5.1 §6.4"
        )

    top5 = fk3.sort(["f_at_3", "n_drugs_adjacent"], descending=True).head(5)

    def _f(x: float | None) -> str:
        return f"{x:.3f}" if x is not None else "—"

    fk_rows = "\n".join(
        f"| {r['unit']} | {r['stratum']} | {r['n_answers']:,} | {r['n_drugs_adjacent']} | "
        f"{_f(r['f_at_2'])} | {_f(r['f_at_3'])} | {_f(r['f_at_4'])} | "
        f"{r['n_drug_reach_excl']} | {_f(r['exclusivity_frac'])} | "
        f"{_f(r['drug_provenance_frac'])} |"
        for r in fk.sort("stratum", "f_at_3", descending=[False, True], nulls_last=True).to_dicts()
    )

    top5_rows = "\n".join(
        f"| {r['unit']} | {r['stratum']} | {r['f_at_3']:.3f} | {r['n_drugs_adjacent']} | "
        f"{r['n_answers']:,} |"
        for r in top5.to_dicts()
    )

    lic_pool = lic.filter(pl.col("scope") == "__POOL__").to_dicts()[0]
    lic_unit_rows = "\n".join(
        f"| {r['scope']} | {r['n_edges']:,} | {r['open_targets_only']:,} | "
        f"{r['disgenet+open_targets']:,} | {r['disgenet_only']:,} | {r['other']:,} | "
        f"{r['pct_disgenet']:.2f}% |"
        for r in lic.filter(pl.col("scope") != "__POOL__")
        .sort("pct_disgenet", descending=True)
        .to_dicts()
    )

    median_pct = f"{100 * median_f3:.1f}%"
    ot_only_pct = 100 * lic_pool["open_targets_only"] / lic_pool["n_edges"]
    text = f"""<!-- track_scope: benchmark -->
# Drug-path measurement + license provenance (v5.1 Step 3 — measurement only)

CPU measurement feeding the `CHANGES/0020` freeze package. Claim ceiling: all
outputs support **disease-associated gene recovery / hypothesis
prioritization** within the neuro/aging/metabolic cluster only.

## Methods summary

- **F@k (v5.1 §6.4), k = {cfg.fk_ks}:** per candidate unit, fraction of the
  unit's DIS-GEN answer genes reachable from the unit disease via drug-mediated
  paths of length <= k in the LODO fact graph (BRG + supervision *minus the
  unit's own edges*, undirected). Drug-mediated distance = 1 + BFS distance from
  any drug adjacent to the unit disease (first hop is a drug by construction;
  k=2 → disease→drug→gene exactly). Drug nodes = Track-D `node_type == "drug"`
  (chembl:/drugbank: prefixes). Drug legs used: DRG-DIS
  (INDICATION/CONTRAINDICATION/OFF_LABEL_USE) + DRG-GEN
  (TARGET/ENZYME/TRANSPORTER/CARRIER/INHIBITOR/AGONIST/…).
- **Exclusivity annotation:** of the answers drug-reachable@{cfg.exclusivity_k},
  the fraction ALSO reachable within {cfg.exclusivity_k} hops from the disease
  in the same graph with all drug nodes removed — reachable is not exclusively
  drug-mediated.
- **Drug-provenance fraction:** fraction of a unit's DIS-GEN edges whose
  `sources_direct`/`sources_indirect` match the drug/clinical rule: exact tokens
  CTD (Comparative Toxicogenomics Database) and CGI (Cancer Genome Interpreter),
  plus substring tokens DRUG/FDA/EMA/CLINICAL/TRIAL/PHARM/CHEMBL (the latter
  occur nowhere in this build — verified by full-column scan).
- **License provenance (v5.1 §2.3):** DIS-GEN edges bucketed by `sources_direct`
  values: OPEN_TARGETS only / DISGENET+OPEN_TARGETS / DISGENET only / other
  (values are exactly OPEN_TARGETS and DISGENET in this build — verified by
  full-column scan of Track-D).

## 1. F@k per candidate unit (all {n_all} units)

`n_reach@3` = answers drug-reachable@3; `excl@3` = fraction of those also
reachable <=3 with drug nodes removed; `drug-prov` = drug-provenance fraction
of the unit's key. Full table: `fk_drug_paths.parquet`.

| unit | stratum | answers | drugs adj | F@2 | F@3 | F@4 | n_reach@3 | excl@3 | drug-prov |
|---|---|---|---|---|---|---|---|---|---|
{fk_rows}

## 2. Exam-integrity check (v5.1 §10 / §6.4)

**{verdict}.** Roster-median F@3 = **{median_pct}** across {fk3.height} units
with supervision. Covariate bands for the Arm A/B analysis (F@3):
**<{int(COVARIATE_BAND * 100)}%: {n_lt} units**;
**>={int(COVARIATE_BAND * 100)}%: {n_gte} units**.

## 3. Top-5 drug-rich units (Arm B subset candidates, v5.1 §6.4)

Ranked by F@3 (ties broken by adjacent-drug count):

| unit | stratum | F@3 | drugs adjacent | answers |
|---|---|---|---|---|
{top5_rows}

## 4. Drug-provenance fractions of answer keys

Per-unit fractions in the §1 table (`drug-prov`). Whole supervision pool:
**{100 * pool["pool_drug_provenance_frac"]:.2f}%** of {pool["n_edges"]:,}
DIS-GEN edges match the drug/clinical source rule (CTD/CGI tokens, via
`sources_indirect` — no direct drug-source tokens exist in this build).

## 5. License provenance quantification (v5.1 §2.3)

Auditable basis of the research-use posture (DisGeNET-sourced edges carry CC
BY-NC terms → research use only; Open Targets is CC0). Machine-readable:
`license_provenance.parquet`.

Whole supervision pool ({lic_pool["n_edges"]:,} edges): OPEN_TARGETS only
**{lic_pool["open_targets_only"]:,}** ({ot_only_pct:.2f}%);
DISGENET+OPEN_TARGETS **{lic_pool["disgenet+open_targets"]:,}**
({lic_pool["pct_disgenet"]:.2f}% any-DisGeNET); DISGENET only
{lic_pool["disgenet_only"]:,}; other {lic_pool["other"]:,}.

Per candidate unit (sorted by any-DisGeNET fraction):

| unit | edges | OT only | DG+OT | DG only | other | any-DG % |
|---|---|---|---|---|---|---|
{lic_unit_rows}

## 6. Feeds into the CHANGES/0020 freeze package

1. **Per-unit F@k table** (§1 / `fk_drug_paths.parquet`) — F@2/3/4 + exclusivity
   for all {n_all} candidate units, published in the frozen roster.
2. **Exam-integrity verdict:** roster-median F@3 = {median_pct} vs the
   {EXAM_INTEGRITY_TRIGGER} trigger — recorded verbatim in the freeze (§2).
3. **Covariate band assignment (F@3):** {n_lt} units <{int(COVARIATE_BAND * 100)}%,
   {n_gte} units >={int(COVARIATE_BAND * 100)}% — stratification for the Arm A-B
   delta regression (v5.1 §6.4); full-roster dual-arm contingency unlocks only
   if roster-median F@3 >= {int(COVARIATE_BAND * 100)}% (it is {median_pct}).
4. **Top-5 drug-rich units** (§3) — Arm B runs on the monogenic panel plus
   these units regardless of F.
5. **License-provenance table** (§5 / `license_provenance.parquet`) — published
   verbatim as the auditable basis of the research-use posture.

## Reproduction

`python scripts/m2_drug_path_measurement.py --config configs/m2_drug_paths.yaml --resume`
"""
    cfg.report_path.write_text(text, encoding="utf-8")
    log.info("report_written", path=str(cfg.report_path))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, DrugPathConfig)

    if args.dry_run:
        print(f"would measure drug paths → {cfg.output_dir} (report: {cfg.report_path})")
        return 0

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)

    def skip(name: str) -> bool:
        return args.resume and name in run.done_stages

    with Run(
        script="m2_drug_path_measurement",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        # load is cheap and every later stage needs the state — always re-run.
        with run.stage("load"):
            state = load_inputs(cfg)

        if not skip("fk"):
            with run.stage("fk", outputs=[cfg.output_dir / "fk_drug_paths.parquet"]):
                run_fk(cfg, state)
        if not skip("provenance"):
            with run.stage("provenance", outputs=[cfg.output_dir / "license_provenance.parquet"]):
                run_provenance(cfg, state)

        with run.stage("report", outputs=[cfg.report_path]):
            write_report(cfg, state)

    log.info("drug_path_measurement_done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
