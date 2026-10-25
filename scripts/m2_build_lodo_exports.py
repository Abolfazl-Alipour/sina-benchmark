#!/usr/bin/env python
"""m2_build_lodo_exports — per-unit LODO BioPathNet exports (v5.1 Phase F).

Implements the CHANGES/0020 freeze package: one self-contained export per
graded disease unit (29 frozen, owner decision D1) under
``data/processed/m2/lodo_exports/<unit>/``::

    train1.txt  BRG = trimmed BRG minus DIS-PHE (freeze §5); identical across
                units — all unit context edges retained (freeze §6)
    train2.txt  supervision pool minus the unit's co-holdout cluster keys
                (freeze §3) minus ancestor-strip (MONDO ancestors of held-out
                members, from Track-D node properties.ancestors) minus
                sibling-strip (overlap coefficient >= 0.3 AND intersection
                >= 20 vs the unit on the frozen matrix, plus the explicit D2
                parkinsonism-plus list, freeze §12)
    valid.txt   15% of the unit's own key, by evidence_group_id (freeze §6)
    test.txt    remaining 85%; reciprocal twins co-assigned by construction
                (splits partition group IDs, never edge rows — §5.4 contract)
    entity_types.txt / entity_names.txt / strip_list.json

Guards (all hard asserts, freeze §6 Phase-F list): no DIS-GEN in train1; the
unit's key size matches the frozen audit; expected ancestor-strips fire
(Werner -> progeria / premature aging syndrome / progeroid syndrome); valid
and test share zero evidence group IDs; zero-vector answer genes (no edge in
train1 + train2) are logged and excluded from valid/test with a coverage note;
the unit's export BRG context count equals the frozen audit value minus its
DIS-PHE edges (the audit counted DIS-PHE; the export drops it).

Usage: python scripts/m2_build_lodo_exports.py --config configs/m2_lodo_exports.yaml
           [--dry-run] [--resume]
"""

import argparse
import json
import sys
from pathlib import Path

import polars as pl
from pydantic import BaseModel

from sina.config import load_config
from sina.logging import configure_logging, get_logger
from sina.manifests import Run
from sina.splits import GROUP_COL, assert_split_disjoint, group_assignment

log = get_logger(__name__)

DIS_GEN = "DIS-GEN"


class UnitEntry(BaseModel):
    curie: str
    name: str
    tier: str
    cluster: str


class LodoExportConfig(BaseModel):
    supervision: Path
    brg: Path
    trackd_nodes: Path
    overlap_matrix: Path
    unit_audit: Path
    output_dir: Path
    report_path: Path
    manifest_path: Path
    split_seed: int = 1337
    valid_fraction: float = 0.15
    drop_brg_labels: list[str] = ["DIS-PHE"]
    sibling_coefficient_min: float = 0.3
    sibling_intersection_min: int = 20
    disclose_coefficient_min: float = 0.15
    clusters: dict[str, list[str]]
    units: list[UnitEntry]
    explicit_sibling_strips: dict[str, list[str]] = {}
    expected_ancestor_strips: dict[str, list[str]] = {}


# ---------------------------------------------------------------------------
# Pure computations (unit-tested on synthetic fixtures)
# ---------------------------------------------------------------------------


def curie_safe(curie: str) -> str:
    """Filesystem-safe unit directory name ('mondo:0010196' -> 'mondo__0010196')."""
    return curie.replace(":", "__")


def normalize_ancestor_id(raw: str) -> str:
    """Track-D properties.ancestors style ('MONDO_0019303') -> CURIE ('mondo:0019303').

    Prefix lowercased, suffix kept verbatim (OBA ids keep their 'VT' uppercase).
    Values without an underscore are returned lowercased as-is.
    """
    prefix, sep, rest = raw.partition("_")
    if not sep:
        return raw.lower()
    return f"{prefix.lower()}:{rest}"


def cluster_ancestors(
    members: list[str] | set[str],
    ancestors_by_node: dict[str, list[str]],
    node_ids: set[str],
) -> tuple[dict[str, set[str]], list[str]]:
    """Mapped MONDO/EFO/... ancestors of all held-out members.

    Returns (ancestor_curie -> source members, sorted unmapped raw values).
    Ancestors that are themselves cluster members are skipped (withheld with
    the cluster anyway).
    """
    member_set = set(members)
    mapped: dict[str, set[str]] = {}
    unmapped: set[str] = set()
    for m in sorted(member_set):
        for raw in ancestors_by_node.get(m, []):
            cur = normalize_ancestor_id(raw)
            if cur in node_ids:
                if cur not in member_set:
                    mapped.setdefault(cur, set()).add(m)
            else:
                unmapped.add(raw)
    return mapped, sorted(unmapped)


def overlap_flagged_pairs(
    unit: str,
    pair_rows: list[dict],
    coefficient_min: float,
    intersection_min: int,
    disclose_min: float,
    exclude: set[str],
) -> tuple[list[dict], list[dict]]:
    """Sibling-strip + disclose-band pairs vs one unit from the frozen matrix.

    Strip: coefficient >= coefficient_min AND intersection >= intersection_min.
    Disclose: disclose_min <= coefficient < coefficient_min AND intersection
    >= intersection_min. Self-pairs, cluster mates and ancestors are excluded
    (handled by their own mechanisms). Returns (strip, disclose), each sorted
    by curie.
    """
    strip: dict[str, dict] = {}
    disclose: dict[str, dict] = {}
    for r in pair_rows:
        a, b = r["disease_a"], r["disease_b"]
        if a != unit and b != unit:
            continue
        other = b if a == unit else a
        if other == unit or other in exclude:
            continue
        coef = float(r["overlap_coefficient"])
        inter = int(r["intersection"])
        if inter < intersection_min:
            continue
        entry = {"curie": other, "overlap_coefficient": coef, "intersection": inter}
        if coef >= coefficient_min:
            strip[other] = entry
        elif coef >= disclose_min:
            disclose[other] = entry
    return ([strip[k] for k in sorted(strip)], [disclose[k] for k in sorted(disclose)])


def split_valid_test(
    unit_edges: pl.DataFrame, valid_fraction: float, seed: int
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """15/85 valid/test split on evidence_group_id (sha1 assignment, §5.4 contract).

    Reciprocal twins share one group ID (sina.evidence_groups), so they can
    never straddle the boundary; disjointness is asserted regardless.
    """
    assignment = group_assignment(
        unit_edges[GROUP_COL], {"valid": valid_fraction, "test": 1 - valid_fraction}, seed
    )
    joined = unit_edges.join(assignment, on=GROUP_COL, how="left")
    split = {
        "valid": joined.filter(pl.col("partition") == "valid").drop("partition"),
        "test": joined.filter(pl.col("partition") == "test").drop("partition"),
    }
    assert_split_disjoint(split)
    return split["valid"], split["test"]


def zero_vector_genes(fact_endpoints: set[str], gene_ids: set[str]) -> set[str]:
    """Answer genes with no edge in the fact graph (train1 + train2)."""
    return gene_ids - fact_endpoints


def write_triples(edges: pl.DataFrame, path: Path) -> None:
    edges.select("src_id", "relation", "dst_id").write_csv(
        path, separator="\t", include_header=False
    )


def write_entity_files(nodes: pl.DataFrame, frames: list[pl.DataFrame], out: Path) -> None:
    entities: set[str] = set()
    for f in frames:
        entities |= set(f["src_id"].to_list()) | set(f["dst_id"].to_list())
    type_ids = {
        "disease": 0,
        "gene": 1,
        "drug": 2,
        "phenotype": 3,
        "biological_process": 4,
        "cellular_component": 5,
        "molecular_function": 6,
        "pathway": 7,
        "anatomy": 8,
        "exposure": 9,
    }
    subset = nodes.filter(pl.col("node_id").is_in(entities))
    with (
        (out / "entity_types.txt").open("w", encoding="utf-8") as ft,
        (out / "entity_names.txt").open("w", encoding="utf-8") as fn,
    ):
        for r in subset.select("node_id", "node_type", "node_name").to_dicts():
            ft.write(f"{r['node_id']}\t{type_ids[r['node_type']]}\n")
            fn.write(f"{r['node_id']}\t{r['node_name']}\n")


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------


def load_inputs(cfg: LodoExportConfig) -> dict:
    sup = pl.read_parquet(cfg.supervision)
    brg = pl.read_parquet(cfg.brg)
    nodes = pl.read_parquet(cfg.trackd_nodes)
    overlap = pl.read_parquet(cfg.overlap_matrix)
    audit = pl.read_parquet(cfg.unit_audit)
    assert brg.filter(pl.col("label") == DIS_GEN).height == 0, "DIS-GEN in BRG input"
    ancestors_by_node = {
        r["node_id"]: (r["properties"] or {}).get("ancestors") or []
        for r in nodes.select("node_id", "properties").to_dicts()
    }
    overlap_rows = overlap.to_dicts()
    state = {
        "sup": sup,
        "brg": brg,
        "nodes": nodes,
        "overlap_rows": overlap_rows,
        "pair_index": {
            (min(r["disease_a"], r["disease_b"]), max(r["disease_a"], r["disease_b"])): r
            for r in overlap_rows
        },
        "audit_ctx": dict(audit.select("unit", "brg_context_edges").iter_rows()),
        "audit_key": dict(audit.select("unit", "key_size").iter_rows()),
        "ancestors_by_node": ancestors_by_node,
        "node_ids": set(nodes["node_id"].to_list()),
        "names": dict(nodes.select("node_id", "node_name").iter_rows()),
        "sup_counts": dict(sup.group_by("src_id").len().iter_rows()),
    }
    log.info(
        "inputs_loaded",
        supervision=sup.height,
        brg=brg.height,
        overlap_pairs=overlap.height,
        units=len(cfg.units),
    )
    return state


def _pair_stats(state: dict, unit: str, other: str) -> dict:
    row = state["pair_index"].get((min(unit, other), max(unit, other)))
    if row is None:
        return {"overlap_coefficient": None, "intersection": None}
    return {
        "overlap_coefficient": float(row["overlap_coefficient"]),
        "intersection": int(row["intersection"]),
    }


def build_exports(cfg: LodoExportConfig, state: dict) -> None:
    sup, brg_full = state["sup"], state["brg"]
    brg = brg_full.filter(~pl.col("label").is_in(cfg.drop_brg_labels))
    assert brg.filter(pl.col("label") == DIS_GEN).height == 0
    assert brg.filter(pl.col("label").is_in(cfg.drop_brg_labels)).height == 0
    # A2: BRG-composition guard (orientation-bug class fails the build itself)
    pge = brg.filter(pl.col("label") == "PHE-GEN").height
    bpo = brg.filter(pl.col("label") == "BPO-GEN").height
    assert pge > 500_000, f"BRG PHE-GEN {pge} <= 500k — annotation layer missing (see 34c1d32)"
    assert bpo > 100_000, f"BRG BPO-GEN {bpo} <= 100k — annotation layer missing"
    brg_endpoints = set(brg["src_id"].to_list()) | set(brg["dst_id"].to_list())

    for u in cfg.units:
        mates = set(cfg.clusters[u.cluster])
        assert u.curie in mates, f"{u.curie} not in its cluster {u.cluster}"

        anc, unmapped_anc = cluster_ancestors(mates, state["ancestors_by_node"], state["node_ids"])
        expected = set(cfg.expected_ancestor_strips.get(u.curie, []))
        missing = expected - set(anc)
        assert not missing, f"{u.curie}: expected ancestor strips missing: {missing}"

        sweep, disclose = overlap_flagged_pairs(
            u.curie,
            state["overlap_rows"],
            cfg.sibling_coefficient_min,
            cfg.sibling_intersection_min,
            cfg.disclose_coefficient_min,
            exclude=mates | set(anc),
        )
        siblings = {s["curie"]: {**s, "source": "overlap_rule"} for s in sweep}
        for c in cfg.explicit_sibling_strips.get(u.curie, []):
            if c == u.curie or c in mates or c in anc:
                continue
            if c in siblings:
                siblings[c]["source"] = "overlap_rule+explicit_D2"
            else:
                siblings[c] = {
                    "curie": c,
                    **_pair_stats(state, u.curie, c),
                    "source": "explicit_D2",
                }

        strip_diseases = set(anc) | set(siblings)
        u_edges = sup.filter(pl.col("src_id") == u.curie)
        assert u_edges.height == state["audit_key"][u.curie], (
            f"{u.curie}: key size {u_edges.height} != frozen audit {state['audit_key'][u.curie]}"
        )
        train2 = sup.filter(~pl.col("src_id").is_in(mates | strip_diseases))
        assert train2.filter(pl.col("src_id").is_in(mates | strip_diseases)).height == 0

        valid, test = split_valid_test(u_edges, cfg.valid_fraction, cfg.split_seed)

        # A2: exact-tuple disjointness (group-level necessary, not sufficient)
        def _tuples(df: pl.DataFrame) -> set[tuple[str, str, str]]:
            return set(
                zip(
                    df["src_id"].to_list(),
                    df["relation"].to_list(),
                    df["dst_id"].to_list(),
                    strict=False,
                )
            )

        train2_tuples = _tuples(train2)
        valid_tuples = _tuples(valid)
        test_tuples = _tuples(test)
        assert train2_tuples.isdisjoint(valid_tuples), f"{u.curie}: train2∩valid tuple overlap"
        assert train2_tuples.isdisjoint(test_tuples), f"{u.curie}: train2∩test tuple overlap"

        # Zero-vector guard: answers with no edge in train1 + train2 cannot be
        # reached by message passing — exclude from valid/test, disclose count.
        answers = set(u_edges["dst_id"].to_list())
        suspect = answers - brg_endpoints
        if suspect:
            t2_small = train2.filter(
                pl.col("src_id").is_in(suspect) | pl.col("dst_id").is_in(suspect)
            )
            t2_ep = set(t2_small["src_id"].to_list()) | set(t2_small["dst_id"].to_list())
            zv = zero_vector_genes(brg_endpoints | t2_ep, answers)
        else:
            zv = set()
        if zv:
            valid = valid.filter(~pl.col("dst_id").is_in(zv))
            test = test.filter(~pl.col("dst_id").is_in(zv))
            log.info("zero_vector_excluded", unit=u.curie, n=len(zv), genes=sorted(zv))

        # Connectivity reconciliation: the frozen audit counted DIS-PHE; the
        # export BRG drops it (freeze §5).
        touches_u = (pl.col("src_id") == u.curie) | (pl.col("dst_id") == u.curie)
        ctx_export = brg.filter(touches_u).height
        disphe_ctx = brg_full.filter(pl.col("label").is_in(cfg.drop_brg_labels) & touches_u).height
        expected_ctx = state["audit_ctx"][u.curie] - disphe_ctx
        assert ctx_export == expected_ctx, (
            f"{u.curie}: export ctx {ctx_export} != audit "
            f"{state['audit_ctx'][u.curie]} - DIS-PHE {disphe_ctx}"
        )

        out = cfg.output_dir / curie_safe(u.curie)
        out.mkdir(parents=True, exist_ok=True)
        write_triples(brg, out / "train1.txt")
        write_triples(train2, out / "train2.txt")
        write_triples(valid, out / "valid.txt")
        write_triples(test, out / "test.txt")
        write_entity_files(state["nodes"], [brg, train2, valid, test], out)

        names, sup_counts = state["names"], state["sup_counts"]
        strip_list = {
            "unit": u.curie,
            "name": u.name,
            "tier": u.tier,
            "cluster": u.cluster,
            "cluster_mates_withheld": sorted(mates),
            "ancestors_stripped": [
                {
                    "curie": c,
                    "name": names.get(c, ""),
                    "supervision_edges": sup_counts.get(c, 0),
                    "ancestor_of": sorted(srcs),
                }
                for c, srcs in sorted(anc.items())
            ],
            "siblings_stripped": [
                {
                    **s,
                    "name": names.get(s["curie"], ""),
                    "supervision_edges": sup_counts.get(s["curie"], 0),
                }
                for s in sorted(siblings.values(), key=lambda s: s["curie"])
            ],
            "disclosed_band": [{**d, "name": names.get(d["curie"], "")} for d in disclose],
            "unmapped_ancestors": unmapped_anc,
            "zero_vector_genes_excluded": sorted(zv),
            "split": {
                "seed": cfg.split_seed,
                "valid_fraction": cfg.valid_fraction,
                "valid_edges": valid.height,
                "test_edges": test.height,
                "valid_groups": valid[GROUP_COL].n_unique(),
                "test_groups": test[GROUP_COL].n_unique(),
            },
            "counts": {
                "train1": brg.height,
                "train2": train2.height,
                "valid": valid.height,
                "test": test.height,
                "key_size": u_edges.height,
                "ctx_edges_export": ctx_export,
                "ctx_edges_audit": state["audit_ctx"][u.curie],
                "ctx_edges_disphe_removed": disphe_ctx,
                "n_ancestors_stripped": len(anc),
                "n_siblings_stripped": len(siblings),
                "n_disclosed_band": len(disclose),
                "n_zero_vector_excluded": len(zv),
            },
        }
        (out / "strip_list.json").write_text(
            json.dumps(strip_list, indent=2) + "\n", encoding="utf-8"
        )
        log.info(
            "unit_exported",
            unit=u.curie,
            tier=u.tier,
            train2=train2.height,
            valid=valid.height,
            test=test.height,
            ctx=ctx_export,
            ancestors=len(anc),
            siblings=len(siblings),
            zv=len(zv),
        )


def build_manifest(cfg: LodoExportConfig) -> None:
    rows = []
    for u in cfg.units:
        sl = json.loads(
            (cfg.output_dir / curie_safe(u.curie) / "strip_list.json").read_text(encoding="utf-8")
        )
        rows.append(
            {
                "unit": u.curie,
                "name": u.name,
                "tier": u.tier,
                "cluster": u.cluster,
                "export_dir": curie_safe(u.curie),
                **sl["counts"],
                "valid_groups": sl["split"]["valid_groups"],
                "test_groups": sl["split"]["test_groups"],
            }
        )
    manifest = pl.DataFrame(rows)
    manifest.write_parquet(cfg.output_dir / "manifest.parquet")
    log.info("manifest_written", units=manifest.height)


def write_report(cfg: LodoExportConfig) -> None:
    manifest = pl.read_parquet(cfg.output_dir / "manifest.parquet")
    strip_lists = {
        u.curie: json.loads(
            (cfg.output_dir / curie_safe(u.curie) / "strip_list.json").read_text(encoding="utf-8")
        )
        for u in cfg.units
    }

    rows = "\n".join(
        f"| {r['unit']} | {r['name']} | {r['tier']} | {r['train2']:,} | "
        f"{r['valid']:,} | {r['test']:,} | {r['ctx_edges_export']} | "
        f"{r['ctx_edges_audit']} | {r['n_ancestors_stripped']} | "
        f"{r['n_siblings_stripped']} | {r['n_zero_vector_excluded']} |"
        for r in manifest.sort("tier", "unit").to_dicts()
    )

    # D2 confirmation: the two named parkinsonism-plus diseases per PD unit.
    d2_lines = []
    for pd_unit, strips in sorted(cfg.explicit_sibling_strips.items()):
        sibs = {s["curie"]: s for s in strip_lists[pd_unit]["siblings_stripped"]}
        for c in strips:
            s = sibs.get(c)
            assert s is not None, f"D2 strip missing for {pd_unit}: {c}"
            coef = s["overlap_coefficient"]
            coef_s = f"{coef:.3f}" if coef is not None else "n/a"
            d2_lines.append(
                f"| {pd_unit} | {c} ({s['name']}) | {coef_s} | "
                f"{s['intersection']} | {s['source']} |"
            )
    d2_rows = "\n".join(d2_lines)

    werner = strip_lists["mondo:0010196"]
    trap = {a["curie"]: a for a in werner["ancestors_stripped"]}
    trap_rows = "\n".join(
        f"| mondo:0010196 | {c} ({trap[c]['name']}) | "
        f"{', '.join(trap[c]['ancestor_of'])} | {trap[c]['supervision_edges']:,} |"
        for c in sorted(cfg.expected_ancestor_strips["mondo:0010196"])
    )

    zv_total = sum(r["n_zero_vector_excluded"] for r in manifest.to_dicts())
    sib_counts = [r["n_siblings_stripped"] for r in manifest.to_dicts()]
    anc_counts = [r["n_ancestors_stripped"] for r in manifest.to_dicts()]
    disc_total = sum(r["n_disclosed_band"] for r in manifest.to_dicts())
    ctx_dropped = manifest.filter(pl.col("ctx_edges_export") < pl.col("ctx_edges_audit"))
    ctx_rows = "\n".join(
        f"| {r['unit']} | {r['name']} | {r['ctx_edges_audit']} | "
        f"{r['ctx_edges_disphe_removed']} | {r['ctx_edges_export']} |"
        for r in ctx_dropped.sort("ctx_edges_export").to_dicts()
    )
    unmapped = sorted({raw for sl in strip_lists.values() for raw in sl["unmapped_ancestors"]})
    if unmapped:
        mapping_note = (
            f"Values with no Track-D node ({len(unmapped)} distinct raws, e.g. "
            f"{', '.join(f'`{x}`' for x in unmapped[:5])}) cannot resolve to a "
            "supervision disease and are listed per unit; none carries DIS-GEN "
            "supervision in the pool."
        )
    else:
        mapping_note = (
            "Every ancestor value resolved to a Track-D node across all 29 "
            "units (0 unmapped). Mapped ancestors without their own DIS-GEN "
            "supervision are listed in `strip_list.json` but strip no edges."
        )
    train1 = manifest["train1"][0]
    keep = not cfg.drop_brg_labels
    train1_desc = (
        "the full trimmed BRG **incl. DIS-PHE**" if keep else "the trimmed BRG minus DIS-PHE"
    )
    train1_rule = "CHANGES/0024 C1 KEEP" if keep else "freeze §5 DROP"
    ctx_recon = (
        "The frozen audit's `brg_context_edges` counted DIS-PHE context and the "
        "KEEP export retains DIS-PHE (CHANGES/0024 C1), so export context equals "
        "the audit count exactly (asserted per unit — all PASS)."
        if keep
        else "The frozen audit's `brg_context_edges` counted DIS-PHE context; the "
        "export BRG drops DIS-PHE (freeze §5), so the export context count equals "
        "audit minus the unit's DIS-PHE edges (asserted per unit — all PASS)."
    )
    disphe_guard = (
        "No DIS-GEN label in train1; DIS-PHE retained (KEEP variant)."
        if keep
        else "No DIS-GEN label in train1; DIS-PHE absent from train1."
    )
    ctx_guard = (
        "Export context count = frozen audit count, per unit (DIS-PHE retained)."
        if keep
        else "Export context count = frozen audit count minus DIS-PHE, per unit."
    )

    text = f"""<!-- track_scope: benchmark -->
# LODO exports (v5.1 Phase F; CHANGES/0020 freeze package)

One self-contained BioPathNet-format export per graded disease unit under
`data/processed/m2/lodo_exports/<unit>/` (unit CURIE with `:` -> `__`).
All roster, cluster, threshold, and strip decisions are the frozen values of
CHANGES/0020; this report confirms the two owner decisions fired as adopted.

- **D1 confirmation:** exactly the {manifest.height} graded diseases of the
  three-panel roster are exported (13 primary confirmatory, 5
  secondary-evidence, 11 secondary-low-context — the last including the two
  double-secondary units). Context-only and supervision-only units have no
  export. Nothing above the 5% key ceiling is graded.
- **D2 confirmation:** Kufor-Rakeb (mondo:0011706) and PSP-parkinsonism
  (mondo:0009839) DIS-GEN edges are stripped from train2 whenever a PD unit
  is graded (table below).

## Export composition

train1.txt is {train1_desc} (**{train1:,} edges**, {train1_rule};
no DIS-GEN label anywhere — asserted per export) and is identical
across units: all unit context edges are retained (freeze §6). train2.txt is
the 191,953-edge supervision pool minus the unit's co-holdout cluster keys,
minus ancestor-strip, minus sibling-strip. valid.txt/test.txt are the unit's
own key split {cfg.valid_fraction:.0%}/{1 - cfg.valid_fraction:.0%} by
evidence group ID (seed {cfg.split_seed}; reciprocal twins share a group ID
and cannot straddle — disjointness asserted per export).

| unit | name | tier | train2 | valid | test | ctx exp | ctx aud | anc | sib | zero-vec |
|---|---|---|---|---|---|---|---|---|---|---|
{rows}

Full per-unit detail: `data/processed/m2/lodo_exports/manifest.parquet` and
each unit's `strip_list.json`.

## D2 — parkinsonism-plus sibling-strip (owner decision, adopted)

Kufor-Rakeb qualifies under the plain overlap rule (coefficient >= 0.3,
intersection >= 20) for all four PD units; PSP-parkinsonism has coefficient
1.0 but a 1-gene intersection, below the 20-gene qualifier, and is stripped
via the explicit D2 list. Both are stripped for every PD graded unit.

| graded unit | stripped disease | coefficient | intersection | source |
|---|---|---|---|---|
{d2_rows}

## Ancestor-strip guard — Werner trap terms

Grading Werner (cluster c1_progeria held out) strips the three named trap
terms via MONDO ancestors of the held-out members (progeria enters through
HGPS's ancestor list — the guard would abort the build otherwise):

| graded unit | stripped ancestor | ancestor of | supervision edges |
|---|---|---|---|
{trap_rows}

## Strip sweep summary

- Ancestor-strip (MONDO/EFO ancestors of held-out cluster members present as
  Track-D supervision diseases): per unit {min(anc_counts)}-{max(anc_counts)}
  diseases.
- Sibling-strip (overlap rule + explicit D2): per unit {min(sib_counts)}-{max(sib_counts)}
  diseases. These lists are long by construction: OptimusKG supervision is
  heavily redundant (the 0.3 coefficient with the 20-gene intersection
  qualifier still means 20+ shared answer genes per flagged disease), and the
  freeze accepts a smaller train2 in exchange for leakage control. Every
  stripped disease is published in the unit's `strip_list.json`.
- Disclose band (0.15-0.3, intersection >= 20; disclosed, NOT stripped —
  freeze §3): {disc_total} unit-disease pairs across the roster, listed per
  unit in `strip_list.json`.

## DIS-PHE connectivity reconciliation

{ctx_recon} Units whose context
drops because of the DIS-PHE removal ({ctx_dropped.height} of {manifest.height}):

| unit | name | ctx (audit) | DIS-PHE removed | ctx (export) |
|---|---|---|---|---|
{ctx_rows}

Tier floors were frozen against the DIS-PHE-inclusive counts (owner decision
D1); units whose export BRG context falls below their panel's floor are
disclosed here, not re-tiered — the roster is frozen.

## Zero-vector guard

Answer genes with no edge in train1 + train2 (isolated in the fact graph,
unreachable by message passing): **{zv_total}** across the roster
(excluded from valid/test with the gene list in the unit's `strip_list.json`).

## Ancestor-ID mapping disclosure

Track-D `properties.ancestors` values are mapped `PREFIX_NNNN` ->
`prefix:NNNN` (prefix lowercased, suffix verbatim). {mapping_note}

## Guards (all PASS at build time)

- {disphe_guard}
- Each unit's key size equals the frozen audit `key_size`.
- valid/test share zero evidence group IDs (assert_split_disjoint).
- train2 contains no edge of any cluster mate, stripped ancestor, or stripped
  sibling.
- Expected ancestor-strips fire (Werner trap terms; hard abort otherwise).
- {ctx_guard}
- D2 strips present in every PD unit's `strip_list.json` (asserted again when
  this report was written).

## Reproduction

`python scripts/m2_build_lodo_exports.py --config configs/m2_lodo_exports.yaml --resume`
"""
    cfg.report_path.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    configure_logging()
    cfg = load_config(args.config, LodoExportConfig)

    if args.dry_run:
        print(
            f"would build {len(cfg.units)} LODO exports → {cfg.output_dir} "
            f"(report: {cfg.report_path})"
        )
        for u in cfg.units:
            print(f"  {u.curie} ({u.tier}, {u.cluster}) → {cfg.output_dir / curie_safe(u.curie)}")
        return 0

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)

    with Run(
        script="m2_build_lodo_exports",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        if not (args.resume and "load" in run.done_stages):
            with run.stage("load"):
                state = load_inputs(cfg)

        if not (args.resume and "exports" in run.done_stages):
            with run.stage("exports", outputs=[cfg.output_dir]):
                build_exports(cfg, state)

        if not (args.resume and "manifest" in run.done_stages):
            with run.stage("manifest", outputs=[cfg.output_dir / "manifest.parquet"]):
                build_manifest(cfg)

        with run.stage("report", outputs=[cfg.report_path]):
            write_report(cfg)

    log.info("lodo_exports_done", units=len(cfg.units))
    return 0


if __name__ == "__main__":
    sys.exit(main())
