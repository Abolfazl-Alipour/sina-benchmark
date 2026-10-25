#!/usr/bin/env python
"""m2_lodo_audits — v5.1 Step 2 freeze-prep audits for cold LODO (plan Phase C).

Measurement only; nothing here is graded. Feeds the CHANGES/0020 freeze package:

  overlap        pairwise answer-key overlap coefficient + Jaccard (v5.1 §5.3)
  bfs_separation answer vs non-answer BFS distance separation on the LODO fact
                 graph (BRG + supervision minus the unit's own edges); the
                 theta-derivation input for the key-size ceiling (v5.1 §5.4)
  connectivity   BRG context-edge census + roster attrition draft (v5.1 §4, §6.1)
  ceiling        key size as % of the gene universe, theta bands (v5.1 §5.4)
  coherence      cluster density/degree vs seeded random same-size subsets (§2.2)
  ppr_disphe     PPR +/- DIS-PHE ablation per candidate unit (v5.1 §4)
  noise_floor    filtered Hits@10 variance vs key size (v5.1 §6.1 evidence floor)
  efo_pinning    EFO/OBA resolution pinning for every roster term (v5.1 §2.1)

Usage: python scripts/m2_lodo_audits.py --config configs/m2_lodo_audits.yaml
           [--dry-run] [--resume]
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl
import scipy.sparse as sp
import yaml
from pydantic import BaseModel
from scipy.sparse import csgraph

from sina.baselines import ppr_baseline
from sina.config import load_config
from sina.logging import configure_logging, get_logger
from sina.manifests import Run
from sina.metrics import link_prediction_metrics
from sina.obo import build_xref_map, descendant_closure, parse_obo

log = get_logger(__name__)

EVIDENCE_FLOOR_PRIMARY = 50  # v5.1 §6.1
EVIDENCE_FLOOR_SECONDARY = 20  # 20-49 → secondary panel
CONTEXT_FLOOR = 10  # v5.1 §4: >=10 BRG context edges
CO_HOLDOUT_THRESHOLD = 0.3  # v5.1 §5.3
DISCLOSE_THRESHOLD = 0.15
BFS_DEPTH_CAP = 5
HITS_K = 10


class LodoAuditConfig(BaseModel):
    supervision: Path
    brg: Path
    trackd_nodes: Path
    trackd_edges: Path
    mondo_obo: Path
    families_config: Path
    trimmed_config: Path
    seed: int = 1337
    n_coherence_nulls: int = 100
    strata_a_primary: list[str]
    strata_a_secondary: list[str]
    strata_b: list[str]
    strata_c: list[str]
    strata_d_root: str
    ppr_restart: float = 0.5
    ppr_sensitivity: list[float] = [0.15, 0.85]
    noise_floor_key_sizes: list[int]
    noise_floor_resamples: int = 1000
    output_dir: Path
    report_path: Path
    manifest_path: Path


# ---------------------------------------------------------------------------
# Pure computations (unit-tested on synthetic fixtures)
# ---------------------------------------------------------------------------


def overlap_coefficient(a: set, b: set) -> float:
    """|A∩B| / min(|A|,|B|) — containment; primary overlap statistic (v5.1 §5.3)."""
    if not a or not b:
        return 0.0
    return len(a & b) / min(len(a), len(b))


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 0.0
    union = len(a | b)
    return len(a & b) / union if union else 0.0


def overlap_band(value: float) -> str:
    """v5.1 §5.3 bands: >=0.3 co-hold-out; 0.15-0.3 disclose; below → ok."""
    if value >= CO_HOLDOUT_THRESHOLD:
        return "co_holdout"
    if value >= DISCLOSE_THRESHOLD:
        return "disclose"
    return "ok"


def key_size_band(pct: float, lo: float = 2.0, hi: float = 5.0) -> str:
    """v5.1 §5.4: key-size-% bands vs the theta working range."""
    if pct < lo:
        return "below_2pct"
    if pct <= hi:
        return "theta_range_2_5pct"
    return "above_5pct"


def evidence_class(key_size: int) -> str:
    """v5.1 §6.1: >=50 primary; 20-49 secondary panel; below → not gradeable."""
    if key_size >= EVIDENCE_FLOOR_PRIMARY:
        return "primary"
    if key_size >= EVIDENCE_FLOOR_SECONDARY:
        return "secondary"
    return "below_floor"


def separation_from_distances(
    answer_dists: np.ndarray,
    nonanswer_dists: np.ndarray,
    ks: tuple[int, ...] = (2, 3),
    cap: int = BFS_DEPTH_CAP,
) -> dict:
    """BFS separation summary: %reachable<=k for answers minus non-answers.

    Distances are hop counts (np.inf = unreachable); the cap bounds "reachable".
    Positive separation = answers systematically closer than non-answers.
    """
    out: dict[str, float] = {
        "n_answers": float(len(answer_dists)),
        "n_nonanswers": float(len(nonanswer_dists)),
        "pct_answers_within_cap": float((answer_dists <= cap).mean() * 100)
        if len(answer_dists)
        else 0.0,
        "pct_nonanswers_within_cap": float((nonanswer_dists <= cap).mean() * 100)
        if len(nonanswer_dists)
        else 0.0,
    }
    for k in ks:
        pa = float((answer_dists <= k).mean() * 100) if len(answer_dists) else 0.0
        pn = float((nonanswer_dists <= k).mean() * 100) if len(nonanswer_dists) else 0.0
        out[f"pct_answers_le{k}"] = pa
        out[f"pct_nonanswers_le{k}"] = pn
        out[f"separation_le{k}"] = pa - pn
    return out


def noise_floor_draws(n_genes: int, key_size: int, resamples: int, seed: int) -> np.ndarray:
    """R seeded resamples of a random ranking: filtered Hits@10 of a random key.

    Rank-one-answer protocol: under a uniform random permutation, one answer is
    equally likely to occupy any of the (n_genes - key_size + 1) filtered slots
    (other answers removed); a hit is <=HITS_K non-answer genes ranked above it.
    Returns one filtered Hits@10 value per resample (shape (resamples,)).
    """
    rng = np.random.default_rng(seed)
    slots = n_genes - key_size + 1
    ranks = 1 + np.floor(rng.random((resamples, key_size)) * slots)
    return (ranks <= HITS_K).mean(axis=1)


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------


def load_inputs(cfg: LodoAuditConfig) -> dict:
    sup = pl.read_parquet(cfg.supervision)
    brg = pl.read_parquet(cfg.brg)
    nodes = pl.read_parquet(cfg.trackd_nodes).select(
        "node_id", "node_type", "node_name", "mondo_id"
    )
    trimmed = yaml.safe_load(cfg.trimmed_config.read_text(encoding="utf-8"))
    fam = yaml.safe_load(cfg.families_config.read_text(encoding="utf-8"))
    terms = list(parse_obo(cfg.mondo_obo))

    trackd_diseases = set(nodes.filter(pl.col("node_type") == "disease")["node_id"].to_list())
    cluster: set[str] = set()
    for root in trimmed["cluster_roots"]:
        cluster |= {m.lower() for m in descendant_closure(terms, root)}
    cluster &= trackd_diseases
    cluster |= (
        set(trimmed["aging_traits"]) | set(trimmed["progeria"]) | set(trimmed.get("trap_terms", []))
    )
    cluster &= trackd_diseases

    genes = sorted(sup["dst_id"].unique().to_list())
    answers = {
        r["src_id"]: set(r["dst_id"])
        for r in sup.group_by("src_id").agg(pl.col("dst_id")).to_dicts()
    }
    key_sizes = {d: len(g) for d, g in answers.items()}

    # Candidate strata enumeration (v5.1 §6.2 — counts re-measured, not asserted).
    units: dict[str, str] = {}
    for u in cfg.strata_a_primary:
        units[u] = "A_primary"
    for u in cfg.strata_a_secondary:
        units[u] = "A_secondary"
    for u in cfg.strata_b:
        units[u] = "B_aging_trait"
    for u in cfg.strata_c:
        units[u] = "C_alzheimer"
    pd50 = sorted(m for m in fam["families"]["parkinson"]["members"] if key_sizes.get(m, 0) >= 50)
    for m in pd50:
        units[m] = "C_parkinson"
    d_desc = {x.lower() for x in descendant_closure(terms, cfg.strata_d_root)}
    d_units = sorted(d for d in (d_desc & cluster) if key_sizes.get(d, 0) >= EVIDENCE_FLOOR_PRIMARY)
    for d in d_units:
        units[d] = "D_diabetes"

    state = {
        "sup": sup,
        "brg": brg,
        "nodes": nodes,
        "terms": terms,
        "trimmed": trimmed,
        "cluster": cluster,
        "genes": genes,
        "answers": answers,
        "key_sizes": key_sizes,
        "units": units,
        "pd50": pd50,
        "d_units": d_units,
        "names": dict(nodes.select("node_id", "node_name").iter_rows()),
        "mondo_ids": dict(nodes.select("node_id", "mondo_id").iter_rows()),
    }
    log.info(
        "inputs_loaded",
        supervision=sup.height,
        brg=brg.height,
        cluster=len(cluster),
        genes=len(genes),
        units=len(units),
        pd50=len(pd50),
        d_units=len(d_units),
    )
    return state


def run_overlap(cfg: LodoAuditConfig, state: dict) -> None:
    diseases = sorted(state["answers"])
    n = len(diseases)
    didx = {d: i for i, d in enumerate(diseases)}
    sizes = np.array([len(state["answers"][d]) for d in diseases], dtype=np.int64)

    inter = np.zeros((n, n), dtype=np.int64)
    gene_to_dis: dict[str, list[int]] = {}
    for d, gs in state["answers"].items():
        for g in gs:
            gene_to_dis.setdefault(g, []).append(didx[d])
    for dis in gene_to_dis.values():
        ix = np.asarray(dis)
        np.add.at(inter, (ix[:, None], ix[None, :]), 1)

    ia, ib = np.triu_indices(n, k=1)
    inter_u = inter[ia, ib].astype(np.float64)
    sz_a, sz_b = sizes[ia].astype(np.float64), sizes[ib].astype(np.float64)
    oc = np.divide(inter_u, np.minimum(sz_a, sz_b))
    jac = np.divide(inter_u, sz_a + sz_b - inter_u)
    pairs = pl.DataFrame(
        {
            "disease_a": [diseases[i] for i in ia],
            "disease_b": [diseases[i] for i in ib],
            "size_a": sizes[ia],
            "size_b": sizes[ib],
            "intersection": inter[ia, ib],
            "overlap_coefficient": oc,
            "jaccard": jac,
        }
    )
    pairs.write_parquet(cfg.output_dir / "overlap_matrix.parquet")

    flags = pairs.filter(pl.col("overlap_coefficient") >= DISCLOSE_THRESHOLD)
    flags.write_parquet(cfg.output_dir / "overlap_flags.parquet")

    # Explicit risk-zone check: every pair involving a candidate unit or trap
    # term at >=0.10, plus each focus term's max-overlap partner.
    focus = set(state["units"]) | set(state["trimmed"].get("trap_terms", []))
    focus = sorted(focus)
    focus_pairs = pairs.filter(
        (pl.col("overlap_coefficient") >= 0.10)
        & (pl.col("disease_a").is_in(focus) | pl.col("disease_b").is_in(focus))
    ).sort("overlap_coefficient", descending=True)
    focus_pairs.write_parquet(cfg.output_dir / "overlap_focus_pairs.parquet")
    log.info(
        "overlap_done",
        pairs=pairs.height,
        co_holdout=flags.filter(pl.col("overlap_coefficient") >= CO_HOLDOUT_THRESHOLD).height,
        disclose=flags.filter(pl.col("overlap_coefficient") < CO_HOLDOUT_THRESHOLD).height,
        focus_pairs=focus_pairs.height,
    )


def run_bfs(cfg: LodoAuditConfig, state: dict) -> None:
    brg = state["brg"].select("src_id", "dst_id")
    sup = state["sup"].select("src_id", "dst_id")
    all_nodes = sorted(
        set(brg["src_id"].to_list())
        | set(brg["dst_id"].to_list())
        | set(sup["src_id"].to_list())
        | set(sup["dst_id"].to_list())
    )
    idx = {nd: i for i, nd in enumerate(all_nodes)}
    n = len(all_nodes)
    bsrc = brg["src_id"].replace_strict(idx).to_numpy()
    bdst = brg["dst_id"].replace_strict(idx).to_numpy()
    ssrc = sup["src_id"].replace_strict(idx).to_numpy()
    sdst = sup["dst_id"].replace_strict(idx).to_numpy()
    dcode = {d: i for i, d in enumerate(sorted(state["answers"]))}
    sup_d = sup["src_id"].replace_strict(dcode).to_numpy()

    genes = state["genes"]
    gene_idx = np.array([idx[g] for g in genes])
    gene_pos = {g: i for i, g in enumerate(genes)}

    rows = []
    for d in sorted(state["answers"]):
        mask = sup_d != dcode[d]  # LODO fact graph: unit's own DIS-GEN removed
        rows_i = np.concatenate([bsrc, ssrc[mask]])
        cols_i = np.concatenate([bdst, sdst[mask]])
        adj = sp.coo_matrix((np.ones(len(rows_i)), (rows_i, cols_i)), shape=(n, n))
        dist = csgraph.shortest_path(
            adj, method="D", directed=False, unweighted=True, indices=idx[d]
        )
        gd = dist[gene_idx]
        ans_mask = np.zeros(len(genes), dtype=bool)
        ans_mask[[gene_pos[g] for g in state["answers"][d]]] = True
        stats = separation_from_distances(gd[ans_mask], gd[~ans_mask])
        rows.append(
            {
                "disease": d,
                "name": state["names"].get(d, ""),
                "key_size": state["key_sizes"][d],
                "key_pct": 100.0 * state["key_sizes"][d] / len(genes),
                **stats,
            }
        )
    pl.DataFrame(rows).write_parquet(cfg.output_dir / "bfs_separation.parquet")
    log.info("bfs_done", diseases=len(rows))


def run_connectivity(cfg: LodoAuditConfig, state: dict) -> None:
    cluster = sorted(state["cluster"])
    brg = state["brg"]
    counts: Counter = Counter()
    for col in ("src_id", "dst_id"):
        counts.update(brg.filter(pl.col(col).is_in(cluster))[col].to_list())
    n_genes = len(state["genes"])
    rows = []
    for d in sorted(cluster):
        ks = state["key_sizes"].get(d, 0)
        ctx = counts.get(d, 0)
        pct = 100.0 * ks / n_genes
        reasons = []
        if ctx < CONTEXT_FLOOR:
            reasons.append("below_context_floor")
        ev = evidence_class(ks)
        if ev == "below_floor":
            reasons.append("below_evidence_floor")
        elif ev == "secondary":
            reasons.append("secondary_panel_20_49")
        band = key_size_band(pct)
        if band == "above_5pct":
            reasons.append("above_key_size_ceiling")
        rows.append(
            {
                "disease": d,
                "name": state["names"].get(d, ""),
                "key_size": ks,
                "key_pct": pct,
                "brg_context_edges": ctx,
                "context_floor_ok": ctx >= CONTEXT_FLOOR,
                "evidence_class": ev,
                "key_size_band": band,
                "is_candidate_unit": d in state["units"],
                "attrition_reasons": ",".join(reasons),
            }
        )
    pl.DataFrame(rows).write_parquet(cfg.output_dir / "connectivity_attrition.parquet")
    log.info("connectivity_done", cluster_diseases=len(rows))


def run_ceiling(cfg: LodoAuditConfig, state: dict) -> None:
    n_genes = len(state["genes"])
    rows = [
        {
            "disease": d,
            "name": state["names"].get(d, ""),
            "key_size": ks,
            "key_pct": 100.0 * ks / n_genes,
            "band": key_size_band(100.0 * ks / n_genes),
        }
        for d, ks in sorted(state["key_sizes"].items())
    ]
    pl.DataFrame(rows).sort("key_size", descending=True).write_parquet(
        cfg.output_dir / "key_size_ceiling.parquet"
    )
    log.info("ceiling_done", diseases=len(rows))


def run_coherence(cfg: LodoAuditConfig, state: dict) -> None:
    disdis = (
        pl.scan_parquet(cfg.trackd_edges)
        .filter(pl.col("label") == "DIS-DIS")
        .select("src_id", "dst_id")
        .collect()
    )
    deg = (
        pl.scan_parquet(cfg.trackd_edges)
        .filter(pl.col("label") == "DIS-GEN")
        .select("src_id")
        .group_by("src_id")
        .len()
        .collect()
    )
    all_diseases = sorted(
        state["nodes"].filter(pl.col("node_type") == "disease")["node_id"].to_list()
    )
    didx = {d: i for i, d in enumerate(all_diseases)}
    dd = disdis.filter(pl.col("src_id").is_in(all_diseases) & pl.col("dst_id").is_in(all_diseases))
    dsrc = dd["src_id"].replace_strict(didx).to_numpy()
    ddst = dd["dst_id"].replace_strict(didx).to_numpy()
    deg_map = dict(deg.iter_rows())

    def subset_stats(subset: list[str]) -> tuple[float, float]:
        mask = np.zeros(len(all_diseases), dtype=bool)
        mask[[didx[d] for d in subset]] = True
        n_edges = int((mask[dsrc] & mask[ddst]).sum())
        pairs = len(subset) * (len(subset) - 1) / 2
        mean_deg = float(np.mean([deg_map.get(d, 0) for d in subset]))
        return n_edges / pairs, mean_deg

    cluster = sorted(state["cluster"])
    cl_density, cl_mean_deg = subset_stats(cluster)
    rng = np.random.default_rng(cfg.seed)
    null_density, null_mean_deg = [], []
    for _ in range(cfg.n_coherence_nulls):
        pick = rng.choice(all_diseases, size=len(cluster), replace=False).tolist()
        d_, m_ = subset_stats(pick)
        null_density.append(d_)
        null_mean_deg.append(m_)
    result = {
        "cluster_size": len(cluster),
        "cluster_disdis_density": cl_density,
        "cluster_mean_supervision_degree": cl_mean_deg,
        "null_density_mean": float(np.mean(null_density)),
        "null_density_sd": float(np.std(null_density, ddof=1)),
        "null_mean_degree_mean": float(np.mean(null_mean_deg)),
        "null_mean_degree_sd": float(np.std(null_mean_deg, ddof=1)),
        "density_percentile": float((np.asarray(null_density) <= cl_density).mean() * 100),
        "mean_degree_percentile": float((np.asarray(null_mean_deg) <= cl_mean_deg).mean() * 100),
        "n_nulls": cfg.n_coherence_nulls,
        "seed": cfg.seed,
    }
    (cfg.output_dir / "coherence.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    log.info("coherence_done", **{k: v for k, v in result.items() if "percentile" in k})


def run_ppr_disphe(cfg: LodoAuditConfig, state: dict) -> None:
    sup = state["sup"]
    brg = state["brg"].select("src_id", "dst_id")
    brg_no_phe = state["brg"].filter(pl.col("label") != "DIS-PHE").select("src_id", "dst_id")
    genes = state["genes"]
    filter_edges = sup.select("src_id", "dst_id")
    alphas = [cfg.ppr_restart, *cfg.ppr_sensitivity]
    rows = []
    for unit, stratum in sorted(state["units"].items()):
        n_answers = state["key_sizes"].get(unit, 0)
        if n_answers == 0:
            rows.append(
                {
                    "unit": unit,
                    "stratum": stratum,
                    "alpha": None,
                    "disphe": "skipped",
                    "filtered_mrr": None,
                    "filtered_hits10": None,
                    "n_answers": 0,
                    "degenerate": None,
                }
            )
            continue
        train = sup.filter(pl.col("src_id") != unit).select("src_id", "dst_id")
        test = sup.filter(pl.col("src_id") == unit).select("src_id", "dst_id")
        for variant, brg_part in (("with", brg), ("without", brg_no_phe)):
            graph = pl.concat([brg_part, train])
            for alpha in alphas:
                scores = ppr_baseline(graph, [unit], genes, alpha=alpha).rename(
                    {"ppr_score": "score"}
                )
                # degenerate: every gene scores 0 → total tie → the ranking
                # carries no information (metric protocol reads ties as rank 1)
                degenerate = bool(scores.height == 0 or scores["score"].max() == 0.0)
                m = link_prediction_metrics(scores, test, filter_edges=filter_edges, ks=(HITS_K,))
                rows.append(
                    {
                        "unit": unit,
                        "stratum": stratum,
                        "alpha": alpha,
                        "disphe": variant,
                        "filtered_mrr": m["filtered"]["mrr"],
                        "filtered_hits10": m["filtered"]["hits@10"],
                        "n_answers": m["n_queries"],
                        "degenerate": degenerate,
                    }
                )
        log.info("ppr_unit_done", unit=unit)
    pl.DataFrame(rows).write_parquet(cfg.output_dir / "ppr_disphe.parquet")
    log.info("ppr_done", rows=len(rows))


def run_noise_floor(cfg: LodoAuditConfig, state: dict) -> None:
    n_genes = len(state["genes"])
    rows = []
    for k in cfg.noise_floor_key_sizes:
        draws = noise_floor_draws(n_genes, k, cfg.noise_floor_resamples, cfg.seed)
        mean = float(draws.mean())
        sd = float(draws.std(ddof=1))
        rows.append(
            {
                "key_size": k,
                "n_genes": n_genes,
                "resamples": cfg.noise_floor_resamples,
                "seed": cfg.seed,
                "mean_hits10": mean,
                "sd_hits10": sd,
                "cv_hits10": sd / mean if mean else float("inf"),
                "p_zero_hits": float((draws == 0).mean()),
            }
        )
    pl.DataFrame(rows).write_parquet(cfg.output_dir / "noise_floor.parquet")
    log.info("noise_floor_done", key_sizes=len(rows))


def run_efo_pinning(cfg: LodoAuditConfig, state: dict) -> None:
    trimmed = state["trimmed"]
    roster: list[tuple[str, str, str]] = []
    for t in trimmed["aging_traits"]:
        roster.append((t, "aging_traits", "explicit_id_list_census"))
    for t in trimmed["progeria"]:
        roster.append((t, "progeria", "explicit_id_list_census"))
    for t in trimmed.get("trap_terms", []):
        roster.append((t, "trap_terms", "explicit_id_list_census"))
    for t in cfg.strata_c:
        roster.append((t, "strata_c_alzheimer", "mondo_closure_family_list"))
    for t in state["pd50"]:
        roster.append((t, "strata_c_parkinson", "mondo_closure_family_list"))
    for t in state["d_units"]:
        roster.append((t, "strata_d_diabetes", "mondo_closure_descendants"))

    xref = build_xref_map(state["terms"], {"EFO", "OBA", "Orphanet"})
    xref_map: dict[str, list[str]] = {}
    for r in xref.select("external_id", "mondo_id").to_dicts():
        xref_map.setdefault(r["external_id"].lower(), []).append(r["mondo_id"].lower())

    node_ids = set(state["names"])
    mondo_in_trackd = {m for m in state["mondo_ids"].values() if m}
    seen: set[str] = set()
    rows = []
    for term, source, method in roster:
        if term in seen:
            continue
        seen.add(term)
        present = term in node_ids
        xrefs = xref_map.get(term, [])
        xref_in_trackd = [m for m in xrefs if m in mondo_in_trackd]
        prefix = term.split(":", 1)[0]
        note = ""
        if prefix in {"efo", "oba", "orphanet"}:
            if xref_in_trackd:
                note = (
                    f"MONDO near-equivalent(s) {','.join(xref_in_trackd)} present in Track-D; "
                    "census name-search could have missed near-synonyms keyed under those "
                    "nodes — bounded by this xref check"
                )
            elif xrefs:
                note = (
                    "MONDO xref exists but not mapped into Track-D; near-synonym miss "
                    "possible, bounded to the unmapped xref"
                )
            else:
                note = "no MONDO xref in crosswalk; near-synonym miss possible, unbounded"
        rows.append(
            {
                "term": term,
                "roster_source": source,
                "resolution_method": method,
                "trackd_node_id": term if present else "",
                "label": state["names"].get(term, ""),
                "mondo_id": state["mondo_ids"].get(term) or "",
                "mondo_xrefs": ",".join(xrefs),
                "xref_mondo_in_trackd": ",".join(xref_in_trackd),
                "present_in_trackd": present,
                "near_synonym_note": note,
            }
        )
    pl.DataFrame(rows).write_parquet(cfg.output_dir / "efo_pinning.parquet")
    log.info("efo_pinning_done", terms=len(rows))


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def _fmt_pct(x: float) -> str:
    return f"{x:.2f}%"


def _name(state: dict, d: str) -> str:
    n = state["names"].get(d, "")
    return f"{d} ({n})" if n else d


def ppr_pivot(ppr: pl.DataFrame, restart: float) -> tuple[pl.DataFrame, int]:
    """Restart-row with/without DIS-PHE pivot; degenerate all-tie pairs excluded.

    A unit enters the delta only if BOTH variants carry ranking information
    (degenerate rows — every gene scoring 0, a total tie — are excluded).
    """
    ppr0 = ppr.filter((pl.col("alpha") == restart) & (pl.col("n_answers") > 0))
    n_degenerate = ppr0.filter(pl.col("degenerate")).height
    nd = ppr0.filter(~pl.col("degenerate"))
    both = nd.group_by("unit").len().filter(pl.col("len") == 2)["unit"].to_list()
    piv = nd.filter(pl.col("unit").is_in(both)).pivot(
        on="disphe", index="unit", values=["filtered_mrr", "filtered_hits10"]
    )
    return piv, n_degenerate


def build_unit_audit(cfg: LodoAuditConfig, state: dict) -> pl.DataFrame:
    ceiling = pl.read_parquet(cfg.output_dir / "key_size_ceiling.parquet")
    conn = pl.read_parquet(cfg.output_dir / "connectivity_attrition.parquet")
    bfs = pl.read_parquet(cfg.output_dir / "bfs_separation.parquet")
    ppr = pl.read_parquet(cfg.output_dir / "ppr_disphe.parquet")
    overlap = pl.read_parquet(cfg.output_dir / "overlap_matrix.parquet")

    units = state["units"]
    unit_ids = sorted(units)
    # max overlap partner per unit
    ov = overlap.filter(pl.col("disease_a").is_in(unit_ids) | pl.col("disease_b").is_in(unit_ids))
    ov_a = ov.filter(pl.col("disease_a").is_in(unit_ids)).select(
        unit=pl.col("disease_a"),
        partner=pl.col("disease_b"),
        oc=pl.col("overlap_coefficient"),
    )
    ov_b = ov.filter(pl.col("disease_b").is_in(unit_ids)).select(
        unit=pl.col("disease_b"),
        partner=pl.col("disease_a"),
        oc=pl.col("overlap_coefficient"),
    )
    ov_long = pl.concat([ov_a, ov_b])
    ov_max = ov_long.group_by("unit").agg(pl.all().sort_by("oc", descending=True).first())
    # max overlap restricted to partners below the >5% ceiling band — the
    # co-holdout-relevant statistic (mega-key partners are never graded)
    below_ceiling = ceiling.filter(pl.col("band") != "above_5pct")["disease"].to_list()
    ov_bc_max = (
        ov_long.filter(pl.col("partner").is_in(below_ceiling))
        .group_by("unit")
        .agg(pl.all().sort_by("oc", descending=True).first())
        .rename({"partner": "partner_bc", "oc": "oc_bc"})
    )

    ppr_piv, _ = ppr_pivot(ppr, cfg.ppr_restart)

    rows = []
    for unit, stratum in sorted(units.items()):
        c = ceiling.filter(pl.col("disease") == unit)
        cn = conn.filter(pl.col("disease") == unit)
        b = bfs.filter(pl.col("disease") == unit)
        o = ov_max.filter(pl.col("unit") == unit)
        obc = ov_bc_max.filter(pl.col("unit") == unit)
        p = ppr_piv.filter(pl.col("unit") == unit)
        row = {
            "unit": unit,
            "stratum": stratum,
            "name": state["names"].get(unit, ""),
            "key_size": c["key_size"][0] if c.height else 0,
            "key_pct": c["key_pct"][0] if c.height else 0.0,
            "key_size_band": c["band"][0] if c.height else "no_supervision",
            "brg_context_edges": cn["brg_context_edges"][0] if cn.height else 0,
            "evidence_class": cn["evidence_class"][0] if cn.height else "below_floor",
            "attrition_reasons": cn["attrition_reasons"][0] if cn.height else "",
            "separation_le2": b["separation_le2"][0] if b.height else None,
            "separation_le3": b["separation_le3"][0] if b.height else None,
            "pct_answers_le3": b["pct_answers_le3"][0] if b.height else None,
            "pct_nonanswers_le3": b["pct_nonanswers_le3"][0] if b.height else None,
            "max_overlap_partner": o["partner"][0] if o.height else "",
            "max_overlap_coefficient": o["oc"][0] if o.height else 0.0,
            "overlap_band": overlap_band(o["oc"][0]) if o.height else "ok",
            "max_overlap_bc_partner": obc["partner_bc"][0] if obc.height else "",
            "max_overlap_bc_coefficient": obc["oc_bc"][0] if obc.height else 0.0,
            "overlap_bc_band": overlap_band(obc["oc_bc"][0]) if obc.height else "ok",
            "ppr_mrr_with": p["filtered_mrr_with"][0] if p.height else None,
            "ppr_mrr_without": p["filtered_mrr_without"][0] if p.height else None,
            "ppr_hits10_with": p["filtered_hits10_with"][0] if p.height else None,
            "ppr_hits10_without": p["filtered_hits10_without"][0] if p.height else None,
        }
        row["ppr_mrr_delta"] = (
            row["ppr_mrr_with"] - row["ppr_mrr_without"]
            if row["ppr_mrr_with"] is not None and row["ppr_mrr_without"] is not None
            else None
        )
        row["ppr_hits10_delta"] = (
            row["ppr_hits10_with"] - row["ppr_hits10_without"]
            if row["ppr_hits10_with"] is not None and row["ppr_hits10_without"] is not None
            else None
        )
        rows.append(row)
    audit = pl.DataFrame(rows)
    audit.write_parquet(cfg.output_dir / "unit_audit.parquet")
    return audit


def write_report(cfg: LodoAuditConfig, state: dict) -> None:
    out = cfg.output_dir
    conn = pl.read_parquet(out / "connectivity_attrition.parquet")
    ceiling = pl.read_parquet(out / "key_size_ceiling.parquet")
    bfs = pl.read_parquet(out / "bfs_separation.parquet")
    overlap = pl.read_parquet(out / "overlap_matrix.parquet")
    focus = pl.read_parquet(out / "overlap_focus_pairs.parquet")
    ppr = pl.read_parquet(out / "ppr_disphe.parquet")
    noise = pl.read_parquet(out / "noise_floor.parquet").sort("key_size")
    efo = pl.read_parquet(out / "efo_pinning.parquet")
    coherence = json.loads((out / "coherence.json").read_text(encoding="utf-8"))
    audit = build_unit_audit(cfg, state)

    n_genes = len(state["genes"])
    co = overlap.filter(pl.col("overlap_coefficient") >= CO_HOLDOUT_THRESHOLD)
    disc = overlap.filter(
        (pl.col("overlap_coefficient") >= DISCLOSE_THRESHOLD)
        & (pl.col("overlap_coefficient") < CO_HOLDOUT_THRESHOLD)
    )
    # co-holdout pairs that matter for grading: both keys below the >5% ceiling
    below_ceiling = ceiling.filter(pl.col("band") != "above_5pct")["disease"].to_list()
    co_bc = co.filter(
        pl.col("disease_a").is_in(below_ceiling) & pl.col("disease_b").is_in(below_ceiling)
    )
    disc_bc = disc.filter(
        pl.col("disease_a").is_in(below_ceiling) & pl.col("disease_b").is_in(below_ceiling)
    )

    # --- attrition summary over candidate units
    units_conn = conn.filter(pl.col("is_candidate_unit"))
    n_fail_ctx = units_conn.filter(~pl.col("context_floor_ok")).height
    n_below_ev = units_conn.filter(pl.col("evidence_class") == "below_floor").height
    n_secondary = units_conn.filter(pl.col("evidence_class") == "secondary").height
    n_above_ceiling = units_conn.filter(pl.col("key_size_band") == "above_5pct").height

    # --- PPR verdict at pinned restart (degenerate all-tie pairs excluded)
    ppr0 = ppr.filter((pl.col("alpha") == cfg.ppr_restart) & (pl.col("n_answers") > 0))
    piv, n_degenerate = ppr_pivot(ppr, cfg.ppr_restart)
    piv = piv.join(ppr0.select("unit", "stratum").unique(), on="unit", how="left")
    d_mrr = float((piv["filtered_mrr_with"] - piv["filtered_mrr_without"]).mean())
    d_h10 = float((piv["filtered_hits10_with"] - piv["filtered_hits10_without"]).mean())
    verdict = (
        "KEEP DIS-PHE in the BRG" if d_mrr >= 0 and d_h10 >= 0 else "DROP DIS-PHE from the BRG"
    )

    # --- tables
    def audit_rows(df: pl.DataFrame) -> str:
        lines = []
        for r in df.sort("stratum", "key_size", descending=[False, True]).to_dicts():
            sep3 = f"{r['separation_le3']:+.1f}" if r["separation_le3"] is not None else "—"
            mrr_d = f"{r['ppr_mrr_delta']:+.4f}" if r["ppr_mrr_delta"] is not None else "—"
            lines.append(
                f"| {r['unit']} | {r['stratum']} | {r['key_size']:,} | "
                f"{r['key_pct']:.2f}% | {r['key_size_band']} | {r['brg_context_edges']} | "
                f"{r['evidence_class']} | {sep3} | "
                f"{r['max_overlap_bc_coefficient']:.3f} | {mrr_d} |"
            )
        return "\n".join(lines)

    focus_rows = "\n".join(
        f"| {_name(state, r['disease_a'])} | {_name(state, r['disease_b'])} | "
        f"{r['overlap_coefficient']:.3f} | {r['jaccard']:.3f} | {r['intersection']:,} |"
        for r in focus.head(30).to_dicts()
    )

    traps = ["mondo:0020732", "mondo:0019303", "mondo:0015333", "efo:0022597", "mondo:0004975"]
    trap_rows = "\n".join(
        f"| {_name(state, r['disease'])} | {r['key_size']:,} | {r['key_pct']:.2f}% | {r['band']} |"
        for r in ceiling.filter(pl.col("disease").is_in(traps)).to_dicts()
    )

    theta = bfs.join(ceiling.select("disease", "band"), on="disease", how="left").sort(
        "key_pct", descending=True
    )
    theta_rows = "\n".join(
        f"| {_name(state, r['disease'])} | {r['key_size']:,} | {r['key_pct']:.2f}% | "
        f"{r['pct_answers_le3']:.1f}% | {r['pct_nonanswers_le3']:.1f}% | "
        f"{r['separation_le3']:+.1f} | {r['separation_le2']:+.1f} |"
        for r in theta.filter(pl.col("key_pct") >= 2.0).head(20).to_dicts()
    )
    sep_by_band = (
        theta.group_by("band")
        .agg(
            pl.len().alias("n"),
            pl.col("separation_le3").mean().alias("mean_sep3"),
            pl.col("separation_le3").median().alias("median_sep3"),
        )
        .sort("band")
    )
    sep_band_rows = "\n".join(
        f"| {r['band']} | {r['n']} | {r['mean_sep3']:+.1f} | {r['median_sep3']:+.1f} |"
        for r in sep_by_band.to_dicts()
    )

    noise_rows = "\n".join(
        f"| {r['key_size']} | {r['mean_hits10']:.5f} | {r['sd_hits10']:.5f} | "
        f"{r['cv_hits10']:.2f} | {r['p_zero_hits']:.3f} |"
        for r in noise.to_dicts()
    )

    efo_rows = "\n".join(
        f"| {r['term']} | {r['label']} | {r['mondo_id'] or '—'} | {r['resolution_method']} | "
        f"{'yes' if r['present_in_trackd'] else 'NO'} | {r['near_synonym_note'] or '—'} |"
        for r in efo.to_dicts()
    )

    def _f4(x: float | None) -> str:
        return f"{x:.4f}" if x is not None else "isolated"

    ppr_rows = "\n".join(
        f"| {r['unit']} | {r['stratum']} | {r['alpha']} | "
        f"{_f4(r['filtered_mrr_with'])} | {_f4(r['filtered_mrr_without'])} | "
        f"{_f4(r['filtered_hits10_with'])} | {_f4(r['filtered_hits10_without'])} |"
        for r in piv.with_columns(pl.lit(cfg.ppr_restart).alias("alpha")).to_dicts()
    )

    n_co_pairs = co.height
    n_disc_pairs = disc.height
    top_co = co.sort("overlap_coefficient", descending=True).head(15)
    co_rows = "\n".join(
        f"| {_name(state, r['disease_a'])} | {_name(state, r['disease_b'])} | "
        f"{r['overlap_coefficient']:.3f} | {r['jaccard']:.3f} |"
        for r in top_co.to_dicts()
    )

    band_counts = "; ".join(
        f"**{r['band']}**: {r['len']}"
        for r in ceiling.group_by("band").len().sort("band").to_dicts()
    )
    density_row = (
        f"| DIS-DIS density (edges/possible pairs) | "
        f"{coherence['cluster_disdis_density']:.6f} | "
        f"{coherence['null_density_mean']:.6f} ± {coherence['null_density_sd']:.6f} | "
        f"{coherence['density_percentile']:.0f}th |"
    )
    degree_row = (
        f"| mean supervision degree | {coherence['cluster_mean_supervision_degree']:.1f} | "
        f"{coherence['null_mean_degree_mean']:.1f} ± "
        f"{coherence['null_mean_degree_sd']:.1f} | "
        f"{coherence['mean_degree_percentile']:.0f}th |"
    )
    verdict_word = "net-positive for" if d_mrr >= 0 and d_h10 >= 0 else "net-negative for"

    text = f"""<!-- track_scope: benchmark -->
# LODO freeze-prep audits (v5.1 Step 2 — measurement only, nothing graded)

CPU audits feeding the `CHANGES/0020` freeze package. Claim ceiling: all
outputs here support **disease-associated gene recovery / hypothesis
prioritization** within the neuro/aging/metabolic cluster only.

**Inputs (re-measured, not asserted):** {state["sup"].height:,} DIS-GEN
supervision edges over {len(state["answers"])} diseases; {state["brg"].height:,}
BRG edges; cluster = {len(state["cluster"])} diseases; candidate gene universe
= **{n_genes:,}** genes; candidate units enumerated = {len(state["units"])}
(A {len([u for u, s in state["units"].items() if s.startswith("A")])},
B {len([u for u, s in state["units"].items() if s.startswith("B")])},
C {len([u for u, s in state["units"].items() if s.startswith("C")])},
D {len([u for u, s in state["units"].items() if s.startswith("D")])}).

## Methods summary

- **Overlap audit (v5.1 §5.3):** per cluster disease, its supervision gene set;
  pairwise overlap coefficient |A∩B|/min(|A|,|B|) (primary — containment;
  Jaccard alongside). Full matrix: `overlap_matrix.parquet`
  ({overlap.height:,} pairs). Bands: ≥{CO_HOLDOUT_THRESHOLD} co-hold-out,
  {DISCLOSE_THRESHOLD}-{CO_HOLDOUT_THRESHOLD} disclose.
- **BFS separation (v5.1 §5.4 θ input):** undirected BFS (depth cap {BFS_DEPTH_CAP})
  from each cluster disease on the LODO fact graph (BRG + supervision *minus the
  unit's own DIS-GEN edges*); % of answer vs non-answer genes reachable ≤2/≤3
  hops; separation = difference. Computed for all {bfs.height} supervised
  cluster diseases: `bfs_separation.parquet`.
- **Connectivity census + attrition (v5.1 §4, §6.1):** BRG context edges per
  cluster disease; floors: ≥{CONTEXT_FLOOR} context edges, ≥{EVIDENCE_FLOOR_PRIMARY}
  supervision edges (primary) / {EVIDENCE_FLOOR_SECONDARY}-{EVIDENCE_FLOOR_PRIMARY - 1}
  (secondary panel). `connectivity_attrition.parquet`.
- **Key-size ceiling (v5.1 §5.4):** key size / {n_genes:,} as %, bands
  <2% / 2-5% (θ working range) / >5%. `key_size_ceiling.parquet`.
- **Coherence control (v5.1 §2.2):** intra-cluster DIS-DIS density (edges among
  cluster diseases / possible pairs) and mean supervision degree vs
  {coherence["n_nulls"]} seeded random same-size Track-D disease subsets
  (seed {coherence["seed"]}). `coherence.json`.
- **PPR ± DIS-PHE (v5.1 §4):** per candidate unit, `sina.baselines.ppr_baseline`
  on BRG + training DIS-GEN (supervision minus the unit's own edges), restart
  {cfg.ppr_restart} pinned, sensitivity {cfg.ppr_sensitivity}, symmetric
  normalization; scored with `sina.metrics.link_prediction_metrics` over the
  {n_genes:,}-gene universe, filtered by the unit's own answers.
  `ppr_disphe.parquet`.
- **Noise-floor simulation (v5.1 §6.1):** {cfg.noise_floor_resamples} seeded
  resamples per key size of a random ranking of {n_genes:,} genes
  (rank-one-answer protocol); distribution of filtered Hits@{HITS_K}.
  `noise_floor.parquet`.
- **EFO/OBA pinning (v5.1 §2.1):** per roster term — Track-D node, label,
  MONDO mapping, resolution method, bounded near-synonym note.
  `efo_pinning.parquet`.

## 1. Overlap audit — risk zones

- Pairs ≥ {CO_HOLDOUT_THRESHOLD} (co-hold-out): **{n_co_pairs:,}** of {overlap.height:,}
- Pairs {DISCLOSE_THRESHOLD}-{CO_HOLDOUT_THRESHOLD} (disclose): **{n_disc_pairs:,}**
- Most ≥{CO_HOLDOUT_THRESHOLD} pairs are containment artifacts of above-ceiling
  mega keys (a small key nested in a 30-70%-of-genome key scores 1.0 by
  construction; those mega keys are never graded). Restricted to pairs where
  **both** keys are below the >5% ceiling band — the co-holdout-relevant set:
  **{co_bc.height:,}** co-hold-out pairs and **{disc_bc.height:,}** disclose pairs.
  Dual artifact, also by containment: pairs whose smaller key is tiny (a
  1-10-gene key sharing one gene scores up to 1.0). The matrix carries
  `size_a`/`size_b`/`intersection` so the freeze can apply a minimum-key-size
  or minimum-intersection qualifier when reading co-holdout pairs.

Top co-hold-out-band pairs:

| disease A | disease B | overlap coeff | Jaccard |
|---|---|---|---|
{co_rows}

Risk-zone pairs involving candidate units / trap terms (overlap ≥ 0.10, top 30;
full table in `overlap_focus_pairs.parquet`):

| disease A | disease B | overlap coeff | Jaccard | shared genes |
|---|---|---|---|---|
{focus_rows}

## 2. BFS separation — the θ-derivation input

Separation vs key-size-% across all supervised cluster diseases
(`bfs_separation.parquet`). Mean/median separation@3 by ceiling band:

| key-size band | n diseases | mean sep@3 | median sep@3 |
|---|---|---|---|
{sep_band_rows}

Diseases at ≥2% of the gene universe (top 20) — separation collapses where the
key covers a large fraction of the genome:

| disease | key size | key % | ans ≤3h | non-ans ≤3h | sep@3 | sep@2 |
|---|---|---|---|---|---|---|
{theta_rows}

## 3. Key-size ceiling classification (θ working range 2-5%)

Named reference terms (never graded — trap parents and ceiling references):

| disease | key size | key % of {n_genes:,} | band |
|---|---|---|---|
{trap_rows}

Band counts across {ceiling.height} supervised cluster diseases: {band_counts}.

## 4. Connectivity census + roster attrition (draft)

Across the {units_conn.height} candidate units: failing the
≥{CONTEXT_FLOOR}-context-edge floor: **{n_fail_ctx}**; below the evidence floor
(<{EVIDENCE_FLOOR_SECONDARY} supervision edges): **{n_below_ev}**; secondary
panel ({EVIDENCE_FLOOR_SECONDARY}-{EVIDENCE_FLOOR_PRIMARY - 1}): **{n_secondary}**;
above the key-size ceiling (>5%): **{n_above_ceiling}**.
Cluster-wide (all {conn.height} diseases incl. non-candidates): failing context
floor {conn.filter(~pl.col("context_floor_ok")).height}, below evidence floor
{conn.filter(pl.col("evidence_class") == "below_floor").height}, secondary
{conn.filter(pl.col("evidence_class") == "secondary").height}, above ceiling
{conn.filter(pl.col("key_size_band") == "above_5pct").height}.

## 5. Per-unit audit table (all candidate units)

Full machine-readable table: `unit_audit.parquet`. Δ MRR = filtered MRR with
DIS-PHE minus without (restart {cfg.ppr_restart}; degenerate all-tie pairs
excluded). "max overlap" is restricted to partners below the >5% ceiling band
(the co-holdout-relevant statistic).

| unit | stratum | key | key % | ceiling band | ctx edges | evidence | sep@3 | max overlap | Δ MRR |
|---|---|---|---|---|---|---|---|---|---|
{audit_rows(audit)}

## 6. Coherence control statistic (v5.1 §2.2)

| statistic | cluster | null mean ± SD (n={coherence["n_nulls"]}) | percentile |
|---|---|---|---|
{density_row}
{degree_row}

## 7. PPR ± DIS-PHE ablation

Per-unit filtered metrics at pinned restart {cfg.ppr_restart} (sensitivity rows
{cfg.ppr_sensitivity} in `ppr_disphe.parquet`):

| unit | stratum | restart | MRR with | MRR without | H@10 with | H@10 without |
|---|---|---|---|---|---|---|
{ppr_rows}

**Mean Δ across {piv.height} informative units: ΔMRR {d_mrr:+.4f}, ΔHits@10 {d_h10:+.4f}.
Verdict: {verdict}** — phenotype context edges are {verdict_word}
topological recovery signal on this graph.
({n_degenerate} unit-variant rows were degenerate — PPR assigns every gene
score 0, a total tie the metric protocol reads as rank 1; excluded from the
delta. These are the same near-isolated units that fail the ≥10-context-edge
floor in §4.)

## 8. Noise-floor simulation — the evidence floor is derived, not heuristic

Filtered Hits@{HITS_K} of a random ranking of {n_genes:,} genes
({cfg.noise_floor_resamples} seeded resamples per key size, seed {cfg.seed}):

| key size | mean H@10 | SD | CV | P(zero hits) |
|---|---|---|---|---|
{noise_rows}

The coefficient of variation inflates sharply at small key sizes — below the
~50-edge region a handful of answers cannot support a stable Hits@{HITS_K}
estimate, which is the derived basis for the ≥{EVIDENCE_FLOOR_PRIMARY}-edge
evidence floor (20-49 → secondary panel, never pooled).

## 9. EFO/OBA resolution pinning (v5.1 §2.1)

| term | label | mondo_id | resolution method | in Track-D | near-synonym bound |
|---|---|---|---|---|---|
{efo_rows}

## 10. Feeds into the CHANGES/0020 freeze package

The freeze entry needs exactly these numbers from this audit:

1. **Co-holdout roster:** the {co_bc.height:,} overlap pairs ≥ {CO_HOLDOUT_THRESHOLD}
   with both keys below the >5% ceiling band (plus {disc_bc.height:,} disclose
   pairs; full matrix in `overlap_matrix.parquet`, flags in
   `overlap_flags.parquet`) + hand-curated seed clusters — each cluster counts
   as ONE statistical unit.
2. **θ (key-size ceiling):** set from the separation-by-band table (§2) — the
   2-5% working range is justified if separation@3 collapses above it; publish
   per-unit key % + separation in the frozen roster.
3. **Roster attrition:** the per-unit attrition table (`unit_audit.parquet`,
   §5) with published reasons — {n_fail_ctx} context-floor, {n_below_ev}
   evidence-floor, {n_secondary} secondary-panel, {n_above_ceiling}
   above-ceiling among candidate units.
4. **BRG composition:** the DIS-PHE keep/drop verdict (§7: ΔMRR {d_mrr:+.4f},
   ΔH@10 {d_h10:+.4f} across {piv.height} informative units, {n_degenerate}
   degenerate all-tie rows excluded) is recorded in the frozen BRG.
5. **Evidence floor derivation:** the noise-floor table (§8) accompanies the
   frozen roster.
6. **Coherence disclosure:** density percentile {coherence["density_percentile"]:.0f}th,
   degree percentile {coherence["mean_degree_percentile"]:.0f}th (§6).
7. **EFO/OBA pinning:** the resolution table (`efo_pinning.parquet`, §9)
   bounds roster incompleteness.

## Reproduction

`python scripts/m2_lodo_audits.py --config configs/m2_lodo_audits.yaml --resume`
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
    cfg = load_config(args.config, LodoAuditConfig)

    if args.dry_run:
        print(f"would run LODO audits → {cfg.output_dir} (report: {cfg.report_path})")
        return 0

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)

    def skip(name: str) -> bool:
        return args.resume and name in run.done_stages

    with Run(
        script="m2_lodo_audits",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        # load is cheap and every later stage needs the state — always re-run.
        with run.stage("load"):
            state = load_inputs(cfg)

        if not skip("overlap"):
            with run.stage("overlap", outputs=[cfg.output_dir / "overlap_matrix.parquet"]):
                run_overlap(cfg, state)
        if not skip("bfs_separation"):
            with run.stage("bfs_separation", outputs=[cfg.output_dir / "bfs_separation.parquet"]):
                run_bfs(cfg, state)
        if not skip("connectivity"):
            with run.stage(
                "connectivity", outputs=[cfg.output_dir / "connectivity_attrition.parquet"]
            ):
                run_connectivity(cfg, state)
        if not skip("ceiling"):
            with run.stage("ceiling", outputs=[cfg.output_dir / "key_size_ceiling.parquet"]):
                run_ceiling(cfg, state)
        if not skip("coherence"):
            with run.stage("coherence", outputs=[cfg.output_dir / "coherence.json"]):
                run_coherence(cfg, state)
        if not skip("ppr_disphe"):
            with run.stage("ppr_disphe", outputs=[cfg.output_dir / "ppr_disphe.parquet"]):
                run_ppr_disphe(cfg, state)
        if not skip("noise_floor"):
            with run.stage("noise_floor", outputs=[cfg.output_dir / "noise_floor.parquet"]):
                run_noise_floor(cfg, state)
        if not skip("efo_pinning"):
            with run.stage("efo_pinning", outputs=[cfg.output_dir / "efo_pinning.parquet"]):
                run_efo_pinning(cfg, state)

        with run.stage("report", outputs=[cfg.report_path]):
            write_report(cfg, state)

    log.info("lodo_audits_done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
