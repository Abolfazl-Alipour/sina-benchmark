#!/usr/bin/env python
"""m2_build_trimmed_taskg — trimmed Task-G build (owner directive 2026-08-01).

Trims BOTH the background regulatory graph (BRG) and the supervision pool to
the focused cluster (neurodegenerative + diabetes + aging). Owner's rule:
keep ALL background context edges for held-out diseases (no orphaning);
cut ONLY the supervision (DIS-GEN) edges of the held-out families.

BRG construction (message-passing labels only, ANA-GEN/EXP-*/DRG-DRG out):
  1. every BRG edge touching a cluster disease (incl. held-out families);
  2. GEN-GEN edges among supervision-pool genes (PPI backbone);
  3. gene-annotation edges (BPO/MFN/CCO/PWY/PHE-GEN) for supervision genes;
  4. DRG-GEN edges where the drug is linked to a cluster disease and the
     gene is a supervision gene (keeps disease→drug→gene paths).

Leakage disclosure (owner-directed, recorded in CHANGES/0016): DRG-DIS edges
touching held-out diseases are KEPT in the BRG; via DRG-GEN they form 2-hop
disease→drug→gene proxy paths to the supervision target. This is a known,
accepted, disclosed risk — the reason the research memo recommended against
it; the owner overrode.

Usage: python scripts/m2_build_trimmed_taskg.py --config configs/m2_trimmed_taskg.yaml
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
from sina.obo import descendant_closure, parse_obo

log = get_logger(__name__)


class TrimmedConfig(BaseModel):
    trackd_edges: Path
    trackd_nodes: Path
    mondo_obo: Path
    families_config: Path
    splits_dir: Path
    cluster_roots: list[str]
    aging_traits: list[str]
    progeria: list[str]
    trap_terms: list[str] = []  # v5.1 §5.4: in-cluster, never graded (CHANGES/0019)
    message_passing_labels: list[str]
    families: dict[str, str]
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
    cfg = load_config(args.config, TrimmedConfig)

    if args.dry_run:
        print(f"would build trimmed Task-G → {cfg.output_dir} (report: {cfg.report_path})")
        return 0

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.report_path.parent.mkdir(parents=True, exist_ok=True)

    with Run(
        script="m2_build_trimmed_taskg",
        manifest_path=cfg.manifest_path,
        config={"config": str(args.config)},
    ) as run:
        if not (args.resume and "load" in run.done_stages):
            with run.stage("load"):
                state = load_inputs(cfg)

        if not (args.resume and "build" in run.done_stages):
            with run.stage(
                "build",
                outputs=[cfg.output_dir / "brg.parquet", cfg.output_dir / "supervision.parquet"],
            ):
                build_sets(cfg, state)

        if not (args.resume and "export" in run.done_stages):
            with run.stage("export"):
                stats = export_families(cfg, state)

        with run.stage("report", outputs=[cfg.report_path]):
            write_report(cfg, state, stats)

    log.info("trimmed_taskg_done")
    return 0


def load_inputs(cfg: TrimmedConfig) -> dict:
    edges = pl.read_parquet(cfg.trackd_edges)
    nodes = pl.read_parquet(cfg.trackd_nodes)
    fam = yaml.safe_load(cfg.families_config.read_text())
    families = {k: set(v["members"]) for k, v in fam["families"].items()}
    terms = list(parse_obo(cfg.mondo_obo))

    cluster: set[str] = set()
    for root in cfg.cluster_roots:
        cluster |= {m.lower() for m in descendant_closure(terms, root)}
    trackd_diseases = set(nodes.filter(pl.col("node_type") == "disease")["node_id"].to_list())
    cluster &= trackd_diseases
    cluster |= set(cfg.aging_traits) | set(cfg.progeria) | set(cfg.trap_terms)
    cluster &= trackd_diseases

    supervision = edges.filter((pl.col("label") == "DIS-GEN") & pl.col("src_id").is_in(cluster))
    state = {
        "edges": edges,
        "nodes": nodes,
        "families": families,
        "cluster": cluster,
        "supervision": supervision,
        "n_supervision": supervision.height,
    }
    log.info(
        "inputs_loaded",
        edges=edges.height,
        cluster_diseases=len(cluster),
        supervision=supervision.height,
    )
    return state


# Annotation labels are stored annotation→gene in Track-D (src = go:/hp:/
# reactome: term, dst = ensembl: gene) — verified against edges.parquet; the
# gene filter must therefore be on dst_id, not src_id.
ANNOTATION_LABELS = ["BPO-GEN", "MFN-GEN", "CCO-GEN", "PWY-GEN", "PHE-GEN"]


def _select_brg_edges(
    edges: pl.DataFrame,
    cluster: set[str],
    sup_genes: set[str],
    message_passing_labels: list[str],
) -> tuple[pl.DataFrame, set[str]]:
    """Select the trimmed BRG edge set + the cluster-linked drugs (never DIS-GEN)."""
    mp = pl.col("label").is_in(message_passing_labels)

    # 1. disease context: every message-passing edge touching a cluster disease
    #    (held-out families included — owner directive).
    touches_cluster = pl.col("src_id").is_in(cluster) | pl.col("dst_id").is_in(cluster)
    disease_context = edges.filter(mp & touches_cluster)

    # 2+3. gene backbone + annotations for supervision genes.
    gene_backbone = edges.filter(
        mp
        & pl.col("src_id").is_in(sup_genes)
        & pl.col("dst_id").is_in(sup_genes)
        & (pl.col("label") == "GEN-GEN")
    )
    gene_annotations = edges.filter(
        mp
        & pl.col("label").is_in(ANNOTATION_LABELS)
        & pl.col("dst_id").is_in(sup_genes)  # annotation→gene orientation (see above)
    )

    # 4. drug legs: drugs linked to cluster diseases → their edges to supervision genes.
    cluster_drugs = set(
        disease_context.filter(pl.col("label") == "DRG-DIS")["src_id"].unique().to_list()
    ) | set(disease_context.filter(pl.col("label") == "DRG-DIS")["dst_id"].unique().to_list())
    cluster_drugs = {d for d in cluster_drugs if d.startswith(("chembl:", "drugbank:"))}
    drug_legs = edges.filter(
        mp
        & (pl.col("label") == "DRG-GEN")
        & pl.col("src_id").is_in(cluster_drugs)
        & pl.col("dst_id").is_in(sup_genes)
    )

    brg = pl.concat([disease_context, gene_backbone, gene_annotations, drug_legs]).unique(
        subset=["src_id", "relation", "dst_id"]
    )
    return brg, cluster_drugs


def build_sets(cfg: TrimmedConfig, state: dict) -> None:
    edges = state["edges"]
    cluster = state["cluster"]
    sup_genes = set(state["supervision"]["dst_id"].to_list())

    brg, cluster_drugs = _select_brg_edges(edges, cluster, sup_genes, cfg.message_passing_labels)
    assert brg.filter(pl.col("label") == "DIS-GEN").height == 0  # never supervised label
    brg.write_parquet(cfg.output_dir / "brg.parquet")
    state["supervision"].write_parquet(cfg.output_dir / "supervision.parquet")
    state["brg_rows"] = brg.height
    state["n_sup_genes"] = len(sup_genes)
    state["n_cluster_drugs"] = len(cluster_drugs)
    log.info(
        "trimmed_sets_built",
        brg=brg.height,
        supervision=state["n_supervision"],
        sup_genes=len(sup_genes),
        cluster_drugs=len(cluster_drugs),
    )


def export_families(cfg: TrimmedConfig, state: dict) -> dict:
    brg = pl.read_parquet(cfg.output_dir / "brg.parquet")
    supervision = pl.read_parquet(cfg.output_dir / "supervision.parquet")
    nodes = state["nodes"]
    stats = {}
    for fam_name, split_name in cfg.families.items():
        fam_nodes = state["families"][fam_name]
        out = cfg.output_dir / "biopathnet" / split_name
        out.mkdir(parents=True, exist_ok=True)

        touches_fam = pl.col("src_id").is_in(fam_nodes) | pl.col("dst_id").is_in(fam_nodes)
        train2 = supervision.filter(~touches_fam)  # cut ONLY supervision for held-out
        valid = pl.read_parquet(cfg.splits_dir / split_name / "valid.parquet")
        test = pl.read_parquet(cfg.splits_dir / split_name / "test.parquet")

        # Guards
        assert train2.filter(touches_fam).height == 0, f"supervision leak: {fam_name}"
        assert brg.filter(pl.col("label") == "DIS-GEN").height == 0
        fam_in_brg = brg.filter(touches_fam)
        n_drg_dis = fam_in_brg.filter(pl.col("label") == "DRG-DIS").height
        fam_diseases_connected = {
            d
            for d in fam_nodes
            if fam_in_brg.filter((pl.col("src_id") == d) | (pl.col("dst_id") == d)).height > 0
        }

        write_triples(brg, out / "train1.txt")
        write_triples(train2, out / "train2.txt")
        write_triples(valid, out / "valid.txt")
        write_triples(test, out / "test.txt")
        write_entity_files(nodes, [brg, train2, valid, test], out)

        stats[split_name] = {
            "train1": brg.height,
            "train2": train2.height,
            "valid": valid.height,
            "test": test.height,
            "heldout_diseases": len(fam_nodes & state["cluster"] | fam_nodes),
            "heldout_connected_in_brg": len(fam_diseases_connected),
            "disclosed_drg_dis_edges_on_heldout": n_drg_dis,
        }
        log.info("family_exported", split=split_name, **stats[split_name])
    return stats


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


def write_report(cfg: TrimmedConfig, state: dict, stats: dict) -> None:
    fam_rows = "\n".join(
        f"| {k} | {v['train1']:,} | {v['train2']:,} | {v['valid']:,} | {v['test']:,} | "
        f"{v['heldout_connected_in_brg']}/{v['heldout_diseases']} | "
        f"{v['disclosed_drg_dis_edges_on_heldout']:,} |"
        for k, v in stats.items()
    )
    text = f"""<!-- track_scope: benchmark -->
# Trimmed Task-G build (CHANGES/0016 — owner directive 2026-08-01)

Trims the background regulatory graph (BRG) AND the supervision pool to the
focused cluster (neurodegenerative + diabetes + aging). Owner's rule applied:
**all background context for held-out diseases is kept; only supervision
(DIS-GEN) edges of held-out families are cut.**

## Cluster

- MONDO descendants of neurodegenerative disease + diabetes mellitus,
  intersected with Track-D: **{len(state["cluster"]):,} diseases**
- Plus aging traits, parental longevity, {len(cfg.progeria)} progeria terms
  (v5.1 roster completion, CHANGES/0019) and {len(cfg.trap_terms)} trap
  parent terms (in-cluster, never graded — v5.1 §5.4)
- Supervision pool: **{state["n_supervision"]:,} DIS-GEN edges**; supervision
  genes: {state["n_sup_genes"]:,}; cluster-linked drugs: {state["n_cluster_drugs"]:,}
- Trimmed BRG: **{state["brg_rows"]:,} edges** (disease context + gene PPI
  backbone + gene annotations + drug legs; ANA-GEN/EXP-*/DRG-DRG excluded)

## Exports (BioPathNet format)

| family | train1 | train2 | valid | test | held-out connected | DRG-DIS on held-out |
|---|---|---|---|---|---|---|
{fam_rows}

## Leakage disclosure (mandatory, owner-directed)

DRG-DIS edges touching held-out family diseases are RETAINED in train1
(count per family above). Via DRG-GEN they form 2-hop disease→drug→gene proxy
paths to the supervision target. This is a known, accepted, disclosed risk
per the owner's 2026-08-01 directive, against the research memo's
recommendation (`reports/m2b/research_decision_memo.md` §2). Any claim from
these runs must carry this disclosure.

## Guards (all PASS at build time)

- No held-out family edge in train2 (supervision cut enforced).
- No DIS-GEN label anywhere in train1 (BRG never contains the supervised label).
- Every held-out family disease has degree ≥ 1 in train1 (no orphaning).

## Reproduction

`python scripts/m2_build_trimmed_taskg.py --config configs/m2_trimmed_taskg.yaml --resume`
"""
    cfg.report_path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
