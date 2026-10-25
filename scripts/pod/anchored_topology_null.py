#!/usr/bin/env python
"""anchored_topology_null — degree/topology-matched null for the anchored table.

Control analysis (CHANGES/0025, addendum 2026-08-25 follow-up): the anchored
EXPLORATORY pathway table exactly tied the degree-null (HGPS) and ORA (Werner)
baselines and survived the within-disease shuffle (Recall@10 0.0069 -> 0.0062),
suggesting the anchored ranking may be pure graph topology rather than learned
signal. This script rebuilds each disease's anchored candidate set CPU-side
(exactly as extract_paths.discover_anchor_candidates + the support>=3
prefilter), ranks candidates by the topology null

    score_null(P) = support(P) / ln(degree(P))

— the anchored aggregation formula with every path weight replaced by 1 — and
compares its Top-10 against the anchored Top-10 (parsed from the funnel
report's EXPLORATORY section) on the same test-genes-only GO BP/Reactome
ground truth as the frozen evaluation (Recall@10/@5 + Top-10 Jaccard).

CPU-only: pandas/polars for parquet, no torch, no GPU. Reuses the pure
functions from scripts/pod/extract_paths.py verbatim (imported module — its
torch-dependent paths are lazy and never touched here).

Usage:
    python scripts/pod/anchored_topology_null.py
"""

import importlib.util
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

EXPORTS_DIR = Path("/opt/sina/data/processed/m2/lodo_exports")
PATHWAYS_DIR = REPO_ROOT / "data" / "processed" / "m3" / "pathways"
FUNNEL_REPORTS = REPO_ROOT / "reports" / "m2b"
REPORT_PATH = FUNNEL_REPORTS / "anchored_topology_null.md"

# Diseases with held-out truth (AD mondo:0004975 is open-world — skipped).
DISEASES = {
    "mondo__0010196": "mondo:0010196",  # Werner syndrome
    "mondo__0008310": "mondo:0008310",  # Hutchinson-Gilford progeria
}

# Frozen pins from configs/m3_pathway_funnel.yaml (read-only reference; do not edit).
TOP_K_GENES = 100  # extraction.top_k_genes
MIN_CANDIDATE_SUPPORT = 3  # anchored.min_candidate_support (addendum 2026-08-25)
GROUND_TRUTH_MIN_GENES = 3  # evaluation.ground_truth_min_genes
SELECT = 10  # aggregation.select

# Anchored Recall@10/@5 as stated in the funnel reports (for cross-checking the
# recomputation from the parsed Top-10).
STATED_ANCHORED_RECALL = {
    "mondo__0010196": (0.0043, 0.0024),
    "mondo__0008310": (0.0069, 0.0034),
}


def _load_extract_paths():
    """Import scripts/pod/extract_paths.py as a module (light module-level code;
    torch is lazily imported inside functions and never loaded here)."""
    spec = importlib.util.spec_from_file_location(
        "extract_paths", REPO_ROOT / "scripts" / "pod" / "extract_paths.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# Pure computations (unit-tested)
# ---------------------------------------------------------------------------


def topology_null_rank(
    candidates: dict[int, set[int]], degrees: dict[int, int], select: int, hub_adjust
) -> list[dict]:
    """Topology null ranking: score_null(P) = support(P) / ln(degree(P)) — the
    anchored aggregation formula with all path weights replaced by 1 (the
    per-seed constant factor cancels in the ranking). Same tie-break as
    aggregate_anchored: (-score, degree, node)."""
    rows = [
        {
            "node": p,
            "support": len(genes),
            "degree": degrees.get(p, 0),
            "score": hub_adjust(float(len(genes)), degrees.get(p, 0)),
        }
        for p, genes in candidates.items()
    ]
    rows.sort(key=lambda r: (-r["score"], r["degree"], r["node"]))
    return rows[:select]


def parse_anchored_top10(report_text: str) -> list[str]:
    """Parse the anchored Top-10 pathway tokens from a funnel report's
    EXPLORATORY section table. Cell format: 'go:0010628 (name)' -> token is the
    first whitespace-separated field."""
    lines = report_text.splitlines()
    start = next(
        (i for i, line in enumerate(lines) if line.startswith("## EXPLORATORY")), None
    )
    if start is None:
        raise ValueError("no EXPLORATORY section in report")
    tokens: list[str] = []
    in_table = False
    for line in lines[start + 1 :]:
        if line.startswith("|"):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if cells[0] in ("rank", "---") or cells[0].startswith("---"):
                in_table = True
                continue  # header / separator
            in_table = True
            if cells[0] == "—":
                continue  # 'no anchored pathway passed min_support' placeholder
            tokens.append(cells[1].split()[0])
        elif in_table:
            break  # first non-table line ends the table
    return tokens


def jaccard(a: list, b: list) -> float:
    sa, sb = set(a), set(b)
    union = sa | sb
    return len(sa & sb) / len(union) if union else 0.0


# ---------------------------------------------------------------------------
# Data loading (CPU-only)
# ---------------------------------------------------------------------------


def load_export(ep, export_dir: Path):
    """Load a LODO keep export as integer-indexed structures matching the
    evaluate stage: entity_vocab (entities.txt order), node_type per index
    (entity_types.txt), train1 edges as (h, t, r) index triples."""
    entity_vocab = (export_dir / "entities.txt").read_text(encoding="utf-8").splitlines()
    token_type = {}
    for line in (export_dir / "entity_types.txt").read_text(encoding="utf-8").splitlines():
        token, _, type_id = line.partition("\t")
        token_type[token] = int(type_id)
    node_type = [token_type[token] for token in entity_vocab]
    inv = {token: i for i, token in enumerate(entity_vocab)}
    rel_inv: dict[str, int] = {}
    train1_edges = []
    for h, r, t in ep.read_tsv_triples(export_dir / "train1.txt"):
        # relations have their own id space (unused by the candidate/degree
        # computations, but kept so edge triples stay well-formed)
        rel_inv.setdefault(r, len(rel_inv))
        train1_edges.append((inv[h], inv[t], rel_inv[r]))
    return entity_vocab, node_type, inv, train1_edges


def analyze_disease(ep, cs: str, curie: str) -> dict:
    """Rebuild candidates + truth CPU-side, rank the topology null, parse the
    anchored Top-10 from the funnel report, and compare on the frozen truth."""
    export_dir = EXPORTS_DIR / f"{cs}_keep_disphe"
    entity_vocab, node_type, inv, train1_edges = load_export(ep, export_dir)
    names = ep.load_entity_names(export_dir)

    # Node classes exactly as the funnel (GO BP split via the hallmark crosswalk).
    _, _, crosswalk_go = ep.load_hallmark_maps()
    classes = ep.classify_nodes(entity_vocab, node_type, crosswalk_go)
    pathway_ids = classes["go_bp"] | classes["reactome"]

    # Consensus Top-100 genes from the pooled paths parquet (gene_rank <= 100).
    parquet = PATHWAYS_DIR / f"{cs}__paths_top100.parquet"
    top_gene_tokens = {
        r["gene_id"]
        for r in ep.read_parquet_dicts(parquet)
        if r["gene_rank"] <= TOP_K_GENES
    }
    top_gene_ids = {inv[t] for t in top_gene_tokens if t in inv}

    # Candidates: direct train1 annotation edge to a consensus Top-K gene,
    # support>=3 prefilter (identical to stage_anchor).
    candidates_all = ep.discover_anchor_candidates(train1_edges, pathway_ids, top_gene_ids)
    candidates = ep.prefilter_candidates(candidates_all, MIN_CANDIDATE_SUPPORT)

    # Undirected train1 degree, same pin as aggregation.normalize.
    degrees: dict[int, int] = defaultdict(int)
    for h, t, _ in train1_edges:
        degrees[h] += 1
        degrees[t] += 1

    null_rows = topology_null_rank(candidates, dict(degrees), SELECT, ep.hub_adjust)
    null_top10 = [entity_vocab[r["node"]] for r in null_rows]

    # Ground truth identical to the evaluate stage: test-genes-only GO/Reactome.
    test_gene_ids = {
        inv[t]
        for s, _, t in ep.read_tsv_triples(export_dir / "test.txt")
        if s == curie and t in inv
    } & classes["gene"]
    neighbors = ep.gene_pathway_neighbors(train1_edges, classes["gene"], pathway_ids)
    truth = ep.ground_truth_pathways(test_gene_ids, neighbors, GROUND_TRUTH_MIN_GENES)

    # Anchored Top-10 from the funnel report's EXPLORATORY section.
    report = (FUNNEL_REPORTS / f"pathway_funnel__{cs}.md").read_text(encoding="utf-8")
    anchored_top10 = parse_anchored_top10(report)

    return {
        "cs": cs,
        "curie": curie,
        "n_candidates_all": len(candidates_all),
        "n_candidates": len(candidates),
        "n_top_genes": len(top_gene_ids),
        "n_test_genes": len(test_gene_ids),
        "n_true_pathways": len(truth),
        "null_rows": null_rows,
        "null_top10": null_top10,
        "anchored_top10": anchored_top10,
        "anchored_hits@10": len({inv[t] for t in anchored_top10 if t in inv} & truth),
        "null_hits@10": len({r["node"] for r in null_rows} & truth),
        "anchored_recall@10": ep.recall_at_k(
            [inv[t] for t in anchored_top10 if t in inv], truth, 10
        ),
        "anchored_recall@5": ep.recall_at_k(
            [inv[t] for t in anchored_top10 if t in inv], truth, 5
        ),
        "null_recall@10": ep.recall_at_k(
            [r["node"] for r in null_rows], truth, 10
        ),
        "null_recall@5": ep.recall_at_k([r["node"] for r in null_rows], truth, 5),
        "jaccard_top10": jaccard(anchored_top10, null_top10),
        "stated": STATED_ANCHORED_RECALL.get(cs),
        "names": names,
    }


def verdict(result: dict) -> str:
    """Plain-language verdict: does the anchored table beat the topology null?"""
    if result["anchored_recall@10"] > result["null_recall@10"]:
        return "anchored adds signal beyond topology"
    return "anchored ≈ topology null — consistent with within-shuffle survival"


def render_report(results: list[dict]) -> str:
    sections = []
    for r in results:
        stated = (
            f"stated in funnel report: R@10 {r['stated'][0]}, R@5 {r['stated'][1]}; "
            if r["stated"]
            else ""
        )
        sections.append(
            f"""## {r['curie']} ({r['cs']})

Candidates: {r['n_candidates_all']} with a direct train1 annotation edge to a
consensus Top-{TOP_K_GENES} gene ({r['n_top_genes']} genes), prefiltered to
support >= {MIN_CANDIDATE_SUPPORT} -> {r['n_candidates']}; truth =
{r['n_true_pathways']} pathways from {r['n_test_genes']} held-out test genes
(>= {GROUND_TRUTH_MIN_GENES} test genes each).

| ranking | Recall@10 | Recall@5 |
|---|---|---|
| anchored (EXPLORATORY table) | {r['anchored_recall@10']:.4f} | {r['anchored_recall@5']:.4f} |
| topology null support/ln(degree) | {r['null_recall@10']:.4f} | {r['null_recall@5']:.4f} |

Top-10 Jaccard (anchored vs null): {r['jaccard_top10']:.2f}
({stated}recomputed here from the parsed Top-10 on the identically
reconstructed truth). Absolute Recall@10 hits: anchored {r['anchored_hits@10']}
of {r['n_true_pathways']} true pathways, null {r['null_hits@10']}.

Null Top-10: {', '.join(r['null_top10'])}

**Verdict: {verdict(r)}.**
"""
        )
    return (
        "<!-- track_scope: benchmark -->\n"
        "# Anchored pathway table vs topology null (CHANGES/0025 control analysis)\n\n"
        "Question: does the anchored EXPLORATORY pathway table carry learned signal "
        "beyond graph topology? Motivation: anchored Recall@10 exactly tied the "
        "degree-null (HGPS 0.0069 vs 0.0069) and ORA (Werner 0.0043 vs 0.0043) "
        "baselines, and the within-disease shuffle control left anchored Recall@10 "
        "nearly intact (0.0069 -> 0.0062) while gene-level MRR collapsed 58-89%. "
        "Null construction: the anchored aggregation formula with every path weight "
        "replaced by 1 — score_null(P) = support(P) / ln(degree(P)) — over the "
        "identically rebuilt candidate set (direct train1 annotation edge to a "
        "consensus Top-100 gene, support >= 3 prefilter), graded on the identical "
        "test-genes-only GO BP/Reactome truth. Anchored per-candidate weights are "
        "not persisted anywhere (checked data/processed/m3/ and reports/m2b/ — only "
        "the Top-10 tables survive in the funnel reports), so the comparison is "
        "metric-level plus Top-10 Jaccard, not per-candidate correlation. This is a "
        "control analysis; it supports no recovery claim either way.\n\n"
        + "\n".join(sections)
    )


def main() -> int:
    ep = _load_extract_paths()
    results = [analyze_disease(ep, cs, curie) for cs, curie in DISEASES.items()]
    REPORT_PATH.write_text(render_report(results), encoding="utf-8")
    for r in results:
        print(
            f"{r['cs']}: anchored R@10 {r['anchored_recall@10']:.4f} vs null "
            f"{r['null_recall@10']:.4f} | Jaccard {r['jaccard_top10']:.2f} | "
            f"{verdict(r)}"
        )
    print(f"report: {REPORT_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
