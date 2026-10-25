#!/usr/bin/env python
"""extract_paths — Pathway Discovery Funnel runner (Job C, CHANGES/0025).

Implements the frozen spec of ``configs/m3_pathway_funnel.yaml``: load a saved
Wave-1 BioPathNet/NBFNet checkpoint (LODO ``*_keep_disphe`` export or the full
graph for open-world AD), score all genes for one disease, extract
gradient-weighted reasoning paths for the Top-K genes (beam search over
per-edge grads of ``sigmoid(score)`` w.r.t. path edge weights — the BioKGC /
BioPathNet interpretation machinery, reused not reinvented), aggregate paths
into hub-penalized GO-BP / Reactome / Hallmark pathway scores, and grade
pathway-level recovery against test-genes-only ground truth with degree-null
and ORA baselines plus a mask-top-paths faithfulness control.

Usage:
    python scripts/pod/extract_paths.py \
        --config configs/m3_pathway_funnel.yaml --disease mondo:0010196 [--resume] [--dry-run]

Environment note: the pod's ``biopathnet`` env (torch/torchdrug) lacks
pydantic/structlog/psutil/polars, so heavy and sina-internal imports are lazy
with small local fallbacks; repo conventions (``sina.config.load_config``,
``sina.logging``, ``sina.manifests.Run``) are used whenever importable.
"""

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict, deque
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

BIOPATHNET_DIR = Path("/opt/sina/third_party/biopathnet")
if not BIOPATHNET_DIR.is_dir():
    BIOPATHNET_DIR = REPO_ROOT / "third_party" / "biopathnet"

# Base torchdrug/biopathnet config template — same one run_wave1_benchmark.py uses.
BASE_CFG_TEMPLATE = REPO_ROOT / "configs" / "m2b_trimmed_ad_long.yaml"

# LODO exports (KEEP graph) live on the pod data disk; repo-relative fallback
# mirrors configs/m2_triviality_keep.yaml.
LODO_EXPORTS_DIRS = [
    Path("/opt/sina/data/processed/m2/lodo_exports"),
    REPO_ROOT / "data" / "processed" / "m2" / "lodo_exports_keep",
]

# Open-world mode (mode: open_world_keep, checkpoint: null) runs de novo
# inference on the full 1,431,945-edge KEEP graph. Assumption: the full-graph
# export is the trackd_family_alzheimer biopathnet export (its train1.txt is
# exactly the 1,431,945-edge KEEP graph; verified on disk).
OPEN_WORLD_EXPORTS = {
    "mondo:0004975": [
        Path("/opt/sina/data/processed/m2/trimmed/biopathnet/trackd_family_alzheimer"),
        REPO_ROOT / "data/processed/m2/trimmed/biopathnet/trackd_family_alzheimer",
    ],
}

# >2-hop strata context (per-unit hop histograms; per-gene hops are computed
# by BFS on train1 — see stage "evaluate").
TRIVIALITY_PARQUETS = [
    Path("/opt/sina/data/processed/m2/lodo_audits_keep/triviality.parquet"),
    REPO_ROOT / "data" / "processed" / "m2" / "lodo_audits_keep" / "triviality.parquet",
]

GENE_PROPERTIES = [
    Path("/opt/sina/data/processed/m1/gene_properties.parquet"),
    REPO_ROOT / "data" / "processed" / "m1" / "gene_properties.parquet",
]

# Positive-control gates (guards.positive_control): gene-level, hard-fail.
# Ensembl IDs verified against entity_names.txt of the LODO exports; the gene
# *name* is the fallback match (e.g. "WRN RecQ like helicase", "lamin A/C").
POSITIVE_CONTROLS = {
    "mondo:0010196": ("WRN", "ensembl:ENSG00000165392"),  # Werner syndrome
    "mondo:0008310": ("LMNA", "ensembl:ENSG00000160789"),  # HGPS
}

STAGES = ["solver", "score", "extract", "aggregate", "anchor", "evaluate", "faithfulness", "report"]


def utcnow_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")  # noqa: UP017


# ---------------------------------------------------------------------------
# sina-internal deps with local fallbacks (pod biopathnet env lacks them)
# ---------------------------------------------------------------------------

try:  # structured logging (repo convention)
    from sina.logging import configure_logging, get_logger

    log = get_logger(__name__)
except Exception:  # ImportError on pod envs without structlog
    import logging

    def configure_logging(log_file=None, level="INFO") -> None:
        logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=level)

    class _StdLogger:
        def __init__(self, name: str) -> None:
            self._log = logging.getLogger(name)

        def _emit(self, lvl: str, event: str, kw: dict) -> None:
            extra = " ".join(f"{k}={v}" for k, v in kw.items())
            getattr(self._log, lvl)("%s%s", event, f" {extra}" if extra else "")

        def info(self, event: str, **kw) -> None:
            self._emit("info", event, kw)

        def warning(self, event: str, **kw) -> None:
            self._emit("warning", event, kw)

        def error(self, event: str, **kw) -> None:
            self._emit("error", event, kw)

    log = _StdLogger(__name__)

try:
    from sina.registry import sha256_file
except Exception:

    def sha256_file(path: Path | str) -> str:
        h = hashlib.sha256()
        with Path(path).open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()


try:
    from sina.manifests import Run
except Exception:  # no psutil/structlog (or py<3.12 syntax) on pod env

    class Run:  # minimal drop-in: stage checkpoints + output sha256, no RSS guard
        """Fallback for sina.manifests.Run when sina/psutil are unavailable."""

        def __init__(self, script: str, manifest_path: Path | str, config: dict | None = None):
            self.manifest_path = Path(manifest_path)
            if self.manifest_path.exists():
                self.manifest = json.loads(self.manifest_path.read_text(encoding="utf-8"))
                log.info("run_resumed", script=script, run_id=self.manifest.get("run_id"))
            else:
                self.manifest = {
                    "run_id": utcnow_iso(),
                    "script": script,
                    "config": config or {},
                    "stages": [],
                    "status": "running",
                }

        @property
        def done_stages(self) -> set[str]:
            return {s["name"] for s in self.manifest["stages"] if s["status"] == "done"}

        @contextmanager
        def stage(self, name: str, outputs: list[Path] | None = None):
            record = {"name": name, "started_at": utcnow_iso(), "status": "running", "outputs": {}}
            self.manifest["stages"].append(record)
            self._save()
            try:
                yield record
            except Exception:
                record["status"] = "failed"
                self.manifest["status"] = "failed"
                self._save()
                raise
            record["status"] = "done"
            record["finished_at"] = utcnow_iso()
            for out in outputs or []:
                out = Path(out)
                if out.is_file():
                    record["outputs"][str(out)] = sha256_file(out)
            self._save()
            log.info("stage_done", stage=name)

        def _save(self) -> None:
            self.manifest_path.parent.mkdir(parents=True, exist_ok=True)
            self.manifest_path.write_text(json.dumps(self.manifest, indent=2), encoding="utf-8")

        def __enter__(self) -> "Run":
            self._save()
            return self

        def __exit__(self, exc_type, *exc) -> None:
            self.manifest["finished_at"] = utcnow_iso()
            self.manifest["status"] = "failed" if exc_type else "done"
            self._save()


# ---------------------------------------------------------------------------
# Config model (pydantic when available; recursive-namespace fallback otherwise)
# ---------------------------------------------------------------------------

try:
    from pydantic import BaseModel

    class DiseaseEntry(BaseModel):
        curie: str
        mode: str
        checkpoint: str | None = None
        note: str = ""

    class ExtractionConfig(BaseModel):
        top_k_genes: int
        top_k_sensitivity: list[int]
        path_topk: int
        max_hops: int
        scoring: str
        scoring_note: str = ""
        min_path_score: float
        dedup: str
        seed: int
        faithfulness_test: str

    class MapEntry(BaseModel):
        source: str
        id_prefix: str
        hallmarks_crosswalk: str | None = None
        top_n_pathways: int | None = None
        weight_in_composite: float | None = None

    class AggregationConfig(BaseModel):
        method: str
        score: str
        normalize: str
        normalize_formula: str = ""
        tie_break: str
        select: int
        genes_per_pathway: int
        min_support: int
        figure_select: int

    class EvaluationConfig(BaseModel):
        ground_truth_source: str
        ground_truth_min_genes: int
        metrics: list[str]
        compare_to: list[str]
        baselines: dict = {}
        significance: str = ""
        strata: str = ""

    class OutputsConfig(BaseModel):
        paths_parquet: str
        report: str
        figure: str
        manifest: str
        benchmark_report: str

    class FunnelConfig(BaseModel):
        status: str = ""
        diseases: list[DiseaseEntry]
        extraction: ExtractionConfig
        maps: dict[str, MapEntry]
        aggregation: AggregationConfig
        evaluation: EvaluationConfig
        consensus: dict = {}
        anchored: dict = {}
        annotation: dict = {}
        guards: dict = {}
        outputs: OutputsConfig

    def load_funnel_config(path: Path) -> FunnelConfig:
        from sina.config import load_config

        return load_config(path, FunnelConfig)

except Exception:  # pydantic/sina absent on the pod env

    class _NS:
        def __init__(self, **kw) -> None:
            for k, v in kw.items():
                setattr(self, k, _wrap(v))

        def get(self, key, default=None):
            return getattr(self, key, default)

    def _wrap(x):
        if isinstance(x, dict):
            return _NS(**x)
        if isinstance(x, list):
            return [_wrap(v) for v in x]
        return x

    def load_funnel_config(path: Path):
        data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            raise ValueError(f"config {path} must be a YAML mapping at top level")
        return _wrap(data)


def config_sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---------------------------------------------------------------------------
# Pure computations (unit-tested on synthetic fixtures)
# ---------------------------------------------------------------------------


def curie_safe(curie: str) -> str:
    """Filesystem-safe unit name ('mondo:0010196' -> 'mondo__0010196')."""
    return curie.replace(":", "__")


def resolve_checkpoint(configured: str | None, curie: str) -> tuple[Path | None, str]:
    """Resolve the config checkpoint path, with the on-disk naming fallback.

    The frozen config names checkpoints ``<curie>_wave1_seed1337.pth`` but the
    actual Wave-1 files are ``<curie>_keep_pass1.pth`` (legacy seed-1337 naming
    in run_wave1_benchmark.py). Config is frozen — fall back, never edit it.
    """
    if not configured:
        return None, "no checkpoint configured (open-world de novo)"
    cfg_path = REPO_ROOT / configured
    if cfg_path.exists():
        return cfg_path, "config path"
    fallback = cfg_path.parent / f"{curie_safe(curie)}_keep_pass1.pth"
    if fallback.exists():
        return fallback, f"fallback {fallback.name} (config name not on disk)"
    raise FileNotFoundError(f"checkpoint missing: neither {cfg_path} nor {fallback} exists")


def resolve_dataset_dir(curie: str, mode: str) -> Path:
    """LODO keep export for lodo_keep; full KEEP graph export for open_world_keep."""
    if mode == "open_world_keep":
        for cand in OPEN_WORLD_EXPORTS.get(curie, []):
            if (cand / "train1.txt").exists():
                return cand
        raise FileNotFoundError(
            f"no open-world export for {curie} in {OPEN_WORLD_EXPORTS.get(curie)}"
        )
    for base in LODO_EXPORTS_DIRS:
        cand = base / f"{curie_safe(curie)}_keep_disphe"
        if (cand / "train1.txt").exists():
            return cand
    raise FileNotFoundError(f"no LODO keep_disphe export for {curie} in {LODO_EXPORTS_DIRS}")


def read_tsv_triples(path: Path) -> list[tuple[str, str, str]]:
    with Path(path).open("r", encoding="utf-8") as fh:
        return [tuple(row) for row in csv.reader(fh, delimiter="\t") if len(row) == 3]


def canonical_path_key(edges: list[tuple[int, int, int]]) -> tuple:
    """Order-insensitive canonical key for a path (extraction.dedup =
    canonicalize_undirected): each edge normalized so (h,t,r) == (t,h,r), the
    path equal to its reverse."""
    norm = tuple(sorted((min(h, t), max(h, t), r) for h, t, r in edges))
    rev = tuple(reversed(norm))
    return min(norm, rev)


def hub_adjust(score: float, degree: int) -> float:
    """aggregation.normalize = hub_penalty_log_degree: score / ln(degree).

    Assumption (spec is silent): degree < 2 makes ln(degree) <= 0, so the
    penalty is undefined — such nodes score 0 (a degree-1 pathway node carries
    no aggregation evidence anyway).
    """
    if degree < 2:
        return 0.0
    return score / math.log(degree)


def seed_checkpoint_name(cs: str, seed: int, primary_seed: int) -> str:
    """Wave-1 checkpoint filename for a consensus seed. The primary seed keeps
    the legacy run_wave1_benchmark name (``{cs}_keep_pass1.pth``); other seeds
    are ``{cs}_keep_seed{N}_pass1.pth``."""
    if seed == primary_seed:
        return f"{cs}_keep_pass1.pth"
    return f"{cs}_keep_seed{seed}_pass1.pth"


def discover_seed_checkpoints(
    ckpt_dir: Path, cs: str, seeds: list[int], primary_seed: int
) -> tuple[list[tuple[int, Path]], list[int]]:
    """Amendment 1 (owner 2026-08-24): auto-discover per-seed checkpoints.

    Returns ([(seed, path)] in config-seed order, [missing seeds]). Missing
    seed checkpoints are skipped with a warning by the caller, so a run works
    with whatever subset is on disk (e.g. 2 seeds today)."""
    found, missing = [], []
    for s in seeds:
        path = ckpt_dir / seed_checkpoint_name(cs, s, primary_seed)
        if path.exists():
            found.append((s, path))
        else:
            missing.append(s)
    return found, missing


def plan_seed_checkpoints(cfg, disease, ckpt, explicit_ckpt: bool = False):
    """Amendment 1 seed plan: [(seed, checkpoint path)]. Open-world
    (checkpoint null) returns [] — the caller trains/saves first; explicit
    --checkpoint or consensus disabled keeps single-checkpoint behavior."""
    ext = cfg.extraction
    seed = ext.seed
    cons = getattr(cfg, "consensus", None) or {}
    if ckpt is None:
        log.info("consensus_not_applicable", reason="open-world de novo single run")
        return []
    if explicit_ckpt or not cons.get("enabled", False):
        return [(seed, ckpt)]
    seeds = [int(s) for s in cons.get("seeds", [seed])]
    ckpt_dir = REPO_ROOT / "reports" / "m2b" / "checkpoints_wave1"
    found, missing = discover_seed_checkpoints(ckpt_dir, curie_safe(disease.curie), seeds, seed)
    if missing:
        log.warning("consensus_seed_ckpts_missing", missing=missing)
    if found:
        log.info("consensus_plan", seeds=[s for s, _ in found])
        return found
    log.warning("consensus_fallback_single", note="no seed checkpoints on disk")
    return [(seed, ckpt)]


def consensus_mean_rank(
    per_seed_ranked: list[list[tuple[int, float]]],
) -> list[tuple[int, float]]:
    """Amendment 1: consensus ranking = mean sigmoid score across seeds,
    tie-break by gene idx. Each seed's list is score_all_genes output over the
    full gene universe; a gene absent from a seed averages over the seeds that
    scored it."""
    sums: dict[int, float] = defaultdict(float)
    counts: dict[int, int] = defaultdict(int)
    for ranked in per_seed_ranked:
        for g, s in ranked:
            sums[g] += s
            counts[g] += 1
    return sorted(
        ((g, sums[g] / counts[g]) for g in sums), key=lambda kv: (-kv[1], kv[0])
    )


def rank_of(ranked: list[tuple[int, float]], gene_idx: int) -> int | None:
    """1-based rank of gene_idx in a ranked [(gene, score)] list (None if absent)."""
    for i, (g, _) in enumerate(ranked, 1):
        if g == gene_idx:
            return i
    return None


def aggregate_pathways(
    gene_paths: list[dict],
    pathway_ids: set[int],
    degrees: dict[int, int],
    min_support: int,
    select: int,
    genes_per_pathway: int,
) -> list[dict]:
    """Pool grad-weighted paths over intermediate pathway nodes.

    gene_paths: [{"gene": idx, "paths": [{"score": float, "nodes": [idx, ...]}]}].
    pathway_score(P) = sum of path scores of paths touching P (each path counts
    once per P), hub-penalized by ln(degree(P)); a pathway needs >= min_support
    distinct genes. Tie-break "specificity": on equal score prefer the lower
    degree (deeper/more specific term — assumption: GO depth is not in the
    export, so degree proxies specificity; root GO:0008150 is the max-degree
    node and loses ties).
    """
    score_sum: dict[int, float] = defaultdict(float)
    support: dict[int, set[int]] = defaultdict(set)
    gene_best: dict[int, dict[int, float]] = defaultdict(dict)
    for gp in gene_paths:
        gene = gp["gene"]
        for path in gp["paths"]:
            for node in set(path["nodes"]) & pathway_ids:
                score_sum[node] += path["score"]
                support[node].add(gene)
                if path["score"] > gene_best[node].get(gene, float("-inf")):
                    gene_best[node][gene] = path["score"]
    rows = []
    for node, total in score_sum.items():
        if len(support[node]) < min_support:
            continue
        rows.append(
            {
                "node": node,
                "raw_score": total,
                "score": hub_adjust(total, degrees.get(node, 0)),
                "degree": degrees.get(node, 0),
                "support": len(support[node]),
                "top_genes": [
                    g
                    for g, _ in sorted(
                        gene_best[node].items(), key=lambda kv: (-kv[1], kv[0])
                    )[:genes_per_pathway]
                ],
            }
        )
    rows.sort(key=lambda r: (-r["score"], r["degree"], r["node"]))
    return rows[:select]


def discover_anchor_candidates(
    train1_edges: list[tuple[int, int, int]],
    pathway_ids: set[int],
    top_genes: set[int],
) -> dict[int, set[int]]:
    """Amendment 2 (owner 2026-08-24): candidate pathways P = GO BP / Reactome
    nodes with a DIRECT annotation edge P->g or g->P in train1 (inverse
    relations folded; any relation token) to at least one consensus Top-K
    gene g. These are the only nodes that can mediate with a direct annotation
    suffix. Returns {P: {g, ...}}."""
    candidates: dict[int, set[int]] = defaultdict(set)
    for h, t, _ in train1_edges:
        if h in pathway_ids and t in top_genes:
            candidates[h].add(t)
        elif t in pathway_ids and h in top_genes:
            candidates[t].add(h)
    return dict(candidates)


def prefilter_candidates(
    candidates: dict[int, set[int]], min_support: int
) -> dict[int, set[int]]:
    """Addendum 2026-08-25: drop candidates with support < min_support BEFORE
    any GPU work. Candidates below the strictest reported min_support
    sensitivity (3) can never appear in any reported table, so the reported
    tables are identical to an exhaustive run (Werner: 3,239 -> 520)."""
    return {p: genes for p, genes in candidates.items() if len(genes) >= min_support}


def aggregate_anchored(
    prefix_weights: dict[int, dict[int, float]],
    candidates: dict[int, set[int]],
    degrees: dict[int, int],
    min_support: int,
    select: int,
    genes_per_pathway: int,
) -> list[dict]:
    """Amendment 2 scoring: anchored_score(P, g, seed) = weight of the best
    grad-weighted D->P prefix path (prefix_weights[seed][P]; the same prefix
    backs every g directly annotated to P, and a (P, g) pair is valid only if
    the direct P—g edge exists in train1 — i.e. g in candidates[P]).
    pathway_score(P) = (sum over top genes g and seeds of anchored_score) /
    ln(degree(P)), degree(P) = train1 edges touching P (same pin as the frozen
    spec). support(P) = distinct consensus Top-K genes with a valid (P, g);
    a pathway needs >= min_support."""
    rows = []
    for node, genes in candidates.items():
        if len(genes) < min_support:
            continue
        per_seed = sum(pw.get(node, 0.0) for pw in prefix_weights.values())
        total = per_seed * len(genes)
        rows.append(
            {
                "node": node,
                "raw_score": total,
                "score": hub_adjust(total, degrees.get(node, 0)),
                "degree": degrees.get(node, 0),
                "support": len(genes),
                "top_genes": sorted(genes)[:genes_per_pathway],
            }
        )
    rows.sort(key=lambda r: (-r["score"], r["degree"], r["node"]))
    return rows[:select]


def cluster_hallmarks(
    pathway_rows: list[dict],
    node_hallmarks: dict[int, set[str]],
    hallmark_names: dict[str, str],
) -> list[dict]:
    """Cluster hub-penalized GO/Reactome pathway rows into Hallmarks of Aging.

    Assumption (spec is silent on hallmark pooling): hallmarks are not graph
    nodes, so a hallmark's score is the sum of its member nodes' *adjusted*
    scores (hub penalty applied once at the member level) and its support the
    union of member genes. Reported separately, never collapsed into the raw
    GO/Reactome ranking (maps.hallmarks_of_aging.weight_in_composite <= 0.15).
    """
    agg: dict[str, dict] = {}
    for row in pathway_rows:
        for hid in node_hallmarks.get(row["node"], ()):
            entry = agg.setdefault(
                hid,
                {"id": hid, "name": hallmark_names.get(hid, hid), "score": 0.0,
                 "members": [], "genes": set()},
            )
            entry["score"] += row["score"]
            entry["members"].append(row["node"])
            entry["genes"] |= set(row["top_genes"])
    rows = [
        {**v, "genes": sorted(v["genes"]), "n_members": len(v["members"])} for v in agg.values()
    ]
    rows.sort(key=lambda r: (-r["score"], r["id"]))
    return rows


def gene_pathway_neighbors(
    train1_edges: list[tuple[int, int, int]],
    gene_ids: set[int],
    pathway_ids: set[int],
) -> dict[int, set[int]]:
    """gene -> adjacent pathway nodes in train1 (any relation; undirected).

    Matches maps.* semantics: GO BP / Reactome membership comes from the
    graph's own edges (the BPO-GEN / PWY-GEN layers of the KEEP graph).
    """
    neighbors: dict[int, set[int]] = defaultdict(set)
    for h, t, _ in train1_edges:
        if h in pathway_ids and t in gene_ids:
            neighbors[t].add(h)
        elif t in pathway_ids and h in gene_ids:
            neighbors[h].add(t)
    return dict(neighbors)


def ground_truth_pathways(
    test_genes: set[int],
    neighbors: dict[int, set[int]],
    min_genes: int,
) -> set[int]:
    """evaluation.ground_truth_source = test_genes_only_go_reactome: a pathway
    is true iff >= min_genes held-out test genes are adjacent to it."""
    counts: dict[int, int] = defaultdict(int)
    for g in test_genes:
        for p in neighbors.get(g, ()):
            counts[p] += 1
    return {p for p, c in counts.items() if c >= min_genes}


def recall_at_k(predicted: list[int], truth: set[int], k: int) -> float:
    if not truth:
        return 0.0
    return len(set(predicted[:k]) & truth) / len(truth)


def precision_at_k(predicted: list[int], truth: set[int], k: int) -> float:
    if k <= 0:
        return 0.0
    return len(set(predicted[:k]) & truth) / k


def mrr(predicted: list[int], truth: set[int]) -> float:
    for rank, node in enumerate(predicted, 1):
        if node in truth:
            return 1.0 / rank
    return 0.0


def ndcg_at_k(predicted: list[int], truth: set[int], k: int) -> float:
    dcg = sum(
        1.0 / math.log2(rank + 1)
        for rank, node in enumerate(predicted[:k], 1)
        if node in truth
    )
    ideal = sum(1.0 / math.log2(rank + 1) for rank in range(1, min(len(truth), k) + 1))
    return dcg / ideal if ideal > 0 else 0.0


def pathway_metrics(predicted: list[int], truth: set[int]) -> dict[str, float]:
    """The evaluation.metrics block: Recall@10/@5, MRR, Precision@3, NDCG@10."""
    return {
        "pathway_recall@10": recall_at_k(predicted, truth, 10),
        "pathway_recall@5": recall_at_k(predicted, truth, 5),
        "pathway_mrr": mrr(predicted, truth),
        "pathway_precision@3": precision_at_k(predicted, truth, 3),
        "pathway_ndcg@10": ndcg_at_k(predicted, truth, 10),
    }


def bh_adjust(pvals: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR adjustment (evaluation.baselines.ora_top100)."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    adjusted = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), 1):
        prev = min(prev, pvals[i] * n / (n - rank + 1))
        adjusted[i] = prev
    return adjusted


def ora_rank(
    top_genes: set[int],
    candidates: set[int],
    neighbors: dict[int, set[int]],
    pathway_genes: dict[int, set[int]],
    n_background: int,
) -> list[dict]:
    """ORA hypergeometric baseline on the flat Top-100 gene list (no paths).

    Ranks candidate pathways by BH-adjusted hypergeometric p-value for the
    overlap between the pathway's train1 gene neighbors and the Top-100 list.
    """
    from scipy.stats import hypergeom

    rows = []
    for p in candidates:
        genes = pathway_genes.get(p, set())
        k = len(genes & top_genes)
        if k == 0:
            continue
        pval = float(hypergeom.sf(k - 1, n_background, len(genes), len(top_genes)))
        rows.append({"node": p, "overlap": k, "set_size": len(genes), "p_value": pval})
    for row, q in zip(rows, bh_adjust([r["p_value"] for r in rows]), strict=False):
        row["q_value"] = q
    rows.sort(key=lambda r: (r["q_value"], r["p_value"], r["node"]))
    return rows


def degree_null_rank(candidates: set[int], degrees: dict[int, int], select: int) -> list[int]:
    """compare_to: degree_null_pathway_recall — rank candidates by hubness only."""
    return sorted(candidates, key=lambda p: (-degrees.get(p, 0), p))[:select]


def bfs_hops(num_nodes: int, edges: list[tuple[int, int, int]], source: int) -> list[int]:
    """Unweighted undirected hop distance from source over train1 (-1 = unreachable)."""
    adj: dict[int, list[int]] = defaultdict(list)
    for h, t, _ in edges:
        adj[h].append(t)
        adj[t].append(h)
    dist = [-1] * num_nodes
    dist[source] = 0
    queue = deque([source])
    while queue:
        u = queue.popleft()
        for v in adj[u]:
            if dist[v] < 0:
                dist[v] = dist[u] + 1
                queue.append(v)
    return dist


def write_parquet(rows: list[dict], path: Path) -> None:
    """Write rows to parquet — polars preferred, pandas/pyarrow fallback."""
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import polars as pl

        pl.DataFrame(rows).write_parquet(path)
        return
    except ImportError:
        pass
    try:
        import pandas as pd

        pd.DataFrame(rows).to_parquet(path, index=False)
        return
    except ImportError as e:
        raise RuntimeError(
            f"cannot write {path}: neither polars nor pandas+pyarrow/fastparquet "
            "is available in this environment"
        ) from e


def read_parquet_dicts(path: Path) -> list[dict]:
    try:
        import polars as pl

        return pl.read_parquet(path).to_dicts()
    except ImportError:
        import pandas as pd

        return pd.read_parquet(path).to_dicts()


# ---------------------------------------------------------------------------
# BioPathNet / torchdrug side (lazy torch — never imported on --dry-run)
# ---------------------------------------------------------------------------


def build_solver(dataset_dir: Path, seed: int, gpus: list[int] | None = None):
    """Mirror scripts/pod/run_wave1_benchmark.py solver construction."""
    import random

    import numpy as np
    import torch
    from torchdrug import core
    from torchdrug.utils import comm

    sys.path.append(str(BIOPATHNET_DIR))
    from biopathnet import dataset, layer, model, task, util  # noqa: F401

    gpus = gpus if gpus is not None else ([0] if torch.cuda.is_available() else [])
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed + int(comm.get_rank()))

    cfg = util.load_config(str(BASE_CFG_TEMPLATE), context={"gpus": gpus})
    cfg.dataset.path = str(dataset_dir)
    cfg.engine.gpus = gpus
    cfg.engine.batch_size = 32  # pod-pinned (run_wave1_benchmark)
    cfg.task.sample_weight = False
    cfg.train.num_epoch = 1

    dataset_obj = core.Configurable.load_config_dict(cfg.dataset)
    solver = util.build_solver(cfg, dataset_obj)
    return solver, dataset_obj


def supervision_relation(dataset_dir: Path, curie: str) -> str:
    """Relation token of the disease->gene supervision edges (read from the
    first matching row of test/valid/train2 rather than hardcoded)."""
    for name in ("test.txt", "valid.txt", "train2.txt"):
        path = dataset_dir / name
        if not path.exists():
            continue
        with path.open("r", encoding="utf-8") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) == 3 and row[0] == curie:
                    return row[1]
    raise ValueError(f"no supervision edges for {curie} in {dataset_dir}")


def load_entity_names(dataset_dir: Path) -> dict[str, str]:
    names = {}
    path = dataset_dir / "entity_names.txt"
    if path.exists():
        with path.open("r", encoding="utf-8") as fh:
            for line in fh:
                token, _, name = line.rstrip("\n").partition("\t")
                names[token] = name
    return names


def classify_nodes(
    entity_vocab: list[str], node_type: list[int], crosswalk_go: set[str]
) -> dict[str, set[int]]:
    """Split node indices into gene / go_bp / reactome / disease classes.

    Genes = ``ensembl:`` prefix, Reactome = ``reactome:`` prefix. GO BP vs
    MF/CC share the ``go:`` prefix, so BP is the node_type id that dominates
    among the hallmark-crosswalk GO terms (all of which are BP) — assumption
    noted; the config's maps only ever use GO Biological Process.
    """
    gene, go, reactome = set(), set(), set()
    for i, token in enumerate(entity_vocab):
        if token.startswith("ensembl:"):
            gene.add(i)
        elif token.startswith("reactome:"):
            reactome.add(i)
        elif token.startswith("go:"):
            go.add(i)
    bp_type_votes: dict[int, int] = defaultdict(int)
    for i in go:
        if entity_vocab[i].upper() in crosswalk_go:
            bp_type_votes[node_type[i]] += 1
    if bp_type_votes:
        bp_type = max(bp_type_votes.items(), key=lambda kv: kv[1])[0]
    else:
        # fallback: modal type among go: nodes (GO CC is a small minority)
        type_counts: dict[int, int] = defaultdict(int)
        for i in go:
            type_counts[node_type[i]] += 1
        bp_type = max(type_counts.items(), key=lambda kv: kv[1])[0]
        log.warning("go_bp_type_fallback", note="no crosswalk GO terms in vocab")
    go_bp = {i for i in go if node_type[i] == bp_type}
    return {"gene": gene, "go_bp": go_bp, "reactome": reactome}


def load_hallmark_maps() -> tuple[dict[str, set[str]], dict[str, str], set[str]]:
    """docs/go_to_hallmarks.tsv + configs/hallmarks_of_aging.yaml.

    Returns (uppercased GO/Reactome id -> hallmark ids, hallmark id -> name,
    set of crosswalk GO ids). The yaml lists the same member ids; the tsv is
    the authoritative crosswalk per maps.go_biological_process.hallmarks_crosswalk.
    """
    node_hallmarks: dict[str, set[str]] = defaultdict(set)
    tsv = REPO_ROOT / "docs" / "go_to_hallmarks.tsv"
    with tsv.open("r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                node_hallmarks[parts[0].upper()].add(parts[1])
    names = {}
    hallmarks_yaml = REPO_ROOT / "configs" / "hallmarks_of_aging.yaml"
    data = yaml.safe_load(hallmarks_yaml.read_text(encoding="utf-8"))
    for entry in data.get("hallmarks_of_aging", []):
        names[entry["id"]] = entry["name"]
        for gid in entry.get("go", []):
            node_hallmarks[gid.upper()].add(entry["id"])
        for rid in entry.get("reactome", []):
            node_hallmarks[rid.upper()].add(entry["id"])
    crosswalk_go = {k for k in node_hallmarks if k.startswith("GO:")}
    return dict(node_hallmarks), names, crosswalk_go


def hallmark_map_for_graph(
    entity_vocab: list[str], node_hallmarks: dict[str, set[str]]
) -> dict[int, set[str]]:
    """Internal node idx -> hallmark ids, normalizing go: -> GO: and
    reactome:R-HSA-n -> R-HSA-n (crosswalk ids are uppercase)."""
    mapping = {}
    for i, token in enumerate(entity_vocab):
        if token.startswith("go:"):
            key = token.upper()
        elif token.startswith("reactome:"):
            key = token.split(":", 1)[1].upper()
        else:
            continue
        if key in node_hallmarks:
            mapping[i] = node_hallmarks[key]
    return mapping


def score_all_genes(solver, dataset_obj, disease_idx: int, rel_idx: int, gene_ids: set[int]):
    """Score every gene as tail of (disease, supervision_relation, ?).

    One NBFNet Bellman-Ford pass over all entities (task.predict eval path),
    then restrict to gene-type nodes. Returns sorted [(gene_idx, sigmoid_score)].
    """
    import torch

    task = solver.model
    task.eval()
    batch = torch.tensor([[disease_idx, 0, rel_idx]], device=solver.device)
    with torch.no_grad():
        pred = task.predict(batch)  # [1, 2, num_entity] on CPU
    scores = pred[0, 0]
    ranked = sorted(
        ((g, float(torch.sigmoid(scores[g]))) for g in gene_ids),
        key=lambda kv: (-kv[1], kv[0]),
    )
    return ranked


def extract_gene_paths(model, fact_graph, h_index, t_index, r_index):
    """Grad-weighted beam-search path extraction (BioKGC/BioPathNet prior art).

    Mirrors NBFNet.visualize (biopathnet/model.py) with one spec-mandated
    change: the score is passed through sigmoid before autograd
    (extraction.scoring_note: "grad = autograd of sigmoid(score) w.r.t. path
    edge weight"). Because NBFNet messages are edge_weight x PNA-aggregated
    (distmult message, pna aggregate), d sigmoid(score)/d edge_weight IS the
    grad x PNA-aggregator weighting the config asks for — NBFNet has no
    transformer attention. Reuses model.bellmanford / beam_search_distance /
    topk_average_length; beam width = model.num_beam, paths = model.path_topk.
    """
    import torch

    assert h_index.numel() == 1 and h_index.ndim == 1
    graph = fact_graph.undirected(add_inverse=True)
    output = model.bellmanford(graph, h_index, r_index, separate_grad=True)
    feature = output["node_feature"]
    step_graphs = output["step_graphs"]

    index = t_index.unsqueeze(0).unsqueeze(-1).expand(-1, -1, feature.shape[-1])
    feature = feature.gather(0, index).squeeze(0)
    score = torch.sigmoid(model.mlp(feature).squeeze(-1))

    edge_weights = [g.edge_weight for g in step_graphs]
    edge_grads = torch.autograd.grad(score, edge_weights)
    for g, edge_grad in zip(step_graphs, edge_grads, strict=False):
        with g.edge():
            g.edge_grad = edge_grad
    distances, back_edges = model.beam_search_distance(
        step_graphs, h_index, t_index, model.num_beam
    )
    return model.topk_average_length(distances, back_edges, t_index, model.path_topk)


def viterbi_forward(edge_lists, edge_grads, num_node: int, source: int):
    """Max-sum DP (Viterbi best path) over the layered step graphs — the
    num_beam=1 core of NBFNet.beam_search_distance, minus the per-target edge
    masking (addendum 2026-08-25; disclosed in the report as an estimator
    difference from the gene-path beam extraction). dist starts at `source`
    (0); per layer, dist[v] = max over edges (u->v) of dist[u] + edge_grad[e].
    Ties break to the lowest edge index (deterministic). Vectorized per layer
    with scatter_reduce (amax values, amin argmax edge index).

    Returns (dists, backptrs): dists[layer][v] = best grad-sum reaching v in
    exactly layer+1 hops (-inf = unreachable); backptrs[layer][v] = argmax edge
    index into that layer's edge list (-1 = unreachable)."""
    import torch

    dist = torch.full((num_node,), float("-inf"), device=edge_grads[0].device)
    dist[source] = 0.0
    dists, backptrs = [], []
    for edges, grads in zip(edge_lists, edge_grads, strict=True):
        node_in, node_out = edges[:, 0], edges[:, 1]
        val = dist[node_in] + grads  # -inf + finite stays -inf (unreachable)
        best = torch.full((num_node,), float("-inf"), device=val.device)
        best.scatter_reduce_(0, node_out, val, reduce="amax")
        # argmax edge per node: val == best[node_out] is exact (amax returns one
        # of the values); among ties the lowest edge index wins via amin
        n_edge = edges.shape[0]
        hit = torch.where(
            val == best[node_out], torch.arange(n_edge, device=val.device), n_edge
        )
        bp = torch.full((num_node,), n_edge, dtype=torch.long, device=val.device)
        bp.scatter_reduce_(0, node_out, hit, reduce="amin")
        dists.append(best)
        backptrs.append(torch.where(bp < n_edge, bp, -1))
        dist = best
    return dists, backptrs


def viterbi_reconstruct(edge_lists, dists, backptrs, source: int, target: int):
    """Best source->target path (<= len(edge_lists) hops) from a viterbi_forward
    pass: pick the layer depth maximizing dists[layer][target] (ties prefer the
    shallowest depth), then walk the backpointers back to `source`. Returns
    (edges, weight) with edges as [(h, t, r)] int tuples; edges is None and
    weight -inf when the target is unreachable."""
    vals = [float(d[target]) for d in dists]
    best = max(vals)
    if best == float("-inf"):
        return None, best
    layer = vals.index(best)
    path = []
    node = target
    while layer >= 0:
        e = int(backptrs[layer][node])
        if e < 0:
            break  # unreachable predecessor (defensive; cannot happen on a finite best)
        h, t, r = (int(x) for x in edge_lists[layer][e])
        path.append((h, t, r))
        node = h
        layer -= 1
    return path[::-1], best


def extract_anchor_paths(model, fact_graph, h_index, t_indices, r_index):
    """Amendment 2: shared-bellmanford variant of extract_gene_paths for many
    targets. One bellmanford pass conditioned on (disease, relation); per
    candidate pathway node P: sigmoid-mlp score, autograd w.r.t. the step-graph
    edge weights (exactly as the gene-path extraction), then the best D->P
    prefix via a max-sum Viterbi DP (viterbi_forward / viterbi_reconstruct)
    instead of beam_search_distance — equivalent to num_beam=1 minus the target
    edge-masking and top-k length averaging; addendum 2026-08-25, disclosed in
    the report as an estimator difference. Relation conditioning stays the
    disease-gene supervision relation — approximation, disclosed in the report.
    Returns {t_idx: (edges, weight)}; edges is None / weight -inf for
    unreachable targets."""
    import torch

    assert h_index.numel() == 1 and h_index.ndim == 1
    graph = fact_graph.undirected(add_inverse=True)
    output = model.bellmanford(graph, h_index, r_index, separate_grad=True)
    feature = output["node_feature"]
    step_graphs = output["step_graphs"]
    edge_weights = [g.edge_weight for g in step_graphs]
    edge_lists = [g.edge_list for g in step_graphs]
    num_node = int(step_graphs[0].num_node)
    h = int(h_index.item())
    results = {}
    for pos, t in enumerate(t_indices):
        t_index = torch.tensor([t], device=feature.device)
        index = t_index.unsqueeze(0).unsqueeze(-1).expand(-1, -1, feature.shape[-1])
        feat = feature.gather(0, index).squeeze(0)
        score = torch.sigmoid(model.mlp(feat).squeeze(-1))
        # all targets share the one bellmanford graph — retain it until the last
        retain = pos < len(t_indices) - 1
        edge_grads = torch.autograd.grad(score, edge_weights, retain_graph=retain)
        dists, backptrs = viterbi_forward(edge_lists, edge_grads, num_node, h)
        results[t] = viterbi_reconstruct(edge_lists, dists, backptrs, h, t)
    return results


def masked_gene_score(model, fact_graph, mask_triplets: set[tuple[int, int, int]],
                      h: int, t: int, r: int) -> float:
    """sigmoid(score) of one (disease, rel, gene) triple with path edges masked.

    Faithfulness control (extraction.faithfulness_test = mask_top_paths_drop):
    path edges live on the undirected graph; here we mask the corresponding
    directed edges of fact_graph (inverse relations folded back)."""
    import torch

    edge_list = fact_graph.edge_list
    keep = torch.ones(fact_graph.num_edge, dtype=torch.bool, device=edge_list.device)
    for eh, et, er in mask_triplets:
        if er >= fact_graph.num_relation:
            eh, et, er = et, eh, er - int(fact_graph.num_relation)
        hit = (edge_list[:, 0] == eh) & (edge_list[:, 1] == et) & (edge_list[:, 2] == er)
        keep &= ~hit
    masked = fact_graph.edge_mask(keep)
    device = edge_list.device
    h_index = torch.tensor([[h]], device=device)
    t_index = torch.tensor([[t]], device=device)
    r_index = torch.tensor([[r]], device=device)
    with torch.no_grad():
        score = model(masked, h_index, t_index, r_index)
    return float(torch.sigmoid(score).item())


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------


def stage_score_and_extract(cfg, cfg_path: Path, disease, dataset_dir, ckpt, cs: str, run,
                            explicit_ckpt: bool = False,
                            skip_positive_control: bool = False) -> dict:
    """Solver + scoring + path extraction. Returns the full extraction state
    (Top-max(K) genes) and writes the configured Top-100 paths parquet (the
    pooled union across seeds) plus per-seed tagged parquets and a sensitivity
    cache parquet (internal; not in the frozen outputs block).

    Amendment 1 (owner 2026-08-24): unless --checkpoint is passed explicitly
    (shuffle-collapse gate / debug) or the disease is open-world de novo, all
    available per-seed checkpoints (consensus.seeds) are discovered, genes are
    scored under EACH separately, and the consensus ranking (mean sigmoid score
    across seeds, tie-break gene idx) defines the Top-K selection and the
    positive-control gate. Paths are then extracted per checkpoint for the
    consensus Top-K and pooled; every row keeps per-seed provenance."""
    import torch

    ext = cfg.extraction
    seed = ext.seed
    base_cs = curie_safe(disease.curie)
    open_world = ckpt is None
    # Amendment 1 seed plan: [(seed, checkpoint path)]; [] while open-world
    seed_plan: list[tuple[int, Path]] = plan_seed_checkpoints(cfg, disease, ckpt, explicit_ckpt)

    ckpt_sha_by_seed: dict[int, str] = {}
    with run.stage(f"solver:{cs}"):
        solver, dataset_obj = build_solver(dataset_dir, seed)
        assert len(solver.model.model.layers) <= ext.max_hops, (
            f"NBFNet depth {len(solver.model.model.layers)} exceeds max_hops {ext.max_hops}"
        )
        if open_world:
            # open_world_keep: checkpoint null -> train 1 pass on the full
            # graph, then save to the keep_pass1 fallback name (assumption:
            # mirrors run_wave1_benchmark's train-when-no-checkpoint branch).
            solver.model.split = "train"
            solver.train(num_epoch=1)
            ckpt = REPO_ROOT / "reports" / "m2b" / "checkpoints_wave1" / f"{base_cs}_keep_pass1.pth"
            ckpt.parent.mkdir(parents=True, exist_ok=True)
            solver.save(str(ckpt))
            seed_plan = [(seed, ckpt)]
            ckpt_sha_by_seed[seed] = sha256_file(ckpt)
            log.info("open_world_trained", checkpoint=str(ckpt))
        else:
            for s, path in seed_plan:
                ckpt_sha_by_seed[s] = sha256_file(path)
            log.info(
                "checkpoints_resolved",
                seeds=[s for s, _ in seed_plan],
                sha256={s: ckpt_sha_by_seed[s][:12] for s, _ in seed_plan},
            )

    entity_vocab = list(dataset_obj.entity_vocab)
    node_type = dataset_obj.graph.node_type.tolist()
    _, hallmark_names, crosswalk_go = load_hallmark_maps()
    classes = classify_nodes(entity_vocab, node_type, crosswalk_go)
    rel_token = supervision_relation(dataset_dir, disease.curie)
    disease_idx = dataset_obj.inv_entity_vocab[disease.curie]
    rel_idx = dataset_obj.inv_relation_vocab[rel_token]
    n_train1 = dataset_obj.num_samples[0]
    log.info(
        "graph_loaded",
        nodes=len(entity_vocab),
        train1_edges=n_train1,
        genes=len(classes["gene"]),
        go_bp=len(classes["go_bp"]),
        reactome=len(classes["reactome"]),
        relation=rel_token,
    )

    control_ranks = None
    with run.stage(f"score:{cs}"):
        # Amendment 1: score all genes under EACH checkpoint separately, then
        # consensus = mean sigmoid score across seeds (tie-break gene idx).
        per_seed_ranked: dict[int, list] = {}
        for s, path in seed_plan:
            solver.load(str(path))
            per_seed_ranked[s] = score_all_genes(
                solver, dataset_obj, disease_idx, rel_idx, classes["gene"]
            )
            log.info("genes_scored", seed=s, total=len(per_seed_ranked[s]))
        ranked = consensus_mean_rank(list(per_seed_ranked.values()))
        k_max = max(ext.top_k_sensitivity)
        top_genes = ranked[:k_max]
        log.info("consensus_ranked", seeds=sorted(per_seed_ranked), top_k=k_max)

        # Positive-control gate (guards.positive_control) — hard-fail loudly.
        # Amendment 1: gate = control gene in the CONSENSUS Top-100; per-seed
        # and consensus ranks are logged and reported.
        if skip_positive_control and disease.curie in POSITIVE_CONTROLS:
            log.warning(
                "positive_control_bypassed",
                reason="--no-positive-control (shuffle-collapse gate runs: a shuffled "
                "model is EXPECTED to fail the gate; bypass is for controls only)",
            )
        elif disease.curie in POSITIVE_CONTROLS:
            symbol, ensembl = POSITIVE_CONTROLS[disease.curie]
            top_ids = {entity_vocab[g] for g, _ in top_genes[: ext.top_k_genes]}
            names = load_entity_names(dataset_dir)
            hit = ensembl in top_ids or any(
                symbol.lower() in names.get(t, "").lower().split()[:1]
                or names.get(t, "").lower().startswith(symbol.lower())
                for t in top_ids
            )
            control_idx = dataset_obj.inv_entity_vocab.get(ensembl)
            control_ranks = {
                "symbol": symbol,
                "ensembl": ensembl,
                "per_seed": (
                    {s: rank_of(r, control_idx) for s, r in per_seed_ranked.items()}
                    if control_idx is not None
                    else {}
                ),
                "consensus": rank_of(ranked, control_idx) if control_idx is not None else None,
            }
            log.info(
                "positive_control_ranks",
                symbol=symbol,
                per_seed=control_ranks["per_seed"],
                consensus=control_ranks["consensus"],
            )
            if not hit:
                raise RuntimeError(
                    f"POSITIVE CONTROL FAILED: {symbol} ({ensembl}) not in "
                    f"{disease.curie} consensus Top-{ext.top_k_genes} — funnel gate violated"
                )
            log.info("positive_control_pass", symbol=symbol, gene=ensembl)

    paths_parquet = REPO_ROOT / cfg.outputs.paths_parquet.replace("{disease}", cs)
    cache_parquet = paths_parquet.with_name(f"{cs}__paths_sensitivity_cache.parquet")
    seed_parquets = [
        (s, paths_parquet.with_name(f"{cs}__paths_top100__seed{s}.parquet"))
        for s, _ in seed_plan
    ]
    with run.stage(
        f"extract:{cs}", outputs=[paths_parquet, cache_parquet, *[p for _, p in seed_parquets]]
    ):
        model = solver.model.model  # unwrap torchdrug task -> NBFNet
        fact_graph = solver.model.fact_graph
        num_relation = int(fact_graph.num_relation)
        cfg_sha = config_sha256(cfg_path)[:12]
        gene_paths_by_gene: dict[int, dict] = {}
        pooled_rows: list[dict] = []
        for s, path in seed_plan:
            solver.load(str(path))
            prov = {
                "disease": disease.curie,
                "seed": s,
                "config_sha256": cfg_sha,
                "checkpoint_sha256": ckpt_sha_by_seed.get(s, "")[:12],
            }
            seed_rows: list[dict] = []
            for gene_rank, (gene, gene_score) in enumerate(top_genes, 1):
                h_t = torch.tensor([disease_idx], device=solver.device)
                t_t = torch.tensor([gene], device=solver.device)
                r_t = torch.tensor([rel_idx], device=solver.device)
                paths, weights = extract_gene_paths(model, fact_graph, h_t, t_t, r_t)
                dedup: dict[tuple, dict] = {}
                for pth, weight in zip(paths, weights, strict=False):
                    if weight < ext.min_path_score:
                        continue  # drop near-zero long tail (extraction.min_path_score)
                    edges = [tuple(e) for e in pth]
                    if len(edges) > ext.max_hops:
                        continue
                    key = canonical_path_key(edges)
                    if key not in dedup or weight > dedup[key]["weight"]:
                        dedup[key] = {"edges": edges, "weight": float(weight)}
                kept = sorted(dedup.values(), key=lambda d: (-d["weight"], d["edges"]))
                gp = gene_paths_by_gene.setdefault(
                    gene,
                    {"gene": gene, "gene_rank": gene_rank, "gene_score": gene_score, "paths": []},
                )
                gp["paths"].extend(
                    {
                        "score": d["weight"],
                        "nodes": sorted({n for e in d["edges"] for n in e[:2]}),
                        "edges": d["edges"],
                        "seed": s,
                    }
                    for d in kept
                )
                if gene_rank <= ext.top_k_genes:  # frozen filename: Top-100 only
                    for path_rank, d in enumerate(kept, 1):
                        for hop, (eh, et, er) in enumerate(d["edges"], 1):
                            inv = er >= num_relation
                            seed_rows.append(
                                {
                                    **prov,
                                    "gene_id": entity_vocab[gene],
                                    "gene_rank": gene_rank,
                                    "gene_score": gene_score,
                                    "path_rank": path_rank,
                                    "path_score": d["weight"],
                                    "hop": hop,
                                    "src_id": entity_vocab[eh],
                                    "relation": dataset_obj.relation_vocab[
                                        er - num_relation if inv else er
                                    ]
                                    + ("^(-1)" if inv else ""),
                                    "dst_id": entity_vocab[et],
                                }
                            )
                if gene_rank % 25 == 0:
                    log.info("extract_progress", seed=s, genes=gene_rank, of=len(top_genes))
            seed_parquet = dict(seed_parquets)[s]
            write_parquet(seed_rows, seed_parquet)
            pooled_rows.extend(seed_rows)
            log.info("seed_paths_extracted", seed=s, rows=len(seed_rows))
        gene_paths = [
            gene_paths_by_gene[g]
            for g in sorted(gene_paths_by_gene, key=lambda g: gene_paths_by_gene[g]["gene_rank"])
        ]
        # primary output: pooled union across seeds (frozen filename)
        write_parquet(pooled_rows, paths_parquet)
        # sensitivity cache: full Top-max(K) extraction for K=50/200 aggregation
        cache_rows = [
            {
                "gene": gp["gene"],
                "gene_rank": gp["gene_rank"],
                "gene_score": gp["gene_score"],
                "seed": p.get("seed"),
                "path_score": p["score"],
                "nodes": p["nodes"],
                "edges": [[int(a), int(b), int(c)] for a, b, c in p["edges"]],
            }
            for gp in gene_paths
            for p in gp["paths"]
        ]
        write_parquet(cache_rows, cache_parquet)
        log.info(
            "paths_extracted",
            genes=len(gene_paths),
            rows=len(pooled_rows),
            seeds=[s for s, _ in seed_plan],
        )

    # restore the primary seed's checkpoint so the faithfulness control is
    # deterministic (its top paths are filtered to that seed below)
    primary_seed, primary_ckpt = next(
        ((s, p) for s, p in seed_plan if s == seed), seed_plan[0]
    )
    solver.load(str(primary_ckpt))

    return {
        "solver": solver,
        "dataset_obj": dataset_obj,
        "entity_vocab": entity_vocab,
        "node_type": node_type,
        "classes": classes,
        "disease_idx": disease_idx,
        "rel_idx": rel_idx,
        "rel_token": rel_token,
        "n_train1": n_train1,
        "gene_paths": gene_paths,
        "ckpt_sha": ckpt_sha_by_seed.get(primary_seed, ""),
        "seed_ckpts": seed_plan,
        "faithfulness_seed": primary_seed,
        "control_ranks": control_ranks,
        "hallmark_names": hallmark_names,
        "crosswalk_go": crosswalk_go,
    }


def load_extraction_cache(cfg, cs: str):
    cache = REPO_ROOT / cfg.outputs.paths_parquet.replace("{disease}", cs)
    cache = cache.with_name(f"{cs}__paths_sensitivity_cache.parquet")
    rows = read_parquet_dicts(cache)
    gene_paths: dict[int, dict] = {}
    for r in rows:
        gp = gene_paths.setdefault(
            r["gene"], {"gene": r["gene"], "gene_rank": r["gene_rank"],
                        "gene_score": r["gene_score"], "paths": []}
        )
        gp["paths"].append(
            {"score": r["path_score"], "nodes": list(r["nodes"]),
             "edges": [tuple(e) for e in r["edges"]], "seed": r.get("seed")}
        )
    return [gene_paths[g] for g in sorted(gene_paths, key=lambda g: gene_paths[g]["gene_rank"])]


def stage_aggregate_evaluate(cfg, disease, dataset_dir, cs: str, state: dict, run) -> dict:
    """Pathway aggregation + evaluation (pure CPU once paths exist)."""
    ext, agg_cfg, ev = cfg.extraction, cfg.aggregation, cfg.evaluation

    if "dataset_obj" not in state:
        # resume from cache: rebuild the dataset (no solver/checkpoint needed)
        _, dataset_obj = _load_dataset_only(dataset_dir, ext.seed)
        entity_vocab = list(dataset_obj.entity_vocab)
        node_type = dataset_obj.graph.node_type.tolist()
        _, hallmark_names, crosswalk_go = load_hallmark_maps()
        classes = classify_nodes(entity_vocab, node_type, crosswalk_go)
        state.update(
            {
                "dataset_obj": dataset_obj,
                "entity_vocab": entity_vocab,
                "node_type": node_type,
                "classes": classes,
                "n_train1": dataset_obj.num_samples[0],
                "hallmark_names": hallmark_names,
                "crosswalk_go": crosswalk_go,
                "disease_idx": dataset_obj.inv_entity_vocab[disease.curie],
            }
        )

    dataset_obj = state["dataset_obj"]
    entity_vocab = state["entity_vocab"]
    classes = state["classes"]
    n_train1 = state["n_train1"]
    pathway_ids = classes["go_bp"] | classes["reactome"]

    # train1 edges + undirected degrees (aggregation.normalize: degree(P) =
    # total train1 edges touching P, undirected BRG degree)
    edge_list = dataset_obj.graph.edge_list[:n_train1].tolist()
    train1_edges = [(int(h), int(t), int(r)) for h, t, r in edge_list]
    degrees: dict[int, int] = defaultdict(int)
    for h, t, _ in train1_edges:
        degrees[h] += 1
        degrees[t] += 1

    gene_paths = state["gene_paths"]
    names = load_entity_names(dataset_dir)
    node_hallmarks, _, _ = load_hallmark_maps()
    node_hm_graph = hallmark_map_for_graph(entity_vocab, node_hallmarks)

    results: dict = {"by_k": {}, "names": names}
    with run.stage(f"aggregate:{cs}"):
        for k in sorted(set(ext.top_k_sensitivity)):
            sub = [gp for gp in gene_paths if gp["gene_rank"] <= k]
            rows = aggregate_pathways(
                sub, pathway_ids, dict(degrees), agg_cfg.min_support,
                agg_cfg.select, agg_cfg.genes_per_pathway,
            )
            hallmarks = cluster_hallmarks(rows, node_hm_graph, state["hallmark_names"])
            results["by_k"][k] = {"pathways": rows, "hallmarks": hallmarks}
        log.info(
            "aggregated",
            top=[
                (entity_vocab[r["node"]], round(r["score"], 4), r["support"])
                for r in results["by_k"][ext.top_k_genes]["pathways"][: agg_cfg.figure_select]
            ],
        )

    # Amendment 2 (owner 2026-08-24): EXPLORATORY ontology-anchored variant,
    # after aggregate and before evaluate; consumed by evaluate + report.
    anchored = None
    anc = getattr(cfg, "anchored", None) or {}
    if anc.get("enabled", False):
        anchored = stage_anchor(
            cfg, disease, cs, state, train1_edges, pathway_ids, dict(degrees), run
        )
    results["anchored"] = anchored

    with run.stage(f"evaluate:{cs}"):
        neighbors = gene_pathway_neighbors(train1_edges, classes["gene"], pathway_ids)
        if disease.mode == "open_world_keep":
            # config note: "no valid/test split — scored by pathway support
            # only, not MRR" — pathway-recall metrics are skipped (assumption:
            # we do NOT fabricate a held-out truth from training supervision).
            results["metrics"] = None
            log.info("evaluate_skipped", reason="open_world_keep: pathway support only")
        else:
            test_genes = {
                dataset_obj.inv_entity_vocab[t]
                for s, _, t in read_tsv_triples(dataset_dir / "test.txt")
                if s == disease.curie and t in dataset_obj.inv_entity_vocab
            }
            test_genes &= classes["gene"]
            truth = ground_truth_pathways(
                test_genes, neighbors, ev.ground_truth_min_genes
            )
            primary = results["by_k"][ext.top_k_genes]["pathways"]
            predicted = [r["node"] for r in primary]
            metrics = pathway_metrics(predicted, truth)

            # degree null: same candidate universe (pathways adjacent to
            # >= min_support Top-100 genes), ranked by degree only.
            top100 = {gp["gene"] for gp in gene_paths if gp["gene_rank"] <= ext.top_k_genes}
            cand_counts: dict[int, int] = defaultdict(int)
            for g in top100:
                for p in neighbors.get(g, ()):
                    cand_counts[p] += 1
            candidates = {p for p, c in cand_counts.items() if c >= agg_cfg.min_support}
            null_pred = degree_null_rank(candidates, dict(degrees), agg_cfg.select)
            null_recall = recall_at_k(null_pred, truth, agg_cfg.select)

            # ORA hypergeometric on the flat Top-100 list (no paths).
            pathway_genes = {p: {g for g in classes["gene"] if p in neighbors.get(g, ())}
                             for p in candidates}
            ora_rows = ora_rank(
                top100, candidates, neighbors, pathway_genes, len(classes["gene"])
            )
            ora_pred = [r["node"] for r in ora_rows[: agg_cfg.select]]
            ora_recall = recall_at_k(ora_pred, truth, agg_cfg.select)

            # K sensitivity: Recall@10 for K = 50/100/200 (pre-registered).
            sensitivity = {
                k: recall_at_k(
                    [r["node"] for r in results["by_k"][k]["pathways"]], truth, agg_cfg.select
                )
                for k in sorted(results["by_k"])
            }

            # >2-hop strata: recompute truth on test genes whose unweighted
            # train1 hop distance from the disease exceeds 2 (triviality
            # parquet carries per-unit histograms only, so per-gene hops come
            # from BFS here; parquet cross-checked for context when present).
            hops = bfs_hops(len(entity_vocab), train1_edges, state["disease_idx"])
            deep_genes = {g for g in test_genes if hops[g] > 2}
            deep_truth = ground_truth_pathways(deep_genes, neighbors, ev.ground_truth_min_genes)
            deep_recall = recall_at_k(predicted, deep_truth, agg_cfg.select) if deep_truth else None

            results["metrics"] = {
                "n_true_pathways": len(truth),
                "n_test_genes": len(test_genes),
                **{k2: round(v, 4) for k2, v in metrics.items()},
                "degree_null_recall@10": round(null_recall, 4),
                "ora_recall@10": round(ora_recall, 4),
                "ora_top10": [
                    {"node": entity_vocab[r["node"]], "q_value": round(r["q_value"], 6),
                     "overlap": r["overlap"]}
                    for r in ora_rows[: agg_cfg.select]
                ],
                "sensitivity_recall@10": {str(k): round(v, 4) for k, v in sensitivity.items()},
                "strata_2hop_plus": {
                    "n_test_genes": len(deep_genes),
                    "n_true_pathways": len(deep_truth),
                    "recall@10": round(deep_recall, 4) if deep_recall is not None else None,
                },
            }
            # Amendment 2: anchored-table Recall@10/5 vs the same test-genes-only
            # truth as the frozen evaluation (exploratory, labeled in the report).
            if anchored:
                a_by_ms = anchored["by_min_support"]
                apred = [r["node"] for r in a_by_ms[agg_cfg.min_support]]
                results["metrics"]["anchored_recall@10"] = round(
                    recall_at_k(apred, truth, 10), 4
                )
                results["metrics"]["anchored_recall@5"] = round(
                    recall_at_k(apred, truth, 5), 4
                )
                results["metrics"]["anchored_sensitivity_recall@10"] = {
                    str(ms): round(
                        recall_at_k([r["node"] for r in rows], truth, 10), 4
                    )
                    for ms, rows in a_by_ms.items()
                }
            log.info("evaluated", **results["metrics"]["sensitivity_recall@10"])

    return results


def stage_anchor(cfg, disease, cs: str, state: dict, train1_edges, pathway_ids,
                 degrees, run) -> dict | None:
    """Amendment 2 (owner 2026-08-24) EXPLORATORY stage: ontology-anchored
    pathway scores — which pathways WOULD mediate if paths were forced through
    them (the frozen unconstrained paths bypass ontology nodes). Candidates =
    GO BP / Reactome nodes with a direct train1 annotation edge to a consensus
    Top-K gene; per seed checkpoint, one shared bellmanford pass conditioned on
    (disease, relation), then the best grad-weighted D->P prefix per candidate
    via the same beam-search primitives (relation conditioning stays the
    disease-gene relation — approximation, disclosed in the report).

    Addendum 2026-08-25 (performance; exhaustive beam was ~40 h/disease):
    candidates are prefiltered to support >= anchored.min_candidate_support
    (default 3) before any GPU work — below the strictest reported min_support
    sensitivity they can never enter a reported table, so the tables are
    identical to exhaustive; the prefix weight is a Viterbi best-path estimator
    (single best path, no beam, no target edge-masking, no length averaging)."""
    if "solver" not in state:
        log.warning(
            "anchor_skipped", reason="resume without solver; run fresh for the anchored table"
        )
        return None
    import torch

    ext, agg_cfg = cfg.extraction, cfg.aggregation
    anc = getattr(cfg, "anchored", None) or {}
    min_cand_support = int(anc.get("min_candidate_support", 3))
    solver = state["solver"]
    model = solver.model.model  # unwrap torchdrug task -> NBFNet
    fact_graph = solver.model.fact_graph
    top_gene_ids = {
        gp["gene"] for gp in state["gene_paths"] if gp["gene_rank"] <= ext.top_k_genes
    }
    all_candidates = discover_anchor_candidates(train1_edges, pathway_ids, top_gene_ids)
    candidates = prefilter_candidates(all_candidates, min_cand_support)
    log.info(
        "anchor_prefilter",
        candidates=len(all_candidates),
        kept=len(candidates),
        min_candidate_support=min_cand_support,
    )
    seed_plan = state.get("seed_ckpts") or [(ext.seed, None)]
    prefix_weights: dict[int, dict[int, float]] = {}
    h_t = torch.tensor([state["disease_idx"]], device=solver.device)
    r_t = torch.tensor([state["rel_idx"]], device=solver.device)
    with run.stage(f"anchor:{cs}"):
        for s, ck in seed_plan:
            if ck is not None:
                solver.load(str(ck))
            found = extract_anchor_paths(model, fact_graph, h_t, sorted(candidates), r_t)
            # unreachable candidates (weight -inf, no path) are excluded so they
            # contribute 0 to the anchored sum via the .get(p, 0.0) default
            prefix_weights[s] = {
                p: w for p, (_, w) in found.items() if w > float("-inf")
            }
            log.info(
                "anchor_seed_done",
                seed=s,
                candidates=len(found),
                unreachable=len(found) - len(prefix_weights[s]),
            )
        by_ms = {
            ms: aggregate_anchored(
                prefix_weights, candidates, degrees, ms,
                agg_cfg.select, agg_cfg.genes_per_pathway,
            )
            for ms in sorted({3, agg_cfg.min_support, 10})
        }
        log.info(
            "anchor_aggregated",
            candidates=len(candidates),
            kept=len(by_ms[agg_cfg.min_support]),
        )
    # restore the primary checkpoint for the faithfulness control
    fck = dict(seed_plan).get(state.get("faithfulness_seed"))
    if fck is not None:
        solver.load(str(fck))
    return {
        "n_candidates": len(all_candidates),
        "n_candidates_prefiltered": len(candidates),
        "min_candidate_support": min_cand_support,
        "seeds": [s for s, _ in seed_plan],
        "by_min_support": by_ms,
    }


def _load_dataset_only(dataset_dir: Path, seed: int):
    """Dataset without solver/training — used by --resume downstream stages."""
    import random

    import numpy as np
    import torch
    from torchdrug import core
    from torchdrug.utils import comm

    sys.path.append(str(BIOPATHNET_DIR))
    from biopathnet import dataset, layer, model, task, util  # noqa: F401

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed + int(comm.get_rank()))
    cfg = util.load_config(str(BASE_CFG_TEMPLATE), context={"gpus": []})
    cfg.dataset.path = str(dataset_dir)
    cfg.engine.gpus = []
    dataset_obj = core.Configurable.load_config_dict(cfg.dataset)
    return None, dataset_obj


def rebuild_solver_state(cfg, disease, dataset_dir, ckpt, cs: str, state: dict,
                         explicit_ckpt: bool = False) -> None:
    """--resume support (addendum 2026-08-25): the extract parquets exist and
    are manifest-consistent, but the anchor / faithfulness stages still need a
    solver. Rebuild it (build_solver + load the primary seed's checkpoint) and
    populate the state keys the downstream stages need, so a killed post-extract
    run can finish anchor + faithfulness without re-running extract."""
    ext = cfg.extraction
    seed = ext.seed
    base_cs = curie_safe(disease.curie)
    seed_plan = plan_seed_checkpoints(cfg, disease, ckpt, explicit_ckpt)
    if not seed_plan:
        # open-world: the de novo pass was saved to the keep_pass1 name by the
        # original run; reload it (never retrain on resume)
        saved = REPO_ROOT / "reports" / "m2b" / "checkpoints_wave1" / f"{base_cs}_keep_pass1.pth"
        if not saved.exists():
            raise FileNotFoundError(
                f"resume needs the open-world checkpoint {saved} (written by the original run)"
            )
        seed_plan = [(seed, saved)]
    primary_seed, primary_ckpt = next(
        ((s, p) for s, p in seed_plan if s == seed), seed_plan[0]
    )
    solver, dataset_obj = build_solver(dataset_dir, seed)
    solver.load(str(primary_ckpt))
    entity_vocab = list(dataset_obj.entity_vocab)
    node_type = dataset_obj.graph.node_type.tolist()
    _, hallmark_names, crosswalk_go = load_hallmark_maps()
    classes = classify_nodes(entity_vocab, node_type, crosswalk_go)
    rel_token = supervision_relation(dataset_dir, disease.curie)
    state.update(
        {
            "solver": solver,
            "dataset_obj": dataset_obj,
            "entity_vocab": entity_vocab,
            "node_type": node_type,
            "classes": classes,
            "disease_idx": dataset_obj.inv_entity_vocab[disease.curie],
            "rel_idx": dataset_obj.inv_relation_vocab[rel_token],
            "rel_token": rel_token,
            "n_train1": dataset_obj.num_samples[0],
            "seed_ckpts": seed_plan,
            "faithfulness_seed": primary_seed,
            "ckpt_sha": sha256_file(primary_ckpt),
            "control_ranks": None,
            "hallmark_names": hallmark_names,
            "crosswalk_go": crosswalk_go,
        }
    )
    log.info("resume_solver_rebuilt", checkpoint=str(primary_ckpt), seed=primary_seed)


def stage_faithfulness(cfg, disease, cs: str, state: dict, run) -> dict | None:
    """mask_top_paths_drop: mask Top-10 scored paths per gene, re-score, the
    prediction must drop (Jain & Wallace control)."""
    if "solver" not in state:
        log.warning("faithfulness_skipped", reason="resume without solver; run fresh to check")
        return None
    import torch

    fact_graph = state["solver"].model.fact_graph
    model = state["solver"].model.model  # unwrap torchdrug task -> NBFNet
    disease_idx, rel_idx = state["disease_idx"], state["rel_idx"]
    # Amendment 1: pooled paths carry per-seed provenance; the control runs on
    # the primary seed's checkpoint, so mask only that seed's paths (paths from
    # legacy caches without a seed tag are kept).
    fseed = state.get("faithfulness_seed")
    drops = []
    with run.stage(f"faithfulness:{cs}"):
        for gp in gene_paths_top(cfg, state):
            top_paths = [
                p for p in gp["paths"] if fseed is None or p.get("seed") in (None, fseed)
            ][:10]  # Top-10 scored paths per gene
            if not top_paths:
                continue
            mask = {e for p in top_paths for e in p["edges"]}
            h_t = torch.tensor([[disease_idx]], device=fact_graph.edge_list.device)
            t_t = torch.tensor([[gp["gene"]]], device=fact_graph.edge_list.device)
            r_t = torch.tensor([[rel_idx]], device=fact_graph.edge_list.device)
            with torch.no_grad():
                base = float(torch.sigmoid(model(fact_graph, h_t, t_t, r_t)).item())
            masked = masked_gene_score(model, fact_graph, mask, disease_idx, gp["gene"], rel_idx)
            drops.append(base - masked)
        if not drops:
            # A model whose extraction yielded zero scored paths (expected for
            # shuffle-collapse controls: the shuffled model's path weights all
            # fall below min_path_score) has nothing to mask — that IS total
            # collapse, not a faithfulness violation. Loud warning, no raise.
            log.warning(
                "faithfulness_vacuous",
                reason="no paths to mask (extraction produced zero scored paths; "
                "expected for shuffle-control models — total collapse)",
            )
            return {"median_drop": None, "frac_dropped": None, "n_genes": 0,
                    "vacuous": "no paths to mask"}
        med = sorted(drops)[len(drops) // 2]
        frac_drop = sum(d > 0 for d in drops) / len(drops)
        log.info("faithfulness", median_drop=round(med, 6), frac_dropped=round(frac_drop, 3))
        if med <= 0:
            raise RuntimeError(
                f"FAITHFULNESS FAILED ({disease.curie}): median score drop after masking "
                f"Top-10 paths = {med:.6f} <= 0 — extracted paths do not drive predictions"
            )
    return {"median_drop": med, "frac_dropped": frac_drop, "n_genes": len(drops)}


def gene_paths_top(cfg, state: dict) -> list[dict]:
    k = cfg.extraction.top_k_genes
    return [gp for gp in state["gene_paths"] if gp["gene_rank"] <= k]


def stage_report(cfg, disease, dataset_dir, cs: str, state: dict, results: dict,
                 faith: dict | None, run) -> None:
    agg_cfg, ext = cfg.aggregation, cfg.extraction
    entity_vocab = state["entity_vocab"]
    names = results.get("names", {})
    primary = results["by_k"][ext.top_k_genes]

    def pname(node: int) -> str:
        token = entity_vocab[node]
        return f"{token} ({names.get(token, '')})" if names.get(token) else token

    def gname(g: int) -> str:
        token = entity_vocab[g]
        return f"{token} ({names.get(token, '')})" if names.get(token) else token

    # optional gene annotation (annotation block; skipped loudly if missing)
    annotation_note = "gene_properties.parquet not available — annotation skipped"
    gene_props: dict[str, dict] = {}
    for cand in GENE_PROPERTIES:
        if cand.exists():
            try:
                props = read_parquet_dicts(cand)
                gene_props = {str(r.get("gene_id")): r for r in props}
                annotation_note = f"annotated from {cand}"
                break
            except Exception as e:
                annotation_note = f"annotation load failed: {e}"

    def gann(g: int) -> str:
        r = gene_props.get(entity_vocab[g], {})
        bits = [str(r.get(k)) for k in ("tractability", "safety_liabilities") if r.get(k)]
        return "; ".join(bits) if bits else "n/a"

    rows_md = "\n".join(
        f"| {i} | {pname(r['node'])} | {r['score']:.4f} | {r['degree']} | {r['support']} | "
        f"{', '.join(gname(g) for g in r['top_genes'])} |"
        for i, r in enumerate(primary["pathways"], 1)
    )
    hm_md = "\n".join(
        f"| {h['id']} | {h['name']} | {h['score']:.4f} | {h['n_members']} | "
        f"{', '.join(gname(g) for g in h['genes'])} |"
        for h in primary["hallmarks"]
    ) or "| — | no hallmark member nodes hit | | | |"

    metrics = results.get("metrics")
    if metrics is None:
        bench_md = (
            f"Open-world mode (`{disease.mode}`): no valid/test split — scored by "
            "pathway support only, per config note.\n"
        )
    else:
        sens = " ".join(f"K={k}: {v}" for k, v in metrics["sensitivity_recall@10"].items())
        strata = metrics["strata_2hop_plus"]
        bench_md = f"""| metric | value |
|---|---|
| true pathways (>= {cfg.evaluation.ground_truth_min_genes} test genes) | \
{metrics['n_true_pathways']} |
| pathway Recall@10 | {metrics['pathway_recall@10']} |
| pathway Recall@5 | {metrics['pathway_recall@5']} |
| pathway MRR | {metrics['pathway_mrr']} |
| pathway Precision@3 | {metrics['pathway_precision@3']} |
| pathway NDCG@10 | {metrics['pathway_ndcg@10']} |
| degree-null Recall@10 | {metrics['degree_null_recall@10']} |
| ORA Top-100 Recall@10 | {metrics['ora_recall@10']} |
| >2-hop strata Recall@10 | {strata['recall@10']} \
({strata['n_test_genes']} genes, {strata['n_true_pathways']} true) |

K sensitivity (pre-registered Recall@10): {sens}
"""

    faith_md = (
        "not run (resume without solver)"
        if faith is None
        else (
            f"VACUOUS — {faith['vacuous']}"
            if faith.get("vacuous")
            else (
                f"median drop {faith['median_drop']:.6f}, "
                f"{faith['frac_dropped']:.1%} of {faith['n_genes']} genes dropped — PASS"
            )
        )
    )

    cr = state.get("control_ranks")
    if cr:
        per_seed_txt = ", ".join(
            f"seed {s}: rank {r}" for s, r in sorted(cr["per_seed"].items())
        )
        control_md = (
            f"Positive control: {cr['symbol']} consensus rank {cr['consensus']} "
            f"(per-seed: {per_seed_txt}) — gate = consensus Top-{ext.top_k_genes} "
            "(owner amendment 2026-08-24).\n"
        )
    else:
        control_md = ""

    # Amendment 2: EXPLORATORY anchored section, separate from the frozen
    # headline numbers above.
    anchor_heading = (
        "## EXPLORATORY — ontology-anchored pathway scores "
        "(post-hoc, not the frozen endpoint)"
    )
    anchored = results.get("anchored")
    if anchored:
        arows = anchored["by_min_support"][agg_cfg.min_support]
        anc_rows_md = "\n".join(
            f"| {i} | {pname(r['node'])} | {r['score']:.4f} | {r['degree']} | {r['support']} | "
            f"{', '.join(gname(g) for g in r['top_genes'])} |"
            for i, r in enumerate(arows, 1)
        ) or "| — | no anchored pathway passed min_support | | | |"
        if metrics is not None and "anchored_recall@10" in metrics:
            a_sens = " ".join(
                f"{ms}: {v}"
                for ms, v in metrics["anchored_sensitivity_recall@10"].items()
            )
            anc_eval_md = (
                f"Anchored Recall@10 {metrics['anchored_recall@10']}, "
                f"Recall@5 {metrics['anchored_recall@5']} vs the same test-genes-only "
                f"GO/Reactome ground truth as the frozen evaluation (exploratory). "
                f"min_support sensitivity Recall@10: {a_sens}."
            )
        else:
            anc_eval_md = (
                "Open-world mode: no held-out truth — scored by pathway support only."
            )
        anchored_md = f"""{anchor_heading}

Unconstrained top-grad paths bypass ontology nodes (frozen result above); this
table answers which pathways WOULD mediate if paths were forced through them.
The D->P prefix uses disease-gene relation conditioning as an approximation, and
the P—gene annotation edges are in train1 (circularity already disclosed below).
Prefix weights are a Viterbi best-path estimator (single best path, no beam, no
target edge-masking, no length averaging) rather than the gene-path beam
extraction (addendum 2026-08-25: exhaustive beam was computationally
infeasible, ~40 h/disease measured).

Candidates: {anchored['n_candidates']} GO BP / Reactome nodes with a direct train1
annotation edge to a consensus Top-{ext.top_k_genes} gene, prefiltered to
support >= {anchored['min_candidate_support']} -> {anchored['n_candidates_prefiltered']}
(candidates below the strictest reported min_support sensitivity can never enter
any reported table, so the tables are identical to exhaustive); seeds {anchored['seeds']}.

| rank | pathway | score | degree | support | top genes |
|---|---|---|---|---|---|
{anc_rows_md}

{anc_eval_md}
"""
    else:
        anchored_md = f"""{anchor_heading}

Not run (anchored.enabled false, or resume without solver).
"""

    report = f"""<!-- track_scope: benchmark -->
# Pathway Discovery Funnel — {disease.curie} (CHANGES/0025)

Mode `{disease.mode}`; checkpoint `{state.get('ckpt_sha', '')[:12] or 'n/a'}`; seed {ext.seed};
K = {ext.top_k_genes} (sensitivity {ext.top_k_sensitivity}); scoring = {ext.scoring}
(grad of sigmoid(score) w.r.t. path edge weights x PNA aggregation, BioKGC/BioPathNet prior art).
{control_md}
## Top-{agg_cfg.select} predicted pathways (hub-penalized: score / ln(degree), \
min_support {agg_cfg.min_support})

| rank | pathway | score | degree | support | top genes |
|---|---|---|---|---|---|
{rows_md}

## Hallmarks of Aging (reported separately, weight <= 0.15, never collapsed)

| id | name | score | member nodes hit | genes |
|---|---|---|---|---|
{hm_md}

## Faithfulness control (mask Top-10 paths per gene, re-score)

{faith_md}

## Pathway benchmark

{bench_md}

{anchored_md}
Gene annotation: {annotation_note}.
Circularity disclosure (evaluation block): BPO-GEN/PWY-GEN edges are in train1
for message passing; the test is held out at disease level, not ontology level.
Raw GO BP and Reactome are reported alongside hallmarks, per the frozen maps.
"""
    report_path = REPO_ROOT / cfg.outputs.report
    bench_path = REPO_ROOT / cfg.outputs.benchmark_report
    figure_path = REPO_ROOT / cfg.outputs.figure
    with run.stage(f"report:{cs}", outputs=[report_path, bench_path, figure_path]):
        # Shared filenames in the frozen outputs block would collide across the
        # three diseases, so each run also writes a disease-suffixed copy; the
        # verbatim configured path mirrors the latest run (manifest records it).
        for path in {report_path, bench_path}:
            path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(report, encoding="utf-8")
        report_path.with_name(f"{report_path.stem}__{cs}{report_path.suffix}").write_text(
            report, encoding="utf-8"
        )
        bench_text = (
            f"<!-- track_scope: benchmark -->\n# Pathway benchmark — {disease.curie}\n\n"
            + bench_md
            + "\nGene MRR stays the pre-registered primary metric (CHANGES/0020); "
            "pathway recall is the secondary funnel metric.\n"
        )
        bench_path.write_text(bench_text, encoding="utf-8")
        bench_path.with_name(f"{bench_path.stem}__{cs}{bench_path.suffix}").write_text(
            bench_text, encoding="utf-8"
        )
        draw_figure(figure_path, disease.curie, primary["pathways"][: agg_cfg.figure_select],
                    pname, gname)
        figure_path.with_name(f"{figure_path.stem}__{cs}{figure_path.suffix}").write_bytes(
            figure_path.read_bytes()
        )
    log.info("report_written", report=str(report_path))


def draw_figure(path: Path, curie: str, top_pathways: list[dict], pname, gname) -> None:
    """aggregation.figure_select pathways radiating to their top genes."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = max(1, len(top_pathways))
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5), squeeze=False)
    for ax, pw in zip(axes[0], top_pathways, strict=False):
        ax.scatter([0], [0], s=400, c="tomato", zorder=3)
        ax.annotate(pname(pw["node"]), (0, 0), ha="center", va="bottom", fontsize=8)
        m = max(1, len(pw["top_genes"]))
        for j, g in enumerate(pw["top_genes"]):
            x, y = math.cos(2 * math.pi * j / m), math.sin(2 * math.pi * j / m)
            ax.plot([0, x], [0, y], c="gray", lw=1, zorder=1)
            ax.scatter([x], [y], s=150, c="steelblue", zorder=2)
            ax.annotate(gname(g), (x, y), ha="center", fontsize=7)
        ax.set_title(f"score {pw['score']:.3f} | support {pw['support']}", fontsize=9)
        ax.axis("off")
    fig.suptitle(f"Pathway funnel — {curie}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--disease", required=True, help="CURIE, e.g. mondo:0010196")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=None,
        help="explicit checkpoint override (e.g. a Wave-2 seed's .pth for the "
        "cross-seed pathway-stability supplement); skips config resolution AND "
        "the consensus multi-seed plan (single-checkpoint behavior, used for "
        "the shuffle-collapse gate and debug)",
    )
    parser.add_argument(
        "--tag",
        default=None,
        help="extra output/stage suffix (e.g. seed1339) so supplement runs never "
        "collide with the pre-registered seed-1337 outputs in the shared manifest",
    )
    parser.add_argument(
        "--no-positive-control",
        action="store_true",
        help="bypass the WRN/LMNA consensus Top-100 hard gate — ONLY for "
        "shuffle-collapse control runs, where the shuffled model is expected "
        "to fail it; never for a quotable funnel run",
    )
    args = parser.parse_args()

    configure_logging()
    cfg = load_funnel_config(args.config)
    matches = [d for d in cfg.diseases if d.curie == args.disease]
    if not matches:
        raise SystemExit(
            f"disease {args.disease} not in config; choices: {[d.curie for d in cfg.diseases]}"
        )
    disease = matches[0]
    cs = curie_safe(disease.curie)
    cfg_hash = config_sha256(args.config)

    dataset_dir = resolve_dataset_dir(disease.curie, disease.mode)
    # config checkpoint: null (open_world_keep) -> (None, note); otherwise the
    # config path with the keep_pass1 on-disk fallback (resolve_checkpoint).
    ckpt, ckpt_note = resolve_checkpoint(disease.checkpoint, disease.curie)
    if args.checkpoint is not None:
        ckpt, ckpt_note = args.checkpoint, f"CLI override ({args.checkpoint.name})"
    if args.tag:
        cs = f"{cs}__{args.tag}"

    outputs = {
        "paths_parquet": cfg.outputs.paths_parquet.replace("{disease}", cs),
        "report": cfg.outputs.report,
        "figure": cfg.outputs.figure,
        "benchmark_report": cfg.outputs.benchmark_report,
        "manifest": cfg.outputs.manifest,
    }

    if args.dry_run:
        print(f"disease: {disease.curie} (mode: {disease.mode})")
        print(f"config: {args.config} (sha256 {cfg_hash[:12]})")
        print(f"dataset: {dataset_dir}")
        print(f"checkpoint: {ckpt} [{ckpt_note}]")
        print(
            f"seed: {cfg.extraction.seed} | top_k: {cfg.extraction.top_k_genes} "
            f"(sensitivity {cfg.extraction.top_k_sensitivity}) | "
            f"max_hops {cfg.extraction.max_hops}"
        )
        print("planned stages:")
        for s in STAGES:
            print(f"  - {s}:{cs}")
        print("outputs:")
        for k, v in outputs.items():
            print(f"  - {k}: {v}")
        return 0

    manifest_path = REPO_ROOT / cfg.outputs.manifest
    with Run(
        script="extract_paths",
        manifest_path=manifest_path,
        config={
            "config": str(args.config),
            "config_sha256": cfg_hash,
            "disease": disease.curie,
            "mode": disease.mode,
            "seed": cfg.extraction.seed,
        },
    ) as run:
        state: dict = {}
        extract_done = args.resume and f"extract:{cs}" in run.done_stages
        if extract_done:
            log.info("resume_from_cache", disease=disease.curie)
            state["gene_paths"] = load_extraction_cache(cfg, cs)
            # Addendum 2026-08-25: a run killed post-extract has no solver in
            # state; rebuild it so the pending anchor / faithfulness stages
            # still run (extract itself is NOT rerun — cache is manifest-checked).
            anc = getattr(cfg, "anchored", None) or {}
            need_anchor = anc.get("enabled", False) and f"anchor:{cs}" not in run.done_stages
            need_faith = f"faithfulness:{cs}" not in run.done_stages
            if need_anchor or need_faith:
                rebuild_solver_state(
                    cfg, disease, dataset_dir, ckpt, cs, state,
                    explicit_ckpt=args.checkpoint is not None,
                )
        else:
            state = stage_score_and_extract(
                cfg, args.config, disease, dataset_dir, ckpt, cs, run,
                explicit_ckpt=args.checkpoint is not None,
                skip_positive_control=args.no_positive_control,
            )

        results = stage_aggregate_evaluate(cfg, disease, dataset_dir, cs, state, run)
        faith = stage_faithfulness(cfg, disease, cs, state, run)
        stage_report(cfg, disease, dataset_dir, cs, state, results, faith, run)

    log.info("funnel_done", disease=disease.curie)
    return 0


if __name__ == "__main__":
    sys.exit(main())
