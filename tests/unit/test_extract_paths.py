"""Unit tests for scripts/pod/extract_paths.py (CHANGES/0025 funnel).

Synthetic fixtures only — no science data, no torch/GPU (spec §0.2.2).
"""

import importlib.util
import math
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "pod" / "extract_paths.py"


@pytest.fixture()
def mod():
    spec = importlib.util.spec_from_file_location("extract_paths", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# path canonicalization / hub penalty
# ---------------------------------------------------------------------------


def test_curie_safe(mod):
    assert mod.curie_safe("mondo:0010196") == "mondo__0010196"


def test_canonical_path_key_undirected(mod):
    # same path walked forward and backward canonicalizes identically
    fwd = [(0, 5, 1), (5, 9, 2)]
    rev = [(9, 5, 2), (5, 0, 1)]
    assert mod.canonical_path_key(fwd) == mod.canonical_path_key(rev)
    # edge direction within a step is normalized too
    flip = [(5, 0, 1), (9, 5, 2)]
    assert mod.canonical_path_key(fwd) == mod.canonical_path_key(flip)


def test_hub_adjust_log_degree(mod):
    assert mod.hub_adjust(10.0, 1) == 0.0  # ln(1) undefined -> 0 (documented guard)
    assert mod.hub_adjust(10.0, math.e**2) == pytest.approx(5.0)
    # hub penalty: same raw score, higher degree -> lower adjusted score
    assert mod.hub_adjust(10.0, 10_000) < mod.hub_adjust(10.0, 10)


# ---------------------------------------------------------------------------
# aggregation
# ---------------------------------------------------------------------------


def _paths(mod):
    # genes 1..6, pathway nodes 100 (specific, deg 8), 101 (hub, deg 10_000)
    degrees = {100: 8, 101: 10_000}
    gene_paths = []
    for g in range(1, 7):
        score = 0.5 if g % 2 else 0.2
        gene_paths.append(
            {"gene": g, "paths": [{"score": score, "nodes": [100, 101, 999]}]}
        )
    return gene_paths, degrees


def test_aggregate_pathways_hub_penalty_and_support(mod):
    gene_paths, degrees = _paths(mod)
    rows = mod.aggregate_pathways(gene_paths, {100, 101}, degrees, 5, 10, 3)
    assert [r["node"] for r in rows] == [100, 101]  # specificity wins via hub penalty
    # raw sums equal (2.1 each); adjusted differs by ln(degree)
    assert rows[0]["raw_score"] == pytest.approx(2.1)
    assert rows[0]["score"] == pytest.approx(2.1 / math.log(8))
    assert rows[1]["score"] == pytest.approx(2.1 / math.log(10_000))
    assert rows[0]["support"] == 6


def test_aggregate_pathways_min_support(mod):
    gene_paths, degrees = _paths(mod)
    rows = mod.aggregate_pathways(gene_paths, {100}, degrees, 7, 10, 3)  # support 6 < 7
    assert rows == []


def test_aggregate_pathways_top_genes(mod):
    gene_paths, degrees = _paths(mod)
    rows = mod.aggregate_pathways(gene_paths, {100}, degrees, 1, 10, 2)
    # genes 1,3,5 score 0.5 > genes 2,4,6 score 0.2 -> top-2 are the two lowest odd ids
    assert rows[0]["top_genes"] == [1, 3]


def test_aggregate_pathways_select_and_tiebreak(mod):
    gene_paths = [{"gene": g, "paths": [{"score": 1.0, "nodes": [100, 101]}]} for g in (1,)]
    degrees = {100: 4, 101: 4}  # identical score and degree -> id tie-break
    rows = mod.aggregate_pathways(gene_paths, {100, 101}, degrees, 1, 1, 1)
    assert rows == [rows[0]] and rows[0]["node"] == 100


# ---------------------------------------------------------------------------
# hallmarks
# ---------------------------------------------------------------------------


def test_cluster_hallmarks(mod):
    rows = [
        {"node": 100, "score": 1.0, "support": 3, "top_genes": [1, 2]},
        {"node": 101, "score": 0.5, "support": 2, "top_genes": [3]},
        {"node": 102, "score": 9.0, "support": 9, "top_genes": [4]},
    ]
    node_hm = {100: {"HALLMARK_A"}, 101: {"HALLMARK_A"}}
    out = mod.cluster_hallmarks(rows, node_hm, {"HALLMARK_A": "test hallmark"})
    assert len(out) == 1
    assert out[0]["score"] == pytest.approx(1.5)
    assert out[0]["n_members"] == 2
    assert out[0]["genes"] == [1, 2, 3]


# ---------------------------------------------------------------------------
# ground truth + metrics
# ---------------------------------------------------------------------------


def test_ground_truth_min_genes(mod):
    neighbors = {1: {100}, 2: {100}, 3: {100}, 4: {101}}
    assert mod.ground_truth_pathways({1, 2, 3, 4}, neighbors, 3) == {100}
    assert mod.ground_truth_pathways({1, 2}, neighbors, 3) == set()


def test_pathway_metrics(mod):
    truth = {10, 11, 12, 13}
    predicted = [10, 99, 11, 98, 97, 12]
    m = mod.pathway_metrics(predicted, truth)
    assert m["pathway_recall@10"] == pytest.approx(3 / 4)
    assert m["pathway_recall@5"] == pytest.approx(2 / 4)
    assert m["pathway_mrr"] == pytest.approx(1.0)
    assert m["pathway_precision@3"] == pytest.approx(2 / 3)
    ideal = sum(1 / math.log2(r + 1) for r in range(1, 5))
    dcg = 1 / math.log2(2) + 1 / math.log2(4) + 1 / math.log2(7)
    assert m["pathway_ndcg@10"] == pytest.approx(dcg / ideal)


def test_mrr_no_hit(mod):
    assert mod.mrr([1, 2, 3], {9}) == 0.0


# ---------------------------------------------------------------------------
# ORA baseline + degree null + BFS strata
# ---------------------------------------------------------------------------


def test_bh_adjust_monotone(mod):
    adj = mod.bh_adjust([0.001, 0.02, 0.5])
    assert adj[0] <= adj[1] <= adj[2]
    assert adj[0] == pytest.approx(0.003)


def test_ora_rank_enriched_first(mod):
    # pathway 100 genes {1..10} fully inside top-100; pathway 101 genes disjoint-ish
    top = set(range(1, 101))
    pathway_genes = {100: set(range(1, 11)), 101: set(range(150, 250))}
    neighbors = {}
    out = mod.ora_rank(top, {100, 101}, neighbors, pathway_genes, n_background=1000)
    assert out[0]["node"] == 100
    assert out[0]["overlap"] == 10
    assert out[0]["q_value"] < 1e-10


def test_degree_null_rank(mod):
    assert mod.degree_null_rank({1, 2, 3}, {1: 5, 2: 50, 3: 10}, 2) == [2, 3]


def test_bfs_hops(mod):
    edges = [(0, 1, 0), (1, 2, 0), (2, 3, 0)]
    assert mod.bfs_hops(5, edges, 0) == [0, 1, 2, 3, -1]


def test_gene_pathway_neighbors_undirected(mod):
    edges = [(100, 1, 0), (2, 101, 0), (1, 2, 0)]  # gene-gene edge ignored
    out = mod.gene_pathway_neighbors(edges, gene_ids={1, 2}, pathway_ids={100, 101})
    assert out == {1: {100}, 2: {101}}


# ---------------------------------------------------------------------------
# checkpoint fallback (config is frozen with wave1_seed1337 names)
# ---------------------------------------------------------------------------


def test_resolve_checkpoint_fallback(mod, tmp_path, monkeypatch):
    ckpt = tmp_path / "reports" / "m2b" / "checkpoints_wave1" / "mondo__0099999_keep_pass1.pth"
    ckpt.parent.mkdir(parents=True)
    ckpt.write_bytes(b"fake")
    monkeypatch.setattr(mod, "REPO_ROOT", tmp_path)
    path, note = mod.resolve_checkpoint(
        "reports/m2b/checkpoints_wave1/mondo__0099999_wave1_seed1337.pth", "mondo:0099999"
    )
    assert path == ckpt
    assert "fallback" in note


def test_resolve_checkpoint_open_world(mod):
    path, note = mod.resolve_checkpoint(None, "mondo:0004975")
    assert path is None and "open-world" in note


# ---------------------------------------------------------------------------
# Amendment 1 (2026-08-24): multi-seed consensus scoring
# ---------------------------------------------------------------------------


def test_seed_checkpoint_name(mod):
    # primary seed keeps the legacy run_wave1_benchmark name
    assert mod.seed_checkpoint_name("mondo__0008310", 1337, 1337) == "mondo__0008310_keep_pass1.pth"
    assert (
        mod.seed_checkpoint_name("mondo__0008310", 1339, 1337)
        == "mondo__0008310_keep_seed1339_pass1.pth"
    )


def test_discover_seed_checkpoints(mod, tmp_path):
    (tmp_path / "mondo__0008310_keep_pass1.pth").write_bytes(b"a")
    (tmp_path / "mondo__0008310_keep_seed5027_pass1.pth").write_bytes(b"b")
    found, missing = mod.discover_seed_checkpoints(
        tmp_path, "mondo__0008310", [1337, 1339, 5027], 1337
    )
    assert [s for s, _ in found] == [1337, 5027]  # config-seed order preserved
    assert missing == [1339]  # skipped with a warning by the caller


def test_consensus_mean_rank(mod):
    seed_a = [(1, 0.9), (2, 0.5), (3, 0.1)]
    seed_b = [(2, 0.7), (1, 0.6), (3, 0.2)]
    ranked = mod.consensus_mean_rank([seed_a, seed_b])
    assert [g for g, _ in ranked] == [1, 2, 3]
    assert ranked[0][1] == pytest.approx(0.75)  # mean of 0.9 and 0.6
    assert ranked[1][1] == pytest.approx(0.6)


def test_consensus_mean_rank_tiebreak_and_partial_coverage(mod):
    # equal means -> lower gene idx wins; a gene scored by one seed averages over one
    ranked = mod.consensus_mean_rank([[(5, 0.5), (2, 0.5)], [(2, 0.5), (9, 0.9)]])
    assert [g for g, _ in ranked] == [9, 2, 5]
    assert ranked[0][1] == pytest.approx(0.9)


def test_rank_of(mod):
    assert mod.rank_of([(7, 0.9), (3, 0.1)], 3) == 2
    assert mod.rank_of([(7, 0.9)], 1) is None


# ---------------------------------------------------------------------------
# Amendment 2 (2026-08-24): ontology-anchored exploratory variant
# ---------------------------------------------------------------------------


def test_discover_anchor_candidates(mod):
    edges = [
        (100, 1, 0),  # P -> g direct annotation
        (2, 101, 7),  # g -> P (inverse direction folded)
        (100, 2, 0),  # second Top-K gene for 100
        (101, 9, 0),  # gene 9 not in the consensus Top-K -> ignored
        (1, 2, 0),  # gene-gene edge ignored
        (300, 1, 0),  # non-pathway node ignored
    ]
    cand = mod.discover_anchor_candidates(edges, pathway_ids={100, 101}, top_genes={1, 2})
    assert cand == {100: {1, 2}, 101: {2}}


def test_aggregate_anchored_scoring(mod):
    candidates = {100: {1, 2, 3, 4, 5}, 101: {1, 2, 3, 4, 5}}
    prefix = {1337: {100: 0.4, 101: 0.4}, 1339: {100: 0.2, 101: 0.2}}
    degrees = {100: 8, 101: 10_000}
    rows = mod.aggregate_anchored(prefix, candidates, degrees, 5, 10, 3)
    # raw = (0.4 + 0.2) per seed-sum x 5 supporting genes = 3.0 each; hub penalty separates
    assert [r["node"] for r in rows] == [100, 101]
    assert rows[0]["raw_score"] == pytest.approx(3.0)
    assert rows[0]["score"] == pytest.approx(3.0 / math.log(8))
    assert rows[1]["score"] == pytest.approx(3.0 / math.log(10_000))
    assert rows[0]["support"] == 5
    assert rows[0]["top_genes"] == [1, 2, 3]


def test_aggregate_anchored_min_support_and_missing_prefix(mod):
    candidates = {100: {1, 2, 3}, 101: {1, 2, 3, 4, 5}}
    prefix = {1337: {101: 0.5}}  # 100 has no prefix weight -> contributes 0
    rows = mod.aggregate_anchored(prefix, candidates, {100: 8, 101: 8}, 5, 10, 5)
    assert [r["node"] for r in rows] == [101]  # support 3 < min_support 5 drops 100
    assert rows[0]["raw_score"] == pytest.approx(2.5)


def test_aggregate_anchored_select_and_tiebreak(mod):
    candidates = {100: {1, 2, 3, 4, 5}, 101: {1, 2, 3, 4, 5}}
    prefix = {1337: {100: 0.3, 101: 0.3}}
    degrees = {100: 4, 101: 4}  # identical score and degree -> id tie-break, select 1
    rows = mod.aggregate_anchored(prefix, candidates, degrees, 5, 1, 3)
    assert [r["node"] for r in rows] == [100]


# ---------------------------------------------------------------------------
# Addendum 2026-08-25: support prefilter + Viterbi best-path estimator
# ---------------------------------------------------------------------------


def test_prefilter_candidates(mod):
    cand = {100: {1, 2, 3}, 101: {1, 2}, 102: {1, 2, 3, 4, 5}}
    assert mod.prefilter_candidates(cand, 3) == {100: {1, 2, 3}, 102: {1, 2, 3, 4, 5}}
    assert mod.prefilter_candidates(cand, 1) == cand  # no-op at the floor


def _layered_graph():
    import torch

    # 6 nodes, 2 layers; edge = (h, t, r)
    layer0 = torch.tensor([[0, 1, 0], [0, 2, 0]])
    layer1 = torch.tensor([[1, 3, 0], [2, 3, 0], [1, 4, 0]])
    return [layer0, layer1], 6


def test_viterbi_best_path_known_winner(mod):
    import torch

    edge_lists, num_node = _layered_graph()
    grads = [torch.tensor([1.0, 0.5]), torch.tensor([2.0, 5.0, 0.1])]
    dists, bps = mod.viterbi_forward(edge_lists, grads, num_node, source=0)
    edges, weight = mod.viterbi_reconstruct(edge_lists, dists, bps, 0, 3)
    assert weight == pytest.approx(5.5)  # 0->2 (0.5) + 2->3 (5.0)
    assert edges == [(0, 2, 0), (2, 3, 0)]


def test_viterbi_tie_breaks_lowest_edge_index(mod):
    import torch

    edge_lists, num_node = _layered_graph()
    # both routes to 3 total 2.5 -> the lower-index edge (1->3) wins
    grads = [torch.tensor([1.0, 0.5]), torch.tensor([1.5, 2.0, 0.1])]
    dists, bps = mod.viterbi_forward(edge_lists, grads, num_node, source=0)
    edges, weight = mod.viterbi_reconstruct(edge_lists, dists, bps, 0, 3)
    assert weight == pytest.approx(2.5)
    assert edges == [(0, 1, 0), (1, 3, 0)]


def test_viterbi_prefers_shallower_depth_on_tie(mod):
    import torch

    # target 1 reachable in 1 hop (0->1 = 2.0) and 2 hops (0->2->1 = 1.0+1.0)
    edge_lists = [torch.tensor([[0, 1, 0], [0, 2, 0]]), torch.tensor([[2, 1, 0]])]
    grads = [torch.tensor([2.0, 1.0]), torch.tensor([1.0])]
    dists, bps = mod.viterbi_forward(edge_lists, grads, 3, source=0)
    edges, weight = mod.viterbi_reconstruct(edge_lists, dists, bps, 0, 1)
    assert weight == pytest.approx(2.0)
    assert edges == [(0, 1, 0)]  # shallowest depth wins the tie


def test_viterbi_disconnected_target(mod):
    import torch

    edge_lists, num_node = _layered_graph()
    grads = [torch.tensor([1.0, 0.5]), torch.tensor([2.0, 5.0, 0.1])]
    dists, bps = mod.viterbi_forward(edge_lists, grads, num_node, source=0)
    edges, weight = mod.viterbi_reconstruct(edge_lists, dists, bps, 0, 5)
    assert edges is None
    assert weight == float("-inf")
