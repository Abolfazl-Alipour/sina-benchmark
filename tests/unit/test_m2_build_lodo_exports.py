"""Unit tests for scripts/m2_build_lodo_exports.py (v5.1 Phase F).

Synthetic fixtures only — no science data in tests (spec §0.2.2).
"""

import importlib.util
import sys
from pathlib import Path

import polars as pl
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "m2_build_lodo_exports.py"


@pytest.fixture()
def mod():
    spec = importlib.util.spec_from_file_location("m2_build_lodo_exports", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


# --- ancestor ID normalization + cluster ancestor mapping --------------------


def test_normalize_ancestor_id(mod):
    assert mod.normalize_ancestor_id("MONDO_0019303") == "mondo:0019303"
    assert mod.normalize_ancestor_id("EFO_0000508") == "efo:0000508"
    assert mod.normalize_ancestor_id("OBA_VT0005372") == "oba:VT0005372"  # suffix verbatim
    assert mod.normalize_ancestor_id("orphanet:75496") == "orphanet:75496"


def test_cluster_ancestors_trap_parents(mod):
    # Werner-like fixture: member d1 has two strip-able trap parents plus one
    # unmapped raw; member d2 adds a third trap parent (the HGPS pattern).
    ancestors_by_node = {
        "mondo:d1": ["MONDO_0000001", "MONDO_0000002", "OTAR_0000018"],
        "mondo:d2": ["MONDO_0000003", "MONDO_0000001", "MONDO_00D2"],  # d2 is its own mate
    }
    node_ids = {"mondo:d1", "mondo:d2", "mondo:0000001", "mondo:0000002", "mondo:0000003"}
    mapped, unmapped = mod.cluster_ancestors(["mondo:d1", "mondo:d2"], ancestors_by_node, node_ids)
    assert set(mapped) == {"mondo:0000001", "mondo:0000002", "mondo:0000003"}
    assert mapped["mondo:0000001"] == {"mondo:d1", "mondo:d2"}  # shared parent, both sources
    assert mapped["mondo:0000003"] == {"mondo:d2"}
    assert unmapped == ["MONDO_00D2", "OTAR_0000018"]  # mate + non-node excluded


def test_cluster_ancestors_excludes_mates(mod):
    # familial-AD pattern: a member's ancestor is itself a cluster mate ->
    # skipped entirely (withheld with the cluster, not stripped, not unmapped)
    ancestors_by_node = {"a": ["MONDO_b"], "mondo:b": []}
    mapped, unmapped = mod.cluster_ancestors(["a", "mondo:b"], ancestors_by_node, {"a", "mondo:b"})
    assert mapped == {}
    assert unmapped == []


# --- overlap-rule sibling strip + disclose band ------------------------------


def _pair(a, b, coef, inter):
    return {"disease_a": a, "disease_b": b, "overlap_coefficient": coef, "intersection": inter}


def test_overlap_flagged_pairs_thresholds(mod):
    rows = [
        _pair("U", "s1", 0.5, 25),  # strip
        _pair("s2", "U", 0.3, 20),  # strip (boundary, reversed orientation)
        _pair("U", "psp", 1.0, 1),  # below intersection qualifier -> nothing
        _pair("U", "d1", 0.29, 30),  # disclose band
        _pair("U", "d2", 0.15, 20),  # disclose band (boundary)
        _pair("U", "ok", 0.149, 50),  # below band -> nothing
        _pair("U", "d3", 0.2, 19),  # band coef but tiny intersection -> nothing
        _pair("U", "mate", 1.0, 100),  # cluster mate -> excluded
        _pair("U", "anc", 1.0, 100),  # ancestor -> excluded
        _pair("U", "U", 1.0, 100),  # diagonal self-pair -> excluded
        _pair("x", "y", 1.0, 100),  # does not touch U -> excluded
    ]
    strip, disclose = mod.overlap_flagged_pairs("U", rows, 0.3, 20, 0.15, exclude={"mate", "anc"})
    assert [s["curie"] for s in strip] == ["s1", "s2"]
    assert [d["curie"] for d in disclose] == ["d1", "d2"]
    assert strip[0] == {"curie": "s1", "overlap_coefficient": 0.5, "intersection": 25}


def test_overlap_flagged_pairs_psp_case(mod):
    # The frozen PSP-parkinsonism row: coefficient 1.0, intersection 1 — the
    # mechanical rule must NOT catch it (hence the explicit D2 list).
    strip, disclose = mod.overlap_flagged_pairs(
        "mondo:0008199",
        [_pair("mondo:0008199", "mondo:0009839", 1.0, 1)],
        0.3,
        20,
        0.15,
        exclude=set(),
    )
    assert strip == [] and disclose == []


# --- valid/test split on evidence groups -------------------------------------


def _sup_frame():
    # 8 edges, 5 groups; g1 twin pair (d->x and x->d, shared group ID),
    # g2 has two duplicate rows (shared group ID).
    rows = [
        ("d", "x", "r", "g1"),
        ("x", "d", "r", "g1"),  # reciprocal twin of row 1
        ("d", "y", "r", "g2"),
        ("d", "y", "r", "g2"),  # multi-source duplicate
        ("d", "z", "r", "g3"),
        ("d", "w", "r", "g4"),
        ("d", "v", "r", "g5"),
        ("d", "u", "r", "g6"),
    ]
    return pl.DataFrame(
        rows, schema=["src_id", "dst_id", "relation", "evidence_group_id"], orient="row"
    )


def test_split_valid_test_twins_coassigned(mod):
    edges = _sup_frame()
    valid, test = mod.split_valid_test(edges, 0.4, 1337)
    # twins and duplicates share a partition
    for frame in (valid, test):
        g1 = frame.filter(pl.col("evidence_group_id") == "g1")
        assert g1.height in (0, 2)
        g2 = frame.filter(pl.col("evidence_group_id") == "g2")
        assert g2.height in (0, 2)
    # partitions cover everything, disjoint by group
    assert valid.height + test.height == edges.height
    vg = set(valid["evidence_group_id"].to_list())
    tg = set(test["evidence_group_id"].to_list())
    assert vg.isdisjoint(tg)
    assert vg | tg == {"g1", "g2", "g3", "g4", "g5", "g6"}


def test_split_valid_test_deterministic(mod):
    edges = _sup_frame()
    v1, t1 = mod.split_valid_test(edges, 0.4, 1337)
    v2, t2 = mod.split_valid_test(edges, 0.4, 1337)
    assert v1.equals(v2) and t1.equals(t2)


def test_split_valid_test_straddle_raises(mod):
    # a group ID forced into both partitions trips the leakage guard
    edges = _sup_frame()
    with pytest.raises(Exception, match="evidence group"):
        mod.assert_split_disjoint({"valid": edges, "test": edges.head(2)})


# --- zero-vector guard --------------------------------------------------------


def test_zero_vector_genes(mod):
    fact_ep = {"g1", "g2", "d1"}
    assert mod.zero_vector_genes(fact_ep, {"g1", "g2", "g3"}) == {"g3"}
    assert mod.zero_vector_genes(fact_ep, {"g1"}) == set()


# --- misc ---------------------------------------------------------------------


def test_curie_safe(mod):
    assert mod.curie_safe("mondo:0010196") == "mondo__0010196"
    assert mod.curie_safe("orphanet:75496") == "orphanet__75496"
