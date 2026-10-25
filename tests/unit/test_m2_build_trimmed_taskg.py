"""Unit tests for scripts/m2_build_trimmed_taskg.py (v5.1 trap_terms key, CHANGES/0019).

Synthetic fixtures only — no science data in tests (spec §0.2.2).
"""

import importlib.util
import sys
from pathlib import Path

import polars as pl
import pytest
import yaml

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "m2_build_trimmed_taskg.py"


@pytest.fixture()
def mod():
    spec = importlib.util.spec_from_file_location("m2_build_trimmed_taskg", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


OBO = """format-version: 1.2
data-version: releases/test

[Term]
id: MONDO:0000001
name: disease

[Term]
id: MONDO:1000001
name: trap parent disease
is_a: MONDO:0000001 ! disease

[Term]
id: MONDO:1000002
name: specific disease
is_a: MONDO:1000001 ! trap parent disease
"""


def _write_fixtures(tmp_path: Path) -> dict:
    nodes = pl.DataFrame(
        {
            "node_id": ["mondo:1000001", "mondo:1000002", "efo:0000001", "ensembl:G1"],
            "node_type": ["disease", "disease", "disease", "gene"],
            "node_name": ["trap parent disease", "specific disease", "aging trait", "gene1"],
            "native_id": [None, None, None, None],
            "properties": [None, None, None, None],
            "mondo_id": ["mondo:1000001", "mondo:1000002", None, None],
        }
    )
    edges = pl.DataFrame(
        {
            "src_id": ["mondo:1000001", "mondo:1000002", "efo:0000001"],
            "dst_id": ["ensembl:G1", "ensembl:G1", "ensembl:G1"],
            "label": ["DIS-GEN", "DIS-GEN", "DIS-GEN"],
            "relation": ["ASSOCIATED_WITH"] * 3,
        }
    )
    nodes.write_parquet(tmp_path / "nodes.parquet")
    edges.write_parquet(tmp_path / "edges.parquet")
    (tmp_path / "mondo.obo").write_text(OBO, encoding="utf-8")
    fam = {"families": {"alzheimer": {"members": ["mondo:1000002"]}}}
    (tmp_path / "families.yaml").write_text(yaml.safe_dump(fam), encoding="utf-8")
    return {
        "trackd_edges": tmp_path / "edges.parquet",
        "trackd_nodes": tmp_path / "nodes.parquet",
        "mondo_obo": tmp_path / "mondo.obo",
        "families_config": tmp_path / "families.yaml",
        "splits_dir": tmp_path,
        "cluster_roots": [],
        "aging_traits": ["efo:0000001"],
        "progeria": [],
        "message_passing_labels": [],
        "families": {},
        "output_dir": tmp_path / "out",
        "report_path": tmp_path / "report.md",
        "manifest_path": tmp_path / "manifest.json",
    }


def test_trap_terms_join_cluster_and_supervision(tmp_path, mod):
    paths = _write_fixtures(tmp_path)
    cfg = mod.TrimmedConfig(**paths, trap_terms=["mondo:1000001"])
    state = mod.load_inputs(cfg)
    assert "mondo:1000001" in state["cluster"]
    assert "efo:0000001" in state["cluster"]
    # the trap term's DIS-GEN edges enter supervision (as-is pending Step 0b)
    assert state["supervision"].filter(pl.col("src_id") == "mondo:1000001").height == 1


def test_trap_terms_default_empty_and_excluded(tmp_path, mod):
    paths = _write_fixtures(tmp_path)
    cfg = mod.TrimmedConfig(**paths)  # trap_terms defaults to []
    assert cfg.trap_terms == []
    state = mod.load_inputs(cfg)
    assert "mondo:1000001" not in state["cluster"]
    assert state["supervision"].filter(pl.col("src_id") == "mondo:1000001").height == 0
    # aging trait still lands via its own list
    assert state["supervision"].filter(pl.col("src_id") == "efo:0000001").height == 1


def test_annotation_edges_orientation_lands_in_brg(mod):
    # Regression: Track-D stores annotation edges annotation→gene
    # (src = go:/hp: term, dst = ensembl: gene); a src-side gene filter
    # silently dropped ALL of them from the BRG.
    edges = pl.DataFrame(
        {
            "src_id": ["mondo:1000002", "go:0000001", "hp:0000002"],
            "dst_id": ["ensembl:G1", "ensembl:G1", "ensembl:G9"],
            "label": ["DIS-GEN", "BPO-GEN", "PHE-GEN"],
            "relation": ["ASSOCIATED_WITH", "ANNOTATES", "ANNOTATES"],
        }
    )
    brg, _drugs = mod._select_brg_edges(
        edges,
        cluster=set(),
        sup_genes={"ensembl:G1"},
        message_passing_labels=["BPO-GEN", "PHE-GEN"],
    )
    # annotation→gene edge for a supervision gene lands in the BRG
    assert (
        brg.filter((pl.col("label") == "BPO-GEN") & (pl.col("dst_id") == "ensembl:G1")).height == 1
    )
    # gene outside the supervision pool → its annotation is dropped
    assert brg.filter(pl.col("label") == "PHE-GEN").height == 0
    # the supervised label never enters the BRG
    assert brg.filter(pl.col("label") == "DIS-GEN").height == 0
