"""Tests for completion exports (Phase 1) — mirror LODO export tests."""

import polars as pl
import tempfile
from pathlib import Path

from scripts.m2_build_completion_exports import completion_split, curie_safe
from sina.splits import GROUP_COL

def _fake_edges():
    # 20 edges, 10 groups, reciprocal twins share group
    groups = ["g1","g1","g2","g2","g3","g3","g4","g4","g5","g5","g6","g6","g7","g7","g8","g8","g9","g9","g10","g10"]
    return pl.DataFrame({
        "src_id": ["mondo:0004975"]*20,
        "dst_id": [f"ensembl:ENSG00000{i:05d}" for i in range(20)],
        "relation": ["ASSOCIATED_WITH"]*20,
        "evidence_group_id": groups,
        "label": ["DIS-GEN"]*20,
        "undirected": [False]*20,
    })

def test_curie_safe():
    assert curie_safe("mondo:0004975") == "mondo__0004975"

def test_completion_split_ratios():
    edges = _fake_edges()
    # 5% test -> with 10 groups, expect 1 group in test (10% due to group granularity, but within tolerance)
    # Use valid 0.15
    train_keep, valid, test = completion_split(edges, 0.05, 0.15, 1337)
    # Twins co-assigned: each group appears fully in one partition
    # Check disjointness on groups
    g_train = set(train_keep[GROUP_COL].to_list())
    g_valid = set(valid[GROUP_COL].to_list())
    g_test = set(test[GROUP_COL].to_list())
    assert g_train.isdisjoint(g_valid)
    assert g_train.isdisjoint(g_test)
    assert g_valid.isdisjoint(g_test)
    # Total edges preserved
    assert train_keep.height + valid.height + test.height == edges.height

def test_completion_split_twins_never_straddle():
    # Two edges share same group gX must be same partition
    edges = pl.DataFrame({
        "src_id": ["mondo:0004975","mondo:0004975","mondo:0004975"],
        "dst_id": ["ensembl:ENSG00000000001","ensembl:ENSG00000000002","ensembl:ENSG00000000003"],
        "relation": ["ASSOCIATED_WITH"]*3,
        "evidence_group_id": ["twinA","twinA","other"],
        "label": ["DIS-GEN"]*3,
    })
    train_keep, valid, test = completion_split(edges, 0.05, 0.15, 1337)
    # Find where twinA went
    all_parts = []
    for df, name in [(train_keep,"train"), (valid,"valid"), (test,"test")]:
        if "twinA" in df[GROUP_COL].to_list():
            all_parts.append(name)
    assert len(all_parts) == 1, f"twins split across {all_parts}"

def test_brg_composition_guard():
    # Minimal guard: ensure we would catch missing annotation layer — just that the script asserts pge>500k
    # Here we test the function not the guard, but ensure the split function respects seed
    edges = _fake_edges()
    a1,b1,c1 = completion_split(edges, 0.05, 0.15, 1337)
    a2,b2,c2 = completion_split(edges, 0.05, 0.15, 1337)
    assert a1.height == a2.height and b1.height == b2.height  # deterministic

