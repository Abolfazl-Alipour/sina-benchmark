"""Unit tests for m2c_build_full_exports (synthetic frames, no trackd needed)."""
import importlib.util
import sys
from pathlib import Path

import polars as pl

SCRIPT = (Path(__file__).resolve().parent.parent.parent
          / "scripts" / "pod" / "m2c_build_full_exports.py")
spec = importlib.util.spec_from_file_location("full_exports", SCRIPT)
mod = importlib.util.module_from_spec(spec)
sys.modules["full_exports"] = mod
spec.loader.exec_module(mod)


def frame():
    return pl.DataFrame({
        "src_id": ["mondo:1", "mondo:1", "mondo:2", "mondo:2", "mondo:1"],
        "relation": ["ASSOCIATED_WITH"] * 5,
        "dst_id": ["ensembl:A", "ensembl:B", "ensembl:A", "ensembl:C",
                   "ensembl:A"],
        "label": ["DIS-GEN", "DIS-GEN", "PHE-GEN", "DIS-GEN", "DIS-GEN"],
        "sources_direct": [["OPEN_TARGETS"], ["DISGENET", "OPEN_TARGETS"],
                           None, ["DISGENET"], ["OPEN_TARGETS"]],
    })


def test_drop_encumbered_tags_disgenet_only():
    tagged, keep = mod.drop_encumbered(frame())
    assert tagged.height == 2
    assert keep.height == 3
    assert "sources_direct" not in keep.columns


def test_drop_encumbered_null_sources_kept():
    _, keep = mod.drop_encumbered(frame())
    assert keep.filter(pl.col("dst_id") == "ensembl:A",
                       pl.col("label") == "PHE-GEN").height == 1


def test_unique_triples_dedupes():
    _, keep = mod.drop_encumbered(frame())
    ded = mod.unique_triples(keep)
    assert ded.height == 2  # (mondo:1,A) x2 -> 1, plus (mondo:2,C)


def test_brg_label_set_excludes_disgen():
    assert "DIS-GEN" not in mod.BRG_LABELS
    assert "ANA-GEN" not in mod.BRG_LABELS
    for lab in ("GEN-GEN", "PHE-GEN", "DRG-DIS", "DIS-PHE", "DIS-DIS"):
        assert lab in mod.BRG_LABELS


def test_type_ids_cover_vocab():
    for t in ("disease", "gene", "drug", "phenotype", "biological_process",
              "cellular_component", "molecular_function", "pathway",
              "anatomy", "exposure"):
        assert t in mod.TYPE_IDS
