"""Unit tests for m4b/tabular_baseline.py (pure functions, no network)."""
import importlib.util
import sys
from pathlib import Path

import pandas as pd

SCRIPT = (Path(__file__).resolve().parent.parent.parent
          / "scripts" / "m4b" / "tabular_baseline.py")
spec = importlib.util.spec_from_file_location("tabular_baseline", SCRIPT)
mod = importlib.util.module_from_spec(spec)
sys.modules["tabular_baseline"] = mod
spec.loader.exec_module(mod)


def test_parse_chrpos_token():
    assert mod.parse_chrpos_token("rs1;2:159516822:A:T:rs2") == ("2", 159516822)
    assert mod.parse_chrpos_token("rs214034") is None
    assert mod.parse_chrpos_token("") is None


def test_first_rsid():
    assert mod.first_rsid("rs214034;rs11295692") == "rs214034"
    assert mod.first_rsid("2:159516822:A:T") is None


def test_chr_folds_cover_disjoint_balanced():
    df = pd.DataFrame({
        "chr": [c for c in ["1", "2", "3", "4", "5"] for _ in range(4)],
        "label": [1, 0, 0, 0] * 5,
        "locus_id": [f"{c}:{i}" for c in "12345" for i in range(4)],
    })
    folds = mod.chr_folds(df, n_folds=5, seed=1337)
    assert len(folds) == len(df)
    assert set(folds) == {0, 1, 2, 3, 4}
    # same chromosome always lands in one fold
    assert (df.assign(f=folds).groupby("chr")["f"].nunique() == 1).all()
    # deterministic
    assert (mod.chr_folds(df, n_folds=5, seed=1337) == folds).all()


def test_locus_metrics_pessimistic_ties():
    frame = pd.DataFrame({
        "locus_id": ["a"] * 4,
        "label": [0, 1, 0, 0],
        "s": [0.9, 0.9, 0.9, 0.1],  # positive tied for first
    })
    m = mod.locus_metrics(frame, "s")
    assert m["ranks"] == [3.0]  # worst tied position, not 1
    assert m["top1"] == [0.0]
    assert m["top5"] == [1.0]


def test_locus_metrics_clear_winner():
    frame = pd.DataFrame({
        "locus_id": ["a"] * 3,
        "label": [0, 1, 0],
        "s": [0.1, 0.9, 0.2],
    })
    m = mod.locus_metrics(frame, "s")
    assert m["ranks"] == [1.0]
    assert m["top1"] == [1.0]


def test_bootstrap_ci_shape():
    out = mod.bootstrap_ci([1.0, 0.0, 1.0, 1.0], n_boot=50, seed=1337)
    assert out["mean"] == 0.75
    assert out["n"] == 4
    lo, hi = out["ci"]
    assert lo <= 0.75 <= hi


def test_auprc_degenerate():
    assert pd.isna(mod.auprc([0, 0], [0.1, 0.2]))
    assert 0.0 <= mod.auprc([0, 1], [0.1, 0.9]) <= 1.0


def test_flames_headtohead(tmp_path):
    tsv = tmp_path / "bench.tsv"
    tsv.write_text(
        "filename\tTP\tXGB_scaled\tOverall L2G score\n"
        "locA\t1\t0.9\t0.2\n"
        "locA\t0\t0.1\t0.8\n"
        "locB\t0\t0.2\t0.1\n"
        "locB\t1\t0.8\t0.9\n", encoding="utf-8")
    out = mod.flames_headtohead(tsv, ["XGB_scaled"])
    assert out["n_loci"] == 2
    assert out["arms"]["XGB_scaled"]["top1"]["mean"] == 1.0
    assert out["arms"]["Overall L2G score"]["top1"]["mean"] == 0.5
    assert out["arms"]["XGB_scaled"]["diff_vs_ot"]["mean"] == 0.5


def test_vep_severity_ladder():
    assert mod.VEP_SEV["stop_gained"] == 3
    assert mod.VEP_SEV["missense_variant"] == 2
    assert mod.VEP_SEV["synonymous_variant"] == 1
    assert mod.VEP_SEV.get("intron_variant", 0) == 0
    assert mod.VEP_SEV.get("intergenic_variant", 0) == 0


def test_store_mapping_picks_grch38():
    cache: dict = {}
    mod._store_mapping("rs1", {"mappings": [
        {"assembly_name": "GRCh37", "seq_region_name": "2", "start": 5},
        {"assembly_name": "GRCh38", "seq_region_name": "2", "start": 9},
    ]}, cache)
    assert cache["rs1"] == {"chr": "2", "pos": 9}
    mod._store_mapping("rs2", {"mappings": []}, cache)
    assert cache["rs2"] is None


def test_featuresets_v3_adds_network():
    assert set(mod.FEATURE_SETS["v3"]) == {
        "dist_kb", "degree", "disgen_degree", "vep_sev",
        "phe_degree", "go_degree", "ppr"}


def test_locus_trait_map_skips_bad_ontology(tmp_path, monkeypatch):
    tsv = tmp_path / "gold.tsv"
    tsv.write_text(
        "trait_info.ontology\tsentinel_variant.locus_GRCh38.chromosome\t"
        "sentinel_variant.locus_GRCh38.position\n"
        "EFO_0000305\t1\t100\n"
        "NOUNDERSCORE\t2\t200\n"
        "EFO_\t3\t300\n", encoding="utf-8")
    monkeypatch.setattr(mod, "OT_GOLD_TSV", tsv)
    out = mod.locus_trait_map()
    assert out == {("1", 100): "efo:0000305"}


def test_featuresets_v2_drops_constraint():
    assert "loeuf" in mod.FEATURE_SETS["v1"]
    assert "loeuf" not in mod.FEATURE_SETS["v2"]
    assert "disgen_degree" in mod.FEATURE_SETS["v2"]
    assert "vep_sev" in mod.FEATURE_SETS["v2"]


def test_paired_top1_sign():
    frame = pd.DataFrame({
        "locus_id": ["a", "a", "b", "b"],
        "label": [1, 0, 1, 0],
        "gbm": [0.9, 0.1, 0.9, 0.1],
        "ot": [0.1, 0.9, 0.9, 0.1],
    })
    d = mod.paired_top1(frame, "gbm", "ot", seed=1337)
    assert d["mean"] == 0.5 and d["n"] == 2
