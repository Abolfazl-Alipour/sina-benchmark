"""Unit tests for sina/lodo_shuffles.py (Directives v2 B3; CHANGES/0022).

Synthetic fixtures only — no science data in tests (spec §0.2.2).
"""

from pathlib import Path

import pytest

from sina.lodo_shuffles import (
    ShuffleError,
    across_disease_shuffle,
    disease_marginal,
    gene_marginal,
    read_triples,
    verify_shuffle,
    within_disease_shuffle,
    write_triples,
)

REL = "DIS-GEN"

# Small fixture: 3 diseases, 4 genes, 7 edges. Gene marginals g1:3, g2:2,
# g3:1, g4:1; disease marginals d1:3, d2:2, d3:2.
TRIPLES = [
    ("d1", REL, "g1"),
    ("d1", REL, "g2"),
    ("d1", REL, "g3"),
    ("d2", REL, "g1"),
    ("d2", REL, "g4"),
    ("d3", REL, "g1"),
    ("d3", REL, "g2"),
]
POOL = ["d1", "d2", "d3", "d4", "d5"]


def test_across_disease_preserves_gene_marginal():
    out = across_disease_shuffle(TRIPLES, POOL, forbidden=set(), seed=1337)
    assert gene_marginal(out) == gene_marginal(TRIPLES)
    assert len(out) == len(TRIPLES)


def test_within_disease_preserves_disease_marginal():
    out = within_disease_shuffle(TRIPLES, forbidden=set(), seed=1337)
    assert disease_marginal(out) == disease_marginal(TRIPLES)
    assert len(out) == len(TRIPLES)


def test_within_disease_uses_only_unit_genes():
    out = within_disease_shuffle(TRIPLES, forbidden=set(), seed=1337)
    assert {t[2] for t in out} <= {t[2] for t in TRIPLES}


def test_determinism_same_seed_identical():
    a = across_disease_shuffle(TRIPLES, POOL, forbidden=set(), seed=1337)
    b = across_disease_shuffle(TRIPLES, POOL, forbidden=set(), seed=1337)
    assert a == b
    c = within_disease_shuffle(TRIPLES, forbidden=set(), seed="1337:u:m")
    d = within_disease_shuffle(TRIPLES, forbidden=set(), seed="1337:u:m")
    assert c == d


def test_determinism_input_order_independent():
    shuffled_input = list(reversed(TRIPLES))
    a = across_disease_shuffle(TRIPLES, POOL, forbidden=set(), seed=7)
    b = across_disease_shuffle(shuffled_input, POOL, forbidden=set(), seed=7)
    assert a == b


def test_different_seed_different_assignment():
    # Large enough fixture that two seeds coincide only with negligible chance;
    # seeds 1 and 2 are fixed, so this is itself deterministic.
    triples = [(f"d{i % 6}", REL, f"g{i}") for i in range(60)]
    pool = [f"d{i}" for i in range(10)]
    a = across_disease_shuffle(triples, pool, forbidden=set(), seed=1)
    b = across_disease_shuffle(triples, pool, forbidden=set(), seed=2)
    assert a != b
    c = within_disease_shuffle(triples, forbidden=set(), seed=1)
    d = within_disease_shuffle(triples, forbidden=set(), seed=2)
    assert c != d


def test_no_duplicate_tuples():
    # Tight pool relative to gene degrees to force resampling pressure.
    triples = [("d1", REL, "g1")] * 1 + [("d2", REL, "g1"), ("d3", REL, "g1")]
    out = across_disease_shuffle(triples, ["d1", "d2", "d3"], forbidden=set(), seed=3)
    assert len(out) == len(set(out))
    out_w = within_disease_shuffle(TRIPLES, forbidden=set(), seed=3)
    assert len(out_w) == len(set(out_w))


def test_forced_resampling_avoids_valid_test_collision():
    # Naive reassignment could keep (d1, REL, g1); it is "held out" here, so
    # the only legal target in the pool is d2 — resampling must find it.
    triples = [("d1", REL, "g1")]
    forbidden = {("d1", REL, "g1")}
    out = across_disease_shuffle(triples, ["d1", "d2"], forbidden, seed=1337)
    assert out == [("d2", REL, "g1")]
    # within_disease analogue: g1 held out, only g2 remains legal.
    triples_w = [("d1", REL, "g1"), ("d2", REL, "g2")]
    forbidden_w = {("d1", REL, "g1")}
    out_w = within_disease_shuffle(triples_w, forbidden_w, seed=1337)
    assert ("d1", REL, "g1") not in out_w
    assert set(out_w).isdisjoint(forbidden_w)


def test_unresolvable_collision_fails_loudly():
    triples = [("d1", REL, "g1")]
    forbidden = {("d1", REL, "g1")}
    with pytest.raises(ShuffleError):
        across_disease_shuffle(triples, ["d1"], forbidden, seed=1, max_attempts=5)
    with pytest.raises(ShuffleError):
        within_disease_shuffle(triples, forbidden, seed=1, max_attempts=5)


def test_output_sorted_and_byte_stable(tmp_path):
    out = across_disease_shuffle(TRIPLES, POOL, forbidden=set(), seed=9)
    assert out == sorted(out)
    p1, p2 = tmp_path / "a.txt", tmp_path / "b.txt"
    write_triples(out, p1)
    write_triples(out, p2)
    assert p1.read_bytes() == p2.read_bytes()
    assert read_triples(p1) == out


def test_verify_shuffle_rows_and_asserts():
    out = across_disease_shuffle(TRIPLES, POOL, forbidden=set(), seed=5)
    row = verify_shuffle(TRIPLES, out, "across_disease", set())
    assert row["marginal_checked"] == "gene_degree"
    assert row["degree_marginal_preserved"] is True
    assert row["n_duplicate_tuples"] == 0
    assert row["valid_test_disjoint"] is True
    assert row["edges_before"] == row["edges_after"] == len(TRIPLES)

    out_w = within_disease_shuffle(TRIPLES, forbidden=set(), seed=5)
    row_w = verify_shuffle(TRIPLES, out_w, "within_disease", set())
    assert row_w["marginal_checked"] == "disease_degree"

    # Broken marginal / duplicates / collisions are hard assertion failures.
    with pytest.raises(AssertionError):
        verify_shuffle(TRIPLES, out[:-1], "across_disease", set())
    with pytest.raises(AssertionError):
        verify_shuffle(TRIPLES, [*out, out[0]], "across_disease", set())
    with pytest.raises(AssertionError):
        verify_shuffle(TRIPLES, out, "across_disease", {out[0]})
    with pytest.raises(ValueError, match="unknown shuffle mode"):
        verify_shuffle(TRIPLES, out, "bogus", set())


# --- script-level dir naming (source/output suffix knobs, CHANGES/0024 KEEP) ---

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "m2_lodo_shuffles.py"


@pytest.fixture()
def script_mod():
    import importlib.util
    import sys

    spec = importlib.util.spec_from_file_location("m2_lodo_shuffles", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_shuffle_dir_suffixes(script_mod, tmp_path):
    cfg = script_mod.ShuffleConfig(
        exports_dir=tmp_path,
        units=["mondo:0010196"],
        manifest_path=tmp_path / "m.json",
        report_path=tmp_path / "r.md",
        source_suffix="_keep_disphe",
        output_suffix="_keep",
    )
    src = script_mod.shuffle_src_dir(cfg, "mondo:0010196")
    out = script_mod.shuffle_out_dir(cfg, "mondo:0010196", "across_disease")
    assert src == tmp_path / "mondo__0010196_keep_disphe"
    assert out == tmp_path / "mondo__0010196_keep_shuffle_across_disease"


def test_shuffle_dir_defaults_unchanged(script_mod, tmp_path):
    cfg = script_mod.ShuffleConfig(
        exports_dir=tmp_path,
        units=["mondo:0010196"],
        manifest_path=tmp_path / "m.json",
        report_path=tmp_path / "r.md",
    )
    assert script_mod.shuffle_src_dir(cfg, "mondo:0010196") == tmp_path / "mondo__0010196"
    assert script_mod.shuffle_out_dir(cfg, "mondo:0010196", "within_disease") == (
        tmp_path / "mondo__0010196_shuffle_within_disease"
    )
