"""Unit tests for scripts/m2_lodo_audits.py (v5.1 Step 2 audit computations).

Synthetic fixtures only — no science data in tests (spec §0.2.2).
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "m2_lodo_audits.py"


@pytest.fixture()
def mod():
    spec = importlib.util.spec_from_file_location("m2_lodo_audits", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


# --- overlap math (v5.1 §5.3) ------------------------------------------------


def test_overlap_coefficient_nested_sets(mod):
    # containment: small key fully nested in a large one → coefficient 1.0
    # even though Jaccard is small (the size-asymmetry case §5.3 calls out)
    a = {"g1", "g2", "g3"}
    b = {"g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9", "g10"}
    assert mod.overlap_coefficient(a, b) == pytest.approx(1.0)
    assert mod.overlap_coefficient(b, a) == pytest.approx(1.0)
    assert mod.jaccard(a, b) == pytest.approx(0.3)


def test_overlap_coefficient_partial_and_empty(mod):
    a = {"g1", "g2", "g3", "g4"}
    b = {"g1", "g2"}
    assert mod.overlap_coefficient(a, b) == pytest.approx(2 / 2)
    assert mod.jaccard(a, b) == pytest.approx(2 / 4)
    assert mod.overlap_coefficient(a, set()) == 0.0
    assert mod.overlap_coefficient(set(), set()) == 0.0
    assert mod.jaccard(set(), set()) == 0.0


def test_overlap_band_thresholds(mod):
    assert mod.overlap_band(0.3) == "co_holdout"
    assert mod.overlap_band(0.5) == "co_holdout"
    assert mod.overlap_band(0.15) == "disclose"
    assert mod.overlap_band(0.299) == "disclose"
    assert mod.overlap_band(0.149) == "ok"
    assert mod.overlap_band(0.0) == "ok"


# --- key-size ceiling classification (v5.1 §5.4) -----------------------------


def test_key_size_band_boundaries(mod):
    assert mod.key_size_band(1.99) == "below_2pct"
    assert mod.key_size_band(2.0) == "theta_range_2_5pct"
    assert mod.key_size_band(5.0) == "theta_range_2_5pct"
    assert mod.key_size_band(5.01) == "above_5pct"
    assert mod.key_size_band(12.3) == "above_5pct"  # aging reference


def test_evidence_class_floors(mod):
    assert mod.evidence_class(50) == "primary"
    assert mod.evidence_class(49) == "secondary"
    assert mod.evidence_class(20) == "secondary"
    assert mod.evidence_class(19) == "below_floor"
    assert mod.evidence_class(0) == "below_floor"


# --- BFS separation statistic ------------------------------------------------


def test_separation_from_distances(mod):
    # answers all within 2 hops; non-answers mostly far/unreachable
    ans = np.array([1.0, 2.0, 2.0, np.inf])
    non = np.array([3.0, 4.0, np.inf, np.inf, 2.0, 5.0])
    out = mod.separation_from_distances(ans, non, ks=(2, 3), cap=5)
    assert out["pct_answers_le2"] == pytest.approx(75.0)
    assert out["pct_nonanswers_le2"] == pytest.approx(100.0 / 6.0)
    assert out["separation_le2"] == pytest.approx(75.0 - 100.0 / 6.0)
    # unreachable (inf) counts as outside the cap
    assert out["pct_answers_within_cap"] == pytest.approx(75.0)
    assert out["pct_nonanswers_within_cap"] == pytest.approx(400.0 / 6.0)


def test_separation_collapsed_key(mod):
    # mega-key behaviour: answers and non-answers nearly identical profiles
    rng = np.random.default_rng(0)
    ans = rng.integers(1, 4, size=1000).astype(float)
    non = rng.integers(1, 4, size=1000).astype(float)
    out = mod.separation_from_distances(ans, non)
    assert abs(out["separation_le3"]) < 10.0  # no topological separation


# --- noise-floor simulation (v5.1 §6.1) --------------------------------------


def test_noise_floor_shape_and_seed(mod):
    n_genes = 1000
    d1 = mod.noise_floor_draws(n_genes, key_size=50, resamples=200, seed=1337)
    d2 = mod.noise_floor_draws(n_genes, key_size=50, resamples=200, seed=1337)
    assert d1.shape == (200,)
    np.testing.assert_array_equal(d1, d2)  # seeded → reproducible
    assert ((d1 >= 0) & (d1 <= 1)).all()
    # mean ≈ theoretical hit probability 10 / (n_genes - key_size + 1)
    assert d1.mean() == pytest.approx(10 / (n_genes - 50 + 1), rel=0.3)


def test_noise_floor_cv_shrinks_with_key_size(mod):
    # the variance-shape claim: relative scatter inflates at small key sizes
    n_genes, resamples = 5000, 500
    cvs = []
    for k in (20, 100, 500):
        draws = mod.noise_floor_draws(n_genes, k, resamples, seed=1337)
        cvs.append(draws.std(ddof=1) / draws.mean())
    assert cvs[0] > cvs[1] > cvs[2]
