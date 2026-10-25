"""Unit tests for scripts/m2_drug_path_measurement.py (v5.1 Step 3).

Synthetic fixtures only — no science data in tests (spec §0.2.2).
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "m2_drug_path_measurement.py"


@pytest.fixture()
def mod():
    spec = importlib.util.spec_from_file_location("m2_drug_path_measurement", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


# --- drug-mediated distance + F@k on a tiny fixture graph --------------------
#
# Fixture: disease d — drug m — gene g1;  m — g2;  d — h (hub) — g3
# BFS-from-adjacent-drugs distances to [g1, g2, g3]: [1, 1, 2]
# drug-mediated distances = 1 + that: [2, 2, 3]


def test_drug_mediated_distance_rule(mod):
    dist_from_drugs = np.array([1.0, 1.0, 2.0, np.inf])
    dmed = mod.drug_mediated_distances(dist_from_drugs)
    np.testing.assert_array_equal(dmed[:3], [2.0, 2.0, 3.0])
    assert np.isinf(dmed[3])  # unreachable stays unreachable


def test_fk_fraction(mod):
    # genes: g1 (answer), g2 (answer), g3 (answer), g4 (non-answer)
    gene_dist = np.array([2.0, 2.0, 3.0, 2.0])  # drug-mediated distances
    ans = np.array([True, True, True, False])
    assert mod.fk_fraction(gene_dist, ans, 2) == pytest.approx(2 / 3)
    assert mod.fk_fraction(gene_dist, ans, 3) == pytest.approx(1.0)
    assert mod.fk_fraction(gene_dist, ans, 4) == pytest.approx(1.0)
    assert mod.fk_fraction(gene_dist, np.zeros(4, dtype=bool), 3) == 0.0


def test_reachable_mask_inf(mod):
    mask = mod.reachable_mask(np.array([0.0, 2.0, np.inf]), 2)
    np.testing.assert_array_equal(mask, [True, True, False])


# --- exclusivity --------------------------------------------------------------


def test_exclusivity_fraction(mod):
    drug_reach = np.array([True, True, True, False])  # 3 answers drug-reachable@3
    nondrug_dist = np.array([2.0, 5.0, np.inf, 1.0])  # only g1 also reachable <=3
    n, frac = mod.exclusivity_fraction(drug_reach, nondrug_dist, 3)
    assert n == 3
    assert frac == pytest.approx(1 / 3)
    n0, frac0 = mod.exclusivity_fraction(np.zeros(4, dtype=bool), nondrug_dist, 3)
    assert (n0, frac0) == (0, 0.0)


# --- source classification ----------------------------------------------------


def test_classify_direct_source(mod):
    assert mod.classify_direct_source(["OPEN_TARGETS"]) == "open_targets_only"
    assert mod.classify_direct_source(["DISGENET", "OPEN_TARGETS"]) == "disgenet+open_targets"
    assert mod.classify_direct_source(["DISGENET"]) == "disgenet_only"
    assert mod.classify_direct_source(["SOMETHING_ELSE"]) == "other"
    assert mod.classify_direct_source([]) == "other"


def test_drug_clinical_source_rule(mod):
    assert mod.is_drug_clinical_source("CTD")
    assert mod.is_drug_clinical_source("cgi")  # case-insensitive
    assert mod.is_drug_clinical_source("DrugBank")
    assert mod.is_drug_clinical_source("FDA_APPROVED")
    assert not mod.is_drug_clinical_source("OPEN_TARGETS")
    assert not mod.is_drug_clinical_source("DISGENET")
    assert not mod.is_drug_clinical_source("GENOMICS_ENGLAND")


def test_has_drug_provenance(mod):
    assert mod.has_drug_provenance(["OPEN_TARGETS"], ["CTD", "UNIPROT"])
    assert not mod.has_drug_provenance(["OPEN_TARGETS"], ["ORPHANET"])
    assert not mod.has_drug_provenance(None, None)


# --- covariate bands ------------------------------------------------------------


def test_covariate_band(mod):
    assert mod.covariate_band(0.099) == "f3_lt_10pct"
    assert mod.covariate_band(0.10) == "f3_gte_10pct"
    assert mod.covariate_band(0.6) == "f3_gte_10pct"
