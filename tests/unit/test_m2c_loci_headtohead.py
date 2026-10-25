"""Unit tests for scripts/pod/m2c_loci_headtohead.py (locus head-to-head).

Synthetic fixtures only — no science data, no network, no torch/GPU.
Covers the pre-registered rules: pessimistic ties, best-known-rank with
multiple knowns, OT-missing-as-miss, coverage classification, seeded chance.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPT = (Path(__file__).resolve().parents[2] / "scripts" / "pod" /
          "m2c_loci_headtohead.py")


@pytest.fixture()
def mod():
    spec = importlib.util.spec_from_file_location("m2c_loci_headtohead",
                                                  SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_pessimistic_ranks_distinct(mod):
    ranks = mod.pessimistic_ranks({"a": 3.0, "b": 1.0, "c": 2.0},
                                  ["a", "b", "c"])
    assert ranks == {"a": 1, "b": 3, "c": 2}


def test_pessimistic_ranks_exact_tie_shares_worst(mod):
    ranks = mod.pessimistic_ranks({"a": 5.0, "b": 5.0, "c": 1.0},
                                  ["a", "b", "c"])
    assert ranks["a"] == 2
    assert ranks["b"] == 2
    assert ranks["c"] == 3


def test_pessimistic_ranks_deterministic(mod):
    s = {"x": 0.5, "y": 0.5, "z": 0.5}
    assert mod.pessimistic_ranks(s, ["x", "y", "z"]) == \
        mod.pessimistic_ranks(s, ["z", "y", "x"])


def test_best_known_rank_min_and_none(mod):
    assert mod.best_known_rank({"a": 1, "b": 4, "c": 2}, {"b", "c"}) == 2
    assert mod.best_known_rank({"a": 1}, {"zzz"}) is None


def test_hit_at_k_boundaries(mod):
    assert mod.hit_at_k(1, 1) == 1.0
    assert mod.hit_at_k(2, 1) == 0.0
    assert mod.hit_at_k(5, 5) == 1.0
    assert mod.hit_at_k(None, 1) == 0.0


def test_locus_arm_metrics_tied_top_is_miss(mod):
    # Known ties for first at equal score -> pessimistic rank 2 -> top-1 miss.
    m = mod.locus_arm_metrics({"g1": 9.0, "g2": 9.0, "g3": 1.0},
                              ["g1", "g2", "g3"], {"g2"})
    assert m["best_rank"] == 2
    assert m["top1"] == 0.0
    assert m["top5"] == 1.0


def test_ot_arm_metrics_missing_is_miss(mod):
    m = mod.ot_arm_metrics({}, ["g1", "g2"], {"g1"})
    assert m == {"best_rank": None, "top1": 0.0, "top5": 0.0}


def test_ot_arm_metrics_present_delegates(mod):
    m = mod.ot_arm_metrics({"g1": 0.9, "g2": 0.1}, ["g1", "g2"], {"g1"})
    assert m["top1"] == 1.0 and m["best_rank"] == 1


def test_classify_locus(mod):
    assert mod.classify_locus(1, 1) == "ungradeable-few-genes"
    assert mod.classify_locus(5, 0) == "no-known"
    assert mod.classify_locus(5, 2) == "gradeable"


def test_bootstrap_ci_deterministic_and_shaped(mod):
    vals = np.array([1.0, 0.0, 1.0, 0.0])
    lo1, hi1 = mod.bootstrap_ci(vals, 200, 1337)
    lo2, hi2 = mod.bootstrap_ci(vals, 200, 1337)
    assert (lo1, hi1) == (lo2, hi2)
    assert lo1 <= 0.5 <= hi1


def test_bootstrap_ci_constant_input(mod):
    lo, hi = mod.bootstrap_ci(np.ones(10), 100, 7)
    assert lo == hi == 1.0


def test_chance_draw_means_matches_analytic(mod):
    probs = [0.25, 0.5, 1.0, 0.1]
    means = mod.chance_draw_means(probs, 4000, 1337)
    assert means.shape == (4000,)
    assert abs(float(np.mean(means)) - float(np.mean(probs))) < 0.02


def test_verdict_for_diff_boundaries(mod):
    assert mod.verdict_for_diff(0.06, 0.01, 0.11) == "WIN"
    assert mod.verdict_for_diff(0.06, -0.01, 0.11) == "TIE"
    assert mod.verdict_for_diff(0.03, 0.01, 0.05) == "TIE"
    assert mod.verdict_for_diff(-0.06, -0.11, -0.01) == "LOSS"
    assert mod.verdict_for_diff(-0.06, -0.11, 0.01) == "TIE"


def test_grade_known_uses_selected_split(mod):
    genes = ["g_train", "g_test", "g_other"]
    assert mod.grade_known(genes, {"g_train"}) == {"g_train"}
    assert mod.grade_known(genes, {"g_test"}) == {"g_test"}
    assert mod.grade_known(genes, set()) == set()


def test_grade_known_test_split_drops_train_only_locus(mod):
    # Locus whose only known gene is a train positive becomes no-known
    # under test-split grading (train recovery must not leak into the
    # unseen-gene grade).
    genes = ["g_train", "g_other"]
    known = mod.grade_known(genes, {"g_unseen"})
    assert mod.classify_locus(len(genes), len(known)) == "no-known"


def _mini_results(unit: str, top1: float, known_split: str) -> dict:
    return {
        "unit": unit, "disease": unit.replace("__", ":"), "checkpoint": "x",
        "pinned_source": "synthetic", "known_split": known_split,
        "n_loci_total": 3, "n_gradeable": 3, "ungradeable": [],
        "no_known": [], "chance_draws": {"mean": 0.5, "ci": [0.4, 0.6]},
        "strata": {"all": {
            "n": 3,
            "model": {"top1": {"mean": top1, "ci": [top1, top1], "n": 3},
                      "top5": {"mean": 1.0, "ci": [1.0, 1.0], "n": 3},
                      "mean_best_rank": {"mean": 1.0, "ci": [1.0, 1.0],
                                         "n": 3}},
            "ot_primary": {"top1": {"mean": 0.5, "ci": [0.0, 1.0], "n": 3},
                           "top5": {"mean": 1.0, "ci": [1.0, 1.0], "n": 3},
                           "mean_best_rank": {"mean": 2.0, "ci": [1.0, 3.0],
                                              "n": 3}},
            "ot_covered_subset": {"n": 2,
                                  "top1": {"mean": 0.5, "ci": [0.0, 1.0],
                                           "n": 2}},
            "degree": {"top1": {"mean": 0.5, "ci": [0.0, 1.0], "n": 3},
                       "top5": {"mean": 1.0, "ci": [1.0, 1.0], "n": 3}},
            "diff_model_minus_ot_top1": {"mean": top1 - 0.5,
                                         "ci": [top1 - 0.5, top1 - 0.5],
                                         "verdict": "TIE"},
        }},
    }


def test_render_report_uses_given_results_paths(mod, tmp_path):
    # Regression: a test-knowns report must never render the full-panel
    # numbers (the 2026-09-25 artifact bug — .md disagreed with its JSON).
    import json

    full = tmp_path / "loci_headtohead_mondo__0005180.json"
    tk = tmp_path / "loci_headtohead_mondo__0005180_testknown.json"
    full.write_text(json.dumps(_mini_results("mondo__0005180", 0.9, "train")))
    tk.write_text(json.dumps(_mini_results("mondo__0005180", 0.3, "test")))
    out = tmp_path / "report.md"
    mod.render_report(out, [tk])
    text = out.read_text()
    assert "0.3000" in text and "0.9000" not in text
    assert "known_split: **test**" in text


def test_render_report_default_glob_unchanged(mod, tmp_path, monkeypatch):
    pd_json = (mod.REPO_ROOT / "reports" / "m2c" /
               "loci_headtohead_mondo__0005180.json")
    if pd_json.exists():  # default path renders whatever full-panel files exist
        out = tmp_path / "report.md"
        mod.render_report(out)
        assert "known_split" in out.read_text()
