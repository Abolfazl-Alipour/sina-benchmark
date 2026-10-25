"""Unit tests for scripts/pod/anchored_topology_null.py (CHANGES/0025 control).

Synthetic fixtures only — no science data, no torch/GPU.
"""

import importlib.util
import math
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "pod" / "anchored_topology_null.py"
EXTRACT_PATHS = Path(__file__).resolve().parents[2] / "scripts" / "pod" / "extract_paths.py"


@pytest.fixture()
def mod():
    spec = importlib.util.spec_from_file_location("anchored_topology_null", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def hub_adjust():
    spec = importlib.util.spec_from_file_location("extract_paths", EXTRACT_PATHS)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.hub_adjust


# ---------------------------------------------------------------------------
# topology null ranking
# ---------------------------------------------------------------------------


def test_topology_null_rank_hub_penalty(mod, hub_adjust):
    # same support (5); the hub (degree 10_000) loses to the specific node (8)
    candidates = {100: {1, 2, 3, 4, 5}, 101: {1, 2, 3, 4, 5}}
    rows = mod.topology_null_rank(candidates, {100: 8, 101: 10_000}, 10, hub_adjust)
    assert [r["node"] for r in rows] == [100, 101]
    assert rows[0]["score"] == pytest.approx(5.0 / math.log(8))
    assert rows[1]["score"] == pytest.approx(5.0 / math.log(10_000))
    assert rows[0]["support"] == 5


def test_topology_null_rank_support_dominates_and_selects(mod, hub_adjust):
    candidates = {100: {1, 2, 3}, 101: {1, 2, 3, 4, 5, 6}}
    rows = mod.topology_null_rank(candidates, {100: 8, 101: 8}, 1, hub_adjust)
    assert [r["node"] for r in rows] == [101]  # higher support wins; select=1
    # degree < 2 -> hub_adjust guard scores 0
    rows = mod.topology_null_rank({102: {1, 2, 3}}, {102: 1}, 10, hub_adjust)
    assert rows[0]["score"] == 0.0


def test_topology_null_rank_tiebreak(mod, hub_adjust):
    # identical score and degree -> degree tie-break, then node id
    candidates = {100: {1, 2, 3}, 101: {4, 5, 6}}
    rows = mod.topology_null_rank(candidates, {100: 4, 101: 4}, 2, hub_adjust)
    assert [r["node"] for r in rows] == [100, 101]


# ---------------------------------------------------------------------------
# report table parsing
# ---------------------------------------------------------------------------

REPORT_SNIPPET = """## EXPLORATORY — ontology-anchored pathway scores (post-hoc)

Disclosure text.

| rank | pathway | score | degree | support | top genes |
|---|---|---|---|---|---|
| 1 | go:0010628 (positive regulation of gene expression) | 0.0262 | 446 | 26 | ensembl:X |
| 2 | reactome:R-HSA-6785807 (some pathway) | 0.0250 | 1679 | 23 | ensembl:Y |

Trailing prose.
"""


def test_parse_anchored_top10(mod):
    assert mod.parse_anchored_top10(REPORT_SNIPPET) == [
        "go:0010628",
        "reactome:R-HSA-6785807",
    ]


def test_parse_anchored_top10_placeholder(mod):
    text = (
        "## EXPLORATORY — x\n\n| rank | pathway | score | degree | support | top genes |\n"
        "|---|---|---|---|---|---|\n| — | no anchored pathway passed min_support | | | |\n"
    )
    assert mod.parse_anchored_top10(text) == []


def test_parse_anchored_top10_missing_section(mod):
    with pytest.raises(ValueError, match="EXPLORATORY"):
        mod.parse_anchored_top10("# some other report\n\nno table here\n")


# ---------------------------------------------------------------------------
# jaccard
# ---------------------------------------------------------------------------


def test_jaccard(mod):
    assert mod.jaccard([1, 2, 3], [2, 3, 4]) == pytest.approx(0.5)
    assert mod.jaccard([], []) == 0.0
    assert mod.jaccard([1], [2]) == 0.0
