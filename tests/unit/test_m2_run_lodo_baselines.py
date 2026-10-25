"""Unit tests for scripts/m2_run_lodo_baselines.py (v5.1 Phase F).

Synthetic fixtures only — no science data in tests (spec §0.2.2).
"""

import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "m2_run_lodo_baselines.py"


@pytest.fixture()
def mod():
    spec = importlib.util.spec_from_file_location("m2_run_lodo_baselines", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_load_export_unit_empty(mod, tmp_path):
    ud = tmp_path / "mondo__0000001"
    ud.mkdir()
    (ud / "strip_list.json").write_text('{"unit": "mondo:0000001", "tier": "primary"}')
    # no train files -> empty frames
    data = mod.load_export_unit(ud)
    assert data["train1"].height == 0
    assert data["test"].height == 0
