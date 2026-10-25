"""Unit tests for the Track-C license gate (scripts/check_trackc_licenses.py).

Synthetic catalog fixtures only.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from check_trackc_licenses import validate_catalog


def _entry(**kw):
    base = {
        "name": "src",
        "license_verbatim": "CC0 1.0",
        "commercial_use": True,
        "status": "candidate",
        "include": False,
    }
    base.update(kw)
    return base


def test_valid_candidate_catalog_passes():
    assert validate_catalog({"sources": [_entry()]}) == []


def test_included_requires_commercial_use():
    v = validate_catalog({"sources": [_entry(include=True, commercial_use=False)]})
    assert any("commercial_use" in x for x in v)


def test_included_requires_owner_approval():
    v = validate_catalog({"sources": [_entry(include=True, status="candidate")]})
    assert any("owner-approved" in x for x in v)


def test_missing_verbatim_license_fails():
    v = validate_catalog({"sources": [_entry(license_verbatim="")]})
    assert any("license_verbatim" in x for x in v)


def test_excluded_requires_reason():
    v = validate_catalog({"sources": [_entry(status="excluded")]})
    assert any("exclusion_reason" in x for x in v)
    ok = validate_catalog({"sources": [_entry(status="excluded", exclusion_reason="paid license")]})
    assert ok == []


def test_empty_catalog_fails():
    assert validate_catalog({}) != []
