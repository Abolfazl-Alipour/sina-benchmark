#!/usr/bin/env python
"""check_trackc_licenses — CI license gate for the Track-C source catalog (§5.3).

Fails (exit 1) if any source marked `include: true` in
configs/m1_trackc_sources.yaml lacks `commercial_use: true` or a non-empty
`license_verbatim`, or if any excluded source lacks a documented reason.
This is the permissive-license gate that must pass BEFORE a source enters the
Track-C pipeline.

Usage: python scripts/check_trackc_licenses.py [--config configs/m1_trackc_sources.yaml]
"""

import argparse
import sys
from pathlib import Path

from sina.config import load_yaml
from sina.logging import configure_logging, get_logger

log = get_logger(__name__)


def validate_catalog(catalog: dict) -> list[str]:
    """Return a list of violations (empty = pass).

    Every entry: {name, license_verbatim, commercial_use, status, include}.
    - include: true requires commercial_use: true AND status: owner-approved
      (the permissive-license gate BEFORE pipeline inclusion, spec §5.3);
    - status: excluded requires a documented exclusion_reason;
    - status: candidate entries are catalogued but not pipeline-included.
    """
    violations = []
    sources = catalog.get("sources")
    if not isinstance(sources, list) or not sources:
        return ["catalog has no 'sources' list"]
    for entry in sources:
        name = entry.get("name", "<unnamed>")
        include = bool(entry.get("include", False))
        status = entry.get("status", "candidate")
        if not entry.get("license_verbatim", "").strip():
            violations.append(f"{name}: missing license_verbatim")
        if not isinstance(entry.get("commercial_use"), bool):
            violations.append(f"{name}: commercial_use must be an explicit boolean")
        if status not in {"candidate", "excluded", "owner-approved"}:
            violations.append(f"{name}: unknown status {status!r}")
        if status == "excluded" and not entry.get("exclusion_reason", "").strip():
            violations.append(f"{name}: excluded without exclusion_reason")
        if include and not entry.get("commercial_use"):
            violations.append(f"{name}: included but commercial_use is not true")
        if include and status != "owner-approved":
            violations.append(f"{name}: included without owner-approved status")
    return violations


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("configs/m1_trackc_sources.yaml"))
    args = parser.parse_args()

    configure_logging()
    catalog = load_yaml(args.config)
    violations = validate_catalog(catalog)
    for v in violations:
        log.error("license_gate_violation", violation=v)
    if violations:
        log.error("license_gate_failed", n_violations=len(violations))
        return 1
    log.info("license_gate_passed", n_sources=len(catalog["sources"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
