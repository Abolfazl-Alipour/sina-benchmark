"""Report language lint — plan §13.11 (claim ceiling) and §13.13 (narrative separation).

Scans Markdown reports for forbidden claim language and for Track-R/D-vs-Track-C
narrative blending. Reports opt into a scope with an HTML comment on any line:
``<!-- track_scope: benchmark -->`` or ``<!-- track_scope: commercial -->``.

Usage:  python -m sina.lint_reports [reports_dir] [--config configs/report_lint.yaml]
Exit code 1 when any violation is found; a missing reports dir passes with a note.
"""

import argparse
import re
import sys
from pathlib import Path

import yaml

SCOPE_RE = re.compile(r"<!--\s*track_scope:\s*(benchmark|commercial)\s*-->")


def _term_re(term: str) -> re.Pattern[str]:
    return re.compile(r"\b" + re.escape(term) + r"\b", re.IGNORECASE)


def lint_file(path: Path, rules: dict) -> list[str]:
    """Return violation strings (file:line — message) for one Markdown file."""
    violations: list[str] = []
    lines = path.read_text(encoding="utf-8").splitlines()
    scope = next((m.group(1) for line in lines if (m := SCOPE_RE.search(line))), None)

    term_sets = [("forbidden", rules["forbidden_terms"])]
    if scope == "benchmark":
        term_sets.append(("benchmark-scope", rules["benchmark_forbidden_terms"]))
    elif scope == "commercial":
        term_sets.append(("commercial-scope", rules["commercial_forbidden_terms"]))

    for label, terms in term_sets:
        for term in terms:
            rx = _term_re(term)
            for lineno, line in enumerate(lines, start=1):
                if rx.search(line):
                    violations.append(
                        f"{path}:{lineno} — {label} term {term!r} in: {line.strip()[:80]!r}"
                    )
    return violations


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports_dir", nargs="?", default="reports")
    parser.add_argument("--config", default="configs/report_lint.yaml")
    args = parser.parse_args(argv)

    reports_dir = Path(args.reports_dir)
    if not reports_dir.is_dir():
        print(f"lint_reports: {reports_dir} does not exist — nothing to lint, PASS")
        return 0

    with Path(args.config).open("r", encoding="utf-8") as fh:
        rules = yaml.safe_load(fh)

    violations: list[str] = []
    md_files = sorted(reports_dir.rglob("*.md"))
    for md in md_files:
        violations.extend(lint_file(md, rules))

    for v in violations:
        print(f"VIOLATION {v}")
    print(f"lint_reports: {len(md_files)} files scanned, {len(violations)} violations")
    return 1 if violations else 0


if __name__ == "__main__":
    sys.exit(main())
