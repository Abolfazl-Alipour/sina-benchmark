"""Unit tests for sina.lint_reports — claim ceiling + narrative separation lints."""

from pathlib import Path

from sina.lint_reports import main

CONFIG = Path(__file__).resolve().parents[2] / "configs" / "report_lint.yaml"


def _write(dirpath: Path, name: str, text: str) -> None:
    (dirpath / name).write_text(text, encoding="utf-8")


def test_clean_reports_pass(tmp_path):
    _write(tmp_path, "ok.md", "# Gene recovery report\n\nDisease-associated gene recovery.\n")
    assert main([str(tmp_path), "--config", str(CONFIG)]) == 0


def test_forbidden_term_fails(tmp_path):
    _write(tmp_path, "bad.md", "# Report\n\nThis is a validated target for AD.\n")
    assert main([str(tmp_path), "--config", str(CONFIG)]) == 1


def test_word_boundary_avoids_false_positives(tmp_path):
    _write(tmp_path, "fine.md", "We safeguard reproducibility; unclear findings are disclosed.\n")
    assert main([str(tmp_path), "--config", str(CONFIG)]) == 0


def test_benchmark_scope_rejects_commercialization(tmp_path):
    _write(tmp_path, "b.md", "<!-- track_scope: benchmark -->\n\nSBIR relevance is discussed.\n")
    assert main([str(tmp_path), "--config", str(CONFIG)]) == 1


def test_commercial_scope_rejects_benchmark_language(tmp_path):
    _write(tmp_path, "c.md", "<!-- track_scope: commercial -->\n\nZero-shot MRR improved.\n")
    assert main([str(tmp_path), "--config", str(CONFIG)]) == 1


def test_missing_dir_passes(tmp_path):
    assert main([str(tmp_path / "nope"), "--config", str(CONFIG)]) == 0
