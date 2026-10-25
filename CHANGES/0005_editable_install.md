# 0005 — `sina` package installed editable

- **change:** The spec does not prescribe how the `sina` package is importable.
  Adopted: `pip install -e . --no-deps` inside the conda env so `python -m sina.*`
  works for scripts and CI; pytest additionally gets `pythonpath=["src"]` from
  `pyproject.toml` so tests run without relying on the install.
- **rationale:** Spec §0.3 requires `python -m`-style script execution and
  reproducibility; an editable install plus pytest `pythonpath` covers both with no
  duplication.
- **requester:** agent (2026-07-24).
- **impact:** Dev ergonomics only; the Docker/CI build performs a non-editable
  install. No spec content changes.
- **date:** 2026-07-24
- **approved:** owner — 2026-07-24
