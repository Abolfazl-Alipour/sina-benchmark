# M5B — beat Open Targets L2G on covered PR AUC

Independent, leakage-controlled attempt to beat OT L2G where v4 lost
(covered PR AUC: OT 0.4767 vs v4 0.22-0.25).

Plan: `.agents/plans/2026-09-26-beat-ot-coloc-hpo.md` (deciding metric: PR AUC).

Layout:
- `scripts/` — pipeline code (imports metric/report helpers from
  `scripts/m4b/`, copies nothing).
- `reports/` — verdicts and manifests (committed).
- `cache/` — downloaded stats and per-triple coloc values (gitignored,
  resume-safe; any session rebuilds from manifests).

Inputs (read-only, never copied): `data/processed/m4b/` panels, features,
OT scores. Tests live in `tests/unit/test_m4b_*.py` per repo convention.

Rules: CPU-only, pre-registered single runs, no OT-derived scores as
features, HPO inside folds only, every number reported as-is.
