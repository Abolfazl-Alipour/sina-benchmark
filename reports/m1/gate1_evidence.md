<!-- track_scope: benchmark -->
# Gate 1 evidence package — M1 three-track evidence graph (spec §5.5)

**Status: Gate 1 ready — owner review**

Gate 1 (spec §5.5): reviewer inspects reciprocal/provenance audits, license
registers, and the Track-C catalog. This file indexes the evidence. All counts
below are measured from the canonical tables by CI-green code; nothing here is
estimated.

## 1. Acceptance criteria (§5.5) → evidence

| Criterion | Result | Evidence |
|---|---|---|
| Characterization reconciles with published counts | Track-R: nodes 129,312 vs 129,375 (0.049% dev, within ±1%); edges 8,100,498 exact. Track-D: **EXACT** 190,531 nodes / 21,813,816 edges | `reports/m1/trackr_characterization.md`, `reports/m1/trackd_characterization.md` |
| FK integrity | 0 dangling endpoints on both tracks (ingest assert + integration tests + DuckDB battery) | `reports/m1/query_battery.md` |
| DuckDB batteries | Identical battery per track, both green | `reports/m1/query_battery.md` |
| Property tables populated | `data/processed/m1/gene_properties.parquet` — 61,306 genes: tractability, safety_liabilities, target_class, constraint_scores, homologues, chemical_probes | `reports/m1/trackd_characterization.md` |
| 1% subsample integration tests | Group-level deterministic subsamples; twin co-assignment guards green on both tracks | `tests/integration/` (49 tests total green) |
| Reciprocal/provenance audits | Track-R: all 30 relations stored reciprocally (fraction 1.0000) → order-insensitive grouping; 4,050,064 groups. Track-D: edges stored once; grouping per §5.2.3 with source-set in key; verified no materialized twins/duplicates | characterization reports + `src/sina/evidence_groups.py` + `tests/unit/test_evidence_groups*.py` |
| License registers | V1 (verbatim licenses) + `data/raw/REGISTRY.json` (38 assets, SHA-256-verified) | `docs/preflight/V1_licenses.md`, REGISTRY |
| Track-C catalog | 9 sources, verbatim licenses, DrugBank/DisGeNET excluded with documented gaps; CI `trackc-license-gate` | `configs/m1_trackc_sources.yaml` |

## 2. Downloads (all SHA-256-verified in `data/raw/REGISTRY.json`)

- **Track-R**: PrimeKG `kg.csv` (981,751,236 B), Dataverse doi:10.7910/DVN/IXA7BM v2.1.
- **Track-D**: OptimusKG 37 stratified per-type Parquet files, Dataverse
  doi:10.7910/DVN/IYNGEV v2.0, via the pinned official `optimuskg` 1.0.0 client.

## 3. Deviations on record (`CHANGES/`)

- 0001–0005: environment/protocol deviations (threads, Docker, deferred deps,
  GitHub workflow, editable install).
- 0006: PrimeKG license verbatim is CC0 1.0 (spec said CC BY-NC-SA 4.0);
  conservative research-only working assumption retained.
- 0007: NIH fall receipt date is Sept 5, 2026 (spec said Sept 8).
- 0008: OptimusKG v2.0 gaps — disease ids not MONDO-preferred (11,724/36,345
  MONDO-native; crosswalk deferred to Phase 2) and pocket/druggable flags
  absent from gene nodes.

## 4. Measured finding with plan impact

- **DIS-GEN temporal coverage: 0.25%** of 9,734,774 DIS-GEN edges carry both
  `year_initial` and `year_final` (threshold 50%, Appendix D assumption 5).
  → The temporal split is **downgraded to secondary** for Phase 2. This is the
  spec-provided decision path, not a deviation; it removes the
  prospective-flavored evaluation as a primary axis on Track-D v2.0.

## 5. Owner sign-off items before/within Gate 1 review

1. Sign V-notes V1–V6 (`docs/preflight/`); author V7 (owner-only eligibility note).
2. Approve the Track-C source catalog (incl. Open Targets per-component terms
   and the DrugBank/DisGeNET coverage gaps).
3. Approve CHANGES/0008 (pending status).

## 6. Test & lint state at this commit

- `pytest tests/unit tests/integration`: 49 green (integration tests run
  against the real canonical tables; they skip in CI where data is absent).
- `ruff check src tests scripts`: clean.
- `python -m sina.lint_reports reports/`: 0 violations (3 report files).
