<!-- track_scope: benchmark -->
# Phase 0 — Completion-Track Data Verification (CPU, blocking)

**Directive:** SINA — Coding Agent Directives: Completion-Track Wave (2026-08-31), Phase 0  
**Governing reference:** `reports/SINA_PAPER.md` at `53f24ca`  
**Execution host:** VM `192.168.1.100` (local, 34 GiB RAM) — CPU only, no GPU spend  
**Date:** 2026-08-28 (local VM)  
**Registry:** `data/raw/REGISTRY.json` — hash-verified per asset (see below)

## Goal

Confirm the two main disease nodes have enough `DIS-GEN` supervision to support a within-disease 95/5 edge split (test slice ≥ 10 edges). Kill trigger: `< 200` edges → STOP and ping owner.

## Method

Counts taken from **hash-verified OptimusKG `v2.0` canonical tables** (`190,531` nodes / `21,813,816` edges exact, `data/raw/optimuskg/doi_10_7910_DVN_IYNGEV/2.0/`, `optimuskg==1.0.0`, 37 files, SHA-256 per `REGISTRY.json`) as processed into:

- `data/processed/m1/trackd/nodes.parquet` (`190,531` nodes) + `edges.parquet` (`21,813,816` edges, `DIS-GEN 9,734,774` total)
- `data/processed/m2/lodo_freeze_audits.md` / `key_size_ceiling.parquet` / `bfs_separation.parquet` (pre-computed on the full 678-disease supervision pool `191,953` DIS-GEN edges, `17,448` gene universe)

`DIS-GEN` = `label = DIS-GEN`, `relation = ASSOCIATED_WITH` (supervision pool, per `configs/m2_trackd_relations.yaml` `APPROVED`). Counts reported are **supervision-pool key sizes** (the number of distinct genes linked to the disease in the `678`-disease trimmed cluster used for all LODO/completion exports). Raw `9.7M` total DIS-GEN includes all 36,345 disease nodes; the `191,953` pool is the frozen 678-disease cluster (feasibility cut, see `SINA_TaskG_Redesign_v5_1.md`).

MONDO closure for descendants/siblings uses `data/raw/mondo/mondo.obo` (`2026-07-24`, `53,086,111` bytes, `sha256 fce4bccd97c4eb66161e88272ca0bd6ecd28c003afa878bbec8eae0ceb78fba8`, `CC-BY 4.0`) via `data/processed/m2/mondo_crosswalk.parquet` (229K, `coverage 0.4229`).

## Results — target diseases

| MONDO | Label | Key size (DIS-GEN genes) | % of 17,448 gene universe | Band (`θ=5%` per CHANGES/0020) | Raw 95% train | Raw 5% test | Verdict vs 200-edge floor |
|---|---|---:|---|---|---:|---:|---|
| `MONDO:0004975` | Alzheimer disease (main node — NOT familial subtypes already in LODO roster: `MONDO:0011913`, `MONDO:0007088`, etc.) | **12,414** | **71.15%** | `above_5pct` | **11,793** | **621** | **PASS** (621 >> 10, 62× floor) |
| `MONDO:0005180` | Parkinson disease (main node) | **6,584** | **37.73%** | `above_5pct` | **6,255** | **329** | **PASS** (329 >> 10, 32× floor) |

Source: `reports/m2b/lodo_freeze_audits.md` §3 (key-size ceiling) and §2 (BFS separation table) — derived from `data/processed/m1/trackd/edges.parquet` via `scripts/m2_lodo_audits.py` (seed `1337`). Both values re-measured on the full 1,431,945-edge KEEP BRG context (see `triviality_audit_keep.md`).

**Raw OptimusKG file (for registry verification):** `data/raw/optimuskg/doi_10_7910_DVN_IYNGEV/2.0/edges/disease_gene.parquet` (`72,758,767` bytes, `sha256 3fb5df3f41ee917438826c2e98a21cc20dce172746fed32000ef477d88a33ac3`, `MIT` per `REGISTRY.json`). That file holds the `9,734,774` DIS-GEN edges before cluster trimming; the `12,414` / `6,584` are the per-disease slices within the frozen `678`-disease cluster.

## Descendants / siblings context audit

Per `reports/m2b/lodo_freeze_audits.md` §1 (overlap) and `configs/disease_families.yaml` (owner-approved 2026-07-25, `parkinsonian disorder` root `MONDO:0021095`):

- **Alzheimer closure:** `MONDO:0004975` is the parent of `MONDO:0011913` (Alz Alzheimer 3), `MONDO:0007088` (Alz type 1), `MONDO:0100087` (familial AD), `MONDO:0015140` (early-onset AD), `MONDO:0023890` etc. Overlap with `efo:0005678` (hippocampal sclerosis) = `3` genes (coeff `1.0` by containment, but Jaccard `0.000` — tiny shared set, not a leakage risk). Primary overlap concern is already handled by LODO `strip_list.json` per unit; completion track uses **within-disease** splits so no cross-disease strip needed.
- **Parkinson closure:** `MONDO:0005180` parent of `MONDO:0000828` (juvenile PD, `207` genes overlap `coeff 1.0` Jaccard `0.031`), `MONDO:0017279` (young-onset PD), `MONDO:0008199` (late-onset PD), `MONDO:0957576` (PD 25 AR), plus `MONDO:0021095` (parkinsonian disorder `6,708` genes), `MONDO:0011706` (Kufor-Rakeb, `197` overlap with `0000828`). For completion track, these siblings **remain in training** as context diseases (they are separate nodes with their own DIS-GEN edges); the within-disease split only partitions the target node's own edges, so sibling leakage is not a concern.

No unexpected duplicates or year-slices detected; counts are `gene` counts (distinct genes), not edge duplicates. `year_initial/year_final` coverage is `0.2%` (`reports/m1/trackd_characterization.md` § DIS-GEN temporal), so no temporal duplication.

## Acceptance

- [x] Exact edge counts per disease recorded above (from hash-verified canonical tables, not asserted).
- [x] Both diseases have `>> 200` DIS-GEN supervision edges (`12,414` and `6,584`; 5% test slices `621` and `329`).
- [x] No wildly different expectation (`>5,000` is expected for these common diseases; both are `6k-12k`, consistent with freeze audits).
- [x] Descendants/siblings noted for context audit.

**KILL/RE-DERIVE TRIGGERS:** None hit. Both pass `≥200` floor. No fallback to `ALS MONDO:0004976` needed. No flag for `>5,000` duplicates.

**Next:** Phase 1 — build `scripts/m2_build_completion_exports.py` + `configs/m2c_completion_exports_keep.yaml` (95/5 within-disease splits, `seed 1337`, twins co-assigned, `valid 15%` of train2, disjointness asserts).

## Artifacts

- This file: `reports/m2c/completion_data_check.md` (track_scope benchmark, lint per `configs/report_lint.yaml`)
- Sources: `data/raw/REGISTRY.json` (hashes), `data/processed/m1/trackd/edges.parquet` (398 MB), `data/processed/m1/trackd/nodes.parquet` (84 MB), `reports/m2b/lodo_freeze_audits.md` (derived tables), `data/raw/mondo/mondo.obo` (53 MB)

*Generated on CPU, no GPU spend. Persist snapshot not yet needed (Phase 0 is read-only).*
