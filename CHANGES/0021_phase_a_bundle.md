# 0021 — CHANGES/0021 bundle (A1-A9, blocking before GPU, 2026-08-10 review round)

- **change:** Phase A CPU bundle (A1-A9) as defined in `SINA_Coding_Agent_Directives_2026-08-10.md` (integrates 7-review round on Phase-F package `abb84f8`+`17bc7a4`, corrected 1.43M BRG). This entry records the rebuild, new guards, audits, floor amendment, statistical pin, theta derivation, near-twin check, pre-Arm wiring, and DisGeNET flag. No GPU spend until this entry is owner-signed.
- **rationale:** Seven-review round (R-A…R-G, 14 points) found: (1) 29 exports on disk predate `splits.py` `unique()` fix — nondeterministic group→partition breaks frozen-export reproducibility (R-D Flaw 3, sharpest catch); (2) missing tuple-level and BRG-composition guards; (3) triviality/degenerate handling absent (PPR 1.0 on 2/29); (4) floor vs D1 inclusive ruling not on record; (5) 13-primary confirmatory population not pinned in code; (6) θ derivation not attempted; (7) 42.3% MONDO crosswalk gaps not spot-checked; (8) DisGeNET secondhand CC BY-NC vs proprietary Subscription Agreement (no redistribution, company use in any capacity requires commercial license — legal teeth, R-C).
- **requester:** Directives 2026-08-10 (owner pasted as amended v5.1 + directives).
- **impact:** Rebuilds 29 LODO exports on patched `splits.py` (hash-verified, guards PASS, baselines re-run on rebuilt exports). New CI guards: `set(train2)∩set(valid)==∅` + `set(train2)∩set(test)==∅` (tuple level, per unit) and BRG-composition `PHE-GEN>500k`/`BPO-GEN>100k` in build script (must fail build, not just test). Triviality audit per unit (shortest-path `train1` disease→answer) feeds degenerate rule: baseline `MRR≥0.99` (query isolation, `98→2` edges) OR `>10% ≤2 hops` → excluded from primary gate, degenerate row with cause; `valid_baseline` flag in parquet; tiny `~10`-node smoke test `MRR<1`. Context-attrition table (per unit `audit→export` `lost`, `below_floor` flag: **all 29 <10 post-DIS-PHE**, exploratory-tier; with/without sensitivity row). Stats: 13 primary confirmatory — exact Wilcoxon (permutation, `2^n` at n=13) + hierarchical BCa `≥1000` (strata→units) + cluster-label permutation; `MDE dz≈0.8` at 13. Theta: null `≥100` shuffles, candidate `35.89%` at 95th-pct vs frozen `5%` (keep `5%`, publish `4%` sensitivity row `+` continuous `pct` bins `<2`/`2–4`/`4–5`/`>5` for Werner `4.91%`). Near-twin: disclose-band `38` pairs among graded units (`0.15–0.3`+`≥20`); `efo:0005422` EFO gap; checklist for manual verify. Pre-Arm wiring: `m2_build_warm_exports.py` (`10%` warm, `remove_pos: no` guard), causal-gene-rank (`WRN`/`LMNA` etc., monogenic only), RotatE/ComplEx `d=256` `GEN`-restricted, `≥100` nulls, `m1_reproduce_biopathnet.py` two-step Gate 0, OT-only filter (DisGeNET no-redistribution). DisGeNET flag: owner/counsel decision whether SBIR research program needs commercial license (Track-C remains only clean commercial path). Hard gate: no GPU until **items 1–4** complete and this entry signed.
- **date:** 2026-08-10
- **approved:** owner — 2026-08-10 (Light Gate 0 plan approved 2026-08-10; hard gate lifted for Light Gate 0 only per CHANGES/0022)

---

## A1 — Rebuild

- Patched `src/sina/splits.py:group_assignment` — single `uniq = group_ids.unique().sort()` then `map_elements` (two `unique()` calls misaligned group→partition due to nondeterministic order).
- Cleared `data/processed/m2/lodo_exports` and rebuilt 29 exports via `m2_build_lodo_exports.py` (run `20260810T002750Z`, `supervision 191,953` `brg 1,431,945` `units 29`, exit 0). Manifest `9.3K` hash-verified; `lodo_baselines.parquet` re-run on rebuilt exports (29 rows) via `m2_run_lodo_baselines.py --resume` (exit 0).

## A2 — Guards

- `scripts/m2_build_lodo_exports.py:build_exports` — BRG-composition `PHE-GEN >500k` `BPO-GEN >100k` (orientation-bug class fails build) + per-unit tuple `train2∩valid` / `train2∩test` disjoint (group-level necessary, not sufficient).

## A3 — Triviality

- `scripts/m2_triviality_audit.py` + `configs/m2_triviality.yaml` + `data/processed/m2/lodo_audits/triviality.parquet` (5.0K) + `reports/m2b/triviality_audit.md` (2.9K): per-unit `dist 1/2/3/4/5/unreachable` + `pct_le2`. Degenerate `2/29` (`mondo:0012853` 201 unreachable + PPR 1.0, `mondo:0014527` 26 unreachable + PPR 1.0) — `pct_le2` all `≤0.015` so degenerate via baseline `≥0.99` only, not triviality `>10%`. Baselines parquet updated with `valid_baseline`/`degenerate` flags.

## A4 — Floor

- `data/processed/m2/lodo_audits/context_attrition.parquet` + `reports/m2b/context_attrition.md` (39 lines): per-unit `audit→export` `lost 0–96`, `below_floor 29/29` (`<10` post-DIS-PHE: max export `9` efo:0005422/mondo:0018911, min `0` orphanet:75496). D1 inclusive ruling on record: **all below-floor are exploratory-tier, with/without sensitivity row**; the freeze already incorporated the `DIS-PHE` drop.

## A5 — Stats

- `src/sina/stats.py` — `exact_wilcoxon` (permutation `2^n`, n=13 → 8192), `bca_bootstrap`/`hierarchical_bca` (strata→units, `≥1000`), `cluster_permutation_test` (permute cluster labels). `compare_models` now uses `exact_wilcoxon`.

## A6 — Theta

- `scripts/m2_theta_derivation.py` + `configs/m2_theta.yaml` + `theta_null.parquet` (515B) + `reports/m2b/theta_derivation.md`: null `35.89%` at 95th-pct vs frozen `5.0%` (keep `5%`, publish `4%` row — Werner drops — + continuous `pct` bins).

## A7 — Near-twin

- `reports/m2b/near_twin_checklist.md` (50 lines): `38` disclose-band pairs among graded units (`0.15–0.3`+`≥20`), plus `efo:0005422` EFO gap; owner/scientist manual verify.

## A8 — Pre-Arm wiring

- `scripts/m2_build_warm_exports.py` + `configs/m2_warm_exports.yaml` (`10%` warm, `remove_pos: no`)
- `scripts/m2_causal_rank.py` + `configs/m2_causal_rank.yaml` (WRN/LMNA, monogenic only)
- `scripts/m2_filter_ot_only.py` (OT-only, DisGeNET no-redistribution)
- `scripts/m1_reproduce_biopathnet.py` + `configs/m1_biopathnet_repro.yaml` (Gate 0 two-step)
- All `--dry-run` `EXIT 0`; `ruff` remaining `SIM117`+`E501` in repro script is cosmetic (nested `with` + docstring length, `EXIT 0`).

## Owner decisions in this entry

- **Company research-use licensing:** whether SBIR-funded research program needs commercial DisGeNET license — **open, counsel-level** (flagged §2.3).
- **Theta:** adopt `35.89%` vs keep `5%` — **keep `5%`** (frozen), publish `4%` row (recorded).
