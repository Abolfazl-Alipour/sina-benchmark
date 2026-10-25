<!-- track_scope: benchmark -->
# SINA Phase F — Comprehensive Review & SitRep (for Codex / Claude / 5 LLMs / 2 Scientists)

**Date:** 2026-08-10 03:10 UTC · **Commit:** `abb84f8` (Phase F: LODO exports verified + baseline battery wired) · **Branch:** `main` · **Repo:** `github.com/Abolfazl-Alipour/sina-platform` (private) · **Recovery:** `refs/tbh/recovery/before-discard/20260810T031314Z-2589738`
**Claim ceiling (binding, CI-linted, per `configs/report_lint.yaml`):** *disease-associated gene recovery / candidate prioritization / computational triage / first benchmark on OptimusKG* — prohibited terms listed in `report_lint.yaml` §13.11 are not used in this report. All reports carry `<!-- track_scope: benchmark -->`.
**Previous commit:** `b171da6` (CHANGES/0020 APPROVED) · **This commit stats:** 10 files, +1414 −18.

---

## 0. TL;DR for busy reviewers

* **Why 1.43M not 295k:** 295k graph was a bug — `src_id.is_in(sup_genes)` on `annotation→gene` edges (`src=GO/PHE`, `dst=gene`) dropped the entire gene-annotation layer (`1.1M` `PHE/BPO/CCO/MFN/PWY-GEN`). Fix `dst_id` → **1,431,945 edges** (`PHE-GEN 768,844` restored) — the **frozen exam** (CHANGES/0020 §1). Old H200 `$21` 4-pass run (`valid 0.0162→0.0109`, test `0.0159/0.0136` vs old gate `0.0695/0.1447`) is **not comparable** — graded on a different graph (trimmed like-for-like degree `0.0094`).
* **Re-evaluated 2026-08-09:** all CPU audits re-run on 1.43M — **DIS-PHE `KEEP(+0.0124) → DROP(-0.0152)`** (295k verdict wrong), `F@3` median `0.000` still **PASS**, `99.36%` OT-only provenance unchanged, BFS separation weakened so `θ=5%` now rests on mega-key pathology (generic AD `71.1%` sep `-2.5`).
* **Phase F status:** LODO exports **verified DONE** (29 units, `train1 1,420,095 =1,431,945-DIS-PHE` identical, `train2` stripped, `valid 2-126`/`test 24-730`, **all guards PASS**). Per-export baseline battery **wired and green** (degree + PPR on `train1+train2`, same-graph v5.1 §8) — `m2_run_lodo_baselines.py` executed headless, `11/11` lodo unit tests + full `unit` suite green, `ruff` clean, `lint_reports 31/0`.
* **Next:** baseline `parquet` already built (29 rows). **Gate 0 only on owner trigger** (two-step PrimeKG + permutation null + **VRAM re-measure on 1.43M** — old `142/144GB` on 295k invalid, cap `$900`).

---

## 1. Context & governing docs (do not re-litigate)

* **Binding spec:** `SINA_Implementation_Plan_v4_FINAL.md` v4.0 for Phases 0–2, plus `SINA_TaskG_Redesign_v5_1.md` v5.1 (CHANGES/0019 adopted 2026-08-07) as **governing Task-G protocol** — supersedes v4 Task-G where they conflict. `CHANGES/0020_v5_1_freeze.md` (approved 2026-08-08, D1 3-panel + D2 sibling-strip) **freezes** roster, clusters, `θ=5%`, `DIS-PHE DROP`, `DRG-BPO EXCLUDE`, seeds `{1337,1339,5027}`, pass-1 only, batch `300→280→250`, power/MDE, arms, license.
* **Data:** OptimusKG `190,531 / 21,813,816 / 27` (repo-confirmed), Dataverse DOI `10.7910/DVN/IYNGEV`, client `optimuskg==1.0.0`. Trimmed cluster `678` diseases / `191,953` DIS-GEN supervision / `17,448` genes / `1,270` drugs. **BRG `1,431,945`** (corrected, commit `34c1d32`).
* **Execution schedule v5.1 §13:** `0 Gate 0 → 1→2→3 → 4a Freeze → 0b pilot → 4b baselines → 5 Arm A → 6 Arm B → 7 stats/release` — cap `$900` (`core $520` + contingencies `+$165`/`+$100`).

---

## 2. Three-pass review (as requested: review 3 times)

### Pass 1 — Lint / tests / claim ceiling (machine-verified)

* `ruff check src tests scripts` — **initial FAIL 15** (`E501` line-too-long ×14 in `m2_run_lodo_baselines.py` + `F401` unused `polars` in `test_m2_run_lodo_baselines.py`) → **fixed** (broke long docstring/CSV/read/filter/PPR/log lines, removed unused import) → **All checks passed! EXIT 0** (re-checked).
* `pytest tests/unit -q` — **EXIT 0**, full unit suite green (68%→100% warnings only `pykeen slcwa shuffle`/`pin_memory` deprecations, expected). Isolated `pytest tests/unit/test_m2_build_lodo_exports.py + test_m2_run_lodo_baselines.py -v` — **11 passed in 0.80s** (10 lodo-export + 1 lodo-baseline fixture).
* `sina.lint_reports reports/ --config configs/report_lint.yaml` — **31 files scanned, 0 violations EXIT 0** — claim ceiling respected (all new reports carry `benchmark` track_scope).

### Pass 2 — Spec compliance (v5.1 §3–§8, CHANGES/0020)

* **BRG fix verified:** `_select_brg_edges` now tests `dst_id.is_in(sup_genes)` for `ANNOTATION_LABELS ["BPO-GEN","MFN-GEN","CCO-GEN","PWY-GEN","PHE-GEN"]` (annotation→gene orientation comment in code), `brg.filter(label==DIS-GEN).height==0` asserted.
* **LODO exports:** `m2_build_lodo_exports.py --help` **EXIT 0**, `--dry-run` **EXIT 0** lists **29** (`mondo:0010196` `0008310` `0012853` `0014157` `0011913` `0008199` `0957576` `0007669` `0010785` `0012192` `0012522` `0100164` `0100165` + 5 secondary-evidence + 11 low-context). Guards hard-asserted: no DIS-GEN in `train1`, `DIS-PHE` dropped (`train1 1,420,095`), key size == frozen audit, `expected_ancestor_strips` Werner→3 traps, `valid/test` `assert_split_disjoint` on `evidence_group_id`, zero-vector `15` excluded, export ctx `== audit ctx - DIS-PHE`.
* **Baseline battery:** `m2_run_lodo_baselines.py --help` **EXIT 0**, `--dry-run` **EXIT 0**, `--resume` **EXIT 0** — **29 units** `degree MRR` `0.0005–0.0157` + `PPR MRR` `0.0–1.0` written to `data/processed/m2/lodo_baselines/lodo_baselines.parquet` (5.1K). Same-graph rule: train is `train1+train2` per export (not global). RotatE deliberately excluded (Step 4b isolated).
* **Config pinning:** `m2_lodo_exports.yaml` inlines **frozen** clusters (5: `c1_progeria` 6, `c2_ad_familial` 5 incl. `mondo:0014265` supervision-only, `c3_pd` 4, `c4_mody` 2, `c5_ndm` 7, `c6_singletons` 6) + 29 units with `tier` + `explicit_sibling_strips` (4 PD units → Kufor-Rakeb + PSP) + `expected_ancestor_strips` Werner→3 traps + `sibling_coefficient_min 0.3` / `intersection 20`.

### Pass 3 — Reproducibility / change control / safety (human-auditable)

* **Manifests & resume:** both new scripts use `sina.config.load_config` + `sina.manifests.Run` with `stage("baselines")` / `stage("build")` etc., `outputs=[...]` hashed, `peak RSS ≤40GB` asserted at run end (spec §0.3.3). `--resume` skips completed stages via `run.done_stages`.
* **Change control:** `HANDOFF.md` updated in **same commit** as `src/`/`scripts/` changes — CI `handoff-freshness` would otherwise fail. `src/sina/splits.py` fix (`group_assignment` now `uniq = group_ids.unique().sort()` then single `map_elements` — two `unique()` calls could misalign partitions due to nondeterministic order) is included and tested.
* **No dummy data:** synthetic fixtures only under `tests/` (spec §0.2.2), verified.
* **Secrets / data:** `data/raw/*` gitignored except `REGISTRY.json`, checkpoints on `model-checkpoints` release, no secrets in commit (`secret-scan` CI).

---

## 3. What was built (Phase F deltas in this commit)

| File | Action | Purpose |
|---|---|---|
| `scripts/m2_build_lodo_exports.py` (679 lines) | **A** new | Per-unit LODO BioPathNet exports per CHANGES/0020 (29 frozen, D1) — `train1=train1.txt` BRG minus `DIS-PHE`, `train2` supervision minus cluster/ancestor/sibling, `valid/test` 15/85 on `evidence_group_id`, `entity_types/names.txt`, `strip_list.json` + all guards |
| `configs/m2_lodo_exports.yaml` (94 lines) | **A** | Frozen roster/clusters/thresholds + D2 strip lists + `expected_ancestor_strips` + `drop_brg_labels: [DIS-PHE]` |
| `reports/m2b/lodo_exports.md` (168 lines) | **A** | 29-unit table (`train2`/`valid`/`test`/`ctx`/`anc`/`sib`/`zv`) + D2 + Werner trap + strip sweep + DIS-PHE reconciliation + zero-vector + guards |
| `scripts/m2_run_lodo_baselines.py` (216 lines) | **A** | Per-export baselines on `train1+train2` (degree, PPR `0.5` symmetric, filtered MRR/Hits via `sina.metrics`) |
| `configs/m2_lodo_baselines.yaml` (14 lines) | **A** | `lodo_exports_dir` `trackd_nodes` `output_dir` `report_path` `null_resamples 100` |
| `reports/m2b/lodo_baselines.md` (39 lines) | **A** | 29-unit degree/PPR table (see §5) |
| `tests/unit/test_m2_build_lodo_exports.py` (164 lines) | **A** | 10 tests: `normalize_ancestor_id`, `cluster_ancestors` trap/ mate exclusion, `overlap_flagged_pairs` thresholds/orientation, `split_valid_test`, `zero_vector`, etc. |
| `tests/unit/test_m2_run_lodo_baselines.py` (30 lines) | **A** | 1 test: `load_export_unit` empty-dir fixture |
| `src/sina/splits.py` | **M** `38d6595→df95223` | `group_assignment` unique-alignment fix |
| `HANDOFF.md` | **M** `6b9aba6→e68fb0d` | Status board M2B-6 + `NEXT ACTION` → Phase F baselines, re-evaluation note 2026-08-09 |

---

## 4. Verification outputs (exact, not paraphrased)

* `m2_build_lodo_exports.py --help` — `EXIT 0` — usage 29 frozen units, guards listed.
* `m2_build_lodo_exports.py --dry-run` — `EXIT 0` — `would build 29 LODO exports → data/processed/m2/lodo_exports` + 29 lines `mondo:0010196 → mondo__0010196` ... `mondo:0019321 → mondo__0019321`.
* `m2_run_lodo_baselines.py --help` — `EXIT 0`
* `m2_run_lodo_baselines.py --dry-run` — `EXIT 0` — `would run LODO baselines ... → data/processed/m2/lodo_baselines`
* `m2_run_lodo_baselines.py --resume` — `EXIT 0` — 29 `unit_baseline` logs, e.g. `efo:0005422 degree_mrr 0.0157 ppr 0.0014`, `mondo:0010196 0.0062/0.0014`, `mondo:0012853 0.0127/1.0`, `mondo:0014527 0.0100/1.0`, `mondo:0008310 0.0070/0.0021` — `baselines_written 29` → `stage_done report` → `lodo_baselines_done`.
* `pytest test_m2_build_lodo_exports.py + test_m2_run_lodo_baselines.py -v` — `11 passed in 0.80s`
* `pytest tests/unit -q` — full unit suite green (warnings only `pykeen` deprecations)
* `ruff check src tests scripts` — `All checks passed! EXIT 0` (after 15→0 fix)
* `lint_reports 31 files scanned, 0 violations EXIT 0`
* **Outputs on disk:** `data/processed/m2/lodo_baselines/lodo_baselines.parquet` (5.1K) + `reports/m2b/lodo_baselines.md`; `data/processed/m2/lodo_exports/manifest.parquet` (9.3K) + 29 `strip_list.json`.

---

## 5. Results (what the numbers mean, plain English)

* **Exports:** `train1` identical `1,420,095` (`1,431,945` trimmed BRG minus `DIS-PHE 11,850` — freeze §5 DROP). `train2` is `191,953` supervision minus held-out cluster + ancestor-strip (`9-45` diseases) + sibling-strip (`5-185`, long by construction because OptimusKG supervision is redundant) — table `reports/m2b/lodo_exports.md` §Export composition.
* **Baselines (same-graph, per export):** degree MRR `0.0005–0.0157`, PPR MRR `0.0–1.0` (two `1.0` on `0012853 Fontaine`/`0014527 progeroid-HCC` — tiny test, likely artefact with `test 201`/`26` — **flag before Gate 0**, not blocking). Median degree ~`0.0015` — popularity barely above random, as expected for ultra-sparse LODO. Every baseline sees **exactly** the model's `train1+train2` — no leakage, no different-graph grading (v5.1 §8).
* **Zero-vector:** `15` answer genes with no edge in `train1+train2` (isolated, unreachable) — excluded from `valid/test` with gene list per `strip_list.json` + coverage note (v5.1 §3 guard).

---

## 6. Risks, flaws found & fixed, and what PPR 1.0 means

* **Fixed pre-commit:** `src/sina/splits.py:group_assignment` used two `group_ids.unique()` calls — `Series.unique()` order nondeterministic across calls → misaligned `GROUP_COL` vs `partition`. Fix: single `uniq = group_ids.unique().sort()` then `map_elements`. Tests still pass.
* **Fixed pre-commit:** `m2_run_lodo_baselines.py` had `15× E501` line-too-long + `F401` unused `polars` — broke into multi-line `pl.read_csv`, `Run(` block, `sorted(d for ...)` generator, `filt`/`deg_metrics`/`ppr_baseline` blocks, `log.info` block, `lines = [...]` block — now `ruff` clean.
* **Flagged not fixed:** PPR `1.0` on two units (`0012853` 201 test, `0014527` 26 test) suggests the query disease is nearly isolated and PPR mass concentrates on its tiny `train2` neighborhood — **do not grade on that number alone**; report alongside degree and (later) RotatE, and check `brg_context_edges` (`98→2` post-DIS-PHE) — the unit is in the low-context tail by design (secondary).
* **Pre-existing:** `reports/m2b/lodo_exports.md` DIS-PHE reconciliation shows `20/29` units drop `1–96` context edges — tier floors were frozen **against DIS-PHE-inclusive counts** (owner D1) — units below floor post-drop are **disclosed not re-tiered** (roster frozen).
* **VRAM debt:** `142/144GB` at `bs300` on `295k` is **invalid** — must re-measure on `1.43M` at Gate 0 (`280→250` fallback pre-registered, extend by CHANGES if `250` OOMs).

---

## 7. Budget & ops

* **Spent:** H200 `~$21` (295k, 4 passes) + earlier A40 smoke/AD runs (all on thin/buggy graphs) — all pre-freeze, not comparable to frozen exam, but gave throughput (`7.39s/batch` steady, `569 batches/epoch ≈70min`, `valid ~8s`, full test `~89s`) and the `pass-1 only` evidence (triple-replicated `0.0162→0.0157→0.0119→0.0109` with falling loss).
* **Remaining cap:** `$900` (v5.1 §13 recomputed at `$2.80/H200-h`: `core $520` = Gate0 `$100`+0b `$25`+Arm A `$245`+Arm B `$147`; contingencies `+$165` 5-seed extension + `+$100` full-roster Arm B; worst `$785` at 25 units, `$835` at 30).
* **Ops:** conda `sina` (py12.13) + `biopathnet` (py10, torch 2.0.1+cu118 patched `strict-negative`+`streaming-eval`), `Run` manifest `data/raw/manifests/*.json` with `done_stages`, `peak RSS ≤40GB` guard, `hourly orphan-instance cron`, auto-terminate, `GITHUB_TOKEN` in `.env` never committed, checkpoints to `model-checkpoints` release + Drive `/workspace`.

---

## 8. Next steps (owner-gated)

1. **Immediate (CPU, no $):** this baseline `parquet` is already built — publish `null` resamples (`≥100` degree-matched) as a follow-up CPU stage and wire `RotatE` isolated (Step 4b, `d=256`, `GEN`-restricted, same filter/tie) — **do not** run `RotatE` on the GPU pool.
2. **Gate 0 (GPU, owner trigger, cap $900):** two-step PrimeKG repro (authors' `torch 2.0.1/py3.10` then ours `2.13.0/py12.13`, per-metric `±10%`), permutation negative control (degree-preserving shuffle, `2–3` units, must collapse to null), `b051a90` + patches reapply, **preemption-resume drill**, **VRAM re-measure on 1.43M**.
3. **Step 0b pilot:** `3–5 units × {supervision as-is vs mega-key-capped} ×1 seed` (`6–10 GPU-h`) → freeze weighting rule + re-confirm `pass-1`/`lr 1e-3`.
4. **Arms A/B:** `~25 units ×3 seeds ×1 pass` + `~15 units ×3 seeds` subset + `warm random 10%` secondary (`remove_pos: no` guard).

No training on any non-frozen exam (v5.1 §10).

---

## 9. Artifacts for reviewers

* **Diff:** `git show abb84f8` — 10 files, `HANDOFF.md` diff shows `M2B-6` bump + `NEXT ACTION` rewrite (old `Phase F LODO export builder + guards + baselines` → new `Phase F — LODO baselines` with re-evaluation note).
* **Reproduce:** `conda run -n sina python scripts/m2_build_lodo_exports.py --config configs/m2_lodo_exports.yaml --dry-run` / `--resume` (already `manifest.parquet` present); `conda run -n sina python scripts/m2_run_lodo_baselines.py --config configs/m2_lodo_baselines.yaml --resume`; `conda run -n sina pytest tests/unit/test_m2_build_lodo_exports.py tests/unit/test_m2_run_lodo_baselines.py -v`; `conda run -n sina ruff check src tests scripts`; `conda run -n sina python -m sina.lint_reports reports/ --config configs/report_lint.yaml`.
* **Raw:** `data/processed/m2/lodo_exports/` (29 `train1.txt`/`train2.txt`/`valid.txt`/`test.txt` + `strip_list.json`), `data/processed/m2/lodo_baselines/lodo_baselines.parquet`, `data/raw/manifests/m2_build_lodo_exports.json` + `m2_lodo_baselines.json` (hashes, `done_stages`).

*Prepared for: Codex, Claude, 5× LLMs (DeepSeek/GEmini/Grok/Mistral/Llama-class), 2 scientists — all primary sources cited by path, no hidden grader material.*
