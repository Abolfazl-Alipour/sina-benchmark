<!-- track_scope: benchmark -->
# M2A-5 Gate 2 evidence package — M2A harness (spec §6.6)

**Status:** DRAFT — reviewer audit of split integrity, temporal-split
construction, nulls, external set, and statistics.

## Acceptance criteria (spec §6.6)

| Criterion | Status | Evidence |
|---|---|---|
| Custom splitter passes all leakage guards in CI | PASS | `src/sina/splits.py`; `assert_split_disjoint` enforced per split; `tests/unit/test_splits.py` + `tests/integration/test_m2_splits.py` |
| Pilot timings | PARTIAL | Degree/PPR/nulls piloted (`reports/m2a/baselines_report.md`); AnyBURL piloted (`reports/m2a/anyburl_report.md`); full-scale KGE pilots not run (local runtime too large) |
| Full baseline+null table | PARTIAL | Degree, PPR, degree-matched null, relation-shuffled null computed (`reports/m2a/metrics_report.md`); KGE and AnyBURL excluded — AnyBURL produced no usable disease→gene scores for zero-shot family holdout (unseen diseases) |
| Regeneration within tolerance | PASS | All scripts are `--resume`/`--dry-run` capable and manifest-checkpointed; canonical tables reconcile with published counts (M1/M2) |

## Split integrity

Six splits were built on `evidence_group_id` (§5.4) with `assert_split_disjoint`
run per split and in CI:

| split | track | protocol | eval relation | test rows | leakage guard |
|---|---|---|---|---|---|
| trackr_random | R | random 90/5/5 | all | 405,988 | PASS |
| trackd_random | D | random 90/5/5 | all | 1,093,156 | PASS |
| trackr_family | R | AD/PD family holdout | disease_protein | 914 | PASS |
| trackd_family | D | AD/PD family holdout | DIS-GEN | 18,625 | PASS |
| trackr_txgnn_neurodegenerative | R | TxGNN area | indication/contraindication | 514 | PASS |
| trackd_temporal | D | temporal (secondary) | DIS-GEN | 252 | PASS |

Track-R family nodes were expanded to `mondo_grouped:*` nodes because individual
MONDO ids in `configs/disease_families.yaml` do not exist as Track-R disease
nodes (CHANGES/0009). Family membership still awaits owner sign-off (§0.3.9).

## Baselines and nulls

- **Degree** — total training degree of each candidate gene.
- **PPR** — personalized PageRank from each query disease (restart α=0.15).
- **Degree-matched null** — negatives sampled from the same degree bin as the
  true answer (popularity-bias control, §0.2.4).
- **Relation-shuffled null** — relation labels permuted (structural null).
- **PyKEEN KGE wrapper** — `src/sina/kge.py` (TransE/ComplEx/RotatE-compatible);
  unit-tested; full-scale pilots not run (local runtime too large).
- **AnyBURL** — `src/sina/anyburl.py` (23-1x, BSD-3); learned rules on both
  family splits; zero-shot family holdout produced no usable disease→gene
  scores because family diseases are unseen in the training graph.

## Metrics and statistics

Filtered and raw MRR/Hits@1/3/10 computed with
`sina.metrics.link_prediction_metrics`; degree-stratified metrics with
`degree_stratified_metrics`. Multi-seed statistics (mean±SD, 95% CI, paired
Wilcoxon) are implemented in `sina.metrics` but require KGE/AnyBURL scores
across seeds and are pending.

Headline results (filtered MRR):

| split | degree | degree-matched null | PPR |
|---|---|---|---|
| trackr_family | 0.0026 | 0.0005 | 0.0000 |
| trackd_family | 0.0620 | 0.0111 | 0.0000 |

The degree baseline beats the degree-matched null on both tracks, indicating
that degree alone captures some disease-association signal. PPR returns no
scores for zero-shot family diseases (unseen in train).

## Limitations and pending items

1. **KGE pilots not run** — PyKEEN wrapper exists; full-scale training on 8M/21M
   edges is a multi-hour job and is deferred to a scheduled long-running job.
2. **AnyBURL unusable for zero-shot family** — rule-based methods cannot score
   unseen diseases; AnyBURL should be evaluated on warm-start or non-family
   splits.
3. **Temporal split is secondary** — DIS-GEN year coverage is 0.25% (Appendix
   D.5), so the temporal split is a coverage-limited probe.
4. **Disease families await owner sign-off** — `configs/disease_families.yaml`
   is a draft (§0.3.9).
5. **External eval set** — pre-registered external set is not yet built; it is
   scheduled for M2B.

## Reproduction

- Splits: `scripts/m2_build_splits.py --config configs/m2_splits.yaml`
- Baselines: `scripts/m2_run_baselines.py --config configs/m2_baselines.yaml`
- AnyBURL: `scripts/m2_run_anyburl.py --config configs/m2_anyburl_run.yaml`
- Metrics: `scripts/m2_compute_metrics.py --config configs/m2_metrics.yaml`
- All scripts support `--resume` and `--dry-run`; manifests in
  `data/raw/manifests/`.
