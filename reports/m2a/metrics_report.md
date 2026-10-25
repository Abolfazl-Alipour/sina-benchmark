<!-- track_scope: benchmark -->
# M2A-4 Metrics + statistics (spec §6.4/§6.5)

**Status:** DRAFT — covers degree/PPR/nulls only; KGE and AnyBURL (unusable on
zero-shot family) are excluded. Multi-seed statistics and paired Wilcoxon tests
require KGE/AnyBURL scores across seeds and are pending.

Metrics are computed with `sina.metrics.link_prediction_metrics` (filtered and
raw MRR/Hits@k; filter set = train+valid+test truths) and
`degree_stratified_metrics` (degree-quantile bins).

**Interpretation note (mandatory, §6.5):** zero-shot family holdout runs are
expected to score 50-80% below transductive settings. The degree baseline
captures popularity bias; the degree-matched null shows what a popularity-matched
random ranking achieves. A model must beat both to be considered useful.

## Headline metrics

| split | model | queries | raw MRR | filt MRR | raw h@1 | filt h@1 | raw h@10 | filt h@10 |
|---|---|---|---|---|---|---|---|---|
| trackr_family | degree_scores | 373 | 0.0027 | 0.0028 | 0.0000 | 0.0000 | 0.0107 | 0.0107 |
| trackr_family | ppr_scores | 0 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| trackr_family | degree_matched_nulls | 7,460 | 0.0004 | 0.0004 | 0.0000 | 0.0000 | 0.0008 | 0.0008 |
| trackd_family | degree_scores | 18,304 | 0.0012 | 0.0629 | 0.0002 | 0.0388 | 0.0014 | 0.1165 |
| trackd_family | ppr_scores | 0 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| trackd_family | degree_matched_nulls | 183,040 | 0.0002 | 0.0113 | 0.0000 | 0.0067 | 0.0002 | 0.0215 |
| trackd_family_alzheimer | degree_scores | 11,766 | 0.0011 | 0.0695 | 0.0001 | 0.0379 | 0.0016 | 0.1447 |
| trackd_family_alzheimer | ppr_scores | 0 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| trackd_family_alzheimer | degree_matched_nulls | 117,660 | 0.0002 | 0.0140 | 0.0000 | 0.0073 | 0.0002 | 0.0297 |
| trackd_family_parkinson | degree_scores | 6,538 | 0.0017 | 0.0523 | 0.0005 | 0.0428 | 0.0018 | 0.0639 |
| trackd_family_parkinson | ppr_scores | 0 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| trackd_family_parkinson | degree_matched_nulls | 65,380 | 0.0002 | 0.0057 | 0.0000 | 0.0045 | 0.0001 | 0.0068 |

## Reproduction

- Config: `configs/m2_metrics.yaml`
- Manifest: `data/raw/manifests/m2_compute_metrics.json`
- Output: `data/processed/m2/metrics`
- Re-run: `python scripts/m2_compute_metrics.py --config configs/m2_metrics.yaml --resume`
