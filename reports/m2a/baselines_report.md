<!-- track_scope: benchmark -->
# M2A-3 Baseline + null-control pilots (spec §6.3)

**Status:** DRAFT — KGE (TransE/ComplEx/RotatE/InductiveNodePiece) and AnyBURL
pilots are not yet run. This report covers degree, PPR, and the two null
controls on the small-query splits (family holdouts + TxGNN neurodegenerative
area) to measure popularity bias and structural baselines.

## Splits scored

| split | train rows | test rows | candidates | queries | degree | PPR | deg-null | rel-null |
|---|---|---|---|---|---|---|---|---|
| trackr_family | 8,098,714 | 746 | 27,610 | 208 | True | True | True | True |
| trackd_family | 21,792,654 | 18,304 | 61,306 | 13 | True | True | True | True |
| trackd_family_alzheimer | 21,800,321 | 11,766 | 61,306 | 6 | True | True | True | True |
| trackd_family_parkinson | 21,806,149 | 6,538 | 61,306 | 8 | True | True | True | True |

## Timing (seconds)

| split | degree | PPR | deg-matched null | rel-shuffled null |
|---|---|---|---|---|
| trackr_family | 3.191 | 75.632 | 7.197 | 0.001 |
| trackd_family | 8.83 | 4.665 | 116.003 | 0.002 |
| trackd_family_alzheimer | 8.983 | 6.135 | 88.366 | 0.002 |
| trackd_family_parkinson | 8.383 | 5.259 | 51.531 | 0.001 |

## Design notes

- **Degree baseline:** ranks candidate genes by total training degree. Expected
  to be a strong popularity-biased baseline on dense disease-gene graphs.
- **PPR baseline:** personalized PageRank from each query disease on the
  symmetric-normalized training graph (restart alpha = 0.15).
- **Degree-matched null:** for each true edge, sample 10
  negative gene(s) from the same degree bin as the true answer.
- **Relation-shuffled null:** permutes relation labels on the test edge set,
  preserving the structural degree distribution but breaking semantics.

## Reproduction

- Config: `data/processed/configs/m2_baselines.yaml`
- Manifest: `data/raw/manifests/m2_run_baselines.json`
- Output: `data/processed/m2/baselines`
- Re-run: `python scripts/m2_run_baselines.py --config configs/m2_baselines.yaml --resume`
