<!-- track_scope: benchmark -->
# M2A-3 AnyBURL rule-mining pilot (spec §6.3)

**Status:** DRAFT — pilot on the family holdout splits only; rule confidence
normalization and ensemble integration happen in M2A-4.

AnyBURL-23-1x.jar (BSD-3) learned rules with snapshots at 10,50,100
seconds per split on the training edges, then applied the 100s
rules to the test edges to produce top-10 rankings. The parsed
confidence scores are stored in
`data/processed/m2/anyburl/<split>/anyburl_scores.parquet`.

## Splits scored

| split | train rows | test rows | queries | scores | learn s | apply s |
|---|---|---|---|---|---|---|
| trackr_family | 8,098,244 | 914 | 55 | 2,280 | 110.247 | 14.954 |
| trackd_family | 21,792,234 | 18,625 | 0 | 0 | 130.981 | 298.389 |

## Reproduction

- Config: `configs/m2_anyburl_run.yaml`
- Manifest: `data/raw/manifests/m2_run_anyburl.json`
- Output: `data/processed/m2/anyburl`
- Re-run: `python scripts/m2_run_anyburl.py --config configs/m2_anyburl_run.yaml --resume`
