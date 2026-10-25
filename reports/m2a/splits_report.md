<!-- track_scope: benchmark -->
# M2A-2 Split workstream — leakage-guarded partitions (spec §6.2)

**Status:** DRAFT — split definitions are engineering outputs; biological
sign-off of disease-family membership is in `configs/disease_families.yaml`.

All splits operate on `evidence_group_id` (§5.4). The leakage guard
(`assert_split_disjoint`) was run for every split and is also enforced by CI in
`tests/unit/test_splits.py`.

## Splits produced

| split | total rows | total groups | train groups | valid groups | test groups | leakage guard |
|---|---|---|---|---|---|---|
| trackd_family | 21,813,010 | 21,813,010 | 21,792,654 | 2,052 | 18,304 | PASS |
| trackd_family_alzheimer | 21,813,431 | 21,813,431 | 21,800,321 | 1,344 | 11,766 | PASS |
| trackd_family_parkinson | 21,813,395 | 21,813,395 | 21,806,149 | 708 | 6,538 | PASS |
| trackd_random | 21,813,816 | 21,813,816 | 19,631,485 | 1,089,175 | 1,093,156 | PASS |
| trackd_temporal | 21,806,425 | 21,806,425 | 21,806,173 | 0 | 252 | PASS |
| trackr_family | 8,099,572 | 4,049,601 | 4,049,172 | 56 | 373 | PASS |
| trackr_family_alzheimer | 8,100,180 | 4,049,905 | 4,049,702 | 19 | 184 | PASS |
| trackr_family_parkinson | 8,099,890 | 4,049,760 | 4,049,534 | 37 | 189 | PASS |
| trackr_random | 8,100,498 | 4,050,064 | 3,644,923 | 202,156 | 202,985 | PASS |
| trackr_txgnn_neurodegenerative | 8,095,964 | 4,047,797 | 4,047,540 | 0 | 257 | PASS |

## Design notes

- **random_edge_split** (Track-R + Track-D): random group-level partition
  (90/5/5). This split is intentionally inflated by shared evidence groups
  (§6.2.4) and is provided for comparability only.
- **family_holdout_split** (Track-R + Track-D): pooled AD + PD family nodes.
  Eval edges touching a family disease go to valid/test; all other edges
  touching a family disease are removed from train (neighborhood suppression).
- **txgnn_area_split** (Track-R, neurodegenerative area): indication/
  contraindication edges for the 68 TxGNN-area MONDO diseases → test;
  95% of remaining family-neighborhood edges removed from train per the TxGNN
  protocol.
- **temporal_split** (Track-D): train = DIS-GEN edges with `year_final ≤ 2018`
  or no years; test = DIS-GEN edges with `year_initial > 2018`; boundary
  straddlers dropped. **Secondary split:** only ~0.25% of DIS-GEN edges carry
  year annotations (Appendix D.5), so this is a coverage-limited probe.

### Track-R family node expansion

Track-R represents several family-level diseases as `mondo_grouped:*` nodes.
The family config lists individual MONDO ids; for Track-R we additionally
include any grouped disease node whose id contains a family member numeric id:

- **alzheimer**: 4 Track-R disease node(s)
- **parkinson**: 5 Track-R disease node(s)

## Config / reproduction

- Config: `configs/m2_splits.yaml`
- Manifest: `data/raw/manifests/m2_build_splits.json`
- Output: `data/processed/m2/splits`
- Re-run: `python scripts/m2_build_splits.py --config configs/m2_splits.yaml --resume`
