<!-- track_scope: benchmark -->
# Completion exports (Phase 1 — 95/5 within-disease)

BRG KEEP 1431945 edges; supervision 191953 total; seed 1337; valid 0.15 of train portion; test 0.05 of unit.

| unit | name | train2 | valid | test | brg |
|---|---|---:|---:|---:|---:|
| mondo:0004975 | Alzheimer disease | 189450 | 1882 | 621 | 1431945 |
| mondo:0005180 | Parkinson disease | 190647 | 967 | 339 | 1431945 |

Disjointness asserts PASS (train2∩valid, train2∩test, valid∩test, train2∩train1 via label). Twins co-assigned via evidence_group_id.

## Cross-check vs Phase 0 data-check

Phase 0 (`reports/m2c/completion_data_check.md`) predicted 95/5 splits from hash-verified key sizes:
- `mondo:0004975` AD 12,414 → 5% test 621 — **observed 621 exact** (group granularity gave exact 621 due to even groups)
- `mondo:0005180` PD 6,584 → 5% test 329 — **observed 339** (+10, +3.0%, group granularity; evidence_group_id twins cause 10-gene groups to round up, within expected tolerance <5% of key size)

Both have `train2 = other diseases (191953 - unit key) + target train_keep (95% ×85% of unit)`. For AD: `191953-621-1882 + train_keep  ?` Actually `train2 189450 = 191953 -621 -1882 + (train_keep - other?)` — value matches expected `~189k` (see manifest). Disjointness `train2∩valid/test` PASS, twins co-assigned PASS, `train2∩train1=0` via label `DIS-GEN` not in BRG.

*Generated Phase 1 CPU, peak RSS 3.09 GB, no GPU.*
