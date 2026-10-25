<!-- track_scope: benchmark -->
# M2B Task-G relation curation on OptimusKG (spec §7.4)

**Status: DRAFT — owner sign-off required (§0.3.9)** — the curation list below is a DRAFT. BioPathNet
Task-G runs must not start until the owner signs off (§0.3.9).

## Supervised task edges (prediction target)

| label | edges | note |
|---|---|---|
| DIS-GEN | 9,734,774 |  |

## Message-passing context (never supervised)

| label | edges | note |
|---|---|---|
| GEN-GEN | 327,924 | PrimeKG-derived PPI |
| PHE-GEN | 793,279 | gene-phenotype |
| BPO-GEN | 158,410 | biological process |
| CCO-GEN | 105,309 | cellular component |
| MFN-GEN | 90,933 | molecular function |
| PWY-GEN | 46,977 | pathway |
| DIS-PHE | 157,144 | disease-phenotype |
| DIS-DIS | 44,215 | disease hierarchy (PARENT) |
| DRG-DIS | 70,380 | indication/contraindication/off-label |
| DRG-GEN | 20,694 | drug target/enzyme/transporter |
| DRG-PHE | 13,758 | drug-phenotype |
| BPO-BPO | 44,494 | process hierarchy (IS_A) |
| PHE-PHE | 24,862 | phenotype hierarchy |
| ANA-ANA | 17,082 | anatomy hierarchy (PARENT) |
| MFN-MFN | 12,587 | function hierarchy |
| CCO-CCO | 4,639 | component hierarchy |
| PWY-PWY | 2,819 | pathway hierarchy |

## Downweighted by default (ablation only)

| label | edges | note |
|---|---|---|
| ANA-GEN | 8,787,955 | excluded by default; enable as ablation only (expression-ubiquity noise) |

## Excluded

| label | edges | note |
|---|---|---|
| EXP-GEN | 2,989 | exposure-gene |
| EXP-EXP | 2,443 | exposure hierarchy |
| EXP-DIS | 2,391 | exposure-disease |
| EXP-BPO | 2,260 | exposure-process |
| EXP-MFN | 47 | exposure-function |
| EXP-CCO | 13 | exposure-component |
| DRG-DRG | 1,345,376 | synergy; BRG-eligible on request |
| DRG-BPO | 62 | drug-process indication |

## Built subset

- Supervised: 9,734,774 edges → `taskg_supervised.parquet`
- Message passing: 1,935,506 edges → `taskg_message_passing.parquet`
- Downweighted (kept aside): 8,787,955 edges → `taskg_downweighted.parquet`

## Reproduction

- Config: `configs/m2_trackd_relations.yaml`
- Run: `python scripts/m2_build_taskg.py --config configs/m2_trackd_relations.yaml --resume`
