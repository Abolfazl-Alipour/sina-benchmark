<!-- track_scope: benchmark -->
# M2B-3 Task-G BioPathNet export + cloud run manifest (spec §7.4, §12)

**Status:** DONE (prep) — data is exported and loader-verified; runs are blocked
on Phase 0.5 (owner cloud approval) and disease-family sign-off.

## Export (BioPathNet `biomedical` format)

`data/processed/m2/biopathnet_taskg/trackd_family/`:

| file | role | rows |
|---|---|---|
| `train1.txt` | BRG message passing (curated per `configs/m2_trackd_relations.yaml`; edges touching the 17 held-out family diseases excluded, spec §7.4(2b)) | 1,934,683 |
| `train2.txt` | supervision DIS-GEN (from the leakage-guarded split's train) | 9,714,066 |
| `valid.txt` | eval DIS-GEN | 2,083 |
| `test.txt` | eval DIS-GEN | 18,625 |
| `entity_types.txt` | type id per entity (10 types, all 155,592 entities) | 155,592 |
| `entity_names.txt` | display names | 155,592 |

**Verified:** BioPathNet's own `biomedical` dataset class loads the export
(155,592 entities, 32 relations, sample counts match).

## Cloud plan (spec §12)

Job manifests: `configs/cloud/m2b_jobs.yaml`. Provider order per §12.1: Lambda
Labs → RunPod → Vast.ai → AWS EC2 (only after GPU quota increase; new AWS
accounts have zero GPU quota).

| job | est GPU-hours | est cost | hard cap |
|---|---|---|---|
| smoke (Phase 0.5) | 0.5 | $1 | $50 |
| Task-D reproduction (PrimeKG, 5 seeds) | 30 | $25 | $300 |
| Task-G OptimusKG AD holdout (5 seeds) | 40 | $35 | $300 |
| Task-G OptimusKG PD holdout (5 seeds) | 40 | $35 | $300 |
| KGE baselines Track-D (3 models x 5 seeds) | 150 | $120 | $300 |
| **total** | ~260 | **~$216** | budget $3,000 |

Note: "Lambda" in the spec means Lambda Labs (GPU rental), **not** AWS Lambda
functions — AWS Lambda's 15-minute limit makes it useless for training. AWS EC2
works but needs a GPU-quota increase first (1-2 days).

## Blockers before first real run

1. Owner cloud approval + billing alert (§12.5).
2. Disease-family sign-off (`configs/disease_families.yaml`).
3. Task-G relation curation sign-off (`configs/m2_trackd_relations.yaml`).
