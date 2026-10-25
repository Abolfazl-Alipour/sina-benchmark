<!-- track_scope: benchmark -->
# AD seed-1337 viability gate — VERDICT: FAIL (spec §7.4, owner gate 2026-07-30)

**Date:** 2026-07-31
**Model:** BioPathNet (NBFNet, hidden 32×6, distmult/pna, num_negative 64,
strict+heterogeneous negatives), seed 1337, OptimusKG Task-G
`trackd_family_alzheimer` export (155,627 entities; train2 9,721,664
supervision edges).
**Protocol:** CHANGES/0012 (owner-approved): batch 4 (46 GB max, verified),
gradient_interval 4, `batch_per_epoch: 1000`, 4 epochs planned.
**Actual training completed:** 2 epochs (4,000 queries/epoch). The run hung
mid-epoch-2 (root cause PROVEN: zero-candidate row in strict negative
sampling → CUDA context fault, `scripts/patches/biopathnet-strict-negative-guard.patch`);
the resume's full test eval was OOM-killed at the finish (46.6 GB container
RAM cap vs ~29 GB prediction accumulation, fix:
`scripts/patches/torchdrug-streaming-eval.patch`). Gate was therefore run on
the best-valid checkpoint (end of epoch 0, valid MRR 0.0179 vs epoch-1's
0.0072 — overfitting trend) via chunked evaluation (3 chunks, weighted mean;
MRR/Hits@k exact under chunking).

## Gate metrics — full AD test set (11,766 triples, filtered ranking vs 155,627 entities)

| metric | BioPathNet seed-1337 (best ckpt) | degree baseline | degree-matched null |
|---|---|---|---|
| filt MRR | **0.0184** | 0.0695 | 0.0140 |
| filt Hits@1 | 0.0116 | 0.0379 | 0.0073 |
| filt Hits@10 | **0.0273** | 0.1447 | 0.0297 |
| filt Hits@100 | 0.0829 | — | — |
| MR | 6,513 | — | — |

**Verdict rule (pre-registered):** PASS iff filt MRR > 0.0695 AND filt
Hits@10 > 0.1447 (beat the degree baseline). **FAIL on both** (0.26× and
0.19× of baseline). The model is above the popularity-matched null on MRR
but not on Hits@10 — i.e. barely better than a degree-matched coin flip.

## Interpretation (honest)

- Training was too light AND unstable: each epoch saw only 4,000 of 9.7M
  supervision edges (0.04%), lr 5e-3 drove train BCE 0.73→0.13 while valid
  MRR FELL 0.0179→0.0072 — classic overfitting to the per-epoch subsample.
- This verdict does NOT condemn the model class: BioPathNet's published
  results use much smaller supervision pools where epochs are real sweeps.
  It condemns this hyperparameter point (lr 5e-3, 1000 batches/epoch, 2
  epochs) on this graph.
- Engineering lessons landed as permanent fixes (both committed and
  auto-applied on rebuild): strict-negative zero-candidate guard; streaming
  eval (no 15 GB prediction accumulation).

## Next step (owner-delegated, in flight)

**tune1** (seed 1338): lr 1.0e-3 (5× lower, anti-overfit), `batch_per_epoch:
2000` (2× supervision per epoch), 4 epochs, batch 4, gradient_interval 4,
streaming-eval patched engine — full valid per epoch + full test at end.
Expected ~21 h. Gate: same rule. If tune1 also fails → one PD pilot (tune2),
then stop and wait for owner.

## Artifacts

- Checkpoints (best = epoch 1 file, plus epoch 2): GitHub release
  `model-checkpoints` — `ad_seed1337_epoch1.pth`, `ad_seed1337_epoch2.pth`
  (SHA-256 recorded in HANDOFF §10).
- Logs: `/workspace/ad_seed1337.log`, `/workspace/ad_seed1337_resume.log`,
  `/workspace/ad_gate_eval.log` (chunked gate eval, 3× ~55 min).
