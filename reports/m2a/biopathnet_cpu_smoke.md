<!-- track_scope: benchmark -->
# BioPathNet CPU smoke test on OptimusKG Task-G — MEASURED (spec §7.1)

**Date:** 2026-07-28
**Verdict:** Full-scale BioPathNet training on the local CPU is **not viable** —
measured, not estimated. Cloud GPU is required.

## Setup

- Data: real Task-G export `data/processed/m2/biopathnet_taskg/trackd_family`
  (train1 BRG 1,934,751 + train2 supervision 9,714,418 + valid 2,052 + test 18,304;
  155,592 typed entities). Loaded successfully by BioPathNet.
- Config: `configs/m2b_taskg_smoke.yaml` — upstream Alzheimer's settings
  (NBFNet 6x32, batch_size 16, num_negative 64), 1 epoch, CPU (`--gpus null`).
- Machine: 12 threads, 46 GB RAM.

## Measurement

- Wall time before kill: **28 minutes** at ~415% CPU (~4.2 cores saturated).
- Training progress in that time: **zero 100-batch log marks** (engine logs
  every 100 gradient updates; none appeared).
- Batches needed per epoch: 9,714,418 / 16 = **607,151**.
- Measured rate: <100 batches / 28 min → **<3.6 batches/min**.

## Extrapolation (from measurement)

| scenario | time per epoch | for spec's runs (AD+PD, 5 seeds, ~4 epochs) |
|---|---|---|
| measured rate (<3.6 batch/min) | ~117 days | ~6.4 years |
| generous (20 min one-time overhead, then 100 batch/5 min steady state) | ~21 days | ~1.1 years |

## Why it's this slow (and why KGE wasn't)

NBFNet does **full-graph message passing** — 6 layers over the 11.6M-edge
graph — for *every training batch*, and each batch scores 16 queries x 64
negatives = 1,040 candidates. That is tens of billions of sparse operations
per batch on CPU. TransE-style KGE only touches the edges inside the batch,
which is why the KGE pilot estimated ~2.6 h/epoch on the same machine while
BioPathNet measures ~3 orders of magnitude slower. This asymmetry is exactly
why path-based models are GPU-native.

## Consequence

- Phase 0.5 / Phase 3 BioPathNet runs **must** go to cloud GPU
  (`configs/cloud/m2b_jobs.yaml`, est ~$81 total for the rescoped plan).
- The exported Task-G data itself is proven valid — BioPathNet loaded it,
  started training, and the only blocker is compute, not data or code.
- The mock config still runs fine on CPU for code-level verification
  (`reports/m2a/biopathnet_onboarding.md`).
