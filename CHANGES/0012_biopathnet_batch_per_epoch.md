# 0012 — BioPathNet Task-G: subsample training batches per epoch (batch_per_epoch)

- **change:** For Task-G BioPathNet runs (spec §7.4), replace full-epoch
  training (all ~9.7M DIS-GEN supervision edges per epoch) with engine-level
  per-epoch subsampling: `train.batch_per_epoch: 1000` (final value per owner,
  sane range 500–2000), `engine.gradient_interval: 4` (effective batch 16 at
  batch-4 memory), and `fast_test: 2000` for per-epoch validation (full test
  ranking once per seed). Model architecture, splits, leakage guards, and the
  evaluation protocol (filtered MRR/Hits@k on the full held-out families) are
  unchanged.
- **rationale:** Measured on the A40 pod 2026-07-30
  (`reports/m2b/gpu_smoke_estimates.md`): 7.9 s/training-batch at batch 4
  (the maximum that fits 46 GB; batch 8 OOM-verified), dominated by PNA
  backward scatter kernels. Full epoch = 2.43M batches ≈ 222 days; the
  committed 4-epoch × 5-seed plan ≈ 12 GPU-years per disease. Both the NBFNet
  authors (ogbl-biokg config, 5.1M edges: `batch_per_epoch: 200`) and the
  BioPathNet authors (`lnctardppi.yaml`: `batch_per_epoch: 150`) use this exact
  mechanism on million-edge graphs. It is the only lever that makes the runs
  fit the $100 budget (measured ~119 GPU-h for AD+PD at 1000 batches/epoch).
- **requester:** measurement-driven (GPU smoke, 2026-07-30); drafted by agent,
  **awaits owner sign-off (§0.3.9)**.
- **impact:** Training sees a random 1000-batch (4,000-query) subset of the
  supervision pool per epoch instead of all 9.7M edges; run manifests must
  record `batch_per_epoch` for reproducibility. Held-out evaluation is
  untouched. Without this change the Task-G BioPathNet runs cannot execute at
  any budget; with it, both family runs fit the $100 credit budget.
- **date:** 2026-07-30
- **approved:** owner — 2026-07-30 (chat: "i confirm the changes")
