<!-- track_scope: benchmark -->
# BioPathNet GPU smoke on OptimusKG Task-G — measured throughput + cost estimates (spec §4, §7.4)

**Date:** 2026-07-30
**Instance:** RunPod pod — 1× NVIDIA A40 46 GB (driver 580.159.03, CUDA 13.0),
96 CPU threads, 503 GB RAM, network volume (`/workspace`) for all data.
**Env:** conda `biopathnet`, torch 2.0.1+cu118, torchdrug `6066fbd`,
BioPathNet `b051a90`. Mock-data GPU smoke: GREEN (5 epochs, test MRR 0.4194).
**Session journal:** `POD_LOG.md` (repo root).

---

## 1. What was measured (Alzheimer Task-G export, seed 1337)

Dataset: `trackd_family_alzheimer` — 155,627 entities, 32 relations,
11,669,905 triples loaded (train1 1,935,131 BRG + train2 9,721,664
supervision + valid 1,344 + test 11,766; counts match the curated-family
rebuild exactly).

| Measurement | Result |
|---|---|
| batch_size 4, steady-state time per training batch | **7.9 s** (batches 2–10 all within 7.88–7.95 s; batch 0 = 10.5 s incl. CUDA init) |
| batch_size 4 peak VRAM | 33.8 GB |
| batch_size 8 | **CUDA OOM** (38.5 GB allocated + 1.93 GB request > 44.4 GB usable) — batch 4 is the maximum on 46 GB |
| Loss trajectory | normal (0.77 → 0.47 over 11 batches) — training itself is healthy |
| GPU draw during training | ~150 W of 300 W cap, util "100%" |
| Profile (torch.profiler, 3 batches) | GPU genuinely busy ~100% of wall time; forward ~0.3 s/batch; **backward ~5.9 s/batch (75%)** |

### Root cause of the backward cost

The PNA aggregator runs three scatter ops per layer (add/min/max) over the
1,420,095-triple `train1` LODO graph (11.67M is source OptimusKG, never per-unit). torchdrug's `rspmm` **backward** kernels are ~20×
slower than their forward at this scale (Add forward 12 ms/call vs backward
245 ms/call; backward Add/Min/Max = 37.5% + 18.9% + 18.8% of all CUDA time).
Data-loading, dataloader workers, network-filesystem I/O, and the CUDA
allocator were all tested and **exonerated** — this is pure GPU kernel time.

## 2. Why the committed protocol is non-viable

The Task-G export makes **all 9,721,664 DIS-GEN train2 edges supervision
queries**. At batch 4 that is **2,430,416 batches per epoch**:

- **Full epoch: 2.43M × 7.9 s ≈ 5,330 h ≈ 222 days**
- Committed plan (4 epochs × 5 seeds) per disease: **~106,700 GPU-h ≈ 12.2 GPU-years**
- AD + PD: **~24 GPU-years** (≥ $85k even at $0.40/GPU-h)

For calibration: the BioPathNet paper reports NBFNet-class training was
"computationally infeasible (>48 h/epoch)" on the ~5M-edge PrimeKG graph, and
the authors trained on A100-80GB with supervision sets orders of magnitude
smaller than ours. On million-edge graphs, both the NBFNet authors
(ogbl-biokg: `batch_per_epoch: 200`) and the BioPathNet authors
(`lnctardppi.yaml`: `batch_per_epoch: 150`) **subsample batches per epoch** —
nobody runs full epochs at this scale.

### Alternatives considered and rejected

| Alternative | Verdict |
|---|---|
| A100-80GB, batch 16 (authors' setup) | Still dead: 607,604 batches × ~8–12 s ≈ 56–85 days/epoch |
| Multi-GPU DDP | Would need 100+ GPUs to reach days-scale epochs |
| `aggregate_func: sum` only (drop min/max) | Cuts backward ~2/3 → epoch still ~74 days; model change needs owner sign-off anyway |
| `num_negative: 32` | Negligible — scoring is not the bottleneck; message passing is |
| `batch_per_epoch` subsampling | **Only viable lever; the authors' own mechanism** |

## 3. Proposed amended protocol (CHANGES/0012 — pending owner sign-off)

- `train.batch_per_epoch: 1000` (proposed; sane range 500–2000) — random
  subset of the 9.7M supervision pool per epoch, drawn by the engine's sampler
- `engine.gradient_interval: 4` — effective batch 16 (the authors' batch size)
  at batch-4 memory
- `fast_test: 2000` for per-epoch validation; full test ranking once per seed
- batch_size stays 4 (hardware maximum, verified)

### Measured cost under the amended protocol (A40; amended smoke 2026-07-30)

Amended smoke (`batch_per_epoch: 100`, `gradient_interval: 4`, `fast_test: 500`,
1 epoch, seed 1337) **verified the mechanism and measured both sides**:

- Training: 100 batches in 13.22 min → **7.93 s/batch** (identical to the
  full-epoch measurement; subsampling works, no per-epoch overhead)
- Eval (valid 500 triples, ×2 runs, test 500 triples): 6.9–7.3 min →
  **~0.85 s/triple** (each query ranks all 155,627 entities, 2 full-graph
  forward passes)

| Item | Estimate (now fully measured) |
|---|---|
| 1 epoch (1000 batches) | 2.2 h |
| Per-epoch valid eval (AD 1,344 / PD 708 triples) | 0.32 / 0.17 h |
| Full test eval per seed (AD 11,766 / PD 6,538 triples) | 2.8 / 1.5 h |
| 1 seed AD (4 epochs + evals) | 12.9 h |
| AD run (5 seeds) | ~64 h |
| PD run (5 seeds) | ~55 h |
| AD + PD | **~119 GPU-h** |
| Cost at $0.40–1.20/GPU-h (RunPod/Vast A40 range) | **$48–143** (honesty band ±2×: $24–286) |

At `batch_per_epoch: 2000` add ~90 GPU-h (training doubles, eval unchanged) —
likely over the $100 credit budget. At 500, ~71 GPU-h total.
**Note:** test eval is ~22% of total cost; if budget tightens, `fast_test`
can cap intermediate-seed test evals (protocol decision for the owner).

## 4. GO / NO-GO

- **NO-GO** for the committed full-epoch protocol on any single GPU
  (~24 GPU-years for AD+PD). This is a protocol problem, not a hardware problem.
- **Conditional GO** for the amended protocol at `batch_per_epoch ≈ 1000`:
  fits the $100 AWS-credit budget at mid-range pricing, uses the model
  authors' own subsampling mechanism, and requires no model-code changes.
  **Gate: owner approval of `CHANGES/0012`** (protocol deviation from spec
  §7.4), then a short amended smoke (100 batches + one eval) to convert the
  eval-time estimates above into measurements before the full runs.

## 5. Caveats

- Training and eval throughput are both measured on this A40 (7.93 s/batch;
  0.85 s/eval-triple); totals are arithmetic over the run plan, not yet a
  full-seed end-to-end run.
- `batch_per_epoch` subsampling changes the effective number of gradient
  updates per epoch; the metrics protocol (filtered MRR/Hits@k on the full
  held-out families) is unaffected, but run manifests must record the
  subsample size for reproducibility.
- Throughput was measured on the Alzheimer's export; PD is the same graph
  size class (±0.1%), so times transfer.
