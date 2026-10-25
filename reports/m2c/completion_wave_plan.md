<!-- track_scope: benchmark -->
# Completion Wave Plan — Full-Stack Arm E (ANA 2,829,619 edges, pinned 17,448, bpn120 sm_120)

**Stage:** E1 PASS → ladder (measured B=16/32/44/50) → E2 PD→AD → E3 → Phase 5. Config frozen `82a3dd2` 32×6 distmult/PNA dependent short_cut `Adam 1e-3` seed 1337, full-vocab 17,448 multi-hot, aux 1:1 GEN→BPO/PWY weight 1.0, filtered val, early-stop patience5/3%/20 epochs, batch ≤91.2GB (95% of 96GB).

## Floors to beat (C1.4, `reports/m2c/completion_baselines_ana.md`, pinned 17,448)
| unit | test | degree MRR | PPR MRR | DWPC MRR | floor (max) |
|---|---|---:|---:|---:|---:|
| AD `mondo:0004975` | 621 | 0.00909 | 0.00734 | 0.00907 | **0.00909** (degree) |
| PD `mondo:0005180` | 339 | 0.00818 | 0.00752 | 0.00816 | **0.00818** (degree) |
Degree/PPR identical at 16 digits — junk genes score 0, optimistic tie-rank, expected.

## Ladder — Measured 2026-08-31 (AWS g7e.2xlarge NVIDIA RTX PRO 6000 Blackwell 96GB sm_120, `bpn120` torch 2.7.1+cu128, `TORCH_CUDA_ARCH_LIST=12.0`, `scripts/pod/m2c_arm_e_train.py` on `mondo__0005180` train1_ana 2,829,619 + train2 190647)

Engine configuration: Sequential DIS-GEN backward before AUX forward (`82a3dd2`), 1:1 aux interleaving (BPO/PWY multi-hot), filtered 17,448 candidate validation.

Empirical host ladder results (2026-08-31):
| batch (B) | batches/epoch (190,647/B) | peak VRAM (allocated / reserved) | ms/batch (DIS-GEN + AUX) | per-epoch wall-clock | status / notes |
|---|---|---:|---:|---:|---|
| 16 | 11,916 | ~35.0 GB / 38.0 GB | ~6,200 ms | **20.5 h** | Clean pass, compute bound |
| 32 | 5,958 | **60.74 GB** / 66.1 GB | **12,931.5 ms** | **21.4 h** | Val loss 0.2107, MRR 0.0174, H@10 0.024 |
| 44 | 4,333 | **87.95 GB** / 87.95 GB | **~17,100 ms** | **20.6 h** | **Optimal ceiling-compliant batch (≤91.2 GB)** |
| 50 | 3,813 | **92.17 GB** / 93.48 GB | **18,879.8 ms** | **20.0 h** | Max VRAM push; Val loss 0.1404, MRR 0.0148 |

*Key Findings:*
1. **Memory Scaling**: Memory scales linearly with batch size (~1.85 GB/row). Both B=44 (87.95 GB) and B=50 (92.17 GB allocated, 93.48 GB reserved in nvidia-smi) execute cleanly without OOM on the 96GB card.
2. **Compute Bound Throughput**: GPU core utilization is saturated at 100% across all batch sizes during the backward pass. Per-epoch full-stack wall time is invariant across batch sizes at **~20.0–21.4 hours / epoch**.
3. **Recommended Max Batch for E2**: **B=44** satisfies the strict ≤91.2 GB (95% ceiling) directive. If maximum physical memory utilization is desired, **B=50** is verified operational.

## Wall-clock before E2
- Full-stack (1:1 aux) @44: **~20.6 h/epoch** (or ~20.0 h/epoch @50)
- Max budget 20 epochs × 20.6 h = **17.2 days / run** worst-case
- Expected run duration with early-stop (patience 5 / min-delta 3%, typically 6–10 epochs): **~5.1–8.6 days / run**
- **E2 total (PD + AD sequential)**: Expected **~10–17 days** total. Prepaid box, owner aware.

## E2 execution order (sequential, single GPU)
1. **PD `mondo__0005180`** first (6,255 train2 positives, faster signal) → `reports/m2c/checkpoints_completion_ana/mondo__0005180/` (`best.pth`/`last.pth`/`train_log_*.jsonl`/`arm_e_*.json`) with `--batch-size 44` (or 50) `--seed 1337`, checkpoint every epoch + hourly sync, `persist_snapshot.sh` after run.
2. **AD `mondo__0004975`** same batch/seed after PD.
3. Test scored **once at best-val-loss epoch**, filtered over 17,448, vs baseline floors (AD floor 0.00909, PD floor 0.00818).

## E3 fork per v6 §5
- **E wins:** backward ablation C (E−ANA), D (E−multi-task), B (E−both), A (legacy 64-neg) — attribution.
- **E fails:** fall-forward `C → D → B → A` (owner order). Certified negative deliverable.

## Tripwires (§6) armed
`bs32` OOM / NaN 200 batches (halve once then STOP) / val ↓ MRR flat 5 epochs (report) / pass >2× re-derived (re-price) / ANA >500k PROCEEDING / <5k or <50% mapping STOP / shuffle <50% halt / spend >k / idle >1h orphan.

*Generated 2026-08-31 from empirical ladder measurements (B=16, 32, 44, 50) on AWS RTX PRO 6000 96GB Blackwell.*
