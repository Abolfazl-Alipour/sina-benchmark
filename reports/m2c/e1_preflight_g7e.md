<!-- track_scope: benchmark -->
# E1 Pre-flight — g7e.2xlarge (RTX PRO 6000 Blackwell, sm_120) — 2026-08-31

**Verdict: PASS** (with recorded environment deviation — CHANGES/0026 "Hardware deviations").

## Environment actually used (deviation from frozen stack)
- Frozen cu118 stack (torch 2.0.1+cu118) **cannot execute on sm_120**: arch list tops at sm_90, no PTX embedded; bare matmul → `CUDA error: no kernel image is available for execution on the device`. The v6 §5 E1 PTX path (`TORCH_CUDA_ARCH_LIST=8.9;9.0+PTX`) is therefore moot — no PTX exists in the wheel to JIT.
- Owner-authorized modern-torch env `bpn120`: python 3.10, torch **2.7.1+cu128** (arch list includes `sm_120`/`compute_120`), torchvision 0.22.1, torch-scatter 2.1.2+pt27cu128, torch-sparse 0.6.18+pt27cu128, torch-cluster 1.6.3+pt27cu128, torchdrug 0.2.1 @ `6066fbd`, BioPathNet `b051a90`. Model config unchanged.
- torchdrug port (2 source fixes, in `third_party/torchdrug`, not version-controlled):
  1. `layers/functional/extension/{spmm,rspmm}.h`: `#include <ATen/SparseTensorUtils.h>` removed upstream → `ATen/native/SparseTensorUtils.h` + `ATen/TensorUtils.h`.
  2. `utils/torch.py` `LazyExtensionLoader`: torch ≥2.7 added positional `extra_sycl_cflags` to `cpp_extension.load` → call made keyword-based.
- CUDA toolchain: conda `cuda-toolkit=12.8` (nvcc 12.8.93) in env + hostbin overlay `~/persist/cuda-12.8-overlay/hostbin` (conda gcc 13.4 — nvcc rejects host gcc 15). `TORCH_CUDA_ARCH_LIST=12.0`.
- **Restore caveat:** the v3 backup tar did NOT contain the `persist/` overlay or conda envs, and its `third_party/` copies had NOT actually absorbed the two SINA patches (`git apply` no-ops on gitignored paths masked this). Both patches (`torchdrug-streaming-eval`, `biopathnet-strict-negative-guard`) re-applied with `patch -p1` and verified in-file.

## Pre-flight run (stock 64-negative path — full-stack runner pending)
- Config: `configs/m2c_arm_e.yaml` (frozen hyperparams: NBFNet 32×6, distmult/PNA, dependent, short_cut, layer_norm, Adam lr 1e-3, seed 1337), unit `mondo__0005180` (PD), `train1_ana` 2,829,619-edge graph.
- 50 training batches @ bs32: **5.02 min → 6.02 s/batch**, average BCE 0.4315 (finite).
- Engine max GPU memory: **60,624.9 MiB (60.6 GB)** @ bs32. (L40S reference: 29.28 GB @ bs32 on the 1.43M graph — ~2.07× for 1.97× edges; linear in edges.)
- 1 validation pass (967 edges, streaming eval): completed in 206.9s. mrr 0.53856, hits@10 0.6070 — **NOT comparable to the Arm E floors**: stock eval ranks over all 60,235 GEN-type nodes (junk loci score ~0 as zero-init non-query nodes, v6 §8 C5), not the pinned 17,448 universe (C1). Bring-up signal only.
- Artifact: `reports/m2c/e1_preflight_g7e.json` + `.manifest.json`.

## Implications for §4 ladder / E2
- VRAM per query-row ≈ 1.9 GB → projected max batch at the 91.2 GB ceiling ≈ **44–48** (ladder will confirm).
- **Epoch time is the binding constraint:** 190,647 DIS-GEN train edges / 32 ≈ 5,958 batches ≈ **9.9 h/epoch**; with 1:1 aux interleave (C4) ≈ **~20 h/epoch**. 2 runs × up to 20 epochs ⇒ potentially **weeks** on this node. Larger batch reduces step count but sublinearly (message passing dominates). Re-derived numbers go into `reports/m2c/completion_wave_plan.md` before E2 per §4.
