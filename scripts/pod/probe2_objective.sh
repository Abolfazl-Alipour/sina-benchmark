#!/usr/bin/env bash
# Probe 2 (owner spec 2026-09-04): side-by-side BCE vs rank loss on IDENTICAL batches.
# Both runs: --limit-batches 200 --seed 1337 --batch-size 32 (ladder-safe ~60GB, NOT B=50).
# Outputs ONLY to reports/m2c/objective_probe/ — no checkpoint writes, no wave interference.
set -euo pipefail
cd "$(dirname "$0")/../.."

P=/home/ubuntu/miniconda3/envs/bpn120
OV=/home/ubuntu/persist/cuda-12.8-overlay
export PATH="$OV/hostbin:$P/bin:$PATH"
export CUDA_HOME=$P
export CXX="$P/bin/x86_64-conda-linux-gnu-g++"
export TORCH_CUDA_ARCH_LIST="12.0"

OUT=reports/m2c/objective_probe
mkdir -p "$OUT"

for variant in bce rank; do
  echo "=== probe variant: $variant — start $(date -u +%FT%TZ) ==="
  $P/bin/python -u scripts/pod/m2c_arm_e_train.py \
    --unit mondo__0005180 \
    --batch-size 32 \
    --seed 1337 \
    --limit-batches 200 \
    --probe-loss "$variant" \
    --probe-out "$OUT" \
    > "$OUT/probe_${variant}_console.log" 2>&1
  echo "=== probe variant: $variant — done $(date -u +%FT%TZ) ==="
done
echo "PROBE2_BOTH_VARIANTS_DONE"
