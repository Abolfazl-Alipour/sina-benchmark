#!/usr/bin/env bash
# Run full-stack batch ladder 36/40/44/48 on bpn120 and gather peak_vram/ms_per_batch
# Usage: bash scripts/pod/run_ladder.sh   (run from ~/SinaProject on host, NOT via muse.bash)
set -uo pipefail
cd "$(dirname "$0")/../.."
P=~/miniconda3/envs/bpn120
OV=~/persist/cuda-12.8-overlay
export PATH="$OV/hostbin:$P/bin:$PATH"
export CUDA_HOME=$P
export CXX="$P/bin/x86_64-conda-linux-gnu-g++"
export TORCH_CUDA_ARCH_LIST="12.0"
# pull latest fix (discover DIS-GEN via train2 + aux strict=False)
git pull --ff-only 2>&1 | tail -5 || true
for B in 16 24 32 36 40; do
  echo "=== B=$B ===" | tee /tmp/ladder_${B}.log
  $P/bin/python -u scripts/pod/m2c_arm_e_train.py --unit mondo__0005180 --batch-size $B --limit-batches 20 2>&1 | tee -a /tmp/ladder_${B}.log || echo "B=$B FAILED" | tee -a /tmp/ladder_${B}.log
done
echo "=== SUMMARY ===" | tee /tmp/ladder_summary.log
grep -E "discovered|aux pools|batches per epoch|peak_vram|ms_per|train loss|val loss|WARNING aux|Traceback|Error" /tmp/ladder_*.log 2>&1 | tee -a /tmp/ladder_summary.log | tail -100
echo "Logs: /tmp/ladder_16.log /tmp/ladder_24.log /tmp/ladder_32.log /tmp/ladder_36.log /tmp/ladder_40.log + /tmp/ladder_summary.log"
ls -lh /tmp/ladder*.log
