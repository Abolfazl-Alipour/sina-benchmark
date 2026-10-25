#!/usr/bin/env bash
# e2_wave_followup — E2 wave sequencer for g7e.2xlarge (2026-08-31).
# Waits for the PD training process to exit; on success persists + launches AD;
# on failure (no summary json) writes a FAILED marker and stops (tripwire semantics:
# do NOT auto-launch AD onto a crashed wave — owner/agent investigates first).
# Launch detached: nohup setsid bash scripts/pod/e2_wave_followup.sh <PD_PID> &
set -uo pipefail

REPO=/home/ubuntu/SinaProject
PD_PID="${1:?usage: e2_wave_followup.sh <pd_training_pid>}"
P=/home/ubuntu/miniconda3/envs/bpn120
OV=/home/ubuntu/persist/cuda-12.8-overlay
export PATH="$OV/hostbin:$P/bin:$PATH"
export CUDA_HOME=$P
export CXX="$P/bin/x86_64-conda-linux-gnu-g++"
export TORCH_CUDA_ARCH_LIST="12.0"

log() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"; }

log "watching PD pid $PD_PID"
while kill -0 "$PD_PID" 2>/dev/null; do sleep 300; done
log "PD process exited"

cd "$REPO"
PD_SUMMARY="$REPO/reports/m2c/arm_e_mondo__0005180.json"
if [[ ! -s "$PD_SUMMARY" ]]; then
    log "FAILED: no PD summary at $PD_SUMMARY — NOT launching AD. Investigate reports/m2c/e2_pd_train.log"
    echo "PD failed $(date -u +%FT%TZ) — no summary json" > "$REPO/reports/m2c/E2_WAVE_FAILED.marker"
    exit 1
fi
log "PD summary found — persisting snapshot"
bash scripts/pod/persist_snapshot.sh || true

log "launching AD mondo__0004975 (batch 50, seed 1337)"
$P/bin/python -u scripts/pod/m2c_arm_e_train.py \
    --unit mondo__0004975 \
    --batch-size 50 \
    --seed 1337 \
    --patience 5 \
    --min-delta 0.03 \
    --max-epochs 20 \
    > "$REPO/reports/m2c/e2_ad_train.log" 2>&1
AD_RC=$?
log "AD exited rc=$AD_RC"

AD_SUMMARY="$REPO/reports/m2c/arm_e_mondo__0004975.json"
if [[ ! -s "$AD_SUMMARY" ]]; then
    log "FAILED: no AD summary at $AD_SUMMARY. Investigate reports/m2c/e2_ad_train.log"
    echo "AD failed $(date -u +%FT%TZ) rc=$AD_RC — no summary json" > "$REPO/reports/m2c/E2_WAVE_FAILED.marker"
    exit 1
fi
log "AD summary found — persisting final snapshot"
bash scripts/pod/persist_snapshot.sh || true
log "E2 wave complete (PD + AD). Summaries: $PD_SUMMARY , $AD_SUMMARY"
