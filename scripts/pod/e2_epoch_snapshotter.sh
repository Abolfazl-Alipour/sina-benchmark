#!/usr/bin/env bash
# e2_epoch_snapshotter — preserve every epoch's checkpoint + auto-commit progress to GitHub.
# Covers the running PD process (started before the runner's per-epoch-save patch) and
# auto-commits artifacts each epoch so progress lands on GitHub without an agent.
# Launch detached: nohup setsid bash scripts/pod/e2_epoch_snapshotter.sh &
set -uo pipefail

REPO=/home/ubuntu/SinaProject
CKPT_ROOT="$REPO/reports/m2c/checkpoints_completion_ana"
LOG="$REPO/reports/m2c/e2_epoch_snapshotter.log"
cd "$REPO"

log() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*" >> "$LOG"; }
log "snapshotter started"

while true; do
    # Wave complete? both summaries exist → final commit and exit
    if [[ -s "$REPO/reports/m2c/arm_e_mondo__0005180.json" && -s "$REPO/reports/m2c/arm_e_mondo__0004975.json" ]]; then
        log "both units complete — final commit and exit"
    fi

    changed=0
    for unit in mondo__0005180 mondo__0004975; do
        dir="$CKPT_ROOT/$unit"
        jsonl="$dir/train_log_$unit.jsonl"
        [[ -f "$jsonl" ]] || continue
        n_epochs=$(wc -l < "$jsonl")
        last="$dir/last.pth"
        [[ -f "$last" ]] || continue
        # Preserve epoch_NNN.pth for any epoch not yet snapshotted.
        # Wait for last.pth to be newer than the jsonl line AND mtime-stable (torch.save finishes),
        # since jsonl is written before last.pth at epoch end.
        target=$(printf "%s/epoch_%03d.pth" "$dir" "$n_epochs")
        if [[ ! -f "$target" ]]; then
            sleep 60
            m1=$(stat -c %Y "$last" 2>/dev/null || echo 0)
            sleep 30
            m2=$(stat -c %Y "$last" 2>/dev/null || echo 0)
            if [[ "$m1" == "$m2" && "$m1" -gt 0 ]]; then
                cp "$last" "$target" && log "snapshot $target" && changed=1
            fi
        fi
    done

    # Auto-commit any artifact drift (checkpoints, jsonl, arm_e summaries, train logs)
    # NOTE: add paths separately — a non-matching glob aborts the whole git add.
    git add reports/m2c/checkpoints_completion_ana/ 2>/dev/null || true
    for f in reports/m2c/arm_e_*.json; do [[ -f "$f" ]] && git add "$f"; done 2>/dev/null
    if ! git diff --cached --quiet; then
        git commit --author='Reza Alipour <reza@alipour.dev>' -m "E2 auto-progress: $(date -u +%Y-%m-%dT%H:%MZ) checkpoints/metrics" >> "$LOG" 2>&1 \
            && git push origin main >> "$LOG" 2>&1 \
            && log "auto-commit pushed"
    elif [[ $changed == 0 ]]; then
        : # nothing new
    fi

    if [[ -s "$REPO/reports/m2c/arm_e_mondo__0005180.json" && -s "$REPO/reports/m2c/arm_e_mondo__0004975.json" ]]; then
        exit 0
    fi
    sleep 600
done
