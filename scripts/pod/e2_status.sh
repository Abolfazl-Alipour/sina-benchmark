#!/usr/bin/env bash
# e2_status — one-shot status summary for the E2 wave on g7e.2xlarge (2026-08-31).
# Usage: bash scripts/pod/e2_status.sh
set -uo pipefail
REPO=/home/ubuntu/SinaProject
LOG_DIR="$REPO/reports/m2c"
CKPT_DIR="$LOG_DIR/checkpoints_completion_ana"

echo "=== E2 wave status $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="

if [[ -f "$LOG_DIR/E2_WAVE_FAILED.marker" ]]; then
    echo "!! FAILURE MARKER: $(cat "$LOG_DIR/E2_WAVE_FAILED.marker")"
fi

echo "--- processes ---"
pgrep -af "m2c_arm_e_train|e2_wave_followup" | grep -v pgrep || echo "  (none running)"

echo "--- GPU ---"
nvidia-smi --query-gpu=memory.used,memory.total,utilization.gpu,temperature.gpu --format=csv,noheader

for unit in mondo__0005180 mondo__0004975; do
    log="$LOG_DIR/e2_$( [[ $unit == mondo__0005180 ]] && echo pd || echo ad )_train.log"
    jsonl="$CKPT_DIR/$unit/train_log_$unit.jsonl"
    summary="$LOG_DIR/arm_e_$unit.json"
    echo "--- $unit ---"
    if [[ -f "$summary" ]]; then
        echo "  DONE: $summary"
        grep -oE '"(test_mrr|test_hits_at_10|best_val_loss|best_epoch)": [^,}]*' "$summary" 2>/dev/null | sed 's/^/    /'
    fi
    if [[ -f "$jsonl" ]]; then
        echo "  per-epoch (latest 3):"
        tail -n 3 "$jsonl" | sed 's/^/    /'
    fi
    if [[ -f "$log" ]]; then
        echo "  log tail:"
        tail -n 3 "$log" | sed 's/^/    /'
        pid=$(pgrep -f "m2c_arm_e_train.py --unit $unit" | head -1 || true)
        if [[ -n "$pid" ]]; then
            el=$(ps -o etimes= -p "$pid" 2>/dev/null || echo 0)
            # 20h/epoch estimate from ladder (bs50): rough ETA to next epoch line
            eph=72000
            rem=$(( (eph - el % eph) / 3600 ))
            echo "  elapsed $(( el / 3600 ))h$(( (el % 3600) / 60 ))m — next epoch report in ~${rem}h (20h/epoch est.)"
        fi
    else
        echo "  (not started)"
    fi
done

echo "--- sequencer ---"
tail -n 3 "$LOG_DIR/e2_wave_followup.log" 2>/dev/null | sed 's/^/  /' || echo "  (no watcher log)"
