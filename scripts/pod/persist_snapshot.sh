#!/usr/bin/env bash
# Snapshot must-survive artifacts to the EBS-backed persist dir (survives instance stop).
set -uo pipefail

P=/home/ubuntu/persist
mkdir -p "$P/agent_sessions" "$P/logs" "$P/checkpoints_shuffle_controls"

echo "[2026-08-17T18:47:55Z] Starting EBS persistent snapshot..."

# 1. Processed graph datasets & raw manifests
rsync -a --delete /opt/sina/data/processed/ "$P/data_processed/"
rsync -a /opt/sina/data/raw/REGISTRY.json "$P/data_raw/" 2>/dev/null || true
rsync -a /opt/sina/data/raw/manifests/ "$P/manifests/" 2>/dev/null || true

# 2. All model checkpoints (Step 0b, Wave 1, Shuffle Controls, Wave 2 Secondaries)
rsync -a /home/ubuntu/sina-platform/reports/m2b/checkpoints_step0b/ "$P/checkpoints_step0b/" 2>/dev/null || true
rsync -a --delete /home/ubuntu/sina-platform/reports/m2b/checkpoints_wave1/ "$P/checkpoints_wave1/" 2>/dev/null || true
rsync -a --delete /home/ubuntu/sina-platform/reports/m2b/checkpoints_shuffle_controls/ "$P/checkpoints_shuffle_controls/" 2>/dev/null || true
rsync -a /home/ubuntu/sina-platform/reports/m2b/checkpoints_wave2_secondary/ "$P/checkpoints_wave2_secondary/" 2>/dev/null || true

# 2b. M3 pathway funnel parquets (live in the repo checkout, not /opt/sina)
rsync -a /home/ubuntu/sina-platform/data/processed/m3/ "$P/data_processed/m3/" 2>/dev/null || true

# 3. Reports & Leaderboards
rsync -a /home/ubuntu/sina-platform/reports/ "$P/reports/" --exclude checkpoints_step0b --exclude checkpoints_wave1 --exclude checkpoints_shuffle_controls

# 4. Agent Sessions & State (Kimi Code, Antigravity IDE, Gemini, Bash history)
rsync -a /home/ubuntu/.kimi-code/ "$P/agent_sessions/kimi-code/" 2>/dev/null || true
rsync -a /home/ubuntu/.gemini/ "$P/agent_sessions/gemini/" 2>/dev/null || true
rsync -a /home/ubuntu/.antigravity-ide-server/ "$P/agent_sessions/antigravity-ide-server/" 2>/dev/null || true
cp -f /home/ubuntu/.bash_history "$P/agent_sessions/bash_history" 2>/dev/null || true

# 5. Live logs & Heartbeats
rsync -a /tmp/wave1_benchmark.log /tmp/sina_*.log /tmp/sina_*.json "$P/logs/" 2>/dev/null || true

# 6. g7e.2xlarge E2 wave artifacts (repo at ~/SinaProject on this box — 2026-08-31)
mkdir -p "$P/checkpoints_completion_ana"
rsync -a /home/ubuntu/SinaProject/reports/m2c/checkpoints_completion_ana/ "$P/checkpoints_completion_ana/" 2>/dev/null || true
rsync -a /home/ubuntu/SinaProject/reports/m2c/*.log /home/ubuntu/SinaProject/reports/m2c/arm_e_*.json "$P/logs/" 2>/dev/null || true
rsync -a /home/ubuntu/SinaProject/reports/ "$P/reports_sinaproject/" --exclude checkpoints_completion_ana 2>/dev/null || true

echo "[2026-08-17T18:47:55Z] Snapshot complete successfully." >> "$P/snapshot.log"
echo "[2026-08-17T18:47:55Z] Snapshot finished. Persist directory summary:"
du -sh "$P" "$P"/*
