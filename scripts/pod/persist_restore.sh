#!/usr/bin/env bash
# Restore all persistent artifacts and agent sessions from EBS /home/ubuntu/persist
set -uo pipefail

P=/home/ubuntu/persist

echo "[2026-08-17T18:48:05Z] Restoring from EBS persist directory..."

# 1. Restore graph datasets
sudo mkdir -p /opt/sina/data/processed /opt/sina/data/raw
sudo chown -R ubuntu:ubuntu /opt/sina
rsync -a "$P/data_processed/" /opt/sina/data/processed/ 2>/dev/null || true
rsync -a "$P/data_raw/" /opt/sina/data/raw/ 2>/dev/null || true
rsync -a "$P/manifests/" /opt/sina/data/raw/manifests/ 2>/dev/null || true

# 2. Restore checkpoints into repository
mkdir -p /home/ubuntu/sina-platform/reports/m2b/checkpoints_step0b
mkdir -p /home/ubuntu/sina-platform/reports/m2b/checkpoints_wave1
mkdir -p /home/ubuntu/sina-platform/reports/m2b/checkpoints_shuffle_controls

rsync -a "$P/checkpoints_step0b/" /home/ubuntu/sina-platform/reports/m2b/checkpoints_step0b/ 2>/dev/null || true
rsync -a "$P/checkpoints_wave1/" /home/ubuntu/sina-platform/reports/m2b/checkpoints_wave1/ 2>/dev/null || true
rsync -a "$P/checkpoints_shuffle_controls/" /home/ubuntu/sina-platform/reports/m2b/checkpoints_shuffle_controls/ 2>/dev/null || true

# 2b. Wave 2 secondary checkpoints + M3 pathway funnel parquets
mkdir -p /home/ubuntu/sina-platform/reports/m2b/checkpoints_wave2_secondary
rsync -a "$P/checkpoints_wave2_secondary/" /home/ubuntu/sina-platform/reports/m2b/checkpoints_wave2_secondary/ 2>/dev/null || true
mkdir -p /home/ubuntu/sina-platform/data/processed
rsync -a "$P/data_processed/m3/" /home/ubuntu/sina-platform/data/processed/m3/ 2>/dev/null || true

# 2c. Repo data-depot symlinks (stats/funnel resolve data/... from the repo root)
ln -sfn /opt/sina/data /home/ubuntu/sina-platform/data/data
mkdir -p /home/ubuntu/sina-platform/data/processed/m2
ln -sfn /opt/sina/data/processed/m2/lodo_baselines_keep /home/ubuntu/sina-platform/data/processed/m2/lodo_baselines_keep 2>/dev/null || true
ln -sfn /opt/sina/data/processed/m2/lodo_audits_keep /home/ubuntu/sina-platform/data/processed/m2/lodo_audits_keep 2>/dev/null || true

# 3. Restore agent sessions (Kimi Code, Antigravity, Gemini)
rsync -a "$P/agent_sessions/kimi-code/" /home/ubuntu/.kimi-code/ 2>/dev/null || true
rsync -a "$P/agent_sessions/gemini/" /home/ubuntu/.gemini/ 2>/dev/null || true
rsync -a "$P/agent_sessions/antigravity-ide-server/" /home/ubuntu/.antigravity-ide-server/ 2>/dev/null || true
[ -f "$P/agent_sessions/bash_history" ] && cat "$P/agent_sessions/bash_history" >> /home/ubuntu/.bash_history

echo "[2026-08-17T18:48:05Z] Restore complete! System is fully primed for any agent to resume."
