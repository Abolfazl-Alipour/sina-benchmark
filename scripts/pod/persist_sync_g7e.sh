#!/usr/bin/env bash
# persist_sync_g7e.sh — mirror everything needed to resume this box's work into
# the persist dir. Called by a 15-min loop; safe to re-run any time (rsync -a).
#
# TRUTH (2026-09-05): this instance has NO separate data volume — ~/persist lives
# on the 100GB ROOT EBS (nvme0n1). This mirror therefore covers: instance STOP
# (root survives stop), root-volume detach+attach to a replacement instance, and
# accidental deletion of ~/SinaProject. It does NOT survive instance TERMINATE
# with DeleteOnTermination — that path is covered by GitHub (code+reports+ckpts
# in-repo, data/processed/m2c as a release asset) + env rebuild per E1.
# If the owner attaches a real data volume, mount it and re-run with:
#   PERSIST_ROOT=/mnt/persist bash scripts/pod/persist_sync_g7e.sh
set -uo pipefail

PERSIST_ROOT="${PERSIST_ROOT:-/home/ubuntu/persist}"
REPO=/home/ubuntu/SinaProject
echo "[$(date -u +%FT%TZ)] persist sync -> $PERSIST_ROOT"

# 1. Repo working tree: code, reports (incl. live checkpoints + logs), .git,
#    data/processed/m2c + data/raw manifests. Excludes m1/m2 processed (8.7G+593M
#    of older-wave data, originals on beast + in the v3 bundle tar).
mkdir -p "$PERSIST_ROOT/sinaproject"
rsync -a --delete \
  --exclude 'data/processed/m1/' --exclude 'data/processed/m2/' \
  --exclude '__pycache__/' --exclude '*.pyc' \
  "$REPO/" "$PERSIST_ROOT/sinaproject/"

# 2. Conda envs (bpn120 GPU stack + sina CPU stack). Set ENVS=skip to skip
#    (they change only on package installs).
if [ "${ENVS:-sync}" != "skip" ]; then
  mkdir -p "$PERSIST_ROOT/miniconda3/envs"
  rsync -a /home/ubuntu/miniconda3/envs/bpn120/ "$PERSIST_ROOT/miniconda3/envs/bpn120/"
  rsync -a /home/ubuntu/miniconda3/envs/sina/ "$PERSIST_ROOT/miniconda3/envs/sina/"
fi

# 3. Live helper scripts + agent state
mkdir -p "$PERSIST_ROOT/tmp_helpers"
cp -f /tmp/width_pusher.sh /tmp/width_monitor.sh /tmp/noana_pusher.sh \
  "$PERSIST_ROOT/tmp_helpers/" 2>/dev/null || true

echo "[$(date -u +%FT%TZ)] persist sync done"
