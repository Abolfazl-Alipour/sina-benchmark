#!/usr/bin/env bash
# cold_start.sh — recreate the full SINA environment on a wiped-clean server
# using GitHub alone (no data disk required). Idempotent; safe to re-run.
#
# Target: Ubuntu 22.04, optional NVIDIA GPU (CUDA-capable). ~1-2 h first time.
#
# Usage:
#   export GITHUB_TOKEN=<fine-grained PAT with repo read access>
#   bash scripts/pod/cold_start.sh [WORKDIR] [--cpu]
#     WORKDIR  repo checkout location (default /opt/sina)
#     --cpu    build the CPU biopathnet env instead of the cu118 GPU env
#
# What it does: base packages -> miniconda -> clone -> `sina` env (locked)
# -> `biopathnet` env (cu118 GPU or CPU) -> BioPathNet/torchdrug pins +
# SINA patches -> hash-verified data downloads -> canonical tables ->
# crosswalk -> splits -> Task-G -> BioPathNet exports -> unit tests.
# Checkpoints are NOT part of this; fetch them from the GitHub release
# `model-checkpoints` if needed (see PERSIST_RESUME.md Scenario C).
# NOTE (2026-09-06): the cu118 GPU env below CANNOT execute on sm_120 Blackwell
# (no PTX in wheel — see reports/m2c/e1_preflight_g7e.md). On sm_120 build the
# bpn120 recipe in PERSIST_RESUME.md Scenario C instead of step [5/8] GPU path.
set -euo pipefail

WORK="${1:-/opt/sina}"
CPU_FLAG="${2:-}"
REPO="Abolfazl-Alipour/sina-platform"
BPN_PIN="b051a9003a1363b5655cdb7dfdad60609c0e4fe3"
TD_PIN="6066fbd82360abb5f270cba1eca560af01b8cc90"

: "${GITHUB_TOKEN:?export GITHUB_TOKEN=<fine-grained PAT with repo read access>}"

echo "=== [1/8] base packages ==="
sudo apt-get update -qq
sudo apt-get install -y -qq git curl tmux rsync build-essential > /dev/null

echo "=== [2/8] miniconda ==="
MC="$HOME/miniconda3"
if [ ! -d "$MC" ]; then
  curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/mc.sh
  bash /tmp/mc.sh -b -p "$MC"
fi
source "$MC/etc/profile.d/conda.sh"
# newer conda requires ToS acceptance for the default channels (idempotent)
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main || true
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r || true

echo "=== [3/8] repo ==="
sudo mkdir -p "$WORK" && sudo chown "$USER" "$WORK"
if [ ! -d "$WORK/.git" ]; then
  git clone "https://${GITHUB_TOKEN}@github.com/${REPO}.git" "$WORK"
fi
cd "$WORK" && git pull --ff-only
git config user.name "Reza Alipour" && git config user.email "reza@alipour.dev"

echo "=== [4/8] sina env (locked) ==="
conda env create -f environment.lock.yml 2>/dev/null || conda env update -f environment.lock.yml --prune
conda run -n sina --no-capture-output pip install -e . --no-deps

echo "=== [5/8] biopathnet env ==="
if [ "$CPU_FLAG" = "--cpu" ]; then
  ENV_YML="environment.biopathnet.yml"
else
  ENV_YML="environment.biopathnet.cu118.yml"   # GPU: torch 2.0.1+cu118 wheel set
fi
conda env create -f "$ENV_YML" 2>/dev/null || conda env update -f "$ENV_YML" --prune

echo "=== [6/8] third_party pins + SINA patches ==="
mkdir -p third_party
[ -d third_party/biopathnet ] || git clone https://github.com/emyyue/BioPathNet.git third_party/biopathnet
git -C third_party/biopathnet fetch --all -q
git -C third_party/biopathnet checkout "$BPN_PIN"
git -C third_party/biopathnet apply --check "$WORK"/scripts/patches/biopathnet-strict-negative-guard.patch 2>/dev/null \
  && git -C third_party/biopathnet apply "$WORK"/scripts/patches/biopathnet-strict-negative-guard.patch \
  || echo "biopathnet guard patch already applied"
[ -d third_party/torchdrug ] || git clone https://github.com/DeepGraphLearning/torchdrug.git third_party/torchdrug
git -C third_party/torchdrug fetch --all -q
git -C third_party/torchdrug checkout "$TD_PIN"
git -C third_party/torchdrug apply --check "$WORK"/scripts/patches/torchdrug-streaming-eval.patch 2>/dev/null \
  && git -C third_party/torchdrug apply "$WORK"/scripts/patches/torchdrug-streaming-eval.patch \
  || echo "torchdrug streaming-eval patch already applied"

echo "=== [7/8] data: hash-verified downloads -> tables -> splits -> Task-G exports ==="
run() { echo "--- $*"; conda run -n sina --no-capture-output "$@"; }
run python scripts/m1_download.py --config configs/m1_trackr.yaml
run python scripts/m1_download_trackd.py --config configs/m1_trackd.yaml
run python scripts/m1_download.py --config configs/m2_mondo.yaml
run python scripts/m1_download.py --config configs/m2_txgnn.yaml
run python scripts/m1_ingest_trackr.py --config configs/m1_trackr_ingest.yaml
run python scripts/m1_ingest_trackd.py --config configs/m1_trackd_ingest.yaml
run python scripts/m2_build_crosswalk.py --config configs/m2_crosswalk.yaml
run python scripts/m2_build_splits.py --config configs/m2_splits.yaml
run python scripts/m2_build_taskg.py --config configs/m2_trackd_relations.yaml
run python scripts/m2_export_taskg_biopathnet.py --config configs/m2b_taskg_export.yaml

echo "=== [8/8] verify ==="
conda run -n sina --no-capture-output python -m pytest tests/unit -q
conda run -n sina --no-capture-output python -m sina.lint_reports reports/ --config configs/report_lint.yaml
if [ "$CPU_FLAG" != "--cpu" ]; then
  conda run -n biopathnet --no-capture-output python -c \
    "import torch; print('torch', torch.__version__, '| cuda?', torch.cuda.is_available())"
  nvidia-smi || echo "WARNING: no GPU visible"
fi
echo "COLD START OK — next: HANDOFF.md §1 (env line) + §4 (run queue), PERSIST_RESUME.md Scenario C"
