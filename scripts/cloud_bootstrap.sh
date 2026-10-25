#!/usr/bin/env bash
# cloud_bootstrap.sh — take a bare cloud Ubuntu 22.04 machine to a runnable
# SINA state (spec §12.3). Idempotent; safe to re-run.
#
# Usage:
#   export GITHUB_TOKEN=<fine-grained PAT, repo read access>
#   bash scripts/cloud_bootstrap.sh [--fast]
#
# --fast: skip data re-download/rebuild; expects data/ to be rsync'd separately:
#   rsync -avP data/ user@host:/opt/sina/data/
#
# After bootstrap, run jobs per configs/cloud/m2b_jobs.yaml, e.g.:
#   cd /opt/sina/third_party/biopathnet
#   conda run -n biopathnet --no-capture-output python script/run.py \
#     -s 1337 -c /opt/sina/configs/m2b_taskg_run.yaml --gpus [0]

set -euo pipefail

REPO="Abolfazl-Alipour/sina-platform"
WORK=/opt/sina
FAST="${1:-}"

echo "=== [1/7] base packages ==="
sudo apt-get update -qq
sudo apt-get install -y -qq git curl tmux rsync build-essential > /dev/null

echo "=== [2/7] miniconda ==="
if [ ! -d "$HOME/miniconda3" ]; then
  curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/mc.sh
  bash /tmp/mc.sh -b -p "$HOME/miniconda3"
fi
source "$HOME/miniconda3/etc/profile.d/conda.sh"

echo "=== [3/7] repo ==="
sudo mkdir -p "$WORK" && sudo chown "$USER" "$WORK"
if [ ! -d "$WORK/.git" ]; then
  git clone "https://${GITHUB_TOKEN}@github.com/${REPO}.git" "$WORK"
fi
cd "$WORK" && git pull --ff-only

echo "=== [4/7] sina env ==="
conda env create -f environment.lock.yml 2>/dev/null || conda env update -f environment.lock.yml --prune
conda run -n sina --no-capture-output pip install -e . --no-deps

echo "=== [5/7] biopathnet env + third_party pins ==="
conda env create -f environment.biopathnet.yml 2>/dev/null || conda env update -f environment.biopathnet.yml --prune
mkdir -p third_party
[ -d third_party/biopathnet ] || git clone https://github.com/emyyue/BioPathNet.git third_party/biopathnet
cd third_party/biopathnet && git fetch --all -q && git checkout b051a9003a1363b5655cdb7dfdad60609c0e4fe3 && cd "$WORK"
[ -d third_party/torchdrug ] || git clone https://github.com/DeepGraphLearning/torchdrug.git third_party/torchdrug
cd third_party/torchdrug && git fetch --all -q && git checkout 6066fbd82360abb5f270cba1eca560af01b8cc90 && cd "$WORK"

if [ "$FAST" != "--fast" ]; then
  echo "=== [6/7] data download (registry-verified) ==="
  conda run -n sina --no-capture-output python scripts/m1_download.py --config configs/m1_trackr.yaml
  conda run -n sina --no-capture-output python scripts/m1_download_trackd.py --config configs/m1_trackd.yaml
  conda run -n sina --no-capture-output python scripts/m1_download.py --config configs/m2_mondo.yaml
  conda run -n sina --no-capture-output python scripts/m1_download.py --config configs/m2_txgnn.yaml
  conda run -n sina --no-capture-output python scripts/m2_download_anyburl.py --config configs/m2_anyburl.yaml

  echo "=== [7/7] rebuild canonical tables → splits → Task-G export ==="
  conda run -n sina --no-capture-output python scripts/m1_ingest_trackr.py --config configs/m1_trackr_ingest.yaml
  conda run -n sina --no-capture-output python scripts/m1_ingest_trackd.py --config configs/m1_trackd_ingest.yaml
  conda run -n sina --no-capture-output python scripts/m2_build_crosswalk.py --config configs/m2_crosswalk.yaml
  conda run -n sina --no-capture-output python scripts/m2_build_splits.py --config configs/m2_splits.yaml
  conda run -n sina --no-capture-output python scripts/m2_build_taskg.py --config configs/m2_trackd_relations.yaml
  conda run -n sina --no-capture-output python scripts/m2_export_taskg_biopathnet.py --config configs/m2b_taskg_export.yaml
else
  echo "=== [6/7]+[7/7] skipped (--fast): expecting rsync'd data/ ==="
fi

echo "=== verify ==="
conda run -n sina --no-capture-output python -m pytest tests/unit -q
conda run -n sina --no-capture-output python -c "import torch" 2>/dev/null || true
nvidia-smi || echo "WARNING: no GPU visible"
echo "BOOTSTRAP OK"
