# GPU KICKSTART — SINA Phase 3 (BioPathNet on OptimusKG)

**Read this first on any new GPU server.** Everything needed is in the repo
(`github.com/Abolfazl-Alipour/sina-platform`, private, branch `main`).

> **2026-07-30 UPDATE (supersedes PART B below):** the smoke + estimates job
> is DONE — see `reports/m2b/gpu_smoke_estimates.md`. Measured: batch 4 max on
> 46 GB (bs8 OOM), 7.9 s/batch, full epochs non-viable (222 days/epoch);
> `CHANGES/0012` (batch_per_epoch: 1000, gradient_interval: 4) OWNER-APPROVED
> and baked into `configs/m2b_taskg_run_ad.yaml`/`_pd.yaml`. On GPU machines
> build the biopathnet env from **`environment.biopathnet.cu118.yml`** (the
> committed `environment.biopathnet.yml` pins CPU wheels). Pod ops journal +
> resume instructions: `HANDOFF.md` §10. PART B is kept below for history.

---

## PART A — Step-by-step for the human (≈20 min of clicking, ≈2 h of waiting)

### 1. Launch the server

- **GPU:** 1× 24 GB (RTX 4090 / A10G / A40 / L40). RunPod/Vast.ai preferred
  (~$0.34-0.69/hr); AWS g5.xlarge works after the GPU quota increase.
- **Disk:** 200 GB. **OS:** Ubuntu 22.04 (Deep Learning AMI if AWS).
- **SSH:** key auth only, port 22 open to your IP only.

### 2. Log in and set the GitHub token

```bash
ssh -i <key>.pem ubuntu@<server-ip>
echo 'export GITHUB_TOKEN=<your fine-grained PAT, repo read access>' > ~/.sina_env
chmod 600 ~/.sina_env && source ~/.sina_env
```

(The PAT lives locally in `/home/reza/projects/SinaProject/.env` and in
`git remote -v`. Do NOT paste it into chat. Rotate it if exposed.)

### 3. Clone + bootstrap

```bash
sudo mkdir -p /opt/sina && sudo chown ubuntu /opt/sina
git clone "https://${GITHUB_TOKEN}@github.com/Abolfazl-Alipour/sina-platform.git" /opt/sina
cd /opt/sina
bash scripts/cloud_bootstrap.sh          # ~2 h full rebuild
# or: rsync -avP data/ from your workstation, then: bash scripts/cloud_bootstrap.sh --fast
```

Bootstrap builds both locked conda envs (`sina`, `biopathnet`), clones
BioPathNet `b051a90` + torchdrug `6066fbd`, downloads all data with hash
verification, rebuilds canonical tables → splits → Task-G exports, and runs
the unit tests. It ends with `BOOTSTRAP OK` and `nvidia-smi`.

### 4. Hand over to the AI agent (Kimi Code)

Install Kimi Code on the server the same way as locally, `cd /opt/sina`,
and paste **PART B** below verbatim as the prompt.

### 5. When jobs finish

```bash
rsync -avP ubuntu@<server-ip>:/opt/sina/third_party/biopathnet/experiments/ ./results/
rsync -avP ubuntu@<server-ip>:/opt/sina/reports/ ./reports_from_gpu/
```

Then **TERMINATE** the instance (stop ≠ terminate — stop still bills the disk).

---

## PART B — Prompt for the next AI agent (paste verbatim into Kimi Code)

```
You are continuing the SINA platform build on a fresh GPU server at /opt/sina
(git repo, branch main, remote = private GitHub Abolfazl-Alipour/sina-platform).

Read first, in this order:
1. HANDOFF.md — status board, constraints (follow strictly).
2. docs/cloud_aws_runbook.md — cloud rules (cost caps, termination checklist).
3. docs/GPU_KICKSTART.md (this file) — you are in PART B.

Environment is already bootstrapped: conda envs `sina` and `biopathnet` exist,
BioPathNet is at /opt/sina/third_party/biopathnet (commit b051a90), and the
Task-G data exports exist under /opt/sina/data/processed/m2/biopathnet_taskg/.

YOUR JOB (smoke test + estimates, nothing else):

1. GPU smoke (pipeline proof): run BioPathNet on mock data on the GPU —
   cd /opt/sina/third_party/biopathnet && conda run -n biopathnet \
     --no-capture-output python script/run.py -s 1024 \
     -c /opt/sina/configs/m2b_server_mock.yaml --gpus [0]
   It must finish 5 epochs and print test metrics. If it fails, root-cause
   before continuing.

2. Real measurement: run ONE epoch on the Alzheimer's Task-G export on GPU —
   conda run -n biopathnet --no-capture-output python script/run.py -s 1337 \
     -c /opt/sina/configs/m2b_smoke_ad_1epoch.yaml --gpus [0]
   Record: wall time, GPU model + VRAM used (nvidia-smi), batches/min from the
   log (engine logs every 100 batches), time for epoch 0, and eval time.

3. Extrapolate honestly: the AD run config is 4 epochs x 5 seeds; PD same.
   From the measured epoch time, compute estimated GPU-hours and cost at the
   instance's hourly rate for: (a) AD full run, (b) PD full run, (c) total.
   Give a range (measured ±2x), not a point estimate.

4. Write reports/m2b/gpu_smoke_estimates.md containing: instance spec,
   measured numbers, extrapolations, and a GO/NO-GO recommendation for the
   full runs at the $100 budget. Respect the claim ceiling (CI-linted):
   "disease-associated gene recovery", never "target discovery/validated/safe".

5. Update HANDOFF.md stage log + NEXT ACTION in the SAME commit (CI enforces),
   commit as "m2b9: GPU smoke measurements + cost estimates", push to main.

Rules: ruff + pytest green via conda run -n sina; no secrets in git; no data
artifacts in git; do not start the full AD/PD runs — that decision comes
after the human reviews your estimates. Do not run anything outside /opt/sina.
```

---

## PART C — What the next agent does AFTER the human approves the estimates

1. Full runs: 5 seeds each on `configs/m2b_taskg_run_ad.yaml` and
   `configs/m2b_taskg_run_pd.yaml` (tmux, one seed at a time, log each).
2. Merge per-seed scores → metrics via `scripts/m2_compute_metrics.py` and
   `src/sina/stats.py` (mean±SD, 95% CI, paired Wilcoxon vs degree/PPR/nulls).
3. KGPFLN arm: contamination audit FIRST (its 57 pretraining benchmarks vs
   the held-out families), then inference on the same splits.
4. Gate 3 evidence package in `reports/m2b/`.

## Gotchas already solved (do not rediscover)

- BioPathNet `create_working_directory` does `os.chdir` → all dataset paths in
  configs are ABSOLUTE (`/opt/sina/...`). Local-only configs with
  removed local-only configs were superseded by
  the workstation — use `m2b_server_mock.yaml` / `m2b_smoke_ad_1epoch.yaml` here.
- torchdrug needs setuptools<81, X11 libs, ninja — handled by
  `environment.biopathnet.yml`.
- CPU training is measured-dead (~117 days/epoch on 12 threads) — GPU only.
- Peak RSS guard in `sina.manifests` is 40 GB; the split builder is streaming.
