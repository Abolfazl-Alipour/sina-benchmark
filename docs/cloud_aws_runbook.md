# AWS runbook — Phase 3 BioPathNet Task-G on OptimusKG (owner + next agent)

**Scope:** run the Phase 0.5 smoke + Task-G AD/PD BioPathNet jobs on AWS within
a $100 credit budget ($5,000 more credits coming later). Binding rules: spec §12
(cost caps, billing alert, auto-terminate, orphan cron) and `configs/cloud/m2b_jobs.yaml`.

**Budget math:** g5.xlarge (A10G 24 GB, x86, 4 vCPU) = $1.006/hr on-demand →
~99 hr on $100. Plan needs ~91 GPU-hr (smoke 1 + AD 40 + PD 40 + KGPFLN 10).
Fits if we don't leave the machine idle. Stop the instance between jobs.

---

## Step 0 — Billing guard (5 min, REQUIRED by spec §12.5)

1. Console → **Billing → Budgets → Create budget**.
2. Cost budget, amount **$50**, email alert at 80% and 100%.
3. Also enable **Billing → Cost alerts** for the credit balance itself.

## Step 1 — GPU quota (the blocker; start today)

New AWS accounts have **0 GPU quota**. Fix:

1. First launch a **t3.micro** (free-tier) and let it bill ~$1 — zero-spend
   accounts often get GPU quota auto-denied.
2. Console → **Service Quotas → Amazon EC2 → "Running On-Demand G and VT
   instances"** → **Request increase to 4 vCPUs** (one g5.xlarge).
3. Wait 1–2 business days. If denied, open a support ticket referencing the
   $100 credit and single-instance usage; the upcoming $5,000 helps.
4. While waiting, everything below Step 3 (rsync/bootstrap on a CPU instance)
   can already run.

## Step 2 — Launch the GPU instance

1. **AMI:** "Deep Learning Base OSS Nvidia Driver GPU AMI (Ubuntu 22.04)" —
   has NVIDIA drivers + CUDA preinstalled (saves an hour of pain).
2. **Type:** `g5.xlarge` (A10G 24 GB). NOT g6 (ARM — our envs are x86_64).
3. **Storage:** 200 GB gp3 root volume.
4. **Network:** default VPC; security group = SSH (22) from **your IP only**.
5. **Key pair:** create/download `sina-aws.pem` (never commit it).
6. SSH in: `ssh -i sina-aws.pem ubuntu@<ip>`.

## Step 3 — Bootstrap (≈2 h full, or ≈20 min with rsync)

```bash
# on the instance
export GITHUB_TOKEN=<fine-grained PAT, repo read>
curl -fsSL https://raw.githubusercontent.com/Abolfazl-Alipour/sina-platform/main/scripts/cloud_bootstrap.sh -o /tmp/bootstrap.sh
bash /tmp/bootstrap.sh            # full re-download (~2 h)
# OR: rsync -avP data/ ubuntu@<ip>:/opt/sina/data/  then: bash /tmp/bootstrap.sh --fast
```

The script builds both conda envs (locked), clones BioPathNet `b051a90` +
torchdrug `6066fbd`, downloads all data with hash verification, rebuilds
canonical tables → splits → Task-G exports, and runs the unit tests.

## Step 4 — Jobs (in tmux, checkpointed)

```bash
tmux new -s sina
cd /opt/sina/third_party/biopathnet

# 4a. Phase 0.5 smoke (~$1): mock config must finish on GPU
conda run -n biopathnet --no-capture-output python script/run.py \
  -s 1024 -c /opt/sina/configs/m2b_server_mock.yaml --gpus [0]

# 4b. Task-G AD, 5 seeds (seed = run id, logged per §12.7)
for SEED in 1337 2024 31415 42 777; do
  conda run -n biopathnet --no-capture-output python script/run.py \
    -s $SEED -c /opt/sina/configs/m2b_taskg_run_ad.yaml --gpus [0] \
    2>&1 | tee /opt/sina/logs/taskg_ad_$SEED.log
done

# 4c. Task-G PD, 5 seeds — same with m2b_taskg_run_pd.yaml
```

Record wall time + `nvidia-smi` hours per run into the run manifest
(GPU-hours + $ are required in every manifest, §12.7).

## Step 5 — After jobs: export results, then KILL THE MACHINE

1. `rsync -avP ubuntu@<ip>:/opt/sina/third_party/biopathnet/experiments/ ./results/`
2. `rsync -avP ubuntu@<ip>:/opt/sina/logs/ ./results/logs/`
3. **Terminate the instance** (not stop — stop still bills EBS). Verify in
   console. The hourly orphan cron (§12.5) pages if anything survives >30 min.

## What the next agent picks up from here

1. Merge per-seed BioPathNet score files → metrics via
   `scripts/m2_compute_metrics.py` + `src/sina/stats.py` (mean±SD, 95% CI,
   paired Wilcoxon vs degree/PPR/nulls — already built).
2. KGPFLN arm: onboard repo, **contamination audit first** (its 57 pretraining
   benchmarks vs held-out AD/PD families), then inference on the same splits.
3. Remaining owner sign-off: `configs/m2_trackd_relations.yaml`.

## Gotchas already solved (do NOT rediscover)

- torchdrug needs `setuptools<81`, X11 libs, `ninja` — handled by
  `environment.biopathnet.yml`.
- BioPathNet's `create_working_directory` does `os.chdir` → all dataset paths
  in configs must be ABSOLUTE (our run configs use `/opt/sina/...`).
- Relation names with spaces break AnyBURL (not BioPathNet) — sanitized in export.
- CPU training is measured-dead (~117 days/epoch) — GPU only.
- Claim ceiling: reports must say "disease-associated gene recovery" — never
  "target discovery/validated/safe/cleared" (CI linted).
