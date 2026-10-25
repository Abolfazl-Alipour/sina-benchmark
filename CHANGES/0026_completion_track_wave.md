# 0026 — Completion-Track Wave (within-disease 95/5) — DRAFT (owner signature required)

- **date:** 2026-08-31
- **status:** DRAFT — pending owner sign-off (this file) + `reports/SINA_PAPER.md` at `53f24ca` remains governing until signed; then this entry governs and the directives file is archived
- **context:** LODO benchmark is complete and is an honest negative (mean MRR 0.00775 p≈0.34; pathway funnel Recall@10 0.0; anchored ties topology null). Owner adds completion track (within-disease 95/5 edge splits) as new headline benchmark on two main disease nodes (AD `MONDO:0004975`, PD `MONDO:0005180`) with three pre-registered arms. LODO remains as hard-mode secondary; nothing deleted.

## Frozen decisions for this wave

- **Model config (frozen, do not change without owner sign-off):** BioPathNet `b051a90`, `input_dim 32`, `hidden_dims [32]*6`, `message_func distmult`, `aggregate_func pna`, `num_mlp_layer 2`, `layer_norm`, `dependent true`, `short_cut`, `batch 32`, `Adam lr 1e-3`, `1 full pass` (no `batch_per_epoch`), **1 seed only per owner 2026-08-31 (1337)**, BRG KEEP `1,431,945` edges, 27 relations.
- **Rejected (do NOT implement):** constrained beam through GO/PWY, hard-constrained message passing through pathway nodes, aging-node supervision (efo:0005422), pathway funnel as grant claim, NeighborLoader, >1 pass, depth ablation (deferred).
- **Spend cap:** $150 total GPU for wave (≈18 passes ×98 min ≈30 GPU-h). If >$150 → STOP, re-derive, ping owner.
- **Every output artifact:** `manifest.json` (SHA-256 + size), `persist_snapshot.sh` after each GPU phase.
- **Gates pre-registered before Arm A:** `configs/m2c_completion_wave_stats.yaml` + `reports/m2c/completion_wave_plan.md` (this entry).

## Phases (cheap and blocking first)

### Phase 0 — Data verification (CPU, blocking, ~1 h) — DONE 2026-08-31

- Query `data/raw` hash-verified per `REGISTRY.json` for DIS-GEN on `0004975` (AD) and `0005180` (PD) + descendants/siblings via `mondo.obo`.
- **Result:** AD `12,414` (71.15% of 17448) → 95% 11793 / 5% 621 PASS; PD `6,584` (37.73%) → 6255 / 329 PASS; both >>200, no kill trigger, no ALS fallback, no >5k duplicate flag. Artifact `reports/m2c/completion_data_check.md` (track_scope benchmark, lint 62 files 0 violations).

### Phase 1 — Completion exports (CPU, blocking) — DONE 2026-08-31

- New `scripts/m2_build_completion_exports.py` + `configs/m2c_completion_exports_keep.yaml` (95/5 within-disease by evidence_group_id seed 1337, twins co-assigned, valid 15% of train portion, disjointness `train2∩train1=0` etc). Output `data/processed/m2c/completion_exports/<unit>/` with `train1 1431945` identical, `train2` 189450 (AD)/190647 (PD), `valid` 1882/967, `test` 621/339, `entity_*` 43k, `manifest.parquet` + `manifest.json` (peak RSS 3.09 GB). Cross-check vs Phase 0: AD 621 exact, PD 339 vs 329 +10 group granularity PASS. Tests `tests/unit/test_m2_build_completion_exports.py` 4 tests PASS (now 166 total).

### Phase 2 — CPU baselines on identical splits (blocking, ~1 day) — DONE (1 seed)

- Extend baselines to completion exports: degree, PPR (exact per-disease `train1+train2`, restart 0.5, `ppr_max_iter 100`), DWPC (degree-weighted path counts damping 0.5 sensitivity 0.3/0.5/0.7, metapaths `DIS-PHE→PHE-GEN`, `DIS-DIS→DIS-GEN`, `DIS-GEN→GEN-GEN`, logistic seed 1337, `strict_negative`).
- **1 seed only (1337)** per owner 2026-08-31 (not 3). Artifacts `reports/m2c/completion_baselines.md` (per-disease filtered MRR/Hits@10 for degree/PPR/DWPC) + `data/raw/manifests/m2c_completion_baselines.json`. All CPU, RSS <40GB.

### Phase D0 — ANA-CNS build (CPU, blocking for Arm E) — DONE 2026-08-31 (amended v6)
- `scripts/m2c_build_ana_cns.py` + `configs/m2c_ana_cns.yaml`; source = hash-verified OptimusKG ANA-GEN 8.79M (anatomy_gene 27c5cba6, anatomy 26e86357, UBERON 2026-04-01).
- Filters: (1) anatomy node ∈ UBERON brain/CNS-subtree (curated CNS keywords + roots UBERON:0000955/0001017/0001018, excludes lens cortex; kept 1325 nodes 0.101, list `data/processed/m2c/completion_exports_ana/kept_anatomy_nodes.txt`); (2) specificity >50% of 13,120 (6560) → 0 removed (max 48 <6560).
- Outputs: `train1_ana` 2,829,619 edges (1,431,945 KEEP + 1,397,674 ANA-GEN, 60077 genes) in `data/processed/m2c/completion_exports_ana/<unit>/` (train1 2,829,619, train2 189450/190647, valid 1882/967, test 621/339), disjointness `train2∩train1_ana=0` etc PASS. Re-ran baselines → `reports/m2c/completion_baselines_ana.md` (AD degree 0.00908 PPR 0.00734 DWPC 0.00987, PD 0.00818/0.00752/0.00816, no lift). Tripwire WARNING >500k per owner 2026-08-31 CNS request (v6 brain-only would be 150/83k <500k), flag for review.

## Hardware deviations (owner-authorized 2026-08-31, g7e.2xlarge)

- **v6 directives on record:** `SINA_Coding_Agent_Directives_v6_2026-08-31.md` is now in the repo (folded verbatim 2026-08-31), including §8 pinned clarifications **C1–C5 (binding)**: C1 pinned 17,448-gene universe binds all losses/rankings/floors (incl. D0 baseline re-run filtered to the pinned list + biotype tripwire); C2 multi-hot BCE with per-row edge+inverse removal; C3 filtered 1-hot val loss; C4 aux GEN-BPO/GEN-PWY 1:1 weight 1.0 multi-hot; C5 junk loci stay in the message graph, never scored.
- **Actual node:** `g7e.2xlarge` — 1× NVIDIA RTX PRO 6000 Blackwell **96GB, sm_120**, 8 vCPU, 62GB RAM, driver 595.91.07 (CUDA 13.2). Supersedes the `g6e.xlarge` L40S/48GB hardware section below and the g7e.8xlarge 91.2GB ladder assumption. Ladder ceiling on this node: **91.2GB (95% of 96GB)**.
- **Frozen cu118 env CANNOT run here:** torch 2.0.1+cu118 arch list tops at `sm_90` with no PTX embedded; bare matmul fails `CUDA error: no kernel image is available for execution on the device`. Driver-JIT fallback is impossible.
- **Owner-authorized divergence (2026-08-31, explicit):** modern-torch env `bpn120` — python 3.10, **torch 2.7.1+cu128 / torchvision 0.22.1** (native sm_120), torch-scatter 2.1.2+pt27cu128, torch-sparse 0.6.18+pt27cu128, torch-cluster 1.6.3+pt27cu128, torchdrug 0.2.1 @ `6066fbd` + SINA patches rebuilt from `third_party/torchdrug`, BioPathNet `b051a90` + patches unchanged. **Science risk (accepted by owner):** env no longer matches Gate-0 reproduction pins (torch 2.0.1/py3.10); model config (b051a90, 32×6, distmult/pna, Adam 1e-3, seed 1337) is unchanged — only the numerical backend differs.
- E1 pre-flight on this node is therefore the sm_120 bring-up test; failure → STOP, report to owner (no fallback instance available from this box).

### Phase E1 — Pre-flight (blocking, ~1 h, on the g7e.2xlarge — see deviation above)
- Build biopathnet env on L40S (proven sm_89, cu118 overlay — no Blackwell PTX needed).
- Run 50 training batches + 1 validation pass of full Arm E config.
- PASS → §4 ladder on same node; FAIL → record deviation, retry with clean env (L40S is proven, no Blackwell fallback needed; if persistent, escalate to p5).

### Phase 4 — Batch ladder to 95% VRAM (blocking, after pre-flight)
- On wave node, full-stack config (train1_ana 2.82M + full-vocab 17,448): binary-search max batch fitting ≤45.6GB (95% of 48GB). bs32 was 29.28GB on L40S; ANA + full-vocab add memory — measure, don't assume. Record VRAM + ms/batch, re-derive epoch time from 20 measured batches (completion train2 ≈190k edges ⇒ old 97.8-min invalid), write per-epoch + total wall-clock (2 runs × up to 20 epochs — potentially days) into `reports/m2c/completion_wave_plan.md` before main wave.

### Phase E2 — Full-stack main wave (exact run specification) — PENDING
2 runs, sequential on single GPU: PD first (`mondo__0005180`, 6,255 train), then AD (`mondo__0004975`, 11,793 train). Seed 1337, completion splits (NOT LODO), full-vocab negatives (17,448), GO/pathway multi-task (GEN-BPO 158k + GEN-PWY 46k, 1:1 interleaved, weight 1.0), Adam lr 1e-3 frozen, batch from ladder, per-epoch train loss / val loss (full-vocab BCE) / val MRR+Hits@10 / VRAM / ms/batch → JSONL, early stop val loss patience5 min_delta 3% relative max20, checkpoint every epoch (best-val + last), test scored ONCE at best-val-loss epoch, `persist_snapshot.sh` after each run, hourly checkpoint stream.

### Phase E3 — fork — PENDING
- E wins → backward ablation C (E-ANA), D (E-multi-task), B (E-both), A (64-neg anchor) 1 seed each → attribution table.
- E fails → fall-forward C → D → B → A (owner order, 1 seed each). Certified negative is deliverable.

### Phase 5 — Stats, report, sign-off — PENDING
Stats per pre-registered config (extended for E-first): primary E vs own-graph best baseline (AD DWPC 0.00987, PD degree 0.00818 on ANA graph), secondary ablations, val-loss+MRR curves. Report `reports/m2c/completion_wave.md` + `§3.7` in `SINA_PAPER.md` (abstract untouched). pytest green, lint clean, final `persist_snapshot.sh` + off-host tarball, amend this file → owner signature.

## 6. Re-derivation triggers (RESTORED v6 — any hit → STOP, ping owner)
| Tripwire | Meaning |
|---|---|
| Pre-flight FAIL on L40S (kernels don't build/load) | → retry clean env, record deviation, escalate if persistent |
| bs32 OOM during ladder / ladder regression | Memory model broken — never silently shrink below 32 |
| Loss NaN/divergence in first 200 batches at ladder batch | Halve batch once; second failure → STOP (lr frozen) |
| Pass time >2× re-derived model | Cost/time model broken → re-price |
| ANA-CNS edges >500k / <5k / anatomy mapping coverage <50% | Filter spec broken → owner decides (current 1,397,674 >500k flagged) |
| Arm E < own-graph floor −20% | Fall-forward C first — owner notified |
| Winning arm shuffle collapse <50% | Leakage → halt claims |
| ANA-graph baselines exceed all ANA arms | "Anatomy helps nulls, not model" — report |
| Val loss still improving at epoch 20 | Cap hit → owner decides extension |
| Val loss ↓ but val MRR flat/↓ for 5+ epochs | Loss/metric divergence → report |
| Instance idle >1 h mid-wave (orphan cron fires) | Silent crash → investigate |
| Projected spend >$3,000 | Runaway detector (not budget cap) |

## Hardware (v6)
- **Node: g6e.xlarge** — 1× NVIDIA L40S 48GB, 4 vCPU, 32GB RAM (owner-selected single GPU, 1 seed × both diseases — back to g6e per 2026-08-31).
- Quota check day 0: G/VT on-demand ≥4 vCPU in us-east-1 (account already ran g6e.2xlarge =8 vCPU and g6e.xlarge =4 vCPU — expected cover; verify literally).
- Known risk: **None for L40S** — sm_89 is proven env (same as prior g6e.2xlarge runs); no Blackwell PTX needed; if pre-flight fails, it's env corruption, not arch.
- Never g6 (ARM) and never g7e Blackwell without PTX — L40S is the proven target.

## Definition of done (v6)
1. D0 + E1 + E2 + E3 + Phase 5 artifacts with manifests; tests/lint green; spend log; backups refreshed.
2. Verdict paragraph: *"Full stack on AD/PD completion: MRR [x] vs own-graph floor [y], shuffle collapse [z]%, attribution [table] → claim / no claim."*
3. This file signed by owner → governs, v5 and owner's edited copy deleted.


*Owner signature: ___________________ Date: ___________*
