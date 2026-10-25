<!-- track_scope: benchmark -->
# Probe 2 — Objective-Function A/B: multi-hot BCE vs masked sampled ranking loss (PD, 200 batches, identical batches)

**Date:** 2026-09-04 · **Base:** runner `scripts/pod/m2c_arm_e_train.py` at `14291bb` (+ final-summary probe guard) · **Spec:** owner Probe-2 refined spec 2026-09-04 (re-armed "then do it") · **Unit:** `mondo__0005180` (PD) · **Seed:** 1337 · **Batch:** 32 (ladder-verified ~60GB, NOT the B=50 wave setting) · **Cap:** `--limit-batches 200`, single pass, no checkpoint writes.

**One line:** the ranking objective separates on training-query MRR (+35% at batch 200, beyond the ±20% band) and is 11% faster per step, but the held-out val MRR gap (+15.5%) sits **inside** the noise band — promising, not decisive at 200 batches from random init.

## 1. What was tested

| | bce (frozen default) | rank (probe variant) |
|---|---|---|
| Loss | multi-hot BCE over full [B,17448] row (all train2 positives = 1) | per row: 1 held-out positive (row's own tail) + K=255 negatives, sampled without replacement from pinned 17,448 **excluding all known positives (train2∪valid∪test) of the row's disease**, gather-down from the same [B,17448] scores (no rescoring), log-softmax temperature 1.0, NLL of the positive |
| Forward / message graph | identical | identical (same forward, cheaper loss) |
| Aux (GEN→BPO/PWY) | untouched, both | untouched, both |
| Validation | filtered 1-hot BCE + filtered MRR over pinned 17,448 (C3) — identical evaluator for both | same |

Identical batches: same seed → same shuffle order, same model init; the rank variant's sampling uses a dedicated CPU `torch.Generator(seed)` and does not perturb the CUDA RNG stream. Commands: `scripts/pod/probe2_objective.sh` (runs bce then rank back-to-back).

## 2. Default path bit-identical (proof)

`--dry-run` config dump without flags vs with `--probe-loss bce` is byte-identical (`probe_loss: "bce", probe_out: null`); all probe code paths are guarded by `probe_mode` / `probe_loss == "rank"`. A run with no probe flags behaves exactly as the wave runner.

## 3. Training curves — identical batches (per-50 log, `probe_{bce,rank}_*_b32_s1337.jsonl`)

| batch | bce loss50 | rank loss50 | bce window MRR | rank window MRR | bce running MRR | rank running MRR |
|---:|---:|---:|---:|---:|---:|---:|
| 50 | 0.15896 | 2.72245 | 0.2894 | 0.3807 | 0.2894 | 0.3807 |
| 100 | 0.00339 | 1.71966 | 0.3677 | 0.4730 | 0.3285 | 0.4269 |
| 150 | 0.00278 | 1.69392 | 0.3889 | 0.5122 | 0.3487 | 0.4553 |
| 200 | 0.00260 | 1.66217 | 0.3790 | **0.5105** | 0.3562 | **0.4691** |

(Window/running MRR = filtered MRR of the queried positive on training batches, computed from the same forward scores, C3 filter: train2∪valid positives of the query disease excluded. Loss scales are not comparable across objectives — BCE mean over 17,448 vs NLL over 256, random-init floor ln(256) = 5.55.)

**Reading:** BCE's loss collapses to ~0.003 almost immediately (memorization regime — the wave saw the same at epoch scale) while its window MRR stalls at ~0.38. Rank's loss keeps descending slowly and its window MRR keeps climbing to 0.51. At batch 200 the rank variant leads by **+34.7%** on window MRR and **+31.7%** on running MRR — outside the ±20% noise band on this metric.

## 4. Held-out validation after 200 batches (existing evaluator, unchanged)

| variant | val loss (filtered 1-hot BCE) | val MRR (filtered) | val H@10 |
|---|---:|---:|---:|
| bce | 0.00096 | 0.01291 | 0.0207 |
| rank | 0.03156 | **0.01492** | 0.0207 |

Val MRR gap: **+15.5%** for rank — inside the ±20% band, single point per variant (967 valid queries, no CI computed). Both are far below the wave's epoch-1 val MRR 0.0195, as expected for 200 batches (~3% of an epoch) from random init.

## 5. Speed & memory (measured, not claimed)

| variant | ms/batch | 200-batch wall | reserved VRAM | projected full epoch (5,958 batches) |
|---|---:|---:|---:|---:|
| bce | 12,745 | 2,549 s | 68.68 GB | ~21.1 h |
| rank | 11,325 | 2,270 s | 68.68 GB | ~18.7 h |

Rank is **11.1% faster per step** (cheaper loss computation; the python multi-hot target loop is skipped). **VRAM is identical** — gather-down keeps the [B,17448] forward, and message passing dominates both memory and time. A variant that scores only 256 candidates (skipping the full-vocab forward) would collapse memory/time but is an architectural change, not a loss swap — out of scope for this probe.

## 6. Contamination audit (mask proof)

Per disease, mask = pinned positions of train2∪valid∪test positives; construction asserted (`allowed == 17448 − known`) and a runtime assert fires if a sampled negative ever equals the row's positive. Unmasked-contamination measured by 200×255 draws from the full 17,448 (seed 1338):

| disease node | known positives masked | expected unmasked contamination | empirical (51,000 draws) |
|---|---:|---:|---:|
| **1430 (PD query disease)** | **6,584** | **0.3773** | **0.3780** |
| 18212 | 5,370 | 0.3078 | 0.3089 |
| 6686 | 5,149 | 0.2951 | 0.2939 |
| 36546 | 4,772 | 0.2735 | 0.2726 |
| 37292 | 2,150 | 0.1232 | 0.1219 |
| (27 more diseases) | 1–846 | 0.0001–0.0485 | all within ±0.001 |

**Without the mask, an unmasked ranking loss would have trained on known positives as negatives 37.7% of the time for PD queries.** The mask is not optional; it is proven working (empirical ≈ expected everywhere, full table in `probe_summary_*.json`).

## 7. Popularity bias (pre-registered pitfall c)

Mean train1 degree of genes in the variants' running top-10s:

| batch | bce top-10 mean degree | rank top-10 mean degree |
|---:|---:|---:|
| 50 | 465 | 427 |
| 100 | 448 | 367 |
| 150 | 488 | 355 |
| 200 | **506** | **347** |

BCE drifts **toward** hubs (+9% over the probe); rank drifts **away** (−19%). The sampled-softmax does not show hub amplification here — if anything it counteracts the degree bias that dominates the frozen stack's rankings (wave top-10s were degree-dominated; cf. `rerank_probe_pd.md`).

## 8. Verdict against the pre-registered bar

- **"Rank separates beyond the ±20% MRR noise band?"** — on training-query MRR: **YES** (+35%). On held-out val MRR: **NO** (+15.5%, inside band, one seed, no CI). Overall: **promising, not decisive at 200 batches from random init.**
- Speedup: measured **11%** per step, zero memory change — the 20 h/epoch wall is message passing, not the loss.
- Any switch of the frozen wave criterion requires a CHANGES entry + owner signature; this probe changes nothing by default.

## 9. Honest limitations

- 200 batches ≈ 3.4% of one epoch, from random init; long-horizon behavior (where the wave's loss/MRR divergence appeared) is exactly what a 200-batch probe cannot see. Rank's slow-descending loss may simply need more steps.
- Train-MRR is measured on training queries (optimistic by construction); the val gap is the decision-grade number and it does not clear the band.
- One seed, one disease, one batch size.
- Bug found + fixed during this probe: probe mode initially fell through to the wave's final-summary write (`reports/m2c/arm_e_<unit>.json`); a stray file was created and deleted, and the write is now guarded by `if not probe_mode`. No wave artifacts existed to be harmed (PD halted at epoch 4).

**Artifacts:** `reports/m2c/objective_probe/` — `probe_{bce,rank}_console.log`, `probe_{bce,rank}_mondo__0005180_b32_s1337.jsonl`, `probe_summary_{bce,rank}_mondo__0005180_b32_s1337.json`, `probe2_wrapper.log`. Flag diff: `git show 14291bb -- scripts/pod/m2c_arm_e_train.py`.
