<!-- track_scope: benchmark -->

> **CORRECTION (2026-08-07):** the baseline/null numbers this run was
> compared against (0.0695/0.1447 gate; 0.0140/0.0297 null) were computed
> on the FULL Task-G export, not on this trimmed export. Re-graded against
> the trimmed export's own baselines, the degree baseline is 0.0094/0.0094
> and the null 0.0020/0.0017; the model's 0.0160/0.0137 is ≈3× the
> like-for-like popularity floor. The pre-registered FAIL-vs-gate verdict
> formally stands; the "null-level, nothing learned" interpretation is
> RETRACTED. Full analysis: `reports/m2b/trimmed_regrade_autopsy.md`.
> The validation-decline trajectory and all raw measurements below are
> unaffected.

# Trimmed-graph learning run — AD holdout, seed 1339 (CHANGES/0015 + CHANGES/0016)

**Status: COMPLETE — gate FAILED, run terminated after epoch 4 by owner
directive (CHANGES/0018). Verdict for the multi-seed program: NO-GO.**

Machine: Vast.ai, 1× NVIDIA H200 NVL (141 GB VRAM), 24 CPU threads, 186 GB RAM,
driver 570.133.20 (CUDA 12.8). Torch stack: 2.0.1+cu118 (verified GPU-visible).

Data: trimmed Task-G export `data/processed/m2/trimmed/biopathnet/trackd_family_alzheimer/`
(rebuilt + hash-verified on this machine 2026-08-06; counts match the recorded
build exactly — train1 291,692 / train2 170,574 / valid 1,344 / test 11,766).

## Leakage disclosure (mandatory, CHANGES/0016 — owner-directed)

DRG-DIS edges touching held-out family diseases are RETAINED in the background
graph (AD: 315 edges). Via DRG-GEN they form 2-hop disease→drug→gene proxy
paths to the supervision target. This is a known, accepted, disclosed risk per
the owner's 2026-08-01 directive, against the research memo's recommendation
(`reports/m2b/research_decision_memo.md` §2). Every claim from this run carries
this disclosure. Note for interpretation: the retained DRG-DIS context can only
*inflate* scores via the 2-hop proxy — so the failure recorded below happened
with the deck stacked in the model's favor, not against it.
Guards held at build time: no held-out AD supervision edge in train2; no
DIS-GEN in train1; every held-out AD disease connected in train1.

## Throughput measurement (2026-08-06, this machine, trimmed AD export)

Batch-size search (probes killed as soon as VRAM reading was clear, per owner):

| batch size | VRAM used | verdict |
|---|---|---|
| 16 | 9.4 GB | fits (early probe) |
| 150 | 71.7 GB | fits |
| **300** | **141.7 GB of 143.8 (98.5%)** | **fits — chosen (maximum that fits)** |

Measured on full epochs (no `batch_per_epoch` — full pass over all 170,574
supervision pairs each epoch):

- **7.39 s/batch**, steady across epochs (12 min 19 s per 100 batches).
- **569 batches/epoch** (170,574 / 300, last batch partial) — full-pass
  coverage verified: 569 × 300 ≥ 170,574. No subsampling anywhere.
- **Epoch wall time ≈ 70 min.**
- Valid scoring (1,344 triples): **~8 s**. Full test scoring (11,766 triples):
  **~89 s**. (On the old A40/full-graph setup the same test scoring took ~2.8 h.)
- GPU util 100% throughout; bs300 survived every epoch + scoring with no OOM.

## Validation trajectory (the learning curve — the decisive measurement)

Config: `configs/m2b_trimmed_ad_long.yaml` (seed 1339, lr 1.0e-3, batch 300,
full passes, valid scored every epoch, checkpoint every epoch;
driver `scripts/pod/run_trimmed_long.py`).

| epoch | train BCE (epoch avg) | valid MRR | vs null (0.0140) | note |
|---|---|---|---|---|
| 1 | 0.2880 | **0.0162** | above | best checkpoint |
| 2 | 0.2484 | 0.0157 | above | no improvement |
| 3 | 0.2395 | 0.0119 | **below** | no improvement |
| 4 | 0.2390 | 0.0109 | **below** | no improvement — terminated |

**The signature is memorization without transfer:** training loss falls
monotonically (0.288 → 0.239) while validation ranking quality *declines*
monotonically, crossing below the degree-matched null by epoch 3. This is the
same failure shape as the two earlier failed runs (seed 1337: 0.0179 → 0.0072;
tune1: 0.0100 → 0.0067) — but this time with every previously suspected cause
removed: held-out diseases fully connected in the background graph
(CHANGES/0016, incl. the disclosed DRG-DIS context), full-pass data exposure
(every epoch sees all 170,574 supervision pairs), gentle learning rate (1e-3),
the strict-negative CUDA guard, and the streaming-eval fix. The model class
fits the training answers and does not transfer to the held-out family under
this protocol.

## Final gate metrics (pre-registered; best-valid checkpoint = epoch 1)

Full held-out AD test (11,766 triples, filtered, heterogeneous evaluation):

| metric | model (epoch 1) | degree-matched null | degree baseline (gate) | verdict |
|---|---|---|---|---|
| MRR | 0.0160 | 0.0140 | **0.0695** | FAIL (0.23× of gate; 1.14× of null) |
| Hits@10 | 0.0137 | 0.0297 | **0.1447** | FAIL (0.09× of gate; 0.46× of null) |

(Auxiliary: test Hits@1 0.0045, Hits@100 0.366, MR 884; valid at the same
checkpoint: MRR 0.0162, Hits@10 0.0141.)

The gate required beating the degree baseline on both metrics. The model does
not beat even the degree-matched null on Hits@10. **Gate: FAILED.**

## Termination (CHANGES/0018, owner directive 2026-08-06)

Owner called the experiment after epoch 4's validation score failed to improve
("go back to the drawing board"). The old in-process early-stop (patience 3,
predating CHANGES/0017) fired at the same moment; per the owner directive the
run was NOT resumed. Epochs 5–20 were not trained. Checkpoints for epochs 1–4
are preserved (run directory + owner-synced copy; sha256 in HANDOFF §11).

## GO/NO-GO for the multi-seed program

**NO-GO.** Three independent training regimes (hot lr + subsampling, gentle lr
+ subsampling, gentle lr + full passes on the trimmed graph) all show the same
memorization-without-transfer signature. Spending multi-seed budget on the
current protocol would measure the same null-level result five times.

Honest caveats on the verdict's scope:

- Single seed (1339), single family (AD), per the CHANGES/0015 design; the
  protocol explicitly used one seed to decide learnability before committing
  multi-seed budget — it answered that question.
- The verdict is about **this model class + this graph + this exam**
  (NBFNet-style propagation on trimmed OptimusKG Task-G, zero-shot family
  holdout). It is not a claim that disease-associated gene recovery on
  OptimusKG is impossible by any method.
- The disclosed DRG-DIS retention biased the run *toward* success; it failed
  anyway.

## What this run established (positive results, reusable)

- Full-pass training on the trimmed graph is practical: **70 min/epoch at
  batch 300 on an H200, scoring in seconds** — the protocol machinery (trimmed
  export, driver, guards, patches, per-epoch scoring) works end to end and is
  committed.
- The trimmed export rebuilds bit-identical counts on a fresh machine from
  GitHub alone (cold-start chain + `scripts/m2_build_trimmed_taskg.py`).
- The failure is now decidable in ~5 GPU-hours and ~$15, not weeks.
