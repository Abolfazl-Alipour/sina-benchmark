<!-- track_scope: benchmark -->
# No-ANA Rank-Loss Full Epoch — PD `mondo__0005180` (K=all, KEEP graph, 2026-09-05)

**One-epoch result of the CHANGES/0027 track:** full filtered-softmax ranking loss
(`--probe-loss rank --probe-k all` — 1 held-out positive vs all 10,864 allowed
negatives per query, knowns masked, gathered from the existing scores tensor),
KEEP-only graph (train1 1,431,945 edges — ANA removed after the inertness proof),
BioPathNet b051a90 32×6 distmult/PNA, Adam lr 1e-3 frozen, seed 1337, B=100,
`--val-every-seconds 3600`. Runner `scripts/pod/m2c_arm_e_train.py` @ `40207d0`.
Log `reports/m2c/noana_rank_pd_epoch1.log`; run JSONL + checkpoints
(`best/epoch_001/last.pth`) under `reports/m2c/checkpoints_completion_noana/mondo__0005180/`.

## Training curve (hourly filtered val over pinned 17,448; 967 valid edges)

| wall (UTC) | batch | val loss | val MRR | val H@10 |
|---|---:|---:|---:|---:|
| 22:16 | 229 | 0.0268 | 0.01497 | 0.0228 |
| 23:20 | 458 | 0.0477 | 0.01459 | 0.0228 |
| 00:24 | 687 | 0.0126 | 0.01473 | 0.0259 |
| 01:28 | 916 | 0.0200 | 0.01649 | 0.0341 |
| 02:31 | 1145 | 0.0805 | 0.01521 | 0.0321 |
| 03:35 | 1374 | 0.2346 | 0.01746 | 0.0279 |
| 04:39 | 1603 | 0.0754 | 0.01573 | 0.0279 |
| 05:44 | 1831 | 0.1919 | **0.02018** | 0.0390 |
| **epoch end** | **1907** | **0.0528** | **0.0165** | **0.0380** |

Epoch wall-clock **8.84 h** (31,815.7 s; 15,752 ms/batch × 1,907 batches), peak VRAM
90.46 GB allocated (94.85 GiB reserved of 95.6 usable), GPU 100% throughout.
Train loss250 9.77 → 2.75 (init ≈ ln(10,865) = 9.29 uniform floor, as expected);
no NaN, no divergence — §6 tripwires never armed. Batch-1 losses of the width
probe runs reproduce the same uniform-init floor.

## Hidden-test result (test scored ONCE, filtered over pinned 17,448; 339 queries)

| model | test MRR [95% CI] | H@10 [95% CI] | median rank |
|---|---:|---:|---:|
| **no-ANA rank epoch (epoch_001.pth)** | **0.01887 [0.01297, 0.02704]** | **0.02950 [0.01180, 0.04720]** | 170 |
| ANA wave ep1 (rank-best, BCE) | 0.01883 [0.01300, 0.02695] | 0.02950 | — |
| ANA wave ep3 (loss-best, BCE) | 0.01863 | 0.02950 | — |
| supervision-only null (no graph, CPU) | 0.01484 | 0.02950 | — |
| degree floor (pinned, C1.4) | 0.00818 | — | — |

Read-out:

1. **No-ANA + rank loss costs nothing.** 0.01887 vs the ANA wave's 0.01883/0.01863 —
   CIs overlap almost completely — on half the graph (1.43M vs 2.83M edges) with a
   2.2× faster epoch (8.8 h vs ~20 h). ANA edges contributed nothing to the scoreboard,
   again (third independent confirmation after the floor and near-constant-feature proofs).
2. **Model > floor (2.31×) and > supervision-only null (+27% MRR), but H@10 equals
   the null exactly** (0.02950 — the same-count easy-hit set). The MRR lift lives in
   mid-rank improvements, not in the top-10 a buyer reads.
3. Val-side hourly MRR peaked at the last hour (0.02018 @ batch 1831) then settled at
   0.0165 on the full epoch-end val — all within the ±20% noise band; no evidence of
   overfitting collapse in-epoch.

## Product-shaped eval (OT panels, this checkpoint)

Re-scored with `epoch_001.pth` (script `scripts/pod/m2c_ot_panel_eval.py`, output
`reports/m2c/ot_panels/score_noana_epoch001.json`; panels re-pulled and committed —
the 2026-09-04 pull was never committed and its coverage numbers reproduce EXACTLY:
PD 87/106/141, AD 52/111/235 at K=100/200/500; ANA-ep3 control reproduces its full
§5 row including the 4/6/6 hit split):

| PD panel (excl. knowns) | K=100 | K=200 | K=500 |
|---|---:|---:|---:|
| OT coverage | 0.2566 | 0.3127 | 0.4159 |
| model two-stage recall@10 | 0.0295 [0.0147, 0.0472] | 0.0295 | 0.0295 |
| OT-ordering null two-stage | 0.0295 | 0.0295 | 0.0295 |
| hit split (both/model-only/null-only) | 4/6/6 | 4/6/6 | 4/6/6 |

Identical totals to the ANA-ep3 checkpoint, but only **6/10 hit genes shared** between
the two models — the "easy set" is ~60% model-independent, not absolute. Still under
half the 0.10 MVP bar; OT coverage remains the only generator signal worth keeping.

## Union ranker + Probe-1 panel — ensemble hope tested and closed (2026-09-05)

The 6/10 hit overlap raised one cheap question: does a two-model union (ANA-ep3 +
no-ANA epoch) lift two-stage recall? Score vectors dumped
(`reports/m2c/ot_panels/score_vectors/`), combined per panel
(`m2c_ot_panel_eval.py union`, output `union_noana_x_anaep3.json`):

| PD panel | union upper bound (A-top10 ∪ B-top10) | RRF fusion two-stage@10 | OT null |
|---|---:|---:|---:|
| K=100 | 0.0413 [14 hits] | 0.0295 | 0.0295 |
| K=200 | 0.0442 [15 hits] | 0.0295 | 0.0295 |
| K=500 | 0.0442 [15 hits] | 0.0295 | 0.0295 |

The union bound shows 14–15 hittable genes exist across the two models, but
rank-fusion (RRF) collapses back to 0.0295 — the model-unique hits do not rank high
enough in the fused list to break into a combined top-10. **Ensembling this stack is
dead without a calibrated fusion rule; certified, not assumed.**

Probe-1 degree panels rebuilt on the KEEP graph (`m2c_rerank_probe.py` +
`--exports-dir/--ckpt-root` flags), no-ANA epoch_001.pth: excl_known coverage
0.0885 (30/339), **recall@10-in-100 0.0295 [0.0147, 0.0472], degree null 0.0088** —
identical to the ANA-stack Probe 1 numbers. Every product-shaped number is invariant
across graph/loss/width-of-stack changes so far: the ceiling sits in the supervision
signal and candidate generation, not in these modeling choices.

## Circularity audit — CLOSED (`reports/m2c/ot_panels/circularity_audit.json`)

Per-gene OT overall-score decomposition (harmonic-sum shares, the OT formula),
evidence classes: literature (europepmc, uniprot_literature), curated_db
(genomics_england, clingen, orphanet, ...), genetics (gwas_credible_sets, eva,
gene_burden, uniprot_variants), functional_other (impc, crispr_screen, ...).

| panel (excl. knowns) | literature | curated_db | genetics | functional | genes >50% lit+curated |
|---|---:|---:|---:|---:|---:|
| PD K=100 | 27.3% | 0.3% | 19.3% | 53.0% | 28% |
| PD K=500 | 17.1% | 0.1% | 6.2% | 76.7% | 16% |
| AD K=100 | 13.1% | 0.0% | 55.2% | 31.7% | 13% |
| AD K=500 | 33.8% | 0.0% | 35.1% | 31.1% | 33% |

The literature-circular channel (same evidence class as the DIS-GEN supervision) is
bounded at ~27% of score mass for the PD top-100 and falls with panel depth; AD panels
are genetics-dominated. Caveat (stated, not fixed): OptimusKG DIS-GEN edges carry no
per-edge provenance in-repo, so overlap is quantified at the evidence-class level.

## Verdict

Rank loss + KEEP-only is the strictly better engineering stack (equal test score, half
the graph, 2.2× faster) and the discovery-track shape is unchanged: genome-wide
candidate prioritization beats floors and the no-graph null on MRR but not on the
top-10 a user reads; the reranker product number stays ~0.03 against a 0.10 bar.
Certified-negative record extended, not reversed. Next lever running: width 32 vs 64
probe (`reports/m2c/width32.log` / `width64.log`, parallel, hourly val — loser killed
on val-curve divergence).

Claim ceiling: disease-associated gene recovery / candidate prioritization only.

## Width probe 32 vs 64 — VERDICT 2026-09-05 11:03Z

Two parallel B=32 K=all rank runs on the KEEP graph, seed 1337, hourly filtered val
(`reports/m2c/ckpt_width32|64/.../train_log_*.jsonl`):

| wall ~h | width32 val MRR (batch) | width64 val MRR (batch) |
|---|---:|---:|
| 1 | 0.01338 (262) | 0.01505 (174) |
| 2 | 0.01346 (526) | 0.01383 (350) |
| 3 | 0.01373 (792) | 0.01457 (524) |
| 4 | **0.01519 (1054)** | 0.01394 (699) |

No val-loss divergence or NaN in either (loss range 0.01–0.10, same noise pattern as
the B=100 run). At hour 4 the gap is +9% for width32 — **inside the ±20% noise band;
no width advantage measurable.** Width64 is strictly dominated on economics
(32.3 h/epoch vs 20.1 h under contention, ~2× activation memory) → **killed at 11:03Z
(pid 138049, plain SIGTERM)**. Width32 continues solo (now full-GPU, ~33.7GB) toward
a full epoch for a B=32-vs-B=100 comparison against this report's 0.01887.

**Answer to the owner's width question: 32. Vector width 64 buys nothing within noise
on this track and costs 1.6× wall-clock; the frozen 32-dim config stands.**

**Post-verdict (18:44Z):** width32 finished its full epoch solo — 12.50 h
(7,072 ms/batch × 5,958 batches, peak 30.74GB), final val MRR 0.0151 / H@10 0.021.
Hidden-test ONCE on its `epoch_001.pth`: **test MRR 0.01887 [0.01299, 0.02760], H@10
0.02950, median rank 170** — the identical point estimate to the B=100 run above.
Batch size 32 vs 100 moves nothing on test; B=100 is 1.43× faster wall-clock and is
the confirmed operating point. **Every knob on this stack (graph ±ANA, loss BCE vs
rank, width 32/64, batch 32/100) lands on test MRR ≈ 0.0189, H@10 0.0295 — the ceiling
is the supervision, full stop.**
Next: AD `mondo__0004975` K=all rank epoch on KEEP, width 32, B=100, seed 1337
(launched 18:46Z, ~1,895 batches ≈ 8.8h) — second-disease datapoint under the
owner's no-idle rule; test-once at completion, then the box awaits owner direction.

