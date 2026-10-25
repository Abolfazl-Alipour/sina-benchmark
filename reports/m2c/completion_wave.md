<!-- track_scope: benchmark -->
# Completion Wave — final report (M2c completion track, AD/PD, 2026-08-31 → 2026-09-06)

**Verdict (v6 §7 format, adapted to the 0027 divergence):** *Full stack on PD/AD
completion: test MRR 0.01887 (PD) vs own-graph floor 0.00818, supervision-only null
0.01484; H@10 0.02950 equal to the null; AD test MRR 0.01061 vs floor 0.00909
(+17% point, CI overlaps floor), H@10 0.01610 equal to the AD null; stack-invariant
across graph (±ANA), loss (BCE/rank), width (32/64), batch (32/100), disease (AD
repeats the PD pattern at lower absolute level). Shuffle-collapse gate not triggered
(no discovery claim to defend — the claim is floors and the certified negative).
Attribution: nothing in the stack carries the signal — the ceiling is the supervision
itself. → no discovery claim; reranker product numbers preserved (PD 0.0295
two-stage, 3.3× degree null; AD 0.0161 two-stage, panel-independent); candidate-generation
ceiling fixable (union coverage 49.3%@500), rerank conversion not.*

## 1. What was frozen and what diverged

Frozen (v6/CHANGES/0020,0026): BioPathNet b051a90, 32-dim ×6 layers, distmult/PNA,
dependent, short_cut; Adam lr 1e-3 (never touched); pinned 17,448-gene universe (C1);
filtered ranking; completion splits 95/5 built at 1d50468; seeds 1337.
Divergences (owner-directed, CHANGES/0027 + amendment): Blackwell sm_120 required
torch 2.7.1+cu128 (0026); ANA removed from the message graph (buried 2026-09-05);
rank loss (full filtered softmax, K=all) replaced multi-hot BCE; single 1-epoch runs
replaced the 20-epoch early-stop wave after the flat-line evidence.
Unintentional deviation found 2026-09-05 (recorded; changes interpretation, not
numbers): on the no-ANA track the aux multi-task was misrouted by the relation/type
fallback (`m2c_arm_e_train.py:315-316`) — actual aux was GO cellular_component
(94,150) + HPO phenotype-gene (768,844), not GO biological_process + reactome
(43,754; present as message edges only). The ANA wave had the intended aux. Aux
never enters validation/test, so all reported metrics stand; the plateau held with
the accidental aux mix as well, which is one more invariance data point.

## 2. Timeline (all artifacts on origin/main)

| date | step | result |
|---|---|---|
| 08-31 | D0 ANA-CNS build (150→1,325 node expansion) | train1_ana 2,829,619 edges; tripwire >500k PROCEEDING per owner |
| 08-31 | C1.4 baselines on pinned 17,448 | PD floor degree 0.00818; AD floor 0.00909 |
| 08-31→09-04 | E2 ANA wave, PD, B=50, 4 epochs | val MRR 0.0195→0.0158→…→0.0166 (flat-to-down); halted at epoch-4 boundary per owner |
| 09-04 | PD test-once-twice (e2_pd_test_once.md) | ep1 (rank-best) 0.01883 [0.01300, 0.02695]; ep3 (loss-best) 0.01863; H@10 0.02950 both |
| 09-04 | CPU nulls (cpu_nulls_and_coverage.md) | supervision-only null PD test MRR 0.01484, H@10 0.02950; ANA inert ×3 proofs; coverage to K=2000 useless (11% of genome) |
| 09-04 | Probe 2 objective A/B (objective_probe.md) | rank loss +35% train-window MRR, −11% step time vs BCE; contamination audit clean |
| 09-04→05 | No-ANA rank epoch PD, K=all, B=100 | 8.84h; test MRR **0.01887** [0.01297, 0.02704], H@10 0.02950, median rank 170 (noana_rank_epoch.md) |
| 09-05 | Width probe 32v64 | tied within ±20% band; width64 killed; width 32 stands |
| 09-05 | width32 B=32 full epoch | 12.5h; test MRR **0.01887** — identical point estimate |
| 09-05 | Union ranker (RRF) | upper bound 0.044, fusion collapses to 0.0295 — ensemble dead |
| 09-05 | Generator wave 2 (generator_panels_wave2.md) | union_best PD coverage 49.3%@500 (record); two-stage 0.0295 on all 7 generators |
| 09-05 18:46Z→09-06 03:53Z | AD no-ANA rank epoch (K=all, B=100, width32) | 9.02h, 1895 batches, 15.73s/batch, 90.46GB; val loss 0.0785, val MRR 0.01214 (flat ~0.011 all epoch); test MRR **0.01061** [0.00735, 0.01506], H@10 0.01610, median 395 |
| 09-06 | AD test-once + OT rescore + prospective/VEP probes | two-stage 0.0161 = test H@10 at all K (panel-independent); Probe A MRR 0.00300; Probe B top1-known 0.954; VEP anecdotes 6/9 top-1 |

## 3. Primary numbers (PD, hidden test n=339, filtered over pinned 17,448)

| model | test MRR [95% CI] | H@10 [95% CI] |
|---|---:|---:|
| no-ANA rank epoch B=100 | 0.01887 [0.01297, 0.02704] | 0.02950 [0.01180, 0.04720] |
| no-ANA rank epoch B=32 (width32) | 0.01887 [0.01299, 0.02760] | 0.02950 |
| ANA wave ep1 (BCE) | 0.01883 [0.01300, 0.02695] | 0.02950 |
| ANA wave ep3 (BCE) | 0.01863 | 0.02950 |
| supervision-only null (no graph) | 0.01484 | 0.02950 |
| degree floor | 0.00818 | — |

## 3b. Primary numbers (AD, hidden test n=621, filtered over pinned 17,448;
5,655 scorable candidates)

| model | test MRR [95% CI] | H@10 [95% CI] |
|---|---:|---:|
| no-ANA rank epoch B=100 (epoch_001) | 0.01061 [0.00735, 0.01506] | 0.01610 [0.00644, 0.02738] |
| supervision-only null (FULL, all 6 variants H@10-identical) | 0.00975 [0.00643, 0.01417] | 0.01610 |
| degree floor | 0.00909 | — |
| OT-null two-stage | — | 0.01450 |

OT rescore (`m2c_ot_panel_eval.py score`, dump
`ot_panels/score_vectors/mondo__0004975__epoch_001.scores.json`, cross-check
reproduces test-once exactly): two-stage **0.0161 at K=100/200/500** vs OT-null
0.0145 — the reranker hit set is panel-independent (same pattern as PD);
coverage 8.37/17.87/37.84%. Reading: AD repeats PD one level down — MRR +13–17%
over floor/null on the point estimate (CI overlaps both), H@10 byte-identical to
counting (10/621). Median rank 395 of 5,655.

## 4. What this wave certifies

1. **The model beats the floor (PD 2.3×; AD 1.17× on the point, CI overlaps floor)
   and the no-graph null on MRR (PD +27%; AD +9% on the point, CIs overlap) but not on
   the top-10 anyone reads** (H@10 identical to counting co-occurrences both diseases:
   10/339 PD, 10/621 AD).
2. **Nothing in the modeling stack matters within the noise band**: anatomy edges
   (inert ×3), loss function (BCE vs rank), width (32/64), batch (32/100), two-model
   ensembles (RRF). Every road leads to test MRR ≈ 0.0189, H@10 0.0295.
3. **Candidate generation is the fixable ceiling** (union coverage 49.3%@500 PD,
   pheno_shared best independent generator, ot_literature strongest but circular);
   **rerank conversion is not fixable on this supervision** (0.0295 across 10+ panels).
4. Circularity audit closed (evidence-class level; DIS-GEN per-edge provenance
   unannotated in OptimusKG — stated assumption).

## 5. Utilization & spend

Prepaid box (no marginal $). GPU utilization 100% during all training windows;
K=all epochs 94.9GiB/96GB. Total GPU time this wave: ANA 4 epochs × 20h (halted) +
ladder (~8h) + probes (~6h) + no-ANA epochs (8.8h + 12.5h) + AD (9.0h).

## 6. Gates & integrity

- Tests 166 pass, ruff 0, report lint 0 at every commit (CI handoff-freshness enforced).
- Pinned universe asserted (17,448) in every eval; biotype tripwire passed (19,999
  protein-coding); contamination audits (rank-loss mask, OT excl-knowns) clean.
- Reproduction controls: OT panels re-pulled → exact coverage match; ANA-ep3 control
  rescore → exact §5 match (0.0295, 4/6/6, MRR 0.01863).
- Persistence: GitHub (code+reports+checkpoints) + release asset data-m2c-2026-09-05
  + on-box persist mirror (15-min loop) + PERSIST_RESUME.md.

Claim ceiling: disease-associated gene recovery / candidate prioritization only.
