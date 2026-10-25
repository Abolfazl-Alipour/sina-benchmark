<!-- track_scope: benchmark -->
# Peer-Review Consultation — SINA Completion Track Status & Open Problems
**Date:** 2026-09-03 · **Repo HEAD at writing:** `0c5064e` · **Author:** SINA GPU agent (for owner Reza Alipour)
**Status 2026-09-06: SUPERSEDED — kept as audit trail only.** The review now lives in `reports/SINA_PAPER.md` §10 (project state, what changed since this snapshot, ordered open questions, feedback protocol), with results in §3.7 and methods in §2.9. Everything below is stale (probes it calls queued are closed; the ANA wave it discusses is buried; the AD epoch it awaits is done: test MRR 0.01061 / H@10 0.01610).
**Original purpose (historical):** everything an external reviewer needs to critique the approach and propose improvements. No conclusions are sought from the reviewer beyond ideas — all binding decisions remain the owner's.

## 0. Context primer (read this first — no repo access needed)
**What SINA is:** a hypothesis-generation and falsification platform built on **OptimusKG**, a biomedical knowledge graph assembled from public sources (PrimeKG-family: genes, diseases, phenotypes, drugs, GO biological processes, pathways, anatomy–gene expression). The scientific goal in this track: given a disease, recover its known associated genes when those associations are hidden from the model — a hard proxy for "can this system prioritize the right candidates."

**The evaluation game (completion track):** hide 5% of known disease–gene associations (test set), train on the rest plus the full background graph, then ask the model to rank all 17,448 protein-coding-relevant genes for that disease. Metric: **filtered MRR** (mean reciprocal rank of the hidden gene, with all other known positives removed from the ranking so they can't inflate or deflate it) and **Hits@10**. "Floors" = the same evaluation done by dumb baselines (degree, personalized PageRank, DWPC) — the GNN must beat them to justify itself.

**Graph schema (26 relations over 86,157 entities):** genes (60,235 nodes, of which the pinned 17,448 are scored), diseases, phenotypes, drugs, GO biological processes (BPO), pathways (PWY), anatomy. Key edge families: DIS-GEN (disease–gene associations, the supervision), DIS-PHE, PHE-GEN, GEN-GEN, GEN-BPO, GEN-PWY, DRUG-*, and ANA-GEN (anatomy→gene expression).

**What Arm E (the thing being trained) adds:** (a) **ANA-CNS edges** — 1.4M anatomy→gene "expressed in" edges filtered to CNS structures (the hypothesis: brain-region expression informs CNS-disease gene ranking); (b) **full-vocab supervision** — every training query scores all 17,448 genes, not 64 sampled negatives; (c) **multi-task aux supervision** — genes are simultaneously trained to predict their GO processes and pathways, 1 aux batch per disease batch.

**Vocabulary used below:** `train1` = background graph (message-passing roads, never labels); `train2` = supervision edges (labels/queries); `KEEP` = the original curated 1.43M-edge background; `LODO` = a separate, stricter leave-one-disease-out track (NOT what runs here — ignore it); Arms A–D = ablations of Arm E (A: legacy 64-neg, B: no ANA no aux, C: no ANA, D: no aux); NBFNet/BioPathNet = a GNN that answers "is (disease, relation, gene) true?" by learning path-based representations from the query node outward — importantly, its scores decompose into contributing paths (explainable).

**Reading guide if reviewers want primary sources:** binding protocol = `SINA_Coding_Agent_Directives_v6_2026-08-31.md` (esp. §8 C1–C5); frozen exam spec = `CHANGES/0020_v5_1_freeze.md`; hardware/env deviation = `CHANGES/0026_completion_track_wave.md`; ladder/wall-clock = `reports/m2c/completion_wave_plan.md`; floors = `reports/m2c/completion_baselines_ana.md`; live metrics = `reports/m2c/checkpoints_completion_ana/mondo__0005180/train_log_mondo__0005180.jsonl`; current state = `HANDOFF.md`.

## 1. The task
Disease-associated **gene recovery / candidate prioritization** on the OptimusKG knowledge graph, completion track. Two diseases:
- Alzheimer's (AD, `MONDO:0004975`): 11,793 train / 621 test disease–gene edges
- Parkinson's (PD, `MONDO:0005180`): 6,255 train / 339 test

Claim ceiling is recovery/prioritization only; stronger downstream-clinical language is forbidden by protocol (see configs/report_lint.yaml).

## 2. Data (all counts verified on disk, hash-manifested)
- Message graph `train1_ana` = **2,829,619 edges** = 1,431,945 KEEP (background KG) + 1,397,674 ANA-GEN (anatomy–gene expression, CNS-filtered: 1,325 CNS anatomy nodes kept via curated keywords + UBERON brain-subtree descendants).
- Supervision `train2` = 190,647 DIS-GEN edges (PD unit) / 189,450 (AD unit); valid 967 / 1,882; test 339 / 621. Disjointness asserts all pass.
- Candidate universe: **pinned 17,448-gene list** (binding; never recomputed). Full GEN node count is 60,235 (19,999 protein-coding) — the ~43k non-pinned loci stay in the message graph but are never scored (directive C5).
- Vocab: 86,157 entities, 26 relations.
- Aux supervision from train1: GEN→BPO 147,477 edges (12,106 GO terms), GEN→PWY 43,754 edges (2,215 pathways).

## 3. Model & training protocol (frozen, v6 directives)
- BioPathNet `b051a90` (NBFNet-family): 32-dim × 6 layers, distmult/PNA, dependent, short_cut, layer_norm. **Adam lr 1e-3 frozen.** 4.8MB checkpoints.
- C2: per-row multi-hot BCE over all 17,448 candidates; only the queried edge + its inverse removed from the message graph per batch row.
- C3: filtered 1-hot val loss (known positives masked out of denominator) + filtered MRR/Hits@10 every epoch.
- C4: aux batches interleaved 1:1 with DIS-GEN batches, weight 1.0, multi-hot over full type-restricted vocab.
- Early stop: val loss, patience 5, min_delta 3% relative, max 20 epochs. Test scored ONCE at best-val epoch. Seed 1337.
- Batch 50 (~92GB VRAM, owner-directed over the 91.2GB ceiling), 3,813 batches/epoch, **20.03h/epoch measured** (compute-bound, 100% GPU).

## 4. Hardware & environment (deviation on record)
- g7e.2xlarge: RTX PRO 6000 Blackwell 96GB, **sm_120**, 8 vCPU, 62GB RAM.
- Frozen cu118 env (torch 2.0.1 + torchdrug) **cannot run on sm_120** (no PTX). Owner authorized modern-torch divergence: `bpn120` env = torch 2.7.1+cu128 + torchdrug 0.2.1 ported (2 source patches). Recorded in `CHANGES/0026`; cedes Gate-0 env-pin reproducibility.

## 5. Results so far
Floors (degree/PPR/DWPC baselines on pinned 17,448, `completion_baselines_ana.md`):
- **AD floor 0.00909 (degree MRR)**, test 621 · **PD floor 0.00818 (degree MRR)**, test 339.

PD full-stack run (E2, in progress) per-epoch validation:
| epoch | train loss | val loss | val MRR | Hits@10 | note |
|---|---:|---:|---:|---:|---|
| 1 | 0.0036 | 0.000686 | **0.0195** | 0.0217 | new best |
| 2 | 0.0012 | 0.000681 | 0.0158 | **0.0372** | flat val (−0.6% < 3%) → patience 1/5 |
| 3 | 0.0012 | **0.000658** | 0.0168 | 0.0259 | −3.4% → new best, patience reset |

## 6. Observed problems (why we're asking for review)
1. **Calibration-vs-rank mismatch.** Val loss improves in tiny absolute steps (0.000686→0.000658) while MRR slides (0.0195→0.0158→0.0168) and H@10 whiplashes (0.022→0.037→0.026). Multi-hot BCE over 17,448 averages 17k independent yes/no calls dominated by already-known answers (AD positives = 68% of the candidate universe). The loss rewards calibration shifts that don't move rank; early-stop watches this noisy signal.
2. **Low absolute utility genome-wide.** ~2× floor but H@10 2–4% over 17,448 = 96%+ miss at top-10. No user ranks 17k genes; the product-shaped question is 50–100-candidate reranking (Probe 1, running).
3. **Weeks-scale wall clock.** 20h/epoch → 5–9 days per disease, 10–17 days for the wave. Only 1 seed affordable; E3 ablations would each cost the same.
4. **Anatomy edge concentration.** Only 51 of 1,324 kept CNS anatomy nodes carry ANA-GEN edges — the "anatomy helps" hypothesis rests on a thin subgraph (E3 ablation Arm C answers this, at full wave cost).
5. **Reproducibility debt.** Modern-torch divergence (torch 2.7.1 vs pinned 2.0.1) accepted for Blackwell; Gate-0 pins ceded. torchdrug-on-2.7 verified only ~50 batches + 3 epochs so far.
6. **>500k ANA-GEN tripwire** flagged (1.39M) — owner chose PROCEED; residual risk measured via Arm C, not argued.

## 7. Probes in flight (peer-review round 1)
- **Probe 1 rerank diagnostic — COMPLETE 2026-09-03 (`7d010db`, `reports/m2c/rerank_probe_pd.md`). VERDICT: FAILS the bar** (needed recall@10-in-100 ≥ 0.40 with CI clear of null). Findings:
  - Pre-registered degree-top-100 shortlist: recall@10-in-100 = 0.0059 = degree null **by construction** — the degree shortlist is 94% saturated with already-known train2/valid positives; coverage of held-out PD genes = **2/339 (0.59%)**.
  - Sensitivity panel (top-100 *excluding* known positives = the realistic "novel candidate" panel): coverage 8.85%; recall@10 = **0.0295, CI [0.0147–0.0472] = 3.3× null**, CI lower bound above the null point; covered positives rank **1–2 of 100**.
  - **Conclusion: the bottleneck is candidate GENERATION, not ranking.** The ranker works when the right gene is in front of it; degree shortlists never contain held-out genes.
  - Machinery proof: probe reproduces the trainer's genome-wide metrics exactly (ep2 0.0158/0.0372, ep3 0.0168/0.0259) — both code paths cross-verified.
  - Correction on record: `best.pth` currently holds epoch-3 weights (overwritten at epoch 3's new best), not epoch 1.
- **Probe 1.5 shortlist-generator comparison — COMPLETE (`b1db56b`, `reports/m2c/shortlist_generator_probe.md`). VERDICT: NO generator + current checkpoint reaches two-stage ≥ 0.10** (best 0.0295, CI [0.0147–0.0472], vs generator null 0.0088). Coverage of 339 held-out PD genes (excl-known): degree 8.9/15.0/27.1% at K=100/200/500; disease-PPR 4.7/11.8/25.7%; PPR+train2 6.2/11.8/24.5% — **no generator beats raw degree**; disease-personalized PPR concentrates on the known-positive neighborhood and covers held-out genes *worse*. The reranker is NOT the bottleneck: given coverage, it ranks hits top-10 in 33–48% of cases vs ~10% for generator ordering (3–5× lift). K=500 two-stage ceiling ≈ 0.09 even with a perfect reranker. **Candidate generation is the open problem.**
- **Probe 1.5 (old running line, superseded):** degree vs PPR-top-K (K=100/200/500) coverage of held-out PD genes on the realistic panel, then rerank within the best generator. The two-stage number that matters: P(gene in shortlist) × recall@10 given coverage; bar ≥ 0.10 with CI clear of null → `reports/m2c/shortlist_generator_probe.md`.
- **Probe 2 ranking-loss — QUEUED at PD epoch-4 boundary (~05:55Z Sep 4, owner-directed pause→probe→resume):** opt-in ranking-loss variant (log-softmax over 1+255 seeded negatives) vs current BCE on 200 identical batches; verdict whether rank separates beyond ±20% MRR noise → `reports/m2c/objective_probe.md`. Secondary goal: if scoring K+1≈256 instead of 17,448 holds up scientifically, epoch cost drops ~20h → ~1–2h, making seeds + E3 ablations affordable. No wave switch without CHANGES entry + owner signature.

## 7a. CPU-side nulls & coverage report (beast CPU agent, commit af53881 `reports/m2c/cpu_nulls_and_coverage.md` — quoted from owner relay; verify on pull)
1. **ANA feature is near-constant**: 99.1% of pinned genes carry it (median 41/48 edges); floors ±ANA move degree 0.00817→0.00831, PPR 0.00752→0.00656, CIs fully overlapping → **Arm C (E−ANA ablation) is dead — it would measure nothing**. The anatomy hypothesis is untestable as an ablation on this graph build.
2. **Supervision-only null** (train2 co-occurrence counts, NO graph, 1.2s CPU): PD test MRR **0.01484** vs the GNN's val 0.01684 — CIs overlap. The 2.8M-edge background graph adds little on this track so far; the model must beat a no-graph null to justify the GPU.
3. **Coverage at K=2000**: degree 0.525, PPR 0.596 — but K=2000 is 11% of the genome; not a shortlist.
4. ppr_train1 was already background-only (train1 holds zero DIS-GEN rows) — no leakage there.

## 7b. Updated strategic picture (post-Probe 1)
The two-stage product shape under test: **generator (PPR/shortlist) → GNN rerank**. Probe 1 said stage 2 works when stage 1 delivers; Probe 1.5 now shows **no graph-statistic generator (degree, disease-PPR, PPR+train2) delivers either** — all cap below the 0.10 two-stage bar at K≤200, and K=500 caps at ≈0.09 even with a perfect reranker. Options reviewers should weigh: (a) external generators (GWAS loci, Open Targets, literature panels) pinned in-repo; (b) larger panels (K=500–1000) with the GNN as reranker — the rerank lift (3–5×) is real and measured; (c) improve the genome-wide model itself (Probe 2 ranking loss) so the GNN becomes the generator; (d) accept the certified-negative on this stack. The measured facts: reranker adds 3–5× within panels; coverage, not ranking, caps the product.

## 8. Questions for reviewers
1. Is the multi-hot BCE + filtered-1-hot-val protocol sound for 68%-positive supervision, or is a listwise/ranking objective strictly better here? Any pitfalls in our per-row edge-removal scheme?
2. Degree-top-100 is now measured as coverage-dead for held-out genes (0.59%; Probe 1). PPR-based generation is being tested (Probe 1.5). Other constructions worth the pinning cost (Open Targets-derived, GWAS-window, pathway co-membership)? Any leakage traps in PPR-on-train1_ana we should check?
3. Given 20h/epoch and one GPU, what is the highest-information-per-hour experiment ordering? (ablation vs seed replication vs loss change vs shortlist eval)
4. The supervision-only co-count null (no graph) already reaches 0.01484 vs model 0.01684 (CIs overlap). Is ANY KG-model claim viable while this gap overlaps the null? What gap size / CI separation would you require before the graph is credited?
5. Any known SOTA benchmark we must match head-to-head before using comparative language (Open Targets Genetics, etc.)?
6. Explainability: NBFNet path-based explanations for ranked candidates — worth building into the MVP, or a distraction at this stage?
7. Red flags we haven't seen: label noise handling, disease-specific degree normalization, leakage vectors in train1_ana we may have missed.

## 9. What is NOT up for debate (frozen by owner/protocol)
Model architecture + hyperparams (b051a90 32×6), Adam lr 1e-3, pinned 17,448 universe, completion splits, claim ceiling, early-stop criterion (unless owner signs a CHANGES amendment).
