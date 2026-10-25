# CLOSED — permanently parked lines of investigation (do not re-spend)

Written 2026-10-07 (T5 of the 2026-10-06 audit handoff). Each entry: what was
hypothesized, the verdict, the evidence pointer. If you are tempted to revive
one of these, read the evidence file first and state in a `CHANGES/` entry what
new information invalidates the listed evidence.

## 1. ANA / CNS anatomy edges as training signal — CLOSED
Hypothesis: adding Bgee anatomy-expression edges (ANA-GEN) to the message-passing
graph improves disease–gene recovery.
Evidence: `reports/m2c/noana_rank_epoch.md` — three inertness proofs; the feature
is near-constant, ablation measures nothing, test scores equal with and without
ANA. Verdict: measures nothing; do not re-add.

## 2. Pathway Discovery Funnel (as frozen) — CLOSED
Hypothesis: pooling gradient-weighted paths into pathway rankings recovers
mechanism-level signal that gene-level MRR hides.
Evidence: `reports/m2b/pathway_funnel.md` — frozen funnel Recall@10 0.0 vs
degree-null 0.0069 / ORA 0.0027 (unconstrained beam bypasses ontology nodes).
`reports/m2b/anchored_topology_null.md` — the exploratory anchored variant ties
the topology null (Jaccard 0.82/0.54) and survives the within-disease shuffle
while gene MRR collapses 58–89%: the anchored table is a topology artifact.
Verdict: NEGATIVE. Any revival requires a constrained beam through BPO/PWY nodes
as a new pre-registered protocol, not a tweak.

## 3. Two-stage GNN rerank over external panels — CLOSED
Hypothesis: model reranking inside a candidate panel beats the panel null.
Evidence: `reports/m2c/completion_wave.md` + `reports/m2c/rerank_probe_pd.md` —
conversion pinned at 0.0295 (= the OT null) across 10+ panels and all K;
the hit set is panel-independent. Generation (coverage) is fixable, conversion
is not. Verdict: NEGATIVE.

## 4. Union / RRF ensembling of candidate generators — CLOSED
Evidence: `reports/m2c/completion_wave.md` (union-ranker + Probe-1 closure,
2026-09-05): union upper bound 0.0413–0.0442 but RRF fusion = 0.0295 = single
models = OT null. Verdict: dead without calibrated fusion; no fusion work is
justified by anything measured here.

## 5. Full-scale KGE (PyKEEN-style embedding models) on this graph — CLOSED
Evidence: `reports/m2a/kge_pilot_1pct_report.md` + `reports/m2a/final_report_m2a.md`
— ~2,600 CPU-hours measured for the full graph; transductive KGE also cannot
score unseen heads (cold exam) by construction. Verdict: cost-prohibitive and
wrong tool for the zero-shot exam.

## 6. GPU GNN training generally (this supervision) — CLOSED
Evidence: `reports/m2c/completion_wave.md` — the plateau is invariant to
±ANA edges, BCE vs rank loss, width 32/64, batch 32/100, and disease (PD/AD).
The ceiling is the supervision, not the model. Verdict: more GPU measures the
same plateau louder. Any new GPU training needs a *supervision change* (new
labels, new task), not a knob change.

## 7. Reweighted-GPU epoch (label-reweighting) — CLOSED (HOLD → closed)
Evidence: `reports/m2c/reweighted_null.md` — CPU gate HOLD both diseases;
PD LITONLY alone = FULL (channels redundant); AD byte-identical H@10 on all
6 variants. The licensed GPU epoch was never launched; nothing in the gate
justifies it.

## 8. Per-locus "rank the knowns" product hypothesis — CLOSED 2026-09-25
Hypothesis: per-locus known-gene recovery (top-1 0.904 PD / 0.954 AD) justifies
workflow software.
Evidence: `reports/m2c/loci_headtohead_mondo__*_testknown.{json,md}` — those
knowns are train positives. Graded on known genes never in training: PD 0.300
vs OT 0.600 (LOSS, n=10), AD 0.1613 vs OT 0.2258 (TIE, n=31), both below
chance-draw means. No unseen-gene signal. Verdict: NEGATIVE.

## 9. Functional-evidence tabular tracks vs Open Targets L2G — CLOSED 2026-09-26
Hypothesis: adding public functional tracks (fine-mapping, ClinVar, regulatory
elements, eQTL) to a tabular model closes the gap to OT L2G.
Evidence: `m5b/reports/verdict.md` + `v5_results.json` — best arm (all five
tracks) 0.300 vs OT 0.477 covered PR AUC; every CI excludes zero on the
negative side. Only ClinVar moved (+0.024). Verdict: LOSS on all arms; the gap
is OT's evidence coverage, not modelable from these features.

## 10. Colocalization engine (M5B Track via disease stats) — CLOSED (never built)
Evidence: `m5b/reports/verdict.md` leakage statement — disease-stats gate failed
at 33.7%; correctly demoted, not faked. Verdict: gate failed; do not build
without new data.
