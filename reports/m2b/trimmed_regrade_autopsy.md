<!-- track_scope: benchmark -->
# Trimmed-run re-grade + exam autopsy (owner-requested, 2026-08-07)

**Why this report exists.** The trimmed long run's verdict
(`trimmed_ad_learning_run.md`: "gate FAILED … model at degree-matched-null
level … memorization without transfer") compared the model's scores against
baseline and null numbers that were computed on a **different graph** — the
full Task-G export (21.8M-edge supervision pool), not the trimmed export
(170,574-edge supervision pool) the model actually trained and was tested on.
This report re-grades the run against baselines computed on the trimmed
export itself, audits the composition of the exam, and measures how much
recoverable signal the graph contains for each held-out disease. No new
training was run; everything below is CPU analysis of committed data
(`data/processed/m2/trimmed/biopathnet/trackd_family_alzheimer/`) with
committed code (`sina.baselines`, `sina.metrics`).

## 1. Protocol replication (validation of the method used here)

Re-computing the pre-registered degree baseline and degree-matched null on
the **full** split (`data/processed/m2/splits/trackd_family_alzheimer/`)
with `sina.baselines` + `sina.metrics` reproduces the official Gate-2
numbers exactly (degree: filt MRR 0.0695 / Hits@10 0.1447; null: 0.0147 /
0.0306 vs official 0.0140 / 0.0297, sampling jitter only). The replication
method is therefore sound.

## 2. The corrected comparison table (the central finding)

Same protocol, same code, **trimmed export's own supervision pool**:

| scorer | protocol | filt MRR | filt Hits@10 |
|---|---|---|---|
| degree baseline, FULL export (the pre-registered gate) | official (gene-direction, 61,306-gene pool) | 0.0695 | 0.1447 |
| degree-matched null, FULL export | official | 0.0140 | 0.0297 |
| **degree baseline, TRIMMED export** | **official (same as gate)** | **0.0094** | **0.0094** |
| **degree-matched null, TRIMMED export** | **official (same as gate)** | **0.0020** | **0.0017** |
| popularity baseline, TRIMMED export | engine protocol (both directions, pessimistic ties, export candidate pools) | 0.0051 | 0.0046 |
| **model (epoch 1), TRIMMED export** | **engine protocol** | **0.0160** | **0.0137** |

Two separate corrections fall out of this table:

1. **The gate was measured on a different graph.** Trimming the supervision
   pool from 21.8M to 170,574 edges removes ~99% of the popularity
   information the degree baseline exploits, so the baseline itself falls
   from 0.0695 to 0.0094 and the null from 0.0140 to 0.0020. Judging the
   trimmed-graph model against 0.0695 / 0.1447 graded it against numbers no
   scorer on that graph could reach. The pre-registered verdict (FAIL vs
   gate) formally stands — the gate was pre-registered — but it does not
   mean what it was read as.
2. **"Model ≈ null, learned nothing" is retracted.** Under the engine's own
   metric protocol, pure popularity scores 0.0051 / 0.0046 on the trimmed
   export; the model scored 0.0160 / 0.0137 — **3.1× / 3.0× above the
   popularity floor**, and ~8× above the trimmed degree-matched null
   (0.0020 / 0.0017). The model did learn transferable signal.

What does NOT change: the validation trajectory (0.0162 → 0.0157 → 0.0119
→ 0.0109 across epochs 1–4) was measured under one fixed protocol, so the
decline is real — after the first pass, additional training degraded
transfer (over-fitting to the training diseases). The best model is the
least-trained one.

## 3. Exam-composition autopsy (what the test actually asks)

The held-out AD family has 6 diseases; the test set has 11,766 triples.

| held-out disease | test answers | share of test |
|---|---|---|
| mondo:0004975 Alzheimer disease | 11,139 | **94.7%** |
| mondo:0100087 familial Alzheimer disease | 368 | 3.1% |
| mondo:0015140 early-onset autosomal dominant Alzheimer disease | 142 | 1.2% |
| mondo:0011913 Alzheimer disease 3 | 69 | 0.6% |
| mondo:0007088 Alzheimer disease type 1 | 29 | 0.2% |
| mondo:0014265 Alzheimer disease 18 | 19 | 0.2% |

- **The exam is one disease.** 94.7% of the grade is generic "Alzheimer
  disease", whose answer key lists **11,139 genes = 64% of all 17,302 genes
  in the graph**. This key is bulk literature-mined association data; the
  five familial/early-onset sub-diseases (the genetically specific keys)
  together carry 5.3% of the weight.
- Connectivity is no longer the problem it was in the orphaned export: in
  the trimmed fact graph the held-out diseases have degree 336 / 23 / 10 /
  4 / 4 / 3, and only 3.0% of test answers are true orphans (no edges of
  any kind); 7.6% have no background-graph edges.
- 88% of answer genes also appear as answers for training diseases, so
  shared-gene paths exist in principle.

## 4. Signal census (is the answer in the graph at all?)

Breadth-first distance from each held-out disease to genes in the fact
graph (train1 + train2, undirected):

- For the **familial/early-onset diseases**, 100% of answers sit within
  2–3 hops while only ~67–85% of non-answers do — answers are
  systematically closer. Signal exists.
- For **generic Alzheimer disease (94.7% of the exam)**, answers and
  non-answers have nearly identical distance profiles (87.9% vs 96.1% at
  3 hops). Topology barely separates a 64%-of-the-genome key.

Personalized PageRank from each held-out disease (random-walk flow from
the disease through the fact graph, genes ranked), official
`sina.metrics` protocol:

| disease | PPR MRR | PPR Hits@10 |
|---|---|---|
| mondo:0011913 (69 answers) | **0.0340** | **0.0870** |
| mondo:0007088 (29) | 0.0176 | 0.0690 |
| mondo:0014265 (19) | 0.0215 | 0.0526 |
| mondo:0015140 (142) | 0.0100 | 0.0141 |
| mondo:0100087 (368) | 0.0041 | 0.0054 |
| mondo:0004975 (11,139) | 0.0048 | 0.0052 |
| **all (edge-weighted)** | **0.0051** | **0.0060** |

A zero-training topological probe already beats the trimmed degree
baseline per-disease on three of the five specific sub-diseases — and is
flat on the mega-key. Conclusion: **recoverable signal is concentrated in
the specific familial/early-onset diseases, and the aggregate metric is
dominated by the one disease where topology carries almost none.**

## 5. Corrected interpretation

1. The model learned real, transferable signal on the trimmed graph
   (≈3× the like-for-like popularity floor), most of it within the first
   pass; further passes over-fit and erode transfer.
2. The absolute scores look tiny chiefly because 94.7% of the exam is one
   disease whose answer key covers 64% of the genome with almost no
   topological footprint. No message-passing model can score high on that
   key from this graph; the gate number (0.0695) reflected the full
   graph's popularity crutch, not learnable biology.
3. The per-disease model breakdown (does the model score well on the
   familial diseases where the signal is?) requires the epoch-1 checkpoint
   (preserved on the owner-synced Drive folder; 31 MB) — pending.
4. The memorization-vs-underfit question is likewise checkpoint-decidable
   (score the epoch-4 checkpoint on its own training triples).

## 6. Options for the redesign (owner decision required)

- **O1 — Re-grade what we already have** (free, CPU): checkpoint autopsy —
  per-disease + per-direction metrics, train-set memorization check.
  Sharpens the diagnosis; changes nothing by itself.
- **O2 — Fix the exam's grading**: pre-registered stratification — report
  per-disease (macro) metrics and evidence-stratified subsets (curated
  vs bulk literature-mined answers) alongside the aggregate. Precedent
  for stratified reporting is already in the research memo. No retraining.
- **O3 — Fix the exam itself**: replace edge-weighted aggregation with
  disease-macro averaging, and/or evaluate on the specific familial keys
  where signal exists; or move to candidate-panel ranking (rank a fixed
  panel of ~200–500 candidates per disease), which matches the triage
  use-case. Requires a CHANGES entry; re-scoring only, no retraining.
- **O4 — Short-training protocol**: if O1 confirms epoch-1 is the peak,
  the correct training regime is 1–2 passes, not 20 — ~70–140 GPU-min per
  run, making multi-seed cheap. Requires a CHANGES entry.
- **O5 — Enrich held-out disease context**: phenotypes are currently dead
  ends (no phenotype→gene edges in the export); adding phenotype→gene and
  disease→disease similarity edges would give the specific diseases more
  than 3–23 bridges. Leakage review needed (CHANGES entry).

**Recommendation:** O1 now (needs the checkpoint from the Drive folder),
then O2+O3 as one exam-redesign package before any further GPU spend; O4
follows from O1's answer. Do not resume multi-seed training on the current
exam: it would measure the mega-key's noise floor five times.

*Correction noted in `trimmed_ad_learning_run.md` (banner at top). The
CHANGES/0016 DRG-DIS leakage disclosure still applies to every number in
this report; it biases scores upward, not downward.*
