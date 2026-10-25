<!-- track_scope: benchmark -->
# AD failure diagnosis — why BioPathNet scored at null level on Task-G (owner-requested, 2026-08-01)

**Context.** AD seed-1337 gate: FAIL (filt MRR 0.0184 / Hits@10 0.0273 vs
degree baseline 0.0695/0.1447, null 0.0140/0.0297). tune1 (lr 1e-3,
batch_per_epoch 2000) tracked the same level (best valid MRR 0.0100). This
report identifies the root cause. Short version: **the exam was impossible by
construction — the held-out diseases are orphaned in the message-passing
graph.**

## Finding 1 (ROOT CAUSE): held-out diseases have ZERO edges in the BRG

Measured on the AD export (`trackd_family_alzheimer`):

- All 6 held-out AD diseases (`mondo:0004975`, `mondo:0007088`,
  `mondo:0011913`, `mondo:0014265`, `mondo:0015140`, `mondo:0100087`) have
  **degree 0 in train1** (the 1.94M-edge BRG). NBFNet's Bellman–Ford
  propagation starts at the query disease — with zero incident edges, **no
  message ever propagates**. The model can only fall back to entity-embedding
  popularity, which is exactly the measured signature: ≈ degree-matched null
  (0.018 vs 0.014), far below the degree baseline (0.0695) that exploits gene
  popularity directly.
- Why they are orphaned: `scripts/m2_export_taskg_biopathnet.py` (lines
  84–95) removes **every BRG edge touching a held-out family disease**
  (src or dst), implementing spec §7.4(2b) ("excludes only BRG edges
  referencing held-out family diseases"). The export faithfully implements
  the spec — **the spec itself amputates the model.**
- The edges exist in full Track-D and were stripped: the 6 diseases have
  **DRG-DIS 315, DIS-PHE 47, DIS-DIS 13, EXP-DIS 10** (385 non-DIS-GEN
  edges). The relation curation (`configs/m2_trackd_relations.yaml`, owner-
  sign-off-pending) explicitly lists DIS-PHE, DIS-DIS, DRG-DIS as
  message-passing context — the intent was there; the export's
  family-exclusion step removed them anyway.
- Contrast with BioPathNet's published PrimeKG setup: held-out diseases stay
  connected via disease–phenotype/disease–disease/drug–disease edges, so
  zero-shot reasoning is possible. Our export removed exactly those bridges.

## Finding 2 (CONTRIBUTING): 18.7% of test triples have BRG-orphan answer genes

2,201 of 11,766 test triples have answer genes with degree 0 in train1
(2,138 distinct genes). These genes only ever appear in DIS-GEN (supervision)
and/or the excluded ANA-GEN expression edges. `train2_in_factgraph` defaults
to True, so supervision edges connect them via train diseases — but they
have no biological context (PPI/pathway/phenotype) in the BRG, capping their
rankability. Report metrics should be stratified by gene reachability.

## Finding 3 (EXONERATED): evaluation protocol fairness

The engine's `heterogeneous_evaluation` masks candidates to the answer's
node type (genes), so its candidate space matches the degree baseline's
gene-only pool. The FAIL is not a candidate-space artifact.

## Finding 4 (CONTRIBUTING, already known): training protocol

lr 5e-3 overfit the 4,000-query/epoch subsample (valid MRR collapse
0.018→0.007); lr 1e-3 delayed but did not cure it (0.0095→0.0100→0.0067).
Even a perfectly tuned optimizer cannot fix Finding 1: there is no signal to
fit for orphaned query nodes.

## Recommendations (owner decision required — protocol deviations)

1. **CHANGES/0013 (draft) — reconnect held-out diseases in the BRG.** Include
   DIS-PHE and DIS-DIS edges (and hierarchies) touching held-out diseases in
   train1. Leakage risk: LOW (phenotype/hierarchy context, not gene
   associations). Optionally include DRG-DIS: risk MEDIUM-HIGH (indication
   edges are close proxies of disease-gene association; if included, this
   must be disclosed in every claim). Then rebuild exports (splits and
   supervision untouched) and rerun the AD seed — ~30 min rebuild + ~13 h run.
2. **Stratify metrics by candidate reachability** (BRG degree 0 vs >0) in
   `scripts/m2_compute_metrics.py` so the 18.7% orphan-gene triples are
   visible separately.
3. If DRG-DIS is deemed too leaky, run Finding-1-fixed with DIS-PHE+DIS-DIS
   only; if metrics still fail, the honest conclusion is that zero-shot
   family holdout on this graph needs a different model class (e.g.
   KGPFLN/ICL arm, already in the roadmap) rather than more NBFNet tuning.

## What NOT to conclude

- Not "BioPathNet doesn't work" — it was never allowed to see the query
  diseases in the graph.
- Not "the GPU/lr/batch is the problem" — those were real but secondary
  (and fixed: batch 4 max verified, lr sweep done, both engineering patches
  committed).
