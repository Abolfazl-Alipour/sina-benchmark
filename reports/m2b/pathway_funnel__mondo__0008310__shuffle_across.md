<!-- track_scope: benchmark -->
# Pathway Discovery Funnel — mondo:0008310 (CHANGES/0025)

Mode `lodo_keep`; checkpoint `797e40d59a29`; seed 1337;
K = 100 (sensitivity [50, 100, 200]); scoring = gradient_weighted
(grad of sigmoid(score) w.r.t. path edge weights x PNA aggregation, BioKGC/BioPathNet prior art).

## Top-10 predicted pathways (hub-penalized: score / ln(degree), min_support 5)

| rank | pathway | score | degree | support | top genes |
|---|---|---|---|---|---|


## Hallmarks of Aging (reported separately, weight <= 0.15, never collapsed)

| id | name | score | member nodes hit | genes |
|---|---|---|---|---|
| — | no hallmark member nodes hit | | | |

## Faithfulness control (mask Top-10 paths per gene, re-score)

VACUOUS — no paths to mask

## Pathway benchmark

| metric | value |
|---|---|
| true pathways (>= 3 test genes) | 1459 |
| pathway Recall@10 | 0.0 |
| pathway Recall@5 | 0.0 |
| pathway MRR | 0.0 |
| pathway Precision@3 | 0.0 |
| pathway NDCG@10 | 0.0 |
| degree-null Recall@10 | 0.0 |
| ORA Top-100 Recall@10 | 0.0 |
| >2-hop strata Recall@10 | 0.0 (142 genes, 258 true) |

K sensitivity (pre-registered Recall@10): K=50: 0.0 K=100: 0.0 K=200: 0.0


## EXPLORATORY — ontology-anchored pathway scores (post-hoc, not the frozen endpoint)

Unconstrained top-grad paths bypass ontology nodes (frozen result above); this
table answers which pathways WOULD mediate if paths were forced through them.
The D->P prefix uses disease-gene relation conditioning as an approximation, and
the P—gene annotation edges are in train1 (circularity already disclosed below).
Prefix weights are a Viterbi best-path estimator (single best path, no beam, no
target edge-masking, no length averaging) rather than the gene-path beam
extraction (addendum 2026-08-25: exhaustive beam was computationally
infeasible, ~40 h/disease measured).

Candidates: 0 GO BP / Reactome nodes with a direct train1
annotation edge to a consensus Top-100 gene, prefiltered to
support >= 3 -> 0
(candidates below the strictest reported min_support sensitivity can never enter
any reported table, so the tables are identical to exhaustive); seeds [1337].

| rank | pathway | score | degree | support | top genes |
|---|---|---|---|---|---|
| — | no anchored pathway passed min_support | | | |

Anchored Recall@10 0.0, Recall@5 0.0 vs the same test-genes-only GO/Reactome ground truth as the frozen evaluation (exploratory). min_support sensitivity Recall@10: 3: 0.0 5: 0.0 10: 0.0.

Gene annotation: annotated from /opt/sina/data/processed/m1/gene_properties.parquet.
Circularity disclosure (evaluation block): BPO-GEN/PWY-GEN edges are in train1
for message passing; the test is held out at disease level, not ontology level.
Raw GO BP and Reactome are reported alongside hallmarks, per the frozen maps.
