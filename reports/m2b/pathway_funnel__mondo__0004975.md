<!-- track_scope: benchmark -->
# Pathway Discovery Funnel — mondo:0004975 (CHANGES/0025)

Mode `open_world_keep`; checkpoint `d23475f36cd1`; seed 1337;
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

median drop 0.018862, 81.8% of 11 genes dropped — PASS

## Pathway benchmark

Open-world mode (`open_world_keep`): no valid/test split — scored by pathway support only, per config note.


## EXPLORATORY — ontology-anchored pathway scores (post-hoc, not the frozen endpoint)

Unconstrained top-grad paths bypass ontology nodes (frozen result above); this
table answers which pathways WOULD mediate if paths were forced through them.
The D->P prefix uses disease-gene relation conditioning as an approximation, and
the P—gene annotation edges are in train1 (circularity already disclosed below).
Prefix weights are a Viterbi best-path estimator (single best path, no beam, no
target edge-masking, no length averaging) rather than the gene-path beam
extraction (addendum 2026-08-25: exhaustive beam was computationally
infeasible, ~40 h/disease measured).

Candidates: 592 GO BP / Reactome nodes with a direct train1
annotation edge to a consensus Top-100 gene, prefiltered to
support >= 3 -> 19
(candidates below the strictest reported min_support sensitivity can never enter
any reported table, so the tables are identical to exhaustive); seeds [1337].

| rank | pathway | score | degree | support | top genes |
|---|---|---|---|---|---|
| 1 | go:0007165 (signal transduction) | 0.2545 | 1679 | 7 | ensembl:ENSG00000149295 (dopamine receptor D2), ensembl:ENSG00000204681 (gamma-aminobutyric acid type B receptor subunit 1), ensembl:ENSG00000136928 (gamma-aminobutyric acid type B receptor subunit 2) |
| 2 | go:0007186 (G protein-coupled receptor signaling pathway) | 0.1505 | 645 | 5 | ensembl:ENSG00000149295 (dopamine receptor D2), ensembl:ENSG00000204681 (gamma-aminobutyric acid type B receptor subunit 1), ensembl:ENSG00000136928 (gamma-aminobutyric acid type B receptor subunit 2) |
| 3 | go:0045944 (positive regulation of transcription by RNA polymerase II) | 0.1150 | 1133 | 5 | ensembl:ENSG00000082701 (glycogen synthase kinase 3 beta), ensembl:ENSG00000149295 (dopamine receptor D2), ensembl:ENSG00000097007 (ABL proto-oncogene 1, non-receptor tyrosine kinase) |
| 4 | go:0009410 (response to xenobiotic stimulus) | 0.0870 | 244 | 5 | ensembl:ENSG00000149295 (dopamine receptor D2), ensembl:ENSG00000142319 (solute carrier family 6 member 3), ensembl:ENSG00000097007 (ABL proto-oncogene 1, non-receptor tyrosine kinase) |

Open-world mode: no held-out truth — scored by pathway support only.

Gene annotation: annotated from /opt/sina/data/processed/m1/gene_properties.parquet.
Circularity disclosure (evaluation block): BPO-GEN/PWY-GEN edges are in train1
for message passing; the test is held out at disease level, not ontology level.
Raw GO BP and Reactome are reported alongside hallmarks, per the frozen maps.
