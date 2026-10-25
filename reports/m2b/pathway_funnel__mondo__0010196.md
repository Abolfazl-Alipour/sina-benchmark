<!-- track_scope: benchmark -->
# Pathway Discovery Funnel — mondo:0010196 (CHANGES/0025)

Mode `lodo_keep`; checkpoint `92ec77b8e372`; seed 1337;
K = 100 (sensitivity [50, 100, 200]); scoring = gradient_weighted
(grad of sigmoid(score) w.r.t. path edge weights x PNA aggregation, BioKGC/BioPathNet prior art).

## Top-10 predicted pathways (hub-penalized: score / ln(degree), min_support 5)

| rank | pathway | score | degree | support | top genes |
|---|---|---|---|---|---|
| 1 | go:0030001 (metal ion transport) | 0.1271 | 30 | 7 | ensembl:ENSG00000110680 (calcitonin related polypeptide alpha), ensembl:ENSG00000126088 (uroporphyrinogen decarboxylase), ensembl:ENSG00000141510 (tumor protein p53) |

## Hallmarks of Aging (reported separately, weight <= 0.15, never collapsed)

| id | name | score | member nodes hit | genes |
|---|---|---|---|---|
| — | no hallmark member nodes hit | | | |

## Faithfulness control (mask Top-10 paths per gene, re-score)

median drop 0.123487, 100.0% of 100 genes dropped — PASS

## Pathway benchmark

| metric | value |
|---|---|
| true pathways (>= 3 test genes) | 2085 |
| pathway Recall@10 | 0.0 |
| pathway Recall@5 | 0.0 |
| pathway MRR | 0.0 |
| pathway Precision@3 | 0.0 |
| pathway NDCG@10 | 0.0 |
| degree-null Recall@10 | 0.0048 |
| ORA Top-100 Recall@10 | 0.0043 |
| >2-hop strata Recall@10 | 0.0 (118 genes, 99 true) |

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

Candidates: 2578 GO BP / Reactome nodes with a direct train1
annotation edge to a consensus Top-100 gene, prefiltered to
support >= 3 -> 384
(candidates below the strictest reported min_support sensitivity can never enter
any reported table, so the tables are identical to exhaustive); seeds [1337, 1339, 5027].

| rank | pathway | score | degree | support | top genes |
|---|---|---|---|---|---|
| 1 | go:0010628 (positive regulation of gene expression) | 0.0262 | 446 | 26 | ensembl:ENSG00000049768 (forkhead box P3), ensembl:ENSG00000232810 (tumor necrosis factor), ensembl:ENSG00000105329 (transforming growth factor beta 1) |
| 2 | go:0007165 (signal transduction) | 0.0250 | 1679 | 23 | ensembl:ENSG00000109320 (nuclear factor kappa B subunit 1), ensembl:ENSG00000112164 (glucagon like peptide 1 receptor), ensembl:ENSG00000010610 (CD4 molecule) |
| 3 | go:0045944 (positive regulation of transcription by RNA polymerase II) | 0.0215 | 1133 | 23 | ensembl:ENSG00000049768 (forkhead box P3), ensembl:ENSG00000232810 (tumor necrosis factor), ensembl:ENSG00000105329 (transforming growth factor beta 1) |
| 4 | go:0008285 (negative regulation of cell population proliferation) | 0.0202 | 395 | 17 | ensembl:ENSG00000049768 (forkhead box P3), ensembl:ENSG00000232810 (tumor necrosis factor), ensembl:ENSG00000105329 (transforming growth factor beta 1) |
| 5 | go:0010629 (negative regulation of gene expression) | 0.0175 | 268 | 18 | ensembl:ENSG00000049768 (forkhead box P3), ensembl:ENSG00000232810 (tumor necrosis factor), ensembl:ENSG00000105329 (transforming growth factor beta 1) |
| 6 | go:0000122 (negative regulation of transcription by RNA polymerase II) | 0.0173 | 902 | 18 | ensembl:ENSG00000049768 (forkhead box P3), ensembl:ENSG00000232810 (tumor necrosis factor), ensembl:ENSG00000105329 (transforming growth factor beta 1) |
| 7 | go:0008284 (positive regulation of cell population proliferation) | 0.0172 | 501 | 14 | ensembl:ENSG00000105329 (transforming growth factor beta 1), ensembl:ENSG00000169083 (androgen receptor), ensembl:ENSG00000138685 (fibroblast growth factor 2) |
| 8 | go:0060736 (prostate gland growth) | 0.0169 | 8 | 5 | ensembl:ENSG00000169083 (androgen receptor), ensembl:ENSG00000017427 (insulin like growth factor 1), ensembl:ENSG00000137869 (cytochrome P450 family 19 subfamily A member 1) |
| 9 | go:0045893 (positive regulation of DNA-templated transcription) | 0.0165 | 730 | 17 | ensembl:ENSG00000049768 (forkhead box P3), ensembl:ENSG00000232810 (tumor necrosis factor), ensembl:ENSG00000105329 (transforming growth factor beta 1) |
| 10 | go:0043410 (positive regulation of MAPK cascade) | 0.0164 | 214 | 13 | ensembl:ENSG00000232810 (tumor necrosis factor), ensembl:ENSG00000105329 (transforming growth factor beta 1), ensembl:ENSG00000010610 (CD4 molecule) |

Anchored Recall@10 0.0043, Recall@5 0.0024 vs the same test-genes-only GO/Reactome ground truth as the frozen evaluation (exploratory). min_support sensitivity Recall@10: 3: 0.0043 5: 0.0043 10: 0.0048.

Gene annotation: annotated from /opt/sina/data/processed/m1/gene_properties.parquet.
Circularity disclosure (evaluation block): BPO-GEN/PWY-GEN edges are in train1
for message passing; the test is held out at disease level, not ontology level.
Raw GO BP and Reactome are reported alongside hallmarks, per the frozen maps.
