<!-- track_scope: benchmark -->
# Pathway Discovery Funnel — mondo:0008310 (CHANGES/0025)

Mode `lodo_keep`; checkpoint `c1041997a6c9`; seed 1337;
K = 100 (sensitivity [50, 100, 200]); scoring = gradient_weighted
(grad of sigmoid(score) w.r.t. path edge weights x PNA aggregation, BioKGC/BioPathNet prior art).
Positive control: LMNA consensus rank 25 (per-seed: seed 1337: rank 16, seed 1339: rank 57, seed 5027: rank 102) — gate = consensus Top-100 (owner amendment 2026-08-24).

## Top-10 predicted pathways (hub-penalized: score / ln(degree), min_support 5)

| rank | pathway | score | degree | support | top genes |
|---|---|---|---|---|---|


## Hallmarks of Aging (reported separately, weight <= 0.15, never collapsed)

| id | name | score | member nodes hit | genes |
|---|---|---|---|---|
| — | no hallmark member nodes hit | | | |

## Faithfulness control (mask Top-10 paths per gene, re-score)

median drop 0.094528, 98.0% of 100 genes dropped — PASS

## Pathway benchmark

| metric | value |
|---|---|
| true pathways (>= 3 test genes) | 1459 |
| pathway Recall@10 | 0.0 |
| pathway Recall@5 | 0.0 |
| pathway MRR | 0.0 |
| pathway Precision@3 | 0.0 |
| pathway NDCG@10 | 0.0 |
| degree-null Recall@10 | 0.0069 |
| ORA Top-100 Recall@10 | 0.0069 |
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

Candidates: 2220 GO BP / Reactome nodes with a direct train1
annotation edge to a consensus Top-100 gene, prefiltered to
support >= 3 -> 301
(candidates below the strictest reported min_support sensitivity can never enter
any reported table, so the tables are identical to exhaustive); seeds [1337, 1339, 5027].

| rank | pathway | score | degree | support | top genes |
|---|---|---|---|---|---|
| 1 | go:0006629 (lipid metabolic process) | 0.0677 | 788 | 21 | ensembl:ENSG00000125730 (complement C3), ensembl:ENSG00000051382 (phosphatidylinositol-4,5-bisphosphate 3-kinase catalytic subunit beta), ensembl:ENSG00000130164 (low density lipoprotein receptor) |
| 2 | go:0010628 (positive regulation of gene expression) | 0.0472 | 446 | 22 | ensembl:ENSG00000232810 (tumor necrosis factor), ensembl:ENSG00000164867 (nitric oxide synthase 3), ensembl:ENSG00000105329 (transforming growth factor beta 1) |
| 3 | go:0008284 (positive regulation of cell population proliferation) | 0.0387 | 501 | 12 | ensembl:ENSG00000112715 (vascular endothelial growth factor A), ensembl:ENSG00000105329 (transforming growth factor beta 1), ensembl:ENSG00000017427 (insulin like growth factor 1) |
| 4 | go:0008202 (steroid metabolic process) | 0.0344 | 129 | 10 | ensembl:ENSG00000130164 (low density lipoprotein receptor), ensembl:ENSG00000160752 (farnesyl diphosphate synthase), ensembl:ENSG00000147852 (very low density lipoprotein receptor) |
| 5 | go:0045944 (positive regulation of transcription by RNA polymerase II) | 0.0342 | 1133 | 20 | ensembl:ENSG00000232810 (tumor necrosis factor), ensembl:ENSG00000112715 (vascular endothelial growth factor A), ensembl:ENSG00000105329 (transforming growth factor beta 1) |
| 6 | go:0008203 (cholesterol metabolic process) | 0.0341 | 104 | 10 | ensembl:ENSG00000130164 (low density lipoprotein receptor), ensembl:ENSG00000160752 (farnesyl diphosphate synthase), ensembl:ENSG00000147852 (very low density lipoprotein receptor) |
| 7 | go:0045893 (positive regulation of DNA-templated transcription) | 0.0327 | 730 | 16 | ensembl:ENSG00000232810 (tumor necrosis factor), ensembl:ENSG00000105329 (transforming growth factor beta 1), ensembl:ENSG00000109320 (nuclear factor kappa B subunit 1) |
| 8 | go:0007165 (signal transduction) | 0.0303 | 1679 | 19 | ensembl:ENSG00000100030 (mitogen-activated protein kinase 1), ensembl:ENSG00000109320 (nuclear factor kappa B subunit 1), ensembl:ENSG00000125730 (complement C3) |
| 9 | go:0006955 (immune response) | 0.0281 | 450 | 15 | ensembl:ENSG00000232810 (tumor necrosis factor), ensembl:ENSG00000125730 (complement C3), ensembl:ENSG00000010610 (CD4 molecule) |
| 10 | go:0010629 (negative regulation of gene expression) | 0.0271 | 268 | 14 | ensembl:ENSG00000232810 (tumor necrosis factor), ensembl:ENSG00000112715 (vascular endothelial growth factor A), ensembl:ENSG00000105329 (transforming growth factor beta 1) |

Anchored Recall@10 0.0069, Recall@5 0.0034 vs the same test-genes-only GO/Reactome ground truth as the frozen evaluation (exploratory). min_support sensitivity Recall@10: 3: 0.0069 5: 0.0069 10: 0.0069.

Gene annotation: annotated from /opt/sina/data/processed/m1/gene_properties.parquet.
Circularity disclosure (evaluation block): BPO-GEN/PWY-GEN edges are in train1
for message passing; the test is held out at disease level, not ontology level.
Raw GO BP and Reactome are reported alongside hallmarks, per the frozen maps.
