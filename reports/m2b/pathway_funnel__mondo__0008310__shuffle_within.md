<!-- track_scope: benchmark -->
# Pathway Discovery Funnel — mondo:0008310 (CHANGES/0025)

Mode `lodo_keep`; checkpoint `a24c10aa3ad8`; seed 1337;
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
| degree-null Recall@10 | 0.0069 |
| ORA Top-100 Recall@10 | 0.0027 |
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

Candidates: 2411 GO BP / Reactome nodes with a direct train1
annotation edge to a consensus Top-100 gene, prefiltered to
support >= 3 -> 402
(candidates below the strictest reported min_support sensitivity can never enter
any reported table, so the tables are identical to exhaustive); seeds [1337].

| rank | pathway | score | degree | support | top genes |
|---|---|---|---|---|---|
| 1 | go:0010628 (positive regulation of gene expression) | 0.0111 | 446 | 17 | ensembl:ENSG00000142192 (amyloid beta precursor protein), ensembl:ENSG00000169083 (androgen receptor), ensembl:ENSG00000141510 (tumor protein p53) |
| 2 | go:0007165 (signal transduction) | 0.0093 | 1679 | 46 | ensembl:ENSG00000100030 (mitogen-activated protein kinase 1), ensembl:ENSG00000120907 (adrenoceptor alpha 1A), ensembl:ENSG00000146648 (epidermal growth factor receptor) |
| 3 | go:0010629 (negative regulation of gene expression) | 0.0083 | 268 | 13 | ensembl:ENSG00000142192 (amyloid beta precursor protein), ensembl:ENSG00000141510 (tumor protein p53), ensembl:ENSG00000179348 (GATA binding protein 2) |
| 4 | go:0045944 (positive regulation of transcription by RNA polymerase II) | 0.0083 | 1133 | 32 | ensembl:ENSG00000142192 (amyloid beta precursor protein), ensembl:ENSG00000146648 (epidermal growth factor receptor), ensembl:ENSG00000100906 (NFKB inhibitor alpha) |
| 5 | reactome:R-HSA-8939211 (ESR-mediated signaling) | 0.0079 | 11 | 6 | ensembl:ENSG00000100030 (mitogen-activated protein kinase 1), ensembl:ENSG00000080824 (heat shock protein 90 alpha family class A member 1), ensembl:ENSG00000102882 (mitogen-activated protein kinase 3) |
| 6 | go:0008284 (positive regulation of cell population proliferation) | 0.0078 | 501 | 13 | ensembl:ENSG00000146648 (epidermal growth factor receptor), ensembl:ENSG00000169083 (androgen receptor), ensembl:ENSG00000178691 (SUZ12 polycomb repressive complex 2 subunit) |
| 7 | go:0071420 (cellular response to histamine) | 0.0074 | 6 | 5 | ensembl:ENSG00000145864 (gamma-aminobutyric acid type A receptor subunit beta2), ensembl:ENSG00000196639 (histamine receptor H1), ensembl:ENSG00000113327 (gamma-aminobutyric acid type A receptor subunit gamma2) |
| 8 | go:0000165 (MAPK cascade) | 0.0072 | 141 | 10 | ensembl:ENSG00000100030 (mitogen-activated protein kinase 1), ensembl:ENSG00000120907 (adrenoceptor alpha 1A), ensembl:ENSG00000146648 (epidermal growth factor receptor) |
| 9 | reactome:R-HSA-4090294 (SUMOylation of intracellular receptors) | 0.0072 | 30 | 9 | ensembl:ENSG00000169083 (androgen receptor), ensembl:ENSG00000188612 (small ubiquitin like modifier 2), ensembl:ENSG00000091831 (estrogen receptor 1) |
| 10 | go:1902895 (positive regulation of miRNA transcription) | 0.0068 | 56 | 10 | ensembl:ENSG00000146648 (epidermal growth factor receptor), ensembl:ENSG00000169083 (androgen receptor), ensembl:ENSG00000141510 (tumor protein p53) |

Anchored Recall@10 0.0062, Recall@5 0.0034 vs the same test-genes-only GO/Reactome ground truth as the frozen evaluation (exploratory). min_support sensitivity Recall@10: 3: 0.0055 5: 0.0062 10: 0.0055.

Gene annotation: annotated from /opt/sina/data/processed/m1/gene_properties.parquet.
Circularity disclosure (evaluation block): BPO-GEN/PWY-GEN edges are in train1
for message passing; the test is held out at disease level, not ontology level.
Raw GO BP and Reactome are reported alongside hallmarks, per the frozen maps.
