<!-- track_scope: benchmark -->
# LODO freeze-prep audits (v5.1 Step 2 — measurement only, nothing graded)

CPU audits feeding the `CHANGES/0020` freeze package. Claim ceiling: all
outputs here support **disease-associated gene recovery / hypothesis
prioritization** within the neuro/aging/metabolic cluster only.

**Inputs (re-measured, not asserted):** 191,953 DIS-GEN
supervision edges over 662 diseases; 1,431,945
BRG edges; cluster = 678 diseases; candidate gene universe
= **17,448** genes; candidate units enumerated = 48
(A 18,
B 6,
C 11,
D 13).

## Methods summary

- **Overlap audit (v5.1 §5.3):** per cluster disease, its supervision gene set;
  pairwise overlap coefficient |A∩B|/min(|A|,|B|) (primary — containment;
  Jaccard alongside). Full matrix: `overlap_matrix.parquet`
  (218,791 pairs). Bands: ≥0.3 co-hold-out,
  0.15-0.3 disclose.
- **BFS separation (v5.1 §5.4 θ input):** undirected BFS (depth cap 5)
  from each cluster disease on the LODO fact graph (BRG + supervision *minus the
  unit's own DIS-GEN edges*); % of answer vs non-answer genes reachable ≤2/≤3
  hops; separation = difference. Computed for all 662 supervised
  cluster diseases: `bfs_separation.parquet`.
- **Connectivity census + attrition (v5.1 §4, §6.1):** BRG context edges per
  cluster disease; floors: ≥10 context edges, ≥50
  supervision edges (primary) / 20-49
  (secondary panel). `connectivity_attrition.parquet`.
- **Key-size ceiling (v5.1 §5.4):** key size / 17,448 as %, bands
  <2% / 2-5% (θ working range) / >5%. `key_size_ceiling.parquet`.
- **Coherence control (v5.1 §2.2):** intra-cluster DIS-DIS density (edges among
  cluster diseases / possible pairs) and mean supervision degree vs
  100 seeded random same-size Track-D disease subsets
  (seed 1337). `coherence.json`.
- **PPR ± DIS-PHE (v5.1 §4):** per candidate unit, `sina.baselines.ppr_baseline`
  on BRG + training DIS-GEN (supervision minus the unit's own edges), restart
  0.5 pinned, sensitivity [0.15, 0.85], symmetric
  normalization; scored with `sina.metrics.link_prediction_metrics` over the
  17,448-gene universe, filtered by the unit's own answers.
  `ppr_disphe.parquet`.
- **Noise-floor simulation (v5.1 §6.1):** 1000 seeded
  resamples per key size of a random ranking of 17,448 genes
  (rank-one-answer protocol); distribution of filtered Hits@10.
  `noise_floor.parquet`.
- **EFO/OBA pinning (v5.1 §2.1):** per roster term — Track-D node, label,
  MONDO mapping, resolution method, bounded near-synonym note.
  `efo_pinning.parquet`.

## 1. Overlap audit — risk zones

- Pairs ≥ 0.3 (co-hold-out): **47,955** of 218,791
- Pairs 0.15-0.3 (disclose): **22,182**
- Most ≥0.3 pairs are containment artifacts of above-ceiling
  mega keys (a small key nested in a 30-70%-of-genome key scores 1.0 by
  construction; those mega keys are never graded). Restricted to pairs where
  **both** keys are below the >5% ceiling band — the co-holdout-relevant set:
  **30,975** co-hold-out pairs and **19,789** disclose pairs.
  Dual artifact, also by containment: pairs whose smaller key is tiny (a
  1-10-gene key sharing one gene scores up to 1.0). The matrix carries
  `size_a`/`size_b`/`intersection` so the freeze can apply a minimum-key-size
  or minimum-intersection qualifier when reading co-holdout pairs.

Top co-hold-out-band pairs:

| disease A | disease B | overlap coeff | Jaccard |
|---|---|---|---|
| efo:0005422 (skin aging) | efo:0022597 (aging) | 1.000 | 0.171 |
| efo:0005422 (skin aging) | mondo:0017948 (ABetaA21G amyloidosis) | 1.000 | 0.003 |
| efo:0005422 (skin aging) | mondo:0017949 (ABeta amyloidosis, Arctic type) | 1.000 | 0.003 |
| efo:0005678 (hippocampal sclerosis of aging) | mondo:0000167 (Huntington disease and related disorders) | 1.000 | 0.001 |
| efo:0005678 (hippocampal sclerosis of aging) | mondo:0004975 (Alzheimer disease) | 1.000 | 0.000 |
| efo:0005678 (hippocampal sclerosis of aging) | mondo:0005148 (type 2 diabetes mellitus) | 1.000 | 0.000 |
| efo:0005678 (hippocampal sclerosis of aging) | mondo:0005180 (Parkinson disease) | 1.000 | 0.000 |
| efo:0005678 (hippocampal sclerosis of aging) | mondo:0007739 (Huntington disease) | 1.000 | 0.001 |
| efo:0005678 (hippocampal sclerosis of aging) | mondo:0018899 (posterior cortical atrophy) | 1.000 | 0.000 |
| efo:0005678 (hippocampal sclerosis of aging) | mondo:0021095 (parkinsonian disorder) | 1.000 | 0.000 |
| efo:0005678 (hippocampal sclerosis of aging) | mondo:0024237 (inherited neurodegenerative disorder) | 1.000 | 0.000 |
| efo:0009885 (frailty measurement) | mondo:0013645 (autosomal recessive spinocerebellar ataxia 11) | 1.000 | 0.091 |
| efo:0022597 (aging) | mondo:0007247 (basal ganglia calcification, idiopathic, childhood-onset) | 1.000 | 0.000 |
| efo:0022597 (aging) | mondo:0009839 (progressive supranuclear palsy-parkinsonism syndrome) | 1.000 | 0.000 |
| efo:0022597 (aging) | mondo:0013179 (hereditary spastic paraplegia 44) | 1.000 | 0.000 |

Risk-zone pairs involving candidate units / trap terms (overlap ≥ 0.10, top 30;
full table in `overlap_focus_pairs.parquet`):

| disease A | disease B | overlap coeff | Jaccard | shared genes |
|---|---|---|---|---|
| efo:0005422 (skin aging) | efo:0022597 (aging) | 1.000 | 0.171 | 368 |
| efo:0005422 (skin aging) | mondo:0017948 (ABetaA21G amyloidosis) | 1.000 | 0.003 | 1 |
| efo:0005422 (skin aging) | mondo:0017949 (ABeta amyloidosis, Arctic type) | 1.000 | 0.003 | 1 |
| efo:0005678 (hippocampal sclerosis of aging) | mondo:0004975 (Alzheimer disease) | 1.000 | 0.000 | 3 |
| efo:0005678 (hippocampal sclerosis of aging) | mondo:0005148 (type 2 diabetes mellitus) | 1.000 | 0.000 | 3 |
| efo:0005678 (hippocampal sclerosis of aging) | mondo:0005180 (Parkinson disease) | 1.000 | 0.000 | 3 |
| efo:0009885 (frailty measurement) | mondo:0013645 (autosomal recessive spinocerebellar ataxia 11) | 1.000 | 0.091 | 1 |
| efo:0022597 (aging) | mondo:0007247 (basal ganglia calcification, idiopathic, childhood-onset) | 1.000 | 0.000 | 1 |
| efo:0022597 (aging) | mondo:0009839 (progressive supranuclear palsy-parkinsonism syndrome) | 1.000 | 0.000 | 1 |
| efo:0022597 (aging) | mondo:0013179 (hereditary spastic paraplegia 44) | 1.000 | 0.000 | 1 |
| efo:0022597 (aging) | mondo:0013905 (autosomal recessive spinocerebellar ataxia 13) | 1.000 | 0.000 | 1 |
| efo:0022597 (aging) | mondo:0015349 (progressive cavitating leukoencephalopathy) | 1.000 | 0.000 | 1 |
| efo:0022597 (aging) | mondo:0017830 (severe Canavan disease) | 1.000 | 0.000 | 1 |
| efo:0022597 (aging) | mondo:0017948 (ABetaA21G amyloidosis) | 1.000 | 0.000 | 1 |
| efo:0022597 (aging) | mondo:0017949 (ABeta amyloidosis, Arctic type) | 1.000 | 0.000 | 1 |
| efo:0022597 (aging) | mondo:0018316 (fatal post-viral neurodegenerative disorder) | 1.000 | 0.000 | 1 |
| efo:0022597 (aging) | mondo:0018416 (autosomal recessive spastic paraplegia type 59) | 1.000 | 0.000 | 1 |
| efo:0022597 (aging) | mondo:0030361 (Aicardi-Goutieres syndrome 8) | 1.000 | 0.000 | 1 |
| efo:0022597 (aging) | mondo:0859241 (neurodegeneration, childhood-onset, with progressive microcephaly) | 1.000 | 0.000 | 1 |
| mondo:0000828 (juvenile-onset Parkinson disease) | mondo:0005180 (Parkinson disease) | 1.000 | 0.031 | 207 |
| mondo:0000828 (juvenile-onset Parkinson disease) | mondo:0009839 (progressive supranuclear palsy-parkinsonism syndrome) | 1.000 | 0.005 | 1 |
| mondo:0000828 (juvenile-onset Parkinson disease) | mondo:0011025 (Cayman type cerebellar ataxia) | 1.000 | 0.005 | 1 |
| mondo:0000828 (juvenile-onset Parkinson disease) | mondo:0011706 (Kufor-Rakeb syndrome) | 1.000 | 0.952 | 197 |
| mondo:0000828 (juvenile-onset Parkinson disease) | mondo:0013060 (autosomal recessive Parkinson disease 14) | 1.000 | 0.005 | 1 |
| mondo:0000828 (juvenile-onset Parkinson disease) | mondo:0013905 (autosomal recessive spinocerebellar ataxia 13) | 1.000 | 0.005 | 1 |
| mondo:0000828 (juvenile-onset Parkinson disease) | mondo:0016089 (infantile Krabbe disease) | 1.000 | 0.010 | 2 |
| mondo:0000828 (juvenile-onset Parkinson disease) | mondo:0017279 (young-onset Parkinson disease) | 1.000 | 0.916 | 207 |
| mondo:0000828 (juvenile-onset Parkinson disease) | mondo:0017721 (Sandhoff disease, infantile form) | 1.000 | 0.010 | 2 |
| mondo:0000828 (juvenile-onset Parkinson disease) | mondo:0017809 (parkinsonism due to ATP13A2 deficiency) | 1.000 | 0.652 | 135 |
| mondo:0000828 (juvenile-onset Parkinson disease) | mondo:0017830 (severe Canavan disease) | 1.000 | 0.005 | 1 |

## 2. BFS separation — the θ-derivation input

Separation vs key-size-% across all supervised cluster diseases
(`bfs_separation.parquet`). Mean/median separation@3 by ceiling band:

| key-size band | n diseases | mean sep@3 | median sep@3 |
|---|---|---|---|
| above_5pct | 35 | +7.2 | +6.7 |
| below_2pct | 598 | +14.1 | +7.2 |
| theta_range_2_5pct | 29 | +11.0 | +9.1 |

Diseases at ≥2% of the gene universe (top 20) — separation collapses where the
key covers a large fraction of the genome:

| disease | key size | key % | ans ≤3h | non-ans ≤3h | sep@3 | sep@2 |
|---|---|---|---|---|---|---|
| mondo:0004975 (Alzheimer disease) | 12,414 | 71.15% | 95.9% | 98.4% | -2.5 | +21.2 |
| mondo:0024237 (inherited neurodegenerative disorder) | 11,636 | 66.69% | 99.7% | 93.1% | +6.7 | +91.3 |
| mondo:0005148 (type 2 diabetes mellitus) | 9,197 | 52.71% | 97.4% | 95.0% | +2.3 | +40.3 |
| mondo:0018899 (posterior cortical atrophy) | 6,835 | 39.17% | 100.0% | 88.3% | +11.7 | +54.8 |
| mondo:0021095 (parkinsonian disorder) | 6,708 | 38.45% | 100.0% | 94.5% | +5.5 | +99.2 |
| mondo:0005180 (Parkinson disease) | 6,584 | 37.73% | 99.9% | 99.4% | +0.5 | +89.7 |
| mondo:0005147 (type 1 diabetes mellitus) | 5,664 | 32.46% | 99.0% | 97.3% | +1.7 | +44.8 |
| mondo:0000167 (Huntington disease and related disorders) | 5,370 | 30.78% | 100.0% | 89.8% | +10.2 | +48.1 |
| mondo:0000437 (cerebellar ataxia) | 5,149 | 29.51% | 99.9% | 93.3% | +6.6 | +94.2 |
| mondo:0007739 (Huntington disease) | 4,935 | 28.28% | 100.0% | 99.5% | +0.5 | +89.1 |
| mondo:0022687 (cerebellar degeneration) | 4,900 | 28.08% | 100.0% | 89.1% | +10.9 | +99.2 |
| mondo:0100310 (hereditary cerebellar ataxia) | 4,854 | 27.82% | 100.0% | 88.7% | +11.3 | +97.5 |
| mondo:0002562 (demyelinating disease) | 4,772 | 27.35% | 99.9% | 88.6% | +11.2 | +96.3 |
| mondo:0005301 (multiple sclerosis) | 4,641 | 26.60% | 99.7% | 97.0% | +2.6 | +92.8 |
| mondo:0020800 (demyelinating disease of central nervous system) | 4,564 | 26.16% | 100.0% | 80.1% | +19.9 | +97.7 |
| mondo:0019303 (premature aging syndrome) | 4,049 | 23.21% | 95.9% | 79.6% | +16.3 | +40.8 |
| mondo:0020380 (autosomal dominant cerebellar ataxia) | 2,932 | 16.80% | 100.0% | 92.2% | +7.8 | +37.9 |
| mondo:0015244 (autosomal recessive cerebellar ataxia) | 2,929 | 16.79% | 100.0% | 87.8% | +12.2 | +86.7 |
| mondo:0019046 (leukodystrophy) | 2,180 | 12.49% | 100.0% | 94.7% | +5.3 | +38.1 |
| efo:0022597 (aging) | 2,150 | 12.32% | 74.8% | 49.5% | +25.3 | +17.1 |

## 3. Key-size ceiling classification (θ working range 2-5%)

Named reference terms (never graded — trap parents and ceiling references):

| disease | key size | key % of 17,448 | band |
|---|---|---|---|
| mondo:0004975 (Alzheimer disease) | 12,414 | 71.15% | above_5pct |
| mondo:0019303 (premature aging syndrome) | 4,049 | 23.21% | above_5pct |
| efo:0022597 (aging) | 2,150 | 12.32% | above_5pct |
| mondo:0015333 (progeroid syndrome) | 1,589 | 9.11% | above_5pct |
| mondo:0020732 (progeria) | 706 | 4.05% | theta_range_2_5pct |

Band counts across 662 supervised cluster diseases: **above_5pct**: 35; **below_2pct**: 598; **theta_range_2_5pct**: 29.

## 4. Connectivity census + roster attrition (draft)

Across the 47 candidate units: failing the
≥10-context-edge floor: **18**; below the evidence floor
(<20 supervision edges): **10**; secondary
panel (20-49): **7**;
above the key-size ceiling (>5%): **6**.
Cluster-wide (all 678 diseases incl. non-candidates): failing context
floor 227, below evidence floor
259, secondary
124, above ceiling
35.

## 5. Per-unit audit table (all candidate units)

Full machine-readable table: `unit_audit.parquet`. Δ MRR = filtered MRR with
DIS-PHE minus without (restart 0.5; degenerate all-tie pairs
excluded). "max overlap" is restricted to partners below the >5% ceiling band
(the co-holdout-relevant statistic).

| unit | stratum | key | key % | ceiling band | ctx edges | evidence | sep@3 | max overlap | Δ MRR |
|---|---|---|---|---|---|---|---|---|---|
| mondo:0010196 | A_primary | 857 | 4.91% | theta_range_2_5pct | 37 | primary | +8.9 | 1.000 | -0.0010 |
| mondo:0008310 | A_primary | 468 | 2.68% | theta_range_2_5pct | 23 | primary | +2.9 | 1.000 | -0.0026 |
| mondo:0012853 | A_primary | 230 | 1.32% | below_2pct | 98 | primary | +2.1 | 1.000 | — |
| mondo:0014157 | A_primary | 60 | 0.34% | below_2pct | 29 | primary | +2.1 | 0.500 | +0.0005 |
| orphanet:75496 | A_secondary | 45 | 0.26% | below_2pct | 88 | secondary | +12.1 | 0.655 | — |
| mondo:0018203 | A_secondary | 41 | 0.23% | below_2pct | 2 | secondary | +15.8 | 1.000 | +0.0000 |
| mondo:0013523 | A_secondary | 35 | 0.20% | below_2pct | 52 | secondary | +11.7 | 1.000 | +0.0295 |
| mondo:0014527 | A_secondary | 29 | 0.17% | below_2pct | 26 | secondary | +11.7 | 0.448 | — |
| mondo:0011150 | A_secondary | 29 | 0.17% | below_2pct | 65 | secondary | +2.1 | 0.724 | -0.0120 |
| mondo:0019321 | A_secondary | 27 | 0.15% | below_2pct | 1 | secondary | +48.7 | 0.704 | +0.0000 |
| mondo:0957953 | A_secondary | 17 | 0.10% | below_2pct | 31 | below_floor | +16.1 | 1.000 | +0.0003 |
| mondo:0030880 | A_secondary | 16 | 0.09% | below_2pct | 61 | below_floor | +11.8 | 0.750 | +0.0282 |
| mondo:0008311 | A_secondary | 12 | 0.07% | below_2pct | 32 | below_floor | +2.1 | 0.500 | +0.0050 |
| orphanet:2963 | A_secondary | 10 | 0.06% | below_2pct | 96 | below_floor | +2.1 | 0.500 | — |
| mondo:0957266 | A_secondary | 8 | 0.05% | below_2pct | 39 | below_floor | +12.6 | 0.500 | -0.0183 |
| mondo:0012590 | A_secondary | 8 | 0.05% | below_2pct | 38 | below_floor | +10.9 | 0.625 | -0.0077 |
| mondo:0008694 | A_secondary | 6 | 0.03% | below_2pct | 11 | below_floor | +13.7 | 0.333 | +0.0017 |
| mondo:0011149 | A_secondary | 0 | 0.00% | no_supervision | 0 | below_floor | — | 0.000 | — |
| efo:0022597 | B_aging_trait | 2,150 | 12.32% | above_5pct | 1 | primary | +25.3 | 1.000 | -0.0000 |
| efo:0005422 | B_aging_trait | 368 | 2.11% | theta_range_2_5pct | 8 | primary | +10.1 | 1.000 | +0.0000 |
| oba:0005494 | B_aging_trait | 150 | 0.86% | below_2pct | 1 | primary | -0.1 | 0.250 | +0.0000 |
| efo:0007796 | B_aging_trait | 68 | 0.39% | below_2pct | 1 | primary | +0.0 | 0.182 | — |
| efo:0009885 | B_aging_trait | 11 | 0.06% | below_2pct | 1 | below_floor | -0.9 | 1.000 | -0.0000 |
| oba:VT0005372 | B_aging_trait | 6 | 0.03% | below_2pct | 2 | below_floor | -0.9 | 0.500 | +0.0000 |
| mondo:0004975 | C_alzheimer | 12,414 | 71.15% | above_5pct | 336 | primary | -2.5 | 1.000 | +0.0005 |
| mondo:0100087 | C_alzheimer | 405 | 2.32% | theta_range_2_5pct | 4 | primary | +3.5 | 1.000 | +0.0007 |
| mondo:0015140 | C_alzheimer | 167 | 0.96% | below_2pct | 4 | primary | +8.7 | 1.000 | +0.0001 |
| mondo:0011913 | C_alzheimer | 75 | 0.43% | below_2pct | 23 | primary | +2.5 | 1.000 | +0.0076 |
| mondo:0007088 | C_alzheimer | 30 | 0.17% | below_2pct | 10 | secondary | +3.4 | 1.000 | +0.0305 |
| mondo:0014265 | C_alzheimer | 19 | 0.11% | below_2pct | 3 | below_floor | +3.4 | 1.000 | -0.0000 |
| mondo:0005180 | C_parkinson | 6,584 | 37.73% | above_5pct | 293 | primary | +0.5 | 1.000 | +0.0003 |
| mondo:0017279 | C_parkinson | 226 | 1.30% | below_2pct | 2 | primary | +12.1 | 1.000 | +0.0000 |
| mondo:0000828 | C_parkinson | 207 | 1.19% | below_2pct | 3 | primary | +36.2 | 1.000 | -0.0796 |
| mondo:0957576 | C_parkinson | 151 | 0.87% | below_2pct | 13 | primary | +4.6 | 1.000 | +0.0070 |
| mondo:0008199 | C_parkinson | 66 | 0.38% | below_2pct | 32 | primary | +2.1 | 1.000 | +0.0070 |
| mondo:0005148 | D_diabetes | 9,197 | 52.71% | above_5pct | 419 | primary | +2.3 | 1.000 | +0.0001 |
| mondo:0005147 | D_diabetes | 5,664 | 32.46% | above_5pct | 242 | primary | +1.7 | 1.000 | +0.0002 |
| mondo:0016391 | D_diabetes | 1,161 | 6.65% | above_5pct | 6 | primary | +16.3 | 1.000 | -0.0019 |
| mondo:0020525 | D_diabetes | 601 | 3.44% | theta_range_2_5pct | 3 | primary | +32.1 | 1.000 | -0.0071 |
| mondo:0100164 | D_diabetes | 532 | 3.05% | theta_range_2_5pct | 15 | primary | +10.8 | 1.000 | -0.3863 |
| mondo:0030089 | D_diabetes | 467 | 2.68% | theta_range_2_5pct | 9 | primary | +2.8 | 1.000 | -0.1003 |
| mondo:0018911 | D_diabetes | 389 | 2.23% | theta_range_2_5pct | 9 | primary | +30.3 | 1.000 | -0.0686 |
| mondo:0100165 | D_diabetes | 381 | 2.18% | theta_range_2_5pct | 10 | primary | +10.9 | 1.000 | -0.0101 |
| mondo:0012192 | D_diabetes | 335 | 1.92% | below_2pct | 33 | primary | +5.0 | 0.778 | -0.0008 |
| mondo:0007669 | D_diabetes | 307 | 1.76% | below_2pct | 38 | primary | +2.2 | 1.000 | -0.0489 |
| mondo:0012522 | D_diabetes | 302 | 1.73% | below_2pct | 12 | primary | +3.5 | 1.000 | -0.0099 |
| mondo:0012480 | D_diabetes | 214 | 1.23% | below_2pct | 4 | primary | +6.8 | 1.000 | -0.0069 |
| mondo:0010785 | D_diabetes | 110 | 0.63% | below_2pct | 17 | primary | +6.2 | 0.800 | +0.0040 |

## 6. Coherence control statistic (v5.1 §2.2)

| statistic | cluster | null mean ± SD (n=100) | percentile |
|---|---|---|---|
| DIS-DIS density (edges/possible pairs) | 0.002749 | 0.000075 ± 0.000122 | 100th |
| mean supervision degree | 283.1 | 271.0 ± 51.5 | 59th |

## 7. PPR ± DIS-PHE ablation

Per-unit filtered metrics at pinned restart 0.5 (sensitivity rows
[0.15, 0.85] in `ppr_disphe.parquet`):

| unit | stratum | restart | MRR with | MRR without | H@10 with | H@10 without |
|---|---|---|---|---|---|---|
| efo:0005422 | B_aging_trait | 0.5 | 0.0017 | 0.0017 | 0.0000 | 0.0000 |
| efo:0009885 | B_aging_trait | 0.5 | 0.0002 | 0.0002 | 0.0000 | 0.0000 |
| efo:0022597 | B_aging_trait | 0.5 | 0.1772 | 0.1773 | 0.1828 | 0.1828 |
| mondo:0000828 | C_parkinson | 0.5 | 0.5771 | 0.6568 | 0.9565 | 0.9662 |
| mondo:0004975 | C_alzheimer | 0.5 | 0.0073 | 0.0068 | 0.0091 | 0.0091 |
| mondo:0005147 | D_diabetes | 0.5 | 0.0033 | 0.0031 | 0.0030 | 0.0030 |
| mondo:0005148 | D_diabetes | 0.5 | 0.0055 | 0.0055 | 0.0084 | 0.0084 |
| mondo:0005180 | C_parkinson | 0.5 | 0.0061 | 0.0058 | 0.0073 | 0.0073 |
| mondo:0007088 | C_alzheimer | 0.5 | 0.0774 | 0.0469 | 0.1333 | 0.0333 |
| mondo:0007669 | D_diabetes | 0.5 | 0.0159 | 0.0648 | 0.0293 | 0.0945 |
| mondo:0008199 | C_parkinson | 0.5 | 0.0251 | 0.0181 | 0.0758 | 0.0152 |
| mondo:0008310 | A_primary | 0.5 | 0.0057 | 0.0083 | 0.0021 | 0.0043 |
| mondo:0008311 | A_secondary | 0.5 | 0.0057 | 0.0007 | 0.0000 | 0.0000 |
| mondo:0008694 | A_secondary | 0.5 | 0.0021 | 0.0004 | 0.0000 | 0.0000 |
| mondo:0010196 | A_primary | 0.5 | 0.0034 | 0.0045 | 0.0047 | 0.0047 |
| mondo:0010785 | D_diabetes | 0.5 | 0.0075 | 0.0036 | 0.0182 | 0.0091 |
| mondo:0011150 | A_secondary | 0.5 | 0.0076 | 0.0196 | 0.0000 | 0.0690 |
| mondo:0011913 | C_alzheimer | 0.5 | 0.0275 | 0.0199 | 0.0667 | 0.0133 |
| mondo:0012192 | D_diabetes | 0.5 | 0.0026 | 0.0034 | 0.0030 | 0.0090 |
| mondo:0012480 | D_diabetes | 0.5 | 0.0100 | 0.0169 | 0.0093 | 0.0280 |
| mondo:0012522 | D_diabetes | 0.5 | 0.0122 | 0.0221 | 0.0066 | 0.0232 |
| mondo:0012590 | A_secondary | 0.5 | 0.0132 | 0.0209 | 0.0000 | 0.1250 |
| mondo:0013523 | A_secondary | 0.5 | 0.0338 | 0.0043 | 0.0286 | 0.0000 |
| mondo:0014157 | A_primary | 0.5 | 0.0019 | 0.0014 | 0.0000 | 0.0000 |
| mondo:0014265 | C_alzheimer | 0.5 | 0.0007 | 0.0007 | 0.0000 | 0.0000 |
| mondo:0015140 | C_alzheimer | 0.5 | 0.5830 | 0.5829 | 0.6407 | 0.6407 |
| mondo:0016391 | D_diabetes | 0.5 | 0.0435 | 0.0454 | 0.0069 | 0.0146 |
| mondo:0017279 | C_parkinson | 0.5 | 0.9127 | 0.9127 | 0.9159 | 0.9159 |
| mondo:0018203 | A_secondary | 0.5 | 0.0004 | 0.0004 | 0.0000 | 0.0000 |
| mondo:0018911 | D_diabetes | 0.5 | 0.4553 | 0.5239 | 0.7172 | 0.7506 |
| mondo:0019321 | A_secondary | 0.5 | 0.0023 | 0.0023 | 0.0000 | 0.0000 |
| mondo:0020525 | D_diabetes | 0.5 | 0.3673 | 0.3744 | 0.7571 | 0.7571 |
| mondo:0030089 | D_diabetes | 0.5 | 0.0269 | 0.1272 | 0.0064 | 0.3833 |
| mondo:0030880 | A_secondary | 0.5 | 0.0342 | 0.0060 | 0.0625 | 0.0000 |
| mondo:0100087 | C_alzheimer | 0.5 | 0.4142 | 0.4135 | 0.4173 | 0.4148 |
| mondo:0100164 | D_diabetes | 0.5 | 0.0857 | 0.4721 | 0.0996 | 0.8327 |
| mondo:0100165 | D_diabetes | 0.5 | 0.0170 | 0.0271 | 0.0079 | 0.0446 |
| mondo:0957266 | A_secondary | 0.5 | 0.0026 | 0.0209 | 0.0000 | 0.1250 |
| mondo:0957576 | C_parkinson | 0.5 | 0.0075 | 0.0005 | 0.0066 | 0.0000 |
| mondo:0957953 | A_secondary | 0.5 | 0.0078 | 0.0075 | 0.0000 | 0.0000 |
| oba:0005494 | B_aging_trait | 0.5 | 0.0004 | 0.0004 | 0.0000 | 0.0000 |
| oba:VT0005372 | B_aging_trait | 0.5 | 0.0003 | 0.0003 | 0.0000 | 0.0000 |

**Mean Δ across 42 informative units: ΔMRR -0.0152, ΔHits@10 -0.0310.
Verdict: DROP DIS-PHE from the BRG** — phenotype context edges are net-negative for
topological recovery signal on this graph.
(4 unit-variant rows were degenerate — PPR assigns every gene
score 0, a total tie the metric protocol reads as rank 1; excluded from the
delta. These are the same near-isolated units that fail the ≥10-context-edge
floor in §4.)

## 8. Noise-floor simulation — the evidence floor is derived, not heuristic

Filtered Hits@10 of a random ranking of 17,448 genes
(1000 seeded resamples per key size, seed 1337):

| key size | mean H@10 | SD | CV | P(zero hits) |
|---|---|---|---|---|
| 20 | 0.00030 | 0.00386 | 12.88 | 0.994 |
| 30 | 0.00030 | 0.00315 | 10.50 | 0.991 |
| 50 | 0.00042 | 0.00287 | 6.83 | 0.979 |
| 75 | 0.00047 | 0.00245 | 5.25 | 0.965 |
| 100 | 0.00048 | 0.00214 | 4.46 | 0.952 |
| 200 | 0.00052 | 0.00158 | 3.03 | 0.898 |
| 500 | 0.00059 | 0.00106 | 1.80 | 0.740 |
| 1000 | 0.00060 | 0.00076 | 1.26 | 0.536 |
| 2000 | 0.00065 | 0.00056 | 0.87 | 0.266 |

The coefficient of variation inflates sharply at small key sizes — below the
~50-edge region a handful of answers cannot support a stable Hits@10
estimate, which is the derived basis for the ≥50-edge
evidence floor (20-49 → secondary panel, never pooled).

## 9. EFO/OBA resolution pinning (v5.1 §2.1)

| term | label | mondo_id | resolution method | in Track-D | near-synonym bound |
|---|---|---|---|---|---|
| efo:0022597 | aging | — | explicit_id_list_census | yes | no MONDO xref in crosswalk; near-synonym miss possible, unbounded |
| efo:0007796 | parental longevity | — | explicit_id_list_census | yes | no MONDO xref in crosswalk; near-synonym miss possible, unbounded |
| oba:0005494 | aging rate | — | explicit_id_list_census | yes | no MONDO xref in crosswalk; near-synonym miss possible, unbounded |
| efo:0005422 | skin aging | — | explicit_id_list_census | yes | no MONDO xref in crosswalk; near-synonym miss possible, unbounded |
| efo:0009885 | frailty measurement | — | explicit_id_list_census | yes | no MONDO xref in crosswalk; near-synonym miss possible, unbounded |
| oba:VT0005372 | life span trait | — | explicit_id_list_census | yes | no MONDO xref in crosswalk; near-synonym miss possible, unbounded |
| efo:0005678 | hippocampal sclerosis of aging | — | explicit_id_list_census | yes | MONDO xref exists but not mapped into Track-D; near-synonym miss possible, bounded to the unmapped xref |
| mondo:0008310 | Hutchinson-Gilford progeria syndrome | mondo:0008310 | explicit_id_list_census | yes | — |
| mondo:0008311 | progeria-short stature-pigmented nevi syndrome | mondo:0008311 | explicit_id_list_census | yes | — |
| mondo:0008694 | pseudoprogeria syndrome | mondo:0008694 | explicit_id_list_census | yes | — |
| mondo:0011149 |  | — | explicit_id_list_census | NO | — |
| mondo:0011150 | acroosteolysis-keloid-like lesions-premature aging syndrome | mondo:0011150 | explicit_id_list_census | yes | — |
| mondo:0010196 | Werner syndrome | mondo:0010196 | explicit_id_list_census | yes | — |
| mondo:0019321 | atypical Werner syndrome | mondo:0019321 | explicit_id_list_census | yes | — |
| mondo:0012853 | Fontaine progeroid syndrome | mondo:0012853 | explicit_id_list_census | yes | — |
| mondo:0018203 | LMNA-related cardiocutaneous progeria syndrome | mondo:0018203 | explicit_id_list_census | yes | — |
| mondo:0013523 | Nestor-Guillermo progeria syndrome | mondo:0013523 | explicit_id_list_census | yes | — |
| mondo:0014157 | mandibular hypoplasia-deafness-progeroid syndrome | mondo:0014157 | explicit_id_list_census | yes | — |
| mondo:0030880 | mandibuloacral dysplasia progeroid syndrome | mondo:0030880 | explicit_id_list_census | yes | — |
| mondo:0012590 | XFE progeroid syndrome | mondo:0012590 | explicit_id_list_census | yes | — |
| mondo:0957266 | RECON progeroid syndrome | mondo:0957266 | explicit_id_list_census | yes | — |
| mondo:0957953 | Garg-Mishra progeroid syndrome | mondo:0957953 | explicit_id_list_census | yes | — |
| orphanet:75496 | Ehlers-Danlos syndrome, progeroid type | — | explicit_id_list_census | yes | MONDO near-equivalent(s) mondo:0007526 present in Track-D; census name-search could have missed near-synonyms keyed under those nodes — bounded by this xref check |
| orphanet:2963 | Progeroid syndrome, Petty type | — | explicit_id_list_census | yes | MONDO near-equivalent(s) mondo:0012853 present in Track-D; census name-search could have missed near-synonyms keyed under those nodes — bounded by this xref check |
| mondo:0014527 | progeroid features-hepatocellular carcinoma predisposition syndrome | mondo:0014527 | explicit_id_list_census | yes | — |
| mondo:0014831 | progeroid and marfanoid aspect-lipodystrophy syndrome | mondo:0014831 | explicit_id_list_census | yes | — |
| mondo:0044634 | retinitis pigmentosa-hearing loss-premature aging-short stature-facial dysmorphism syndrome | mondo:0044634 | explicit_id_list_census | yes | — |
| orphanet:276432 | Premature aging appearance-developmental delay-cardiac arrhythmia syndrome | — | explicit_id_list_census | yes | MONDO near-equivalent(s) mondo:0010457 present in Track-D; census name-search could have missed near-synonyms keyed under those nodes — bounded by this xref check |
| mondo:0020477 | progeria-associated arthropathy | mondo:0020477 | explicit_id_list_census | yes | — |
| mondo:0020732 | progeria | mondo:0020732 | explicit_id_list_census | yes | — |
| mondo:0019303 | premature aging syndrome | mondo:0019303 | explicit_id_list_census | yes | — |
| mondo:0015333 | progeroid syndrome | mondo:0015333 | explicit_id_list_census | yes | — |
| mondo:0100087 | familial Alzheimer disease | mondo:0100087 | mondo_closure_family_list | yes | — |
| mondo:0015140 | early-onset autosomal dominant Alzheimer disease | mondo:0015140 | mondo_closure_family_list | yes | — |
| mondo:0011913 | Alzheimer disease 3 | mondo:0011913 | mondo_closure_family_list | yes | — |
| mondo:0007088 | Alzheimer disease type 1 | mondo:0007088 | mondo_closure_family_list | yes | — |
| mondo:0014265 | Alzheimer disease 18 | mondo:0014265 | mondo_closure_family_list | yes | — |
| mondo:0004975 | Alzheimer disease | mondo:0004975 | mondo_closure_family_list | yes | — |
| mondo:0000828 | juvenile-onset Parkinson disease | mondo:0000828 | mondo_closure_family_list | yes | — |
| mondo:0005180 | Parkinson disease | mondo:0005180 | mondo_closure_family_list | yes | — |
| mondo:0008199 | late-onset Parkinson disease | mondo:0008199 | mondo_closure_family_list | yes | — |
| mondo:0017279 | young-onset Parkinson disease | mondo:0017279 | mondo_closure_family_list | yes | — |
| mondo:0957576 | parkinson disease 25, autosomal recessive early-onset, with impaired intellectual development | mondo:0957576 | mondo_closure_family_list | yes | — |
| mondo:0005147 | type 1 diabetes mellitus | mondo:0005147 | mondo_closure_descendants | yes | — |
| mondo:0005148 | type 2 diabetes mellitus | mondo:0005148 | mondo_closure_descendants | yes | — |
| mondo:0007669 | renal cysts and diabetes syndrome | mondo:0007669 | mondo_closure_descendants | yes | — |
| mondo:0010785 | maternally-inherited diabetes and deafness | mondo:0010785 | mondo_closure_descendants | yes | — |
| mondo:0012192 | permanent neonatal diabetes mellitus-pancreatic and cerebellar agenesis syndrome | mondo:0012192 | mondo_closure_descendants | yes | — |
| mondo:0012480 | diabetes mellitus, transient neonatal, 2 | mondo:0012480 | mondo_closure_descendants | yes | — |
| mondo:0012522 | diabetes mellitus, transient neonatal, 3 | mondo:0012522 | mondo_closure_descendants | yes | — |
| mondo:0016391 | neonatal diabetes mellitus | mondo:0016391 | mondo_closure_descendants | yes | — |
| mondo:0018911 | maturity-onset diabetes of the young | mondo:0018911 | mondo_closure_descendants | yes | — |
| mondo:0020525 | transient neonatal diabetes mellitus | mondo:0020525 | mondo_closure_descendants | yes | — |
| mondo:0030089 | diabetes mellitus, permanent neonatal 4 | mondo:0030089 | mondo_closure_descendants | yes | — |
| mondo:0100164 | permanent neonatal diabetes mellitus | mondo:0100164 | mondo_closure_descendants | yes | — |
| mondo:0100165 | permanent neonatal diabetes mellitus 1 | mondo:0100165 | mondo_closure_descendants | yes | — |

## 10. Feeds into the CHANGES/0020 freeze package

The freeze entry needs exactly these numbers from this audit:

1. **Co-holdout roster:** the 30,975 overlap pairs ≥ 0.3
   with both keys below the >5% ceiling band (plus 19,789 disclose
   pairs; full matrix in `overlap_matrix.parquet`, flags in
   `overlap_flags.parquet`) + hand-curated seed clusters — each cluster counts
   as ONE statistical unit.
2. **θ (key-size ceiling):** set from the separation-by-band table (§2) — the
   2-5% working range is justified if separation@3 collapses above it; publish
   per-unit key % + separation in the frozen roster.
3. **Roster attrition:** the per-unit attrition table (`unit_audit.parquet`,
   §5) with published reasons — 18 context-floor, 10
   evidence-floor, 7 secondary-panel, 6
   above-ceiling among candidate units.
4. **BRG composition:** the DIS-PHE keep/drop verdict (§7: ΔMRR -0.0152,
   ΔH@10 -0.0310 across 42 informative units, 4
   degenerate all-tie rows excluded) is recorded in the frozen BRG.
5. **Evidence floor derivation:** the noise-floor table (§8) accompanies the
   frozen roster.
6. **Coherence disclosure:** density percentile 100th,
   degree percentile 59th (§6).
7. **EFO/OBA pinning:** the resolution table (`efo_pinning.parquet`, §9)
   bounds roster incompleteness.

## Reproduction

`python scripts/m2_lodo_audits.py --config configs/m2_lodo_audits.yaml --resume`
