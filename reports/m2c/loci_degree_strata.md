<!-- track_scope: benchmark -->
# Locus panels: degree-stratified rescore

Generated 2026-10-26T03:49:28Z. Starting point `scripts/pod/m2c_loci_headtohead.py` (imported verbatim grading: same score dumps, loci, known definitions, pessimistic ranks). Degree = total endpoint occurrences over the completion-track training graph (train1 + train2, all relations), the `src/sina/baselines.py::degree_baseline` definition. Supervision count = disease-headed (`mondo:`/`efo:`/`orphanet:`) `ASSOCIATED_WITH` edges with the gene as tail, train1 + train2. Strata split gradeable loci at the median of the locus max-known degree. Bootstrap CIs (n=1000, seed 1337) match the head-to-head.

## mondo__0005180 / train2 (n=52)

Verification: rescore top-1 0.9038 vs committed 0.9038 (EXACT); checkpoint `reports/m2c/checkpoints_completion_noana/mondo__0005180/epoch_001.pth`; chance-draw mean 0.6767.
Known-gene degree: median 110.5, Q1 55.0, Q3 170.2 (universe median 76.0, bottom-decile cut 15.0). Fraction of graded known genes in universe bottom decile: 0.0424; with zero disease supervision: 0.0.

| stratum | n loci | top-1 | 95% CI | median locus degree |
|---|---|---|---|---|
| bottom_half | 26 | 0.9231 | [0.8077, 1.0] | 95.5 |
| top_half | 26 | 0.8846 | [0.7308, 1.0] | 236.0 |

Reading: the below-chance test-knowns rates (n/a — train panel) in the high-degree stratum: top-half top-1 0.8846 [0.7308, 1.0] vs chance 0.6767.

## mondo__0005180 / test (n=10)

Verification: rescore top-1 0.3000 vs committed 0.3000 (EXACT); checkpoint `reports/m2c/checkpoints_completion_noana/mondo__0005180/epoch_001.pth`; chance-draw mean 0.33.
Known-gene degree: median 131.0, Q1 90.5, Q3 274.5 (universe median 76.0, bottom-decile cut 15.0). Fraction of graded known genes in universe bottom decile: 0.0; with zero disease supervision: 0.0.

| stratum | n loci | top-1 | 95% CI | median locus degree |
|---|---|---|---|---|
| bottom_half | 5 | 0.2000 | [0.0, 0.6] | 93.0 |
| top_half | 5 | 0.4000 | [0.0, 0.8] | 516.0 |

Reading: the below-chance test-knowns rates (persist) in the high-degree stratum: top-half top-1 0.4000 [0.0, 0.8] vs chance 0.33.

## mondo__0004975 / train2 (n=239)

Verification: rescore top-1 0.9540 vs committed 0.9540 (EXACT); checkpoint `reports/m2c/checkpoints_completion_noana/mondo__0004975/epoch_001.pth`; chance-draw mean 0.7391.
Known-gene degree: median 103.0, Q1 54.0, Q3 172.8 (universe median 75.0, bottom-decile cut 15.0). Fraction of graded known genes in universe bottom decile: 0.0631; with zero disease supervision: 0.0.

| stratum | n loci | top-1 | 95% CI | median locus degree |
|---|---|---|---|---|
| bottom_half | 119 | 0.9328 | [0.8908, 0.9748] | 97.0 |
| top_half | 120 | 0.9750 | [0.9417, 1.0] | 254.0 |

Reading: the below-chance test-knowns rates (n/a — train panel) in the high-degree stratum: top-half top-1 0.9750 [0.9417, 1.0] vs chance 0.7391.

## mondo__0004975 / test (n=31)

Verification: rescore top-1 0.1613 vs committed 0.1613 (EXACT); checkpoint `reports/m2c/checkpoints_completion_noana/mondo__0004975/epoch_001.pth`; chance-draw mean 0.3277.
Known-gene degree: median 96.5, Q1 40.2, Q3 168.2 (universe median 75.0, bottom-decile cut 15.0). Fraction of graded known genes in universe bottom decile: 0.1579; with zero disease supervision: 0.1842.

| stratum | n loci | top-1 | 95% CI | median locus degree |
|---|---|---|---|---|
| bottom_half | 15 | 0.0667 | [0.0, 0.2] | 39.0 |
| top_half | 16 | 0.2500 | [0.0625, 0.5] | 188.0 |

Reading: the below-chance test-knowns rates (persist) in the high-degree stratum: top-half top-1 0.2500 [0.0625, 0.5] vs chance 0.3277.

## Limitations

1. Train-known genes are train positives — high top-1 there is recovery, not discovery.
2. Degree is computed on the training graph the model saw; high-degree test-knowns remain on-distribution by construction.
3. One seed (1337), frozen checkpoints, no retraining — CIs are over loci, not seeds.
