<!-- track_scope: benchmark -->
# Phase D wave stats — pre-registered pins (CHANGES/0020 §13, CHANGES/0024)

Model vs best same-graph null (degree/PPR) on per-unit test MRR, 3 seed(s) [1337, 1339, 5027], 13 primary LODO units, KEEP graph (train1 = 1,431,945). Macro test MRR over seeds: 0.00775 ± 0.00008 (95% CI [0.00754, 0.00796]). Disease-associated gene recovery benchmark; MRR exact.

**Degenerate-flag conflict (CHANGES/0024 §(b), owner decision pending):** the frozen A3 ≤2-hop prong flags all 13 primaries on the KEEP graph. Both variants are reported; neither preempts the owner call.

## Primary gate — all 13 primaries
- units: 13
- exact Wilcoxon (2^n enumeration): stat 31.0, p = 0.3396
- hierarchical BCa (n=1000): mean 0.00775, 95% CI [0.00548, 0.01101]
- cluster-label permutation: p = 0.9381

## Primary gate — degenerate-excluded (frozen A3, literal)
- <2 non-degenerate units

## Per-unit

| unit | cluster | model MRR (seed mean) | best null MRR | A3 degenerate (KEEP) |
|---|---|---|---|---|
| mondo:0008199 | c3_pd | 0.01639 | 0.0170 | True |
| mondo:0100164 | c5_ndm | 0.01624 | 0.0144 | True |
| mondo:0011913 | c2_ad_familial | 0.01239 | 0.0157 | True |
| mondo:0008310 | c1_progeria | 0.01007 | 0.0070 | True |
| mondo:0007669 | c4_mody | 0.00994 | 0.0137 | True |
| mondo:0012853 | c1_progeria | 0.00641 | 0.0128 | True |
| mondo:0010196 | c1_progeria | 0.00620 | 0.0062 | True |
| mondo:0012192 | c5_ndm | 0.00542 | 0.0075 | True |
| mondo:0012522 | c5_ndm | 0.00447 | 0.0019 | True |
| mondo:0010785 | c6_singletons | 0.00398 | 0.0082 | True |
| mondo:0014157 | c6_singletons | 0.00354 | 0.0018 | True |
| mondo:0100165 | c5_ndm | 0.00348 | 0.0028 | True |
| mondo:0957576 | c3_pd | 0.00222 | 0.0031 | True |

Reproduction: `python scripts/m2_wave_stats.py --config configs/m2_wave_stats.yaml`. With 1 seed (Wave 1 only) all stats are directional; the pre-registered read requires seeds {1337,1339,5027}.
