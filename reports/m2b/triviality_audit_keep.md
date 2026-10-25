<!-- track_scope: benchmark -->
# Triviality audit (A3, CHANGES/0021)

Per-unit shortest-path distance in `train1` (BRG, 1,431,945 (incl. DIS-PHE, KEEP graph) edges) from held-out disease to test answers. Degenerate if baseline MRR>=0.99 OR >10% <=2 hops.

| unit | tier | test | ≤1 | ≤2 | ≤2% | 3 | 4 | 5 | unreachable | degenerate |
|---|---|---|---|---|---|---|---|---|---|---|
| mondo:0014157 | primary | 52 | 0 | 49 | 0.942 | 3 | 0 | 0 | 0 | True |
| mondo:0012853 | primary | 201 | 0 | 183 | 0.910 | 18 | 0 | 0 | 0 | True |
| mondo:0007669 | primary | 258 | 0 | 234 | 0.907 | 21 | 3 | 0 | 0 | True |
| orphanet:75496 | secondary_evidence | 41 | 0 | 36 | 0.878 | 5 | 0 | 0 | 0 | True |
| mondo:0013523 | secondary_evidence | 32 | 0 | 28 | 0.875 | 4 | 0 | 0 | 0 | True |
| mondo:0010196 | primary | 730 | 0 | 612 | 0.838 | 115 | 3 | 0 | 0 | True |
| mondo:0012192 | primary | 281 | 0 | 221 | 0.786 | 51 | 9 | 0 | 0 | True |
| mondo:0010785 | primary | 96 | 0 | 74 | 0.771 | 21 | 1 | 0 | 0 | True |
| mondo:0014527 | secondary_evidence | 26 | 0 | 20 | 0.769 | 5 | 1 | 0 | 0 | True |
| mondo:0011150 | secondary_evidence | 24 | 0 | 18 | 0.750 | 6 | 0 | 0 | 0 | True |
| mondo:0008199 | primary | 57 | 0 | 42 | 0.737 | 13 | 2 | 0 | 0 | True |
| mondo:0007088 | secondary_evidence | 27 | 0 | 19 | 0.704 | 6 | 2 | 0 | 0 | True |
| mondo:0011913 | primary | 68 | 0 | 44 | 0.647 | 21 | 3 | 0 | 0 | True |
| mondo:0008310 | primary | 394 | 0 | 252 | 0.640 | 132 | 10 | 0 | 0 | True |
| mondo:0100165 | primary | 329 | 0 | 162 | 0.492 | 149 | 18 | 0 | 0 | True |
| mondo:0100164 | primary | 459 | 0 | 220 | 0.479 | 216 | 23 | 0 | 0 | True |
| mondo:0012522 | primary | 251 | 0 | 113 | 0.450 | 127 | 11 | 0 | 0 | True |
| mondo:0030089 | secondary_low_context | 394 | 0 | 166 | 0.421 | 203 | 25 | 0 | 0 | True |
| mondo:0957576 | primary | 127 | 0 | 53 | 0.417 | 61 | 13 | 0 | 0 | True |
| efo:0005422 | secondary_low_context | 325 | 0 | 5 | 0.015 | 93 | 227 | 0 | 0 | False |
| mondo:0000828 | secondary_low_context | 170 | 0 | 0 | 0.000 | 130 | 32 | 8 | 0 | False |
| mondo:0012480 | secondary_low_context | 178 | 0 | 0 | 0.000 | 0 | 165 | 8 | 0 | False |
| mondo:0015140 | secondary_low_context | 144 | 0 | 0 | 0.000 | 100 | 40 | 4 | 0 | False |
| mondo:0017279 | secondary_low_context | 188 | 0 | 0 | 0.000 | 148 | 34 | 6 | 0 | False |
| mondo:0018203 | secondary_low_context | 34 | 0 | 0 | 0.000 | 0 | 33 | 1 | 0 | False |
| mondo:0018911 | secondary_low_context | 319 | 0 | 0 | 0.000 | 296 | 20 | 3 | 0 | False |
| mondo:0019321 | secondary_low_context | 25 | 0 | 0 | 0.000 | 0 | 21 | 2 | 0 | False |
| mondo:0020525 | secondary_low_context | 498 | 0 | 0 | 0.000 | 221 | 246 | 31 | 0 | False |
| mondo:0100087 | secondary_low_context | 336 | 0 | 0 | 0.000 | 218 | 118 | 0 | 0 | False |

Degenerate units: 19 / 29 (excluded from primary gate, reported separately).

Tiny smoke test: ~10-node graph baseline MRR <1 asserted in unit tests.

