<!-- track_scope: benchmark -->
# Triviality audit (A3, CHANGES/0021)

Per-unit shortest-path distance in `train1` (BRG, 1,420,095 edges) from held-out disease to test answers. Degenerate if baseline MRR≥0.99 OR >10% ≤2 hops.

| unit | tier | test | ≤1 | ≤2 | ≤2% | 3 | 4 | 5 | unreachable | degenerate |
|---|---|---|---|---|---|---|---|---|---|---|
| efo:0005422 | secondary_low_context | 325 | 0 | 5 | 0.015 | 93 | 227 | 0 | 0 | False |
| mondo:0008310 | primary | 394 | 0 | 3 | 0.008 | 57 | 334 | 0 | 0 | False |
| mondo:0000828 | secondary_low_context | 170 | 0 | 0 | 0.000 | 0 | 19 | 82 | 0 | False |
| mondo:0007088 | secondary_evidence | 27 | 0 | 0 | 0.000 | 0 | 0 | 5 | 0 | False |
| mondo:0007669 | primary | 258 | 0 | 0 | 0.000 | 0 | 0 | 4 | 0 | False |
| mondo:0008199 | primary | 57 | 0 | 0 | 0.000 | 6 | 36 | 15 | 0 | False |
| mondo:0010196 | primary | 730 | 0 | 0 | 0.000 | 0 | 3 | 79 | 0 | False |
| mondo:0010785 | primary | 96 | 0 | 0 | 0.000 | 0 | 0 | 15 | 0 | False |
| mondo:0011150 | secondary_evidence | 24 | 0 | 0 | 0.000 | 0 | 0 | 0 | 0 | False |
| mondo:0011913 | primary | 68 | 0 | 0 | 0.000 | 0 | 0 | 11 | 0 | False |
| mondo:0012192 | primary | 281 | 0 | 0 | 0.000 | 3 | 41 | 236 | 0 | False |
| mondo:0012480 | secondary_low_context | 178 | 0 | 0 | 0.000 | 0 | 1 | 8 | 0 | False |
| mondo:0012522 | primary | 251 | 0 | 0 | 0.000 | 0 | 1 | 16 | 0 | False |
| mondo:0012853 | primary | 201 | 0 | 0 | 0.000 | 0 | 0 | 0 | 201 | True |
| mondo:0013523 | secondary_evidence | 32 | 0 | 0 | 0.000 | 0 | 8 | 24 | 0 | False |
| mondo:0014157 | primary | 52 | 0 | 0 | 0.000 | 0 | 0 | 0 | 0 | False |
| mondo:0014527 | secondary_evidence | 26 | 0 | 0 | 0.000 | 0 | 0 | 0 | 26 | True |
| mondo:0015140 | secondary_low_context | 144 | 0 | 0 | 0.000 | 0 | 19 | 71 | 0 | False |
| mondo:0017279 | secondary_low_context | 188 | 0 | 0 | 0.000 | 22 | 91 | 75 | 0 | False |
| mondo:0018203 | secondary_low_context | 34 | 0 | 0 | 0.000 | 0 | 0 | 0 | 0 | False |
| mondo:0018911 | secondary_low_context | 319 | 0 | 0 | 0.000 | 0 | 3 | 20 | 0 | False |
| mondo:0019321 | secondary_low_context | 25 | 0 | 0 | 0.000 | 0 | 0 | 2 | 0 | False |
| mondo:0020525 | secondary_low_context | 498 | 0 | 0 | 0.000 | 3 | 28 | 464 | 0 | False |
| mondo:0030089 | secondary_low_context | 394 | 0 | 0 | 0.000 | 0 | 3 | 20 | 0 | False |
| mondo:0100087 | secondary_low_context | 336 | 0 | 0 | 0.000 | 35 | 201 | 100 | 0 | False |
| mondo:0100164 | primary | 459 | 0 | 0 | 0.000 | 2 | 23 | 433 | 0 | False |
| mondo:0100165 | primary | 329 | 0 | 0 | 0.000 | 0 | 3 | 28 | 0 | False |
| mondo:0957576 | primary | 127 | 0 | 0 | 0.000 | 8 | 59 | 60 | 0 | False |
| orphanet:75496 | secondary_evidence | 41 | 0 | 0 | 0.000 | 0 | 0 | 0 | 41 | False |

Degenerate units: 2 / 29 (excluded from primary gate, reported separately).

Tiny smoke test: ~10-node graph baseline MRR <1 asserted in unit tests.

