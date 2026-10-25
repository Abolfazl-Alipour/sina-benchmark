<!-- track_scope: benchmark -->
# Phase D Wave 1 — 13 primary LODO units (seed 1337, 1 pass, KEEP graph)

BioPathNet (NBFNet full-graph message passing, CHANGES/0020 pinned config) trained per unit on the keep_disphe exports (train1 = 1,431,945 edges incl. DIS-PHE), batch 32, seed 1337. Run window: 2026-08-16T04:43:59Z -> 2026-08-17T04:50:13Z.

Framing (frozen, CHANGES/0024): Wave 1 is **directional only** — 1 seed, no significance claims. Island-band expectation: MRR 0.01-0.03; beating the same-graph degree/PPR nulls IS the result. This benchmark measures disease-associated gene recovery under leave-one-disease-out.

Caveats: AUROC/AP below are batch-averaged approximations (streaming-eval patch) — MRR/Hits@k are exact and are the only gate metrics. The frozen A3 ≤2-hop triviality prong flags all 13 units on the KEEP graph (triviality_audit_keep.md); stats treatment of that conflict is pending owner decision (CHANGES/0024 §(b)).

| unit | name | test | model MRR | H@10 | H@100 | degree MRR | PPR MRR | lift | ≤2-hop% |
|---|---|---|---|---|---|---|---|---|---|
| mondo:0100164 | Permanent neonatal diabetes mellitus | 459 | 0.01662 | 1.63% | 3.38% | 0.0007 | 0.0144 | 1.16x | 47.9% |
| mondo:0008199 | Parkinson disease 1 | 57 | 0.01425 | 1.75% | 8.77% | 0.0016 | 0.0170 | 0.84x | 73.7% |
| mondo:0011913 | Mandibuloacral dysplasia type A | 68 | 0.01292 | 2.21% | 13.24% | 0.0057 | 0.0157 | 0.82x | 64.7% |
| mondo:0008310 | HGPS | 394 | 0.00918 | 1.14% | 5.20% | 0.0070 | 0.0022 | 1.31x | 64.0% |
| mondo:0007669 | Renal cysts and diabetes syndrome | 258 | 0.00811 | 0.78% | 4.46% | 0.0015 | 0.0137 | 0.59x | 90.7% |
| mondo:0010196 | Werner syndrome | 730 | 0.00788 | 0.55% | 5.07% | 0.0062 | 0.0031 | 1.26x | 83.8% |
| mondo:0012853 | Mandibuloacral dysplasia type B | 201 | 0.00681 | 0.50% | 6.22% | 0.0128 | 0.0017 | 0.53x | 91.0% |
| mondo:0012192 | Permanent neonatal diabetes mellitus with pancreatic and cerebellar agenesis | 281 | 0.00582 | 0.18% | 5.69% | 0.0075 | 0.0009 | 0.78x | 78.6% |
| mondo:0010785 | Maternally-inherited diabetes and deafness | 96 | 0.00511 | 0.52% | 3.12% | 0.0014 | 0.0082 | 0.62x | 77.1% |
| mondo:0014157 | Mandibuloacral dysplasia | 52 | 0.00470 | 0.96% | 4.81% | 0.0018 | 0.0018 | 2.54x | 94.2% |
| mondo:0012522 | Transient neonatal diabetes mellitus 3 | 251 | 0.00455 | 0.20% | 1.00% | 0.0013 | 0.0019 | 2.40x | 45.0% |
| mondo:0100165 | Permanent neonatal diabetes mellitus 1 | 329 | 0.00313 | 0.30% | 1.37% | 0.0007 | 0.0028 | 1.13x | 49.2% |
| mondo:0957576 | Wiedemann-Rautenstrauch syndrome | 127 | 0.00243 | 0.39% | 2.36% | 0.0006 | 0.0031 | 0.78x | 41.7% |

Units evaluated: 13/13 — 6/13 beat the best same-graph null (degree/PPR); mean test MRR 0.00781.

Reproduction: `python scripts/pod/run_wave1_benchmark.py --batch-size 32 --seed 1337` (results + manifest: `reports/m2b/wave1_benchmark_results.json` / `.manifest.json`; checkpoints gitignored, EBS-persisted). Nulls: `reports/m2b/lodo_baselines_keep.md`; negative controls: `reports/m2b/lodo_shuffles_keep.md` + GPU collapse runs (`--unit-set shuffle`).
