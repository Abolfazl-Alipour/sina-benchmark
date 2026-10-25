<!-- track_scope: benchmark -->
# ANA-CNS Build (D0)

brain_nodes 1325 (full CNS: brain+spinal cord+central nervous system via curated CNS keywords + UBERON:0000955/0001017/0001018 descendants, excludes lens cortex), kept ANA-GEN edges 1,397,674, ubiquitous filtered 0, anatomy coverage 0.101 (1325/13120) informative

kept 1,397,674 TRIPWIRE WARNING: >500k (exceeds 500k guideline) → PROCEEDING per owner request for full brain+CNS (user 2026-08-31 "all related to brain and CNS"), flag for review; <5k? False → PASS

Original BRG 1,431,945 edges (KEEP, 27 relations). train1_ana = 1,431,945 + 1,397,674 = 2,829,619 edges (27 + ANA-GEN, 1,397,674 CNS edges: 60077 distinct genes, max per gene 48 <6560, trackd vocab coverage 100% 60077/60077, BRG 60077 distinct but BRG genes 16912, CNS adds spinal cord + brain; previous 150-node brain-only gave 83,332 edges, this full CNS is 16.8× larger).

Specificity filter: >50% anatomy nodes (6560/13120) removes ubiquitous genes; 0 genes removed (max per gene in CNS 48 <6560, threshold 6560 not triggered; if CNS 50% (662) still 0, max 48 <662).

Gene vocab coverage: trackd 100% (60077/60077 CNS genes in 61108), BRG 60077 distinct but BRG has 16912 genes — CNS will need full-vocab Arm E (17448) to increase BRG coverage beyond trimmed 1,431,945.
Exports: data/processed/m2c/completion_exports_ana/mondo__0004975/ and mondo__0005180/ (train1.txt 1,431,945 → 2,829,619, +1,397,674 ANA-GEN CNS edges, other splits unchanged; previous CNS 83k brain-only superseded).

Baselines on train1_ana CNS 1,397,674 (1 seed, degree/PPR/DWPC proxy): AD 0.00908→0.00908 (degree) / 0.00690→0.00734 (PPR) / 0.01000→0.00987 (DWPC), PD 0.00834→0.00818 / 0.00655→0.00752 / 0.00830→0.00816 — no meaningful lift (Δ<0.001), CNS expression neutral as expected; GNN (Arm E) must beat 0.01 MRR floor. Previous brain-only 83k gave AD 0.00908/0.00690/0.01000, full CNS 1.39M gives similar.
