<!-- track_scope: benchmark -->
# CPU nulls and coverage — ANA saturation, background-only PPR, extended-K curves, supervision-only null (PD, zero training)

Date 2026-09-04. Unit `mondo__0005180` (PD). All four probes ran on CPU on beast;
the live GPU training run was never touched; no repo-tracked file was modified to
produce these numbers. Scripts (scratch, uncommitted): `/tmp/ana_saturation.py`,
`/tmp/gen_coverage.py`, `/tmp/supervision_null.py`, `/tmp/floors_ana.py`.
Pinned universe derived from `entity_types.txt` type token 1 (17,448 genes, verified
count; the `pinned_genes_17448.txt` artifact is absent on this host).

## Input hashes (sha256)

| input | sha256 |
|---|---|
| train1.txt (ANA) | `dff240c6087b304283838a8e9276fc1374c32fb7e342ee929d2de76d6cdab39f` |
| train2.txt (ANA) | `18ce464837ee1e5c43a4a84c3d01f155c341bb5e5d3e36ee81a50d99d12bfeb3` |
| valid.txt (ANA) | `6151e3e9a7200b5a4532aa889078fdc7c0dd9005d49d9f768c72c622794b1ba1` |
| test.txt (ANA) | `18bbbcff754bf60db85620cf5b9d72cc81a3ce1a5e0cb4ac4985190590b8040e` |
| entity_types.txt (base) | `0f887e29588c4cdc249aab680e99af06e8b951b7022803f8ba64f4f83701dae7` |

All metrics below are filtered over the pinned 17,448 (train2+valid positives of PD
excluded from ranking; ties broken by gene token, mirroring the runner and Probe 1).
Bootstrap CIs: 1000 resamples, seed 1337. Test n = 339.

## 1. ANA saturation stat — the anatomy feature is near-constant

ANA edges in `train1.txt` (UBERON-head rows): 1,397,674 over **51** active heads.
Pinned genes with ≥1 ANA edge: **17,290/17,448 = 99.1%**; per-gene median 41,
mean 37.51, max 48. Buckets over pinned: 0: 158 · 1: 12 · 2–5: 43 · 6–20: 2855 ·
21–48: 14380. A feature 99% of candidates carry at near-identical strength cannot
discriminate rank positions; it is background hum the message passing must average out.

## 2. Floors ± ANA — the 1.4M edges move nothing, not even the baselines

Graphs mirror C1.4 (train1+train2), PD test:

| graph | degree MRR [95% CI] | degree H@10 | PPR MRR [95% CI] | PPR H@10 |
|---|---|---:|---:|---:|---:|
| with ANA | 0.00817 [0.00385, 0.01601] | 0.0118 | 0.00752 [0.00268, 0.01476] | 0.00885 |
| without ANA | 0.00831 [0.00374, 0.01499] | 0.0118 | 0.00656 [0.00217, 0.01422] | 0.00885 |

Recomputed with-ANA degree (0.00817) matches the published floor (0.00818),
validating the eval. CIs fully overlap on both metrics: the ANA addition carries no
measurable signal at the floor level. An Arm C ablation (strip ANA, ~20h/epoch)
cannot be justified on these numbers.

## 3. Extended-K coverage — degree still leads to K=1000; PPR passes it only at K=2000

Test coverage (train2+valid excluded), reproduced Probe 1.5 values exactly at
K=100/200/500 before extending:

| generator | K=100 | K=200 | K=500 | K=1000 | K=2000 |
|---|---:|---:|---:|---:|---:|
| degree_train1 | 0.0885 (30) | 0.1504 (51) | 0.2714 (92) | 0.3687 (125) | 0.5251 (178) |
| ppr_train1 | 0.0472 (16) | 0.1180 (40) | 0.2566 (87) | 0.3953 (134) | 0.5959 (202) |
| ppr_train1train2 | 0.0619 (21) | 0.1180 (40) | 0.2448 (83) | 0.3776 (128) | 0.5664 (192) |
| ppr_background_only | 0.0472 (16) | 0.1180 (40) | 0.2566 (87) | 0.3953 (134) | 0.5959 (202) |

Background-only PPR (train1 minus DIS-GEN rows) is **identical** to `ppr_train1` at
every K: train1 contains zero DIS-GEN rows (stripped 0 of 2,829,619), so PPR-on-train1
was already background-only — there is no leakage to strip and no rescue there.
A K=2000 panel with the measured 33–48% rerank conversion implies two-stage recovery
~0.18–0.29, but a 2000-gene panel is 11% of the universe, not a usable shortlist.

## 4. Supervision-only null — counting ties the GPU model in 1.2 seconds

Score(PD,g) = Σ over d′ in D(g) of |G_pd ∩ genes(d′)|, train2 only, no background
graph, no learned weights. PD test, filtered:

| method | MRR [95% CI] | H@10 [95% CI] | cost |
|---|---:|---:|---|
| supervision-only null | 0.01484 [0.00884, 0.02289] | 0.02950 [0.01180, 0.05015] | 1.2 s CPU |
| GNN val ep3 (reference) | 0.01684 | 0.02585 | ~60 GPU-h |
| degree floor | 0.00818 | 0.0118 | seconds |

The 2.8M-edge message passing, aux tasks, and three 20-hour epochs buy ~0.002 MRR
over co-occurrence counting, with CIs overlapping each other and the floor. (Caveat:
null scored on test, GNN numbers are validation — the GNN test is unscored; the
comparison is directional, and its direction favors the null's cheapness.)

## Verdict (graph-proximity part)

Candidate generation from graph proximity is exhausted (best K=100 coverage 8.85%;
K=2000 is not a shortlist); the anatomy feature is near-constant (99.1%) and moves
no floor; the background graph adds ~0.002 over supervision counting on this track.
The remaining direction is outside-graph candidate generation (below), plus
the Probe 2 ranking-loss decision on economics. No further multi-day GPU training is
justified on the current stack.

## 5. Outside-graph generators — Open Targets overall-score panels (PD + AD)

Pulled 2026-09-04 from `https://api.platform.opentargets.org/api/v4/graphql`
(`associatedTargets`, overall score): PD (`MONDO:0005180`) 8,923 targets, AD
(`MONDO:0004975`) 13,367 targets; 8,274 / 12,579 map to the pinned universe.
Panels = top-K by overall score among eligible (train2+valid positives excluded),
tiebreak gene token. Null = OT's own score ordering on the same panel. PD panels
rescored with the saved ep3 checkpoint (one CPU forward; genome-wide cross-check
reproduces the runner log: valid 0.01685 vs logged 0.01684, test 0.0186).

| PD (339 hidden) | K=100 | K=200 | K=500 |
|---|---|---|---|
| OT coverage | 0.2566 (87) | 0.3127 (106) | 0.4159 (141) |
| model two-stage | 0.0295 [0.0147, 0.0472] | 0.0295 [0.0118, 0.0501] | 0.0295 [0.0118, 0.0501] |
| OT-ordering two-stage | 0.0295 [0.0118, 0.0473] | 0.0295 [0.0118, 0.0473] | 0.0295 [0.0118, 0.0473] |
| hit split (both/model-only/null-only) | 4 / 6 / 6 | 4 / 6 / 6 | 4 / 6 / 6 |

| AD (621 hidden, coverage only — no AD checkpoint exists) | K=100 | K=200 | K=500 |
|---|---|---|---|
| OT coverage | 0.0837 (52) | 0.1787 (111) | 0.3784 (235) |

OT covers 3× degree on PD (25.7% vs 8.85%) with the same recall the model gets
(10 hits each, only 4 shared — union 16/339 = 0.047). The reranker adds 6 OT
misses, but the total stays under half the 0.10 bar, and OT top-100s are 92%
already-known positives (92/100 PD; 184/200; 464/500) — the novel-candidate panel
is thin. Circularity asterisk stands: OT shares evidence sources with DIS-GEN
supervision, so part of the 25.7% is the same literature counted twice; the
per-gene source-overlap audit is still open. Panel shas: K=100
`9d89095c…`, K=200 `48eb7860…`, K=500 `c1fd0494…` (full values in probe output).

## 6. Outside-graph generators — GWAS ±500kb window panels (PD + AD)

Lead SNPs from GWAS Catalog REST (`efoTraits/{MONDO_0004975,MONDO_0005180}/associations`,
2026-09-04): PD 814 associations → 602 unique rsIDs (all positioned, 0 failed);
AD 6,716 → 4,033 rsIDs (4,031 positioned). SNP positions via
`singleNucleotidePolymorphisms/{rsID}`. Windows ±500kb overlapped against
OptimusKG `gene.parquet` `genomic_location` (59,757 genes with coordinates);
ranked by distance to nearest lead SNP. PD rescored with ep3 checkpoint as in §5.

| PD (339 hidden) | K=100 | K=200 | K=500 |
|---|---|---|---|
| GWAS coverage | 0.0590 (20) | 0.0708 (24) | 0.0914 (31) |
| model two-stage | 0.0295 [0.0147, 0.0472] | 0.0295 [0.0118, 0.0473] | 0.0295 [0.0118, 0.0501] |
| distance-ordering two-stage | 0.0147 [0.0029, 0.0295] | 0.0147 [0.0029, 0.0295] | 0.0147 [0.0029, 0.0295] |
| hit split (both/model-only/null-only) | 3 / 7 / 2 | 3 / 7 / 2 | 2 / 8 / 3 |

| AD (621 hidden, coverage only) | K=100 | K=200 | K=500 |
|---|---|---|---|
| GWAS coverage | 0.0403 (25) | 0.0853 (53) | 0.1514 (94) |

Positional windows underperform degree on both diseases (PD 5.9% vs 8.85% at
K=100; AD 4.0% vs degree-equivalent). AD's 4,031 loci sprawl over 11,896 pinned
genes (68% of the universe) yet catch 4% of hidden answers at K=100 — proximity
to a lead SNP barely overlaps the literature-curated answer set. The model doubles
the distance null where it matters (both/model-only/null-only 3/7/2: 10 vs 5 hits),
and its 10 hits are the same easy set seen on every panel — the reranker lifts
what is liftable and nothing else. Panel shas: K=100 `f7285811…`, K=200
`cde0b0df…`, K=500 `89809a03…`.

## Verdict (outside-graph part)

OT is the better generator by 3–5× coverage and the only one worth keeping;
GWAS windows go in the bin. But OT + reranker tops out at ~0.03–0.05 two-stage,
still under half the bar, with 92%-known panels and an open circularity audit.
No generator tested yields a product number; the reranker's 10 hits are panel-
independent (same genes everywhere). Combined with §1–4: the discovery track on
this stack is a certified negative with a working reranker attached.

Claim ceiling: disease-associated gene recovery / candidate prioritization only.
