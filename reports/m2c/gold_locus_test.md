# Final per-locus test — functional gold standards (VEP coding culprits)

Owner directive 2026-09-06 ("go"): the pre-registered final test — per-locus
top-1 on loci where the culprit is named by lab-grade function, not by our
training labels. Script `scripts/pod/m2c_vep_culprits.py` (ruff-clean).

## Why VEP instead of the curated lists

- OT genetics-gold-standards repo: 1 PD gene + 6 drug-class AD genes. Too
  thin, and drug-class is not causal proof.
- Bulk L2G predictions: ~20GB hash-partitioned parquet (17k parts). No
  study-level slicing without downloading all of it. Dead end.
- VEP coding culprits: lead SNP per GWAS locus (min-p genome-wide row; every
  row carries exactly one rsID — verified) through Ensembl REST VEP; culprit
  loci = EXACTLY ONE pinned gene with a strict-coding consequence
  (missense/LoF). Mirrors L2G's pathogenicity feature and the gold-standard
  "functional" class. 1,006 unique rsIDs queried, responses cached.

## Yield (power rule: pooled novel n>=15 to run the test)

| unit | loci | single-coding loci | novel culprits |
|---|---|---|---|
| PD | 137 | 5 | 2 |
| AD | 869 | 9 | 0 |

Pooled novel n=2 (PD locus022 ENSG00000188906 / rs34637584; PD locus094
ENSG00000127419 / rs34311866). **Power rule FAILS — no test run, no
noise sold.** Per-locus discovery stays UNDECIDED (neither alive nor dead).

Pipeline positive controls (VEP attribution works): rs429358 -> APOE
recovered as the single-coding culprit in BOTH diseases; PD recovers 5/137
loci at coding singletons, the expected rarity for GWAS hits.

## Descriptive anecdotes, PD (n=5, NO statistics)

In-locus rank of the functional culprit, model (no-ANA epoch_001) vs degree:

| locus | n | culprit | novel | model rank | degree rank |
|---|---|---|---|---|---|
| locus022 | 4 | ENSG00000188906 | yes | 1 | 1 |
| locus046 | 3 | ENSG00000177455 | no | 3 (LAST) | 1 |
| locus064 | 2 | APOE | no | 2 (LAST) | 1 |
| locus094 | 4 | ENSG00000127419 | yes | 3 | 3 |
| locus120 | 2 | ENSG00000146409 | no | 2 (LAST) | 2 |

Two notes, both honest. (1) The single model top-1 is a NOVEL culprit
(locus022 — also a Probe-1 genome-wide hit at rank 6 and a hidden-test
gene). Anecdote, not evidence. (2) The misses cut against Probe-B
optimism with precision: Probe B showed the model top-ranks A known gene
90% of the time — but here, twice, the top-ranked known gene is NOT the
functional culprit. The model retrieves the already-known; the causal can
be the gene it ranks last. Known-recovery and causal-discovery are
different tasks, and only the first is proven.

## Descriptive anecdotes, AD (n=9, NO statistics, 2026-09-06)

Same construction with the AD no-ANA epoch_001 score vector (pinned-order
dump asserted; model sort by (-score, gene); degree = train1 edge count).
All 9 AD culprits are already-known genes (0 novel — consistent with the
power-rule failure; this table cannot speak to discovery either way).

| locus | n | culprit | novel | model rank | degree rank |
|---|---|---|---|---|---|
| locus056 | 3 | ENSG00000170382 | no | 3 (LAST) | 3 |
| locus136 | 6 | ENSG00000110077 | no | 2 | 4 |
| locus154 | 2 | ENSG00000087053 | no | 1 | 2 |
| locus277 | 7 | ENSG00000138613 | no | 1 | 5 |
| locus348 | 2 | ENSG00000108798 | no | 1 | 1 |
| locus405 | 49 | ENSG00000130203 (APOE) | no | 2 | 1 |
| locus551 | 2 | ENSG00000144648 | no | 1 | 1 |
| locus676 | 6 | ENSG00000253910 | no | 1 | 6 |
| locus708 | 11 | ENSG00000095970 | no | 1 | 1 |

Model top-1 in 6/9. Two notes. (1) locus405 is the APOE positive control
(rs429358): in a 49-gene locus the model puts APOE 2nd, degree 1st — the
pipeline finds the famous gene but popularity finds it first. (2) The two
clean model-beats-degree cases (locus277: 1 vs 5; locus676: 1 vs 6) are the
same shape as Probe-B recovery, not discovery — every culprit here is a
train positive. The PD warning stands unrebutted: top-ranked known ≠ causal.

## Verdict

- The final test as designed cannot run: functional single-coding novel
  loci are too rare in PD/AD GWAS (2 pooled vs 15 required).
- Per-locus discovery: UNDECIDED. Recovery: alive (0.904). The gap between
  them is now measured, not hypothesized.
- Cheapest next gold (if wanted): QTL-colocalization single-gene loci
  (needs coloc tables — heavy); or accept recovery-grade per-locus product
  scoped to "rank the knowns," which is workflow software, not discovery.

Raw: `vep_culprits_{mondo__0005180,mondo__0004975}.json` (culprit lists,
VEP cache notes); VEP response cache ephemeral in /tmp/vep_cache.json
(re-runnable: cache hits skip the API).
