# Prospective + per-locus probes — temporal track, first results (2026-09-06)

Owner directive 2026-09-06: temporal-split + per-locus reframing. Data: EBI
GWAS Catalog bulk `gwas-catalog-associations_ontology-annotated-split.zip`
(FTP releases/latest; DisGeNET bulk is login-walled). Year-split parts,
genome-wide only (p <= 5e-8). PD: 814 rows -> 137 loci, 265 pinned genes with
first-reported years. AD: 6,716 rows (proxy/biomarker traits dominate) ->
869 loci, 1,330 pinned genes with years. Cutoff 2020 throughout.
Scripts: `m2c_gwas_loci_build.py`, `m2c_prospective_null.py`,
`m2c_prospective_model.py` (all ruff-clean). Inputs committed:
`gwas_loci_{mondo__0005180,mondo__0004975}.json`.

## Probe 0 — prospective null (counting under cutoff semantics)

Post-cutoff genes as test queries (all-traits genetics clock, cutoff 2020).

| unit | n post | FULL (full knowledge) | PRE (strict prospective) |
|---|---|---|---|
| PD | 82 | MRR 0.01238, H@10 0.04878 (4/82) | MRR 0.00009, H@10 0/82 |
| AD | 1086 | MRR 0.00509, H@10 0.00737 (8/1086) | MRR 0.00019, H@10 0/1086 |

PRE_ALL == PRE_DIRECT to 5 decimals on both diseases. Reading: with only
pre-cutoff knowledge, co-occurrence scores every future gene ~0 — counting
has no discovery mechanism, as predicted structurally. Note FULL still
recovers a few (4/82, 8/1086): future genes are partially countable *after*
the fact through pre-cutoff disease links. Raw: `prospective_null_{pd,ad}.json`.

## Probe A — prospective MODEL (PD, no-ANA epoch_001 score vector)

Clean queries = 29 post-cutoff pinned genes in NEITHER train2 NOR valid (the
checkpoint never saw them as positives; 6 are hidden-test genes, 23 unseen
entirely). Standard filtered genome-wide eval:

- MRR 0.00373, H@10 0/29, median rank 5,594 of 11,203 (uniform-random region).

The model also fails prospective discovery — graph features (GO, phenotypes)
do not lift future genes above noise. Discovery is dead on this stack for
both null and model. (AD model probe landed 2026-09-06: Probe A2, same
verdict at n=239.)

## Probe B — per-locus known-gene recovery (PD, 52 multi-gene loci w/ known)

In each GWAS locus (>=2 pinned genes, >=1 known PD positive): is the model's
in-locus top-1 a known gene?

| method | top-1-known rate |
|---|---|
| model | **0.904** (47/52) |
| degree | 0.615 |
| nearest-gene (Catalog distances) | 0.595 (42 loci w/ distances) |
| chance | 0.677 |

Model-degree agreement 0.462, model-nearest agreement 0.286,
MRR-of-known-in-locus 0.613. Reading: the FIRST task where the model beats
degree and chance by a wide margin, and not by copying either (low
agreements). Caveat, stated plainly: the known genes are train positives —
this is recovery, not discovery. It proves in-locus discrimination among
seen genes, which is the machinery a per-locus product needs; the
novel-gene version needs fine-mapping gold standards (not in-repo).

## Probe A2 — prospective MODEL (AD, no-ANA epoch_001 score vector, 2026-09-06)

Same construction, AD loci (cutoff 2020, seed 1337, bootstrap 1000). Raw:
`prospective_model_ad.json`.

Clean queries = 239 post-cutoff pinned genes in NEITHER train2 NOR valid (847
post-cutoff genes were seen as positives):

- MRR 0.00300 [0.00096, 0.00661], H@10 0.00837 [0, 0.02092] (2/239), median
  rank 2,283 of 5,655.

Same verdict as PD at 8× the query count: no prospective discovery. (AD's
median looks less bad than PD's 5,594 only because the filtered universe is
smaller — 5,655 vs 11,203 — not because ranks are better; MRR is identical
within noise: 0.00300 vs 0.00373.)

## Probe B2 — per-locus known-gene recovery (AD, 239 multi-gene loci w/ known)

| method | top-1-known rate |
|---|---|
| model | **0.954** (228/239) |
| degree | 0.816 |
| nearest-gene (Catalog distances) | 0.780 (177 loci w/ distances) |
| chance | 0.734 |

Model-degree agreement 0.695, model-nearest agreement 0.446,
MRR-of-known-in-locus 0.582. Recovery replicates STRONGER than PD (0.954 vs
0.904) — but so does every baseline (degree 0.816 vs 0.615, chance 0.734 vs
0.677: AD loci are smaller/sparser, so popularity and luck go further). Margin
over the best baseline is actually thinner than PD (+0.138 over degree vs
PD's +0.289). Same caveat: knowns are train positives — recovery, not
discovery.

## What this decides

- Temporal: null collapse under cutoff semantics is total (0/82, 0/1086);
  model collapse on clean future genes is total (0/29). Genome-wide
  discovery from this supervision is certified dead from a second,
  independent direction. No retrain can fix it — the information is absent,
  not unlearned.
- Per-locus: alive as recovery (PD 0.904 vs 0.677 chance; AD 0.954 vs 0.734
  chance, thinner margin over degree). Next steps: (1) AD model probes — DONE
  2026-09-06 (Probe A2/B2 above); (2) novel-gene locus validation against
  fine-mapping/L2G gold standards (data pull, CPU); (3) product spec for the
  per-locus "of these 5, which one" task if (2) holds.
