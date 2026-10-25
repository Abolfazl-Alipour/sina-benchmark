<!-- track_scope: benchmark -->
# M4B tabular baseline — SINA-free LightGBM vs OT L2G

Generated 2026-09-12T05:31:10Z. LightGBM 4.7.0 binary logloss, seeds [1337, 1339, 5027], 5-fold chromosome-disjoint CV. Features: dist_kb, degree, loeuf, pli. Panel: gene TSS within ±500kb of lead SNP (OT convention). Ties pessimistic everywhere.

## Claim corrections (verified 2026-09-12)
FLAMES is Schipper/de Leeuw/Posthuma (VU Amsterdam), not Mountjoy. OT Genetics GraphQL is DNS-dead; OT scores via the Platform credibleSets API. HF model exists (opentargets/locus_to_gene_25.12) with an MIT tag, not Apache-2.0. GTEx v8 = Aug 2019 (not 2020). GWAS Catalog is EBI Terms of Use (some datasets CC-BY-NC) — verify terms before any use beyond research. FLAMES comparison is vs-OT-only in v1 (annotation pipeline deferred, per brief fallback).

## Feature freshness (temporal rule)
TMP view allows ONLY pre-2021-datable features: dist_kb, gnomAD v2 LOEUF/pLI (2018). Excluded: graph degree (current graph), FLAMES/PoPS scores (2023+ sources), eQTL (not yet built). FLAMES loci assumed 2018 (UK Biobank release).

## View CV (all loci, n=1817)
GBM top-1 0.4150 [0.3919,0.4386] (n=1817); distance-only 0.4078 [0.3853,0.4293] (n=1817); OT L2G 0.1723 [0.1552,0.1904] (n=1817). Paired diff GBM-OT {'mean': 0.2427, 'ci': [0.219, 0.2686], 'n': 1817}.
OT auto-misses 1,326 FLAMES (UK Biobank, non-Catalog) loci — the full-set OT number is a coverage artifact, not a fair fight. Fair comparison is the covered subset:

## View CV covered subset (OT-scored loci, n=491, all OT-gold)
GBM 0.5275 [0.4847,0.5723] (n=491); distance 0.5356 [0.4888,0.5784] (n=491); OT 0.5866 [0.5397,0.6273] (n=491). Paired diff GBM-OT {'mean': -0.0591, 'ci': [-0.1079, -0.0041], 'n': 491}. Structural note: OT trained on ALL these gold loci (home advantage); our GBM sees 4/5 folds. OT wins here by ~6pts, CI clear.

## View TMP
cutoff_2020: EMPTY — zero post-2020 loci exist in the positives (gold curation ends 2020; FLAMES loci assumed 2018). The pre-registered temporal comparison cannot fire; verdict N/A. This is itself a finding: neither OT nor FLAMES published one because the test set does not publicly exist.
cutoff_2015 fallback (train<=2015, test>2015, n=1157, features dist_kb, loeuf, pli): GBM 0.3967 [0.3682,0.4270] (n=1157) vs distance 0.3570 [0.3302,0.3855] (n=1157) vs OT 0.1582 [0.1366,0.1798] (n=1157).

## Ablations (CV top-1 deltas vs full 0.4150)
no_dist_kb 0.2658 [0.2455,0.2862] (n=1817) (-0.149); no_degree 0.3368 [0.3148,0.3583] (n=1817) (-0.078); no_constraint 0.4540 [0.4304,0.4766] (n=1817) (+0.039 — constraint features hurt causal-gene ranking; biologically consistent: GWAS causal genes are often unconstrained).

## Verdicts vs frozen rules
REPLICATION FAILED on v1 features: 4-feature GBM ties distance and trails OT by 5.9pts (CI clear) on OT's own gold loci. The FLAMES precedent used rich functional annotations; v1 lacks eQTL/functional features, so this is a features gap, not a modeling gap. Do not proceed to lift claims on v1. UPSTREAM comparison vs FLAMES: deferred to annotation build. Ceiling: gene prioritization at GWAS loci (benchmark track) — no target-discovery or clinical language.

## FLAMES replication (their model, their loci, our metric)
L2G set (n=103): XGB 0.6019 vs OT 0.6893, diff {'mean': -0.0874, 'ci': [-0.165, -0.0097]}; full FLAMES 0.7476, diff {'mean': 0.0583, 'ci': [-0.0194, 0.1359]}. 9-pheno set (n=140): full FLAMES 0.6071 vs OT 0.5429, diff {'mean': 0.0643, 'ci': [0.0, 0.1286]}. Verdict: raw XGB loses to OT clearly; full FLAMES leads by +2 to +6pts, never significant at 95% on top-1. The brief's 'FLAMES beat OT' precedent does not survive their own shipped benchmark tables on this metric.

## Reproduce
`scripts/m4b/tabular_baseline.py positives panels features train --view both` (OT cache reused); this report regenerates via `report`. External sources: OT gold repo, FLAMES repos (Marijn-Schipper), GENCODE v46, gnomAD v2.1.1, Ensembl REST (rsID cache), PubMed esummary (dates), OT Platform API.
