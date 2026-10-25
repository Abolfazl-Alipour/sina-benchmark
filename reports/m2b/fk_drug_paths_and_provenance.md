<!-- track_scope: benchmark -->
# Drug-path measurement + license provenance (v5.1 Step 3 — measurement only)

CPU measurement feeding the `CHANGES/0020` freeze package. Claim ceiling: all
outputs support **disease-associated gene recovery / hypothesis
prioritization** within the neuro/aging/metabolic cluster only.

## Methods summary

- **F@k (v5.1 §6.4), k = [2, 3, 4]:** per candidate unit, fraction of the
  unit's DIS-GEN answer genes reachable from the unit disease via drug-mediated
  paths of length <= k in the LODO fact graph (BRG + supervision *minus the
  unit's own edges*, undirected). Drug-mediated distance = 1 + BFS distance from
  any drug adjacent to the unit disease (first hop is a drug by construction;
  k=2 → disease→drug→gene exactly). Drug nodes = Track-D `node_type == "drug"`
  (chembl:/drugbank: prefixes). Drug legs used: DRG-DIS
  (INDICATION/CONTRAINDICATION/OFF_LABEL_USE) + DRG-GEN
  (TARGET/ENZYME/TRANSPORTER/CARRIER/INHIBITOR/AGONIST/…).
- **Exclusivity annotation:** of the answers drug-reachable@3,
  the fraction ALSO reachable within 3 hops from the disease
  in the same graph with all drug nodes removed — reachable is not exclusively
  drug-mediated.
- **Drug-provenance fraction:** fraction of a unit's DIS-GEN edges whose
  `sources_direct`/`sources_indirect` match the drug/clinical rule: exact tokens
  CTD (Comparative Toxicogenomics Database) and CGI (Cancer Genome Interpreter),
  plus substring tokens DRUG/FDA/EMA/CLINICAL/TRIAL/PHARM/CHEMBL (the latter
  occur nowhere in this build — verified by full-column scan).
- **License provenance (v5.1 §2.3):** DIS-GEN edges bucketed by `sources_direct`
  values: OPEN_TARGETS only / DISGENET+OPEN_TARGETS / DISGENET only / other
  (values are exactly OPEN_TARGETS and DISGENET in this build — verified by
  full-column scan of Track-D).

## 1. F@k per candidate unit (all 48 units)

`n_reach@3` = answers drug-reachable@3; `excl@3` = fraction of those also
reachable <=3 with drug nodes removed; `drug-prov` = drug-provenance fraction
of the unit's key. Full table: `fk_drug_paths.parquet`.

| unit | stratum | answers | drugs adj | F@2 | F@3 | F@4 | n_reach@3 | excl@3 | drug-prov |
|---|---|---|---|---|---|---|---|---|---|
| mondo:0008310 | A_primary | 468 | 4 | 0.011 | 1.000 | 1.000 | 468 | 1.000 | 0.009 |
| mondo:0010196 | A_primary | 857 | 1 | 0.000 | 0.000 | 1.000 | 0 | 0.000 | 0.001 |
| mondo:0012853 | A_primary | 230 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0014157 | A_primary | 60 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.017 |
| mondo:0008311 | A_secondary | 12 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0008694 | A_secondary | 6 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0011150 | A_secondary | 29 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.034 |
| mondo:0012590 | A_secondary | 8 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.125 |
| mondo:0013523 | A_secondary | 35 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.029 |
| mondo:0014527 | A_secondary | 29 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.034 |
| mondo:0018203 | A_secondary | 41 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0019321 | A_secondary | 27 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0030880 | A_secondary | 16 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0957266 | A_secondary | 8 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0957953 | A_secondary | 17 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| orphanet:2963 | A_secondary | 10 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| orphanet:75496 | A_secondary | 45 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0011149 | A_secondary | 0 | 0 | — | — | — | 0 | — | — |
| efo:0005422 | B_aging_trait | 368 | 7 | 0.014 | 0.976 | 1.000 | 359 | 1.000 | 0.000 |
| efo:0007796 | B_aging_trait | 68 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| efo:0009885 | B_aging_trait | 11 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| efo:0022597 | B_aging_trait | 2,150 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| oba:0005494 | B_aging_trait | 150 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| oba:VT0005372 | B_aging_trait | 6 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0004975 | C_alzheimer | 12,414 | 315 | 0.039 | 0.895 | 0.988 | 11109 | 0.997 | 0.007 |
| mondo:0007088 | C_alzheimer | 30 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0011913 | C_alzheimer | 75 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0014265 | C_alzheimer | 19 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.053 |
| mondo:0015140 | C_alzheimer | 167 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0100087 | C_alzheimer | 405 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.074 |
| mondo:0005180 | C_parkinson | 6,584 | 247 | 0.056 | 0.999 | 1.000 | 6579 | 1.000 | 0.012 |
| mondo:0000828 | C_parkinson | 207 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.048 |
| mondo:0008199 | C_parkinson | 66 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.015 |
| mondo:0017279 | C_parkinson | 226 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0957576 | C_parkinson | 151 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0005147 | D_diabetes | 5,664 | 221 | 0.045 | 0.963 | 0.996 | 5454 | 0.999 | 0.007 |
| mondo:0005148 | D_diabetes | 9,197 | 411 | 0.060 | 0.928 | 0.994 | 8536 | 0.997 | 0.018 |
| mondo:0016391 | D_diabetes | 1,161 | 1 | 0.008 | 0.828 | 1.000 | 961 | 0.998 | 0.000 |
| mondo:0007669 | D_diabetes | 307 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.003 |
| mondo:0010785 | D_diabetes | 110 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0012192 | D_diabetes | 335 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.003 |
| mondo:0012480 | D_diabetes | 214 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.005 |
| mondo:0012522 | D_diabetes | 302 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.003 |
| mondo:0018911 | D_diabetes | 389 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0020525 | D_diabetes | 601 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0030089 | D_diabetes | 467 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |
| mondo:0100164 | D_diabetes | 532 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.009 |
| mondo:0100165 | D_diabetes | 381 | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 |

## 2. Exam-integrity check (v5.1 §10 / §6.4)

**roster-median F@3 = 0.000 < 0.5 — below the pharmacological-reconstruction trigger; F@3 enters as a covariate per v5.1 §6.4.** Roster-median F@3 = **0.0%** across 47 units
with supervision. Covariate bands for the Arm A/B analysis (F@3):
**<10%: 40 units**;
**>=10%: 7 units**.

## 3. Top-5 drug-rich units (Arm B subset candidates, v5.1 §6.4)

Ranked by F@3 (ties broken by adjacent-drug count):

| unit | stratum | F@3 | drugs adjacent | answers |
|---|---|---|---|---|
| mondo:0008310 | A_primary | 1.000 | 4 | 468 |
| mondo:0005180 | C_parkinson | 0.999 | 247 | 6,584 |
| efo:0005422 | B_aging_trait | 0.976 | 7 | 368 |
| mondo:0005147 | D_diabetes | 0.963 | 221 | 5,664 |
| mondo:0005148 | D_diabetes | 0.928 | 411 | 9,197 |

## 4. Drug-provenance fractions of answer keys

Per-unit fractions in the §1 table (`drug-prov`). Whole supervision pool:
**0.49%** of 191,953
DIS-GEN edges match the drug/clinical source rule (CTD/CGI tokens, via
`sources_indirect` — no direct drug-source tokens exist in this build).

## 5. License provenance quantification (v5.1 §2.3)

Auditable basis of the research-use posture (DisGeNET-sourced edges carry CC
BY-NC terms → research use only; Open Targets is CC0). Machine-readable:
`license_provenance.parquet`.

Whole supervision pool (191,953 edges): OPEN_TARGETS only
**190,734** (99.36%);
DISGENET+OPEN_TARGETS **1,219**
(0.64% any-DisGeNET); DISGENET only
0; other 0.

Per candidate unit (sorted by any-DisGeNET fraction):

| unit | edges | OT only | DG+OT | DG only | other | any-DG % |
|---|---|---|---|---|---|---|
| mondo:0012590 | 8 | 7 | 1 | 0 | 0 | 12.50% |
| mondo:0100087 | 405 | 374 | 31 | 0 | 0 | 7.65% |
| mondo:0000828 | 207 | 196 | 11 | 0 | 0 | 5.31% |
| mondo:0014265 | 19 | 18 | 1 | 0 | 0 | 5.26% |
| mondo:0019321 | 27 | 26 | 1 | 0 | 0 | 3.70% |
| mondo:0018911 | 389 | 375 | 14 | 0 | 0 | 3.60% |
| mondo:0011150 | 29 | 28 | 1 | 0 | 0 | 3.45% |
| mondo:0014527 | 29 | 28 | 1 | 0 | 0 | 3.45% |
| mondo:0008199 | 66 | 64 | 2 | 0 | 0 | 3.03% |
| mondo:0013523 | 35 | 34 | 1 | 0 | 0 | 2.86% |
| mondo:0018203 | 41 | 40 | 1 | 0 | 0 | 2.44% |
| orphanet:75496 | 45 | 44 | 1 | 0 | 0 | 2.22% |
| mondo:0016391 | 1,161 | 1,136 | 25 | 0 | 0 | 2.15% |
| mondo:0005148 | 9,197 | 9,033 | 164 | 0 | 0 | 1.78% |
| mondo:0014157 | 60 | 59 | 1 | 0 | 0 | 1.67% |
| mondo:0011913 | 75 | 74 | 1 | 0 | 0 | 1.33% |
| mondo:0100164 | 532 | 525 | 7 | 0 | 0 | 1.32% |
| mondo:0005180 | 6,584 | 6,505 | 79 | 0 | 0 | 1.20% |
| mondo:0008310 | 468 | 464 | 4 | 0 | 0 | 0.85% |
| mondo:0005147 | 5,664 | 5,621 | 43 | 0 | 0 | 0.76% |
| mondo:0004975 | 12,414 | 12,320 | 94 | 0 | 0 | 0.76% |
| mondo:0007669 | 307 | 305 | 2 | 0 | 0 | 0.65% |
| mondo:0012480 | 214 | 213 | 1 | 0 | 0 | 0.47% |
| mondo:0012522 | 302 | 301 | 1 | 0 | 0 | 0.33% |
| mondo:0012192 | 335 | 334 | 1 | 0 | 0 | 0.30% |
| mondo:0010196 | 857 | 856 | 1 | 0 | 0 | 0.12% |
| efo:0005422 | 368 | 368 | 0 | 0 | 0 | 0.00% |
| efo:0007796 | 68 | 68 | 0 | 0 | 0 | 0.00% |
| efo:0009885 | 11 | 11 | 0 | 0 | 0 | 0.00% |
| efo:0022597 | 2,150 | 2,150 | 0 | 0 | 0 | 0.00% |
| mondo:0007088 | 30 | 30 | 0 | 0 | 0 | 0.00% |
| mondo:0008311 | 12 | 12 | 0 | 0 | 0 | 0.00% |
| mondo:0008694 | 6 | 6 | 0 | 0 | 0 | 0.00% |
| mondo:0010785 | 110 | 110 | 0 | 0 | 0 | 0.00% |
| mondo:0012853 | 230 | 230 | 0 | 0 | 0 | 0.00% |
| mondo:0015140 | 167 | 167 | 0 | 0 | 0 | 0.00% |
| mondo:0017279 | 226 | 226 | 0 | 0 | 0 | 0.00% |
| mondo:0020525 | 601 | 601 | 0 | 0 | 0 | 0.00% |
| mondo:0030089 | 467 | 467 | 0 | 0 | 0 | 0.00% |
| mondo:0030880 | 16 | 16 | 0 | 0 | 0 | 0.00% |
| mondo:0100165 | 381 | 381 | 0 | 0 | 0 | 0.00% |
| mondo:0957266 | 8 | 8 | 0 | 0 | 0 | 0.00% |
| mondo:0957576 | 151 | 151 | 0 | 0 | 0 | 0.00% |
| mondo:0957953 | 17 | 17 | 0 | 0 | 0 | 0.00% |
| oba:0005494 | 150 | 150 | 0 | 0 | 0 | 0.00% |
| oba:VT0005372 | 6 | 6 | 0 | 0 | 0 | 0.00% |
| orphanet:2963 | 10 | 10 | 0 | 0 | 0 | 0.00% |

## 6. Feeds into the CHANGES/0020 freeze package

1. **Per-unit F@k table** (§1 / `fk_drug_paths.parquet`) — F@2/3/4 + exclusivity
   for all 48 candidate units, published in the frozen roster.
2. **Exam-integrity verdict:** roster-median F@3 = 0.0% vs the
   0.5 trigger — recorded verbatim in the freeze (§2).
3. **Covariate band assignment (F@3):** 40 units <10%,
   7 units >=10% — stratification for the Arm A-B
   delta regression (v5.1 §6.4); full-roster dual-arm contingency unlocks only
   if roster-median F@3 >= 10% (it is 0.0%).
4. **Top-5 drug-rich units** (§3) — Arm B runs on the monogenic panel plus
   these units regardless of F.
5. **License-provenance table** (§5 / `license_provenance.parquet`) — published
   verbatim as the auditable basis of the research-use posture.

## Reproduction

`python scripts/m2_drug_path_measurement.py --config configs/m2_drug_paths.yaml --resume`
