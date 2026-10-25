<!-- track_scope: benchmark -->
# LODO exports (v5.1 Phase F; CHANGES/0020 freeze package)

One self-contained BioPathNet-format export per graded disease unit under
`data/processed/m2/lodo_exports/<unit>/` (unit CURIE with `:` -> `__`).
All roster, cluster, threshold, and strip decisions are the frozen values of
CHANGES/0020; this report confirms the two owner decisions fired as adopted.

- **D1 confirmation:** exactly the 29 graded diseases of the
  three-panel roster are exported (13 primary confirmatory, 5
  secondary-evidence, 11 secondary-low-context — the last including the two
  double-secondary units). Context-only and supervision-only units have no
  export. Nothing above the 5% key ceiling is graded.
- **D2 confirmation:** Kufor-Rakeb (mondo:0011706) and PSP-parkinsonism
  (mondo:0009839) DIS-GEN edges are stripped from train2 whenever a PD unit
  is graded (table below).

## Export composition

train1.txt is the full trimmed BRG **incl. DIS-PHE** (**1,431,945 edges**, CHANGES/0024 C1 KEEP;
no DIS-GEN label anywhere — asserted per export) and is identical
across units: all unit context edges are retained (freeze §6). train2.txt is
the 191,953-edge supervision pool minus the unit's co-holdout cluster keys,
minus ancestor-strip, minus sibling-strip. valid.txt/test.txt are the unit's
own key split 15%/85% by
evidence group ID (seed 1337; reciprocal twins share a group ID
and cannot straddle — disjointness asserted per export).

| unit | name | tier | train2 | valid | test | ctx exp | ctx aud | anc | sib | zero-vec |
|---|---|---|---|---|---|---|---|---|---|---|
| mondo:0007669 | renal cysts and diabetes syndrome | primary | 81,624 | 49 | 258 | 38 | 38 | 21 | 24 | 0 |
| mondo:0008199 | late-onset Parkinson disease | primary | 39,274 | 9 | 57 | 32 | 32 | 9 | 108 | 0 |
| mondo:0008310 | Hutchinson-Gilford progeria syndrome | primary | 58,021 | 73 | 394 | 23 | 23 | 28 | 28 | 1 |
| mondo:0010196 | Werner syndrome | primary | 59,454 | 126 | 730 | 37 | 37 | 28 | 36 | 1 |
| mondo:0010785 | maternally-inherited diabetes and deafness | primary | 66,714 | 14 | 96 | 17 | 17 | 45 | 27 | 0 |
| mondo:0011913 | Alzheimer disease 3 | primary | 41,627 | 7 | 68 | 23 | 23 | 17 | 76 | 0 |
| mondo:0012192 | permanent neonatal diabetes mellitus-pancreatic and cerebellar agenesis | primary | 61,995 | 53 | 281 | 33 | 33 | 29 | 28 | 1 |
| mondo:0012522 | diabetes mellitus transient neonatal 3 | primary | 85,648 | 49 | 251 | 12 | 12 | 29 | 17 | 2 |
| mondo:0012853 | Fontaine progeroid syndrome | primary | 59,970 | 29 | 201 | 98 | 98 | 28 | 27 | 0 |
| mondo:0014157 | mandibular hypoplasia-deafness-progeroid syndrome | primary | 78,480 | 8 | 52 | 29 | 29 | 45 | 20 | 0 |
| mondo:0100164 | permanent neonatal diabetes mellitus | primary | 85,648 | 70 | 459 | 15 | 15 | 29 | 17 | 3 |
| mondo:0100165 | permanent neonatal diabetes mellitus 1 | primary | 85,648 | 50 | 329 | 10 | 10 | 29 | 17 | 2 |
| mondo:0957576 | parkinson disease 25 autosomal recessive | primary | 52,982 | 24 | 127 | 13 | 13 | 9 | 109 | 0 |
| mondo:0007088 | Alzheimer disease type 1 | secondary_evidence | 80,944 | 3 | 27 | 10 | 10 | 17 | 17 | 0 |
| mondo:0011150 | acroosteolysis-keloid-like lesions-premature aging syndrome | secondary_evidence | 143,039 | 5 | 24 | 65 | 65 | 45 | 5 | 0 |
| mondo:0013523 | Nestor-Guillermo progeria syndrome | secondary_evidence | 84,898 | 3 | 32 | 52 | 52 | 45 | 16 | 0 |
| mondo:0014527 | progeroid features-hepatocellular carcinoma predisposition | secondary_evidence | 138,762 | 3 | 26 | 26 | 26 | 45 | 5 | 0 |
| orphanet:75496 | Ehlers-Danlos syndrome progeroid type | secondary_evidence | 81,545 | 4 | 41 | 88 | 88 | 45 | 19 | 0 |
| efo:0005422 | skin aging | secondary_low_context | 58,021 | 43 | 325 | 8 | 8 | 28 | 28 | 0 |
| mondo:0000828 | juvenile-onset Parkinson disease | secondary_low_context | 34,283 | 37 | 170 | 3 | 3 | 9 | 181 | 0 |
| mondo:0012480 | diabetes mellitus transient neonatal 2 | secondary_low_context | 84,952 | 36 | 178 | 4 | 4 | 29 | 19 | 0 |
| mondo:0015140 | early-onset autosomal dominant Alzheimer disease | secondary_low_context | 52,672 | 23 | 144 | 4 | 4 | 17 | 40 | 0 |
| mondo:0017279 | young-onset Parkinson disease | secondary_low_context | 31,958 | 38 | 188 | 2 | 2 | 9 | 185 | 0 |
| mondo:0018203 | LMNA-related cardiocutaneous progeria syndrome | secondary_low_context | 74,633 | 7 | 34 | 2 | 2 | 28 | 19 | 0 |
| mondo:0018911 | maturity-onset diabetes of the young | secondary_low_context | 77,101 | 70 | 319 | 9 | 9 | 21 | 26 | 0 |
| mondo:0019321 | atypical Werner syndrome | secondary_low_context | 117,789 | 2 | 25 | 1 | 1 | 28 | 8 | 0 |
| mondo:0020525 | transient neonatal diabetes mellitus | secondary_low_context | 85,341 | 101 | 498 | 3 | 3 | 29 | 18 | 2 |
| mondo:0030089 | diabetes mellitus permanent neonatal 4 | secondary_low_context | 85,648 | 70 | 394 | 9 | 9 | 29 | 17 | 3 |
| mondo:0100087 | familial Alzheimer disease | secondary_low_context | 57,798 | 69 | 336 | 4 | 4 | 17 | 38 | 0 |

Full per-unit detail: `data/processed/m2/lodo_exports/manifest.parquet` and
each unit's `strip_list.json`.

## D2 — parkinsonism-plus sibling-strip (owner decision, adopted)

Kufor-Rakeb qualifies under the plain overlap rule (coefficient >= 0.3,
intersection >= 20) for all four PD units; PSP-parkinsonism has coefficient
1.0 but a 1-gene intersection, below the 20-gene qualifier, and is stripped
via the explicit D2 list. Both are stripped for every PD graded unit.

| graded unit | stripped disease | coefficient | intersection | source |
|---|---|---|---|---|
| mondo:0000828 | mondo:0011706 (Kufor-Rakeb syndrome) | 1.000 | 197 | overlap_rule+explicit_D2 |
| mondo:0000828 | mondo:0009839 (progressive supranuclear palsy-parkinsonism syndrome) | 1.000 | 1 | explicit_D2 |
| mondo:0008199 | mondo:0011706 (Kufor-Rakeb syndrome) | 0.591 | 39 | overlap_rule+explicit_D2 |
| mondo:0008199 | mondo:0009839 (progressive supranuclear palsy-parkinsonism syndrome) | 1.000 | 1 | explicit_D2 |
| mondo:0017279 | mondo:0011706 (Kufor-Rakeb syndrome) | 1.000 | 197 | overlap_rule+explicit_D2 |
| mondo:0017279 | mondo:0009839 (progressive supranuclear palsy-parkinsonism syndrome) | 1.000 | 1 | explicit_D2 |
| mondo:0957576 | mondo:0011706 (Kufor-Rakeb syndrome) | 0.556 | 84 | overlap_rule+explicit_D2 |
| mondo:0957576 | mondo:0009839 (progressive supranuclear palsy-parkinsonism syndrome) | 0.000 | 0 | explicit_D2 |

## Ancestor-strip guard — Werner trap terms

Grading Werner (cluster c1_progeria held out) strips the three named trap
terms via MONDO ancestors of the held-out members (progeria enters through
HGPS's ancestor list — the guard would abort the build otherwise):

| graded unit | stripped ancestor | ancestor of | supervision edges |
|---|---|---|---|
| mondo:0010196 | mondo:0015333 (progeroid syndrome) | mondo:0008310, mondo:0010196 | 1,589 |
| mondo:0010196 | mondo:0019303 (premature aging syndrome) | mondo:0008310, mondo:0010196, mondo:0018203 | 4,049 |
| mondo:0010196 | mondo:0020732 (progeria) | mondo:0008310 | 706 |

## Strip sweep summary

- Ancestor-strip (MONDO/EFO ancestors of held-out cluster members present as
  Track-D supervision diseases): per unit 9-45
  diseases.
- Sibling-strip (overlap rule + explicit D2): per unit 5-185
  diseases. These lists are long by construction: OptimusKG supervision is
  heavily redundant (the 0.3 coefficient with the 20-gene intersection
  qualifier still means 20+ shared answer genes per flagged disease), and the
  freeze accepts a smaller train2 in exchange for leakage control. Every
  stripped disease is published in the unit's `strip_list.json`.
- Disclose band (0.15-0.3, intersection >= 20; disclosed, NOT stripped —
  freeze §3): 628 unit-disease pairs across the roster, listed per
  unit in `strip_list.json`.

## DIS-PHE connectivity reconciliation

The frozen audit's `brg_context_edges` counted DIS-PHE context and the KEEP export retains DIS-PHE (CHANGES/0024 C1), so export context equals the audit count exactly (asserted per unit — all PASS). Units whose context
drops because of the DIS-PHE removal (0 of 29):

| unit | name | ctx (audit) | DIS-PHE removed | ctx (export) |
|---|---|---|---|---|


Tier floors were frozen against the DIS-PHE-inclusive counts (owner decision
D1); units whose export BRG context falls below their panel's floor are
disclosed here, not re-tiered — the roster is frozen.

## Zero-vector guard

Answer genes with no edge in train1 + train2 (isolated in the fact graph,
unreachable by message passing): **15** across the roster
(excluded from valid/test with the gene list in the unit's `strip_list.json`).

## Ancestor-ID mapping disclosure

Track-D `properties.ancestors` values are mapped `PREFIX_NNNN` ->
`prefix:NNNN` (prefix lowercased, suffix verbatim). Every ancestor value resolved to a Track-D node across all 29 units (0 unmapped). Mapped ancestors without their own DIS-GEN supervision are listed in `strip_list.json` but strip no edges.

## Guards (all PASS at build time)

- No DIS-GEN label in train1; DIS-PHE retained (KEEP variant).
- Each unit's key size equals the frozen audit `key_size`.
- valid/test share zero evidence group IDs (assert_split_disjoint).
- train2 contains no edge of any cluster mate, stripped ancestor, or stripped
  sibling.
- Expected ancestor-strips fire (Werner trap terms; hard abort otherwise).
- Export context count = frozen audit count, per unit (DIS-PHE retained).
- D2 strips present in every PD unit's `strip_list.json` (asserted again when
  this report was written).

## Reproduction

`python scripts/m2_build_lodo_exports.py --config configs/m2_lodo_exports.yaml --resume`
