<!-- track_scope: benchmark -->
# Track-R (PrimeKG) characterization — M1 (spec §5.1/§5.5)

**Status: DRAFT — reviewer-auditable evidence (Gate 1)**

Source: `data/raw/primekg/kg.csv` (Dataverse doi:10.7910/DVN/IXA7BM v2.1,
sha256 in `data/raw/REGISTRY.json`).

## Empirical facts (from the file, not from the paper)

- Delimiter sniffed: `','`; rows: **8,100,498**
- Node types (10): `anatomy`, `biological_process`, `cellular_component`, `disease`, `drug`, `effect/phenotype`, `exposure`, `gene/protein`, `molecular_function`, `pathway`
- Relations (30): `anatomy_anatomy`, `anatomy_protein_absent`, `anatomy_protein_present`, `bioprocess_bioprocess`, `bioprocess_protein`, `cellcomp_cellcomp`, `cellcomp_protein`, `contraindication`, `disease_disease`, `disease_phenotype_negative`, `disease_phenotype_positive`, `disease_protein`, `drug_drug`, `drug_effect`, `drug_protein`, `exposure_bioprocess`, `exposure_cellcomp`, `exposure_disease`, `exposure_exposure`, `exposure_molfunc`, `exposure_protein`, `indication`, `molfunc_molfunc`, `molfunc_protein`, `off-label use`, `pathway_pathway`, `pathway_protein`, `phenotype_phenotype`, `phenotype_protein`, `protein_protein`
- Symmetric relations (reverse-fraction ≥ 0.99):
  `anatomy_anatomy`, `anatomy_protein_absent`, `anatomy_protein_present`, `bioprocess_bioprocess`, `bioprocess_protein`, `cellcomp_cellcomp`, `cellcomp_protein`, `contraindication`, `disease_disease`, `disease_phenotype_negative`, `disease_phenotype_positive`, `disease_protein`, `drug_drug`, `drug_effect`, `drug_protein`, `exposure_bioprocess`, `exposure_cellcomp`, `exposure_disease`, `exposure_exposure`, `exposure_molfunc`, `exposure_protein`, `indication`, `molfunc_molfunc`, `molfunc_protein`, `off-label use`, `pathway_pathway`, `pathway_protein`, `phenotype_phenotype`, `phenotype_protein`, `protein_protein`

## Reciprocal audit

| relation | edges | with reverse twin | fraction |
|---|---|---|---|
| anatomy_anatomy | 28064 | 28064 | 1.0000 |
| anatomy_protein_absent | 39774 | 39774 | 1.0000 |
| anatomy_protein_present | 3036406 | 3036406 | 1.0000 |
| bioprocess_bioprocess | 105772 | 105772 | 1.0000 |
| bioprocess_protein | 289610 | 289610 | 1.0000 |
| cellcomp_cellcomp | 9690 | 9690 | 1.0000 |
| cellcomp_protein | 166804 | 166804 | 1.0000 |
| contraindication | 61350 | 61350 | 1.0000 |
| disease_disease | 64388 | 64388 | 1.0000 |
| disease_phenotype_negative | 2386 | 2386 | 1.0000 |
| disease_phenotype_positive | 300634 | 300634 | 1.0000 |
| disease_protein | 160822 | 160822 | 1.0000 |
| drug_drug | 2672628 | 2672628 | 1.0000 |
| drug_effect | 129568 | 129568 | 1.0000 |
| drug_protein | 50936 | 50936 | 1.0000 |
| exposure_bioprocess | 3250 | 3250 | 1.0000 |
| exposure_cellcomp | 20 | 20 | 1.0000 |
| exposure_disease | 4608 | 4608 | 1.0000 |
| exposure_exposure | 4140 | 4140 | 1.0000 |
| exposure_molfunc | 90 | 90 | 1.0000 |
| exposure_protein | 2424 | 2424 | 1.0000 |
| indication | 18776 | 18776 | 1.0000 |
| molfunc_molfunc | 27148 | 27148 | 1.0000 |
| molfunc_protein | 139060 | 139060 | 1.0000 |
| off-label use | 5136 | 5136 | 1.0000 |
| pathway_pathway | 5070 | 5070 | 1.0000 |
| pathway_protein | 85292 | 85292 | 1.0000 |
| phenotype_phenotype | 37472 | 37472 | 1.0000 |
| phenotype_protein | 6660 | 6660 | 1.0000 |
| protein_protein | 642150 | 642150 | 1.0000 |

## Canonical tables

- `data/processed/m1/trackr/nodes.parquet` — **129,312** nodes
  (`node_id` CURIEs: `ncbigene:` / `mondo:` / `drugbank:` / …)
- `data/processed/m1/trackr/edges.parquet` — **8,100,498** edges with
  `evidence_group_id` (reciprocal twins + duplicates share one ID, spec §5.4)

## Count reconciliation (published: 129,375 nodes / 8,100,498
directed edges; tolerance ±1%)

- nodes: 129,312 (deviation 0.0487%)
- edges: 8,100,498 (deviation 0.0000%)
- FK integrity: all edge endpoints present in `nodes.parquet` — **PASS**
