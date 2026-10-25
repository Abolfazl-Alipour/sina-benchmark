<!-- track_scope: benchmark -->
# Track-D (OptimusKG) characterization — M1 (spec §5.2/§5.5)

**Status: DRAFT — reviewer-auditable evidence (Gate 1)**

Source: stratified per-type Parquet, Dataverse doi:10.7910/DVN/IYNGEV v2.0
(37 files, sha256 each in `data/raw/REGISTRY.json`; downloaded via pinned
`optimuskg` 1.0.0 client).

## Canonical tables

- `data/processed/m1/trackd/nodes.parquet` — **190,531** nodes,
  10 types, CURIE `node_id`; MONDO-native diseases: 11,724
  (crosswalk gap documented in `CHANGES/0008`)
- `data/processed/m1/trackd/edges.parquet` — **21,813,816** edges,
  `label`/`relation`/`undirected`, provenance `sources_direct`/`sources_indirect`,
  properties superset struct preserved; **21,813,816** evidence groups
  (§5.2.3: order-insensitive only when `undirected=True`)
- `data/processed/m1/gene_properties.parquet` — per-gene property table
  (§5.2.4): tractability, safety_liabilities, target_class, constraint_scores,
  homologues, chemical_probes. NOTE: structure/pocket/druggable flags are NOT
  present in OptimusKG v2.0 — see `CHANGES/0008`.

## Node types

| node_type | nodes |
|---|---|
| anatomy | 13,120 |
| biological_process | 25,754 |
| cellular_component | 4,052 |
| disease | 36,345 |
| drug | 16,766 |
| exposure | 881 |
| gene | 61,306 |
| molecular_function | 10,161 |
| pathway | 2,805 |
| phenotype | 19,341 |

## Edge labels

| label | edges |
|---|---|
| ANA-ANA | 17,082 |
| ANA-GEN | 8,787,955 |
| BPO-BPO | 44,494 |
| BPO-GEN | 158,410 |
| CCO-CCO | 4,639 |
| CCO-GEN | 105,309 |
| DIS-DIS | 44,215 |
| DIS-GEN | 9,734,774 |
| DIS-PHE | 157,144 |
| DRG-BPO | 62 |
| DRG-DIS | 70,380 |
| DRG-DRG | 1,345,376 |
| DRG-GEN | 20,694 |
| DRG-PHE | 13,758 |
| EXP-BPO | 2,260 |
| EXP-CCO | 13 |
| EXP-DIS | 2,391 |
| EXP-EXP | 2,443 |
| EXP-GEN | 2,989 |
| EXP-MFN | 47 |
| GEN-GEN | 327,924 |
| MFN-GEN | 90,933 |
| MFN-MFN | 12,587 |
| PHE-GEN | 793,279 |
| PHE-PHE | 24,862 |
| PWY-GEN | 46,977 |
| PWY-PWY | 2,819 |

## Exact count reconciliation (§5.5 — must match exactly)

- nodes: 190,531 vs published 190,531 — **EXACT**
- edges: 21,813,816 vs published 21,813,816 — **EXACT**
- FK integrity: all edge endpoints present in nodes.parquet — **PASS**

## DIS-GEN temporal coverage (Appendix D assumption 5)

- DIS-GEN edges: 9,734,774
- with BOTH `year_initial` and `year_final`: 23,987 (0.2%)
- with either year field: 23,987 (0.2%)
- threshold for a primary temporal split: 50%
  → temporal split is **DOWNGRADED to secondary**
