# V4 — TxGNN split inventory vs SINA's AD/PD family (verified 2026-07-24)

**Status: DRAFT — awaiting owner signature**

## Inventory (repo `mims-harvard/TxGNN`, branch `main`, `data/disease_files/`)

Nine disease-area split CSVs, verbatim filenames:

| File | Rows (excl. header) |
|---|---|
| `adrenal_gland.csv` | ~8 |
| `anemia.csv` | ~37 |
| `autoimmune.csv` | ~29 |
| `cardiovascular.csv` | ~150 |
| `cell_proliferation.csv` | ~301 |
| `diabetes.csv` | ~4 |
| `mental_health.csv` | ~137 |
| `metabolic_disorder.csv` | ~156 |
| **`neurodigenerative.csv`** | **68** |

⚠️ The neurodegenerative file is spelled **`neurodigenerative.csv`** (upstream
typo). Any code referencing "neurodegenerative" verbatim will not match; record
this in the Phase-2 splitter config.

Note on the spec's pinned commit: `1000aac` ("add disease files", 2023-04-13)
added only 5 of the 9 area files (adrenal_gland, anemia, cardiovascular,
cell_proliferation, mental_health). The other four — including the
neurodegenerative split — landed later; the nine-area inventory above is from
`main` (per §3.3-V4 the nine-area claim is about the current repo state).
Phase-2 pinning will record the exact commit used for the split files.

## Overlap with SINA's AD/PD family (spot check, `neurodigenerative.csv`)

The file holds 68 MONDO disease nodes. String-matched AD/PD rows:

- `juvenile onset Parkinson disease 19A` (MONDO:0014231-ish id 14231)
- `parkinsonian-pyramidal syndrome` (9830)
- `postencephalitic Parkinson disease` (1945)
- `hereditary late onset Parkinson disease` (18466)
- `Alzheimer disease without neurofibrillary tangles` (11401)

**The canonical Alzheimer disease (MONDO:0004975) and Parkinson disease
(MONDO:0005180) nodes are NOT in TxGNN's neurodegenerative held-out set** — only
atypical variants/subtypes are. Implication for the novelty claim (§0.0-A.2):
TxGNN's zero-shot neurodegenerative evaluation does not hold out the canonical
AD/PD nodes; SINA's AD/PD-family holdout (defined via MONDO descendants in
`configs/disease_families.yaml`, Phase 2, owner sign-off) is a distinct, harder
evaluation for these specific diseases. This does not change the novelty
relocation already in the spec (first disease→gene task evaluation + first
OptimusKG benchmark); it sharpens the overlap documentation §6.2.1 requires.

Full family-membership quantification (descendants of MONDO:0004975 / 0005180
vs all nine TxGNN areas) is a Phase-2 deliverable with owner sign-off; this note
records the inventory and the spot-check only.

Owner signature: ______________________  Date: __________
