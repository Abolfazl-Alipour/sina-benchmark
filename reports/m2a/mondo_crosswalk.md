<!-- track_scope: benchmark -->
# MONDO crosswalk + disease families — M2A-1 (spec §6.2.2)

**Status: DRAFT — disease families await owner sign-off (§0.3.9)**

MONDO release: `releases/2026-07-06` (`data/raw/mondo/mondo.obo`, sha256 in
REGISTRY). Parser: `src/sina/obo.py` (unit-tested on fixture OBO). All MONDO
ids below use the canonical node_id form (`mondo:…`).

## Crosswalk coverage (closes the CHANGES/0008 gap)

- Track-D disease nodes: **36,345**
- MONDO-native (no mapping needed): 11,724
- MONDO after crosswalk: **15,372
  (42.3%)**
- Unmapped: 20,973 — dominated by **OBA** (Ontology of
  Biological Attributes; these are attribute terms, not diseases — they can
  never map) and **EFO** terms without a MONDO xref in this release. The
  mappable disease universe is therefore ~25k nodes, not 36k.

Per-prefix coverage:

| prefix | Track-D nodes | MONDO-mapped | coverage |
|---|---|---|---|
| EFO | 11,960 | 2,384 | 19.9% |
| MONDO | 11,724 | 11,724 | 100.0% |
| OBA | 10,997 | 0 | 0.0% |
| Orphanet | 1,637 | 1,258 | 76.8% |
| OTAR | 9 | 0 | 0.0% |
| DOID | 6 | 6 | 100.0% |
| OBI | 4 | 0 | 0.0% |
| NCIT | 3 | 0 | 0.0% |
| MP | 2 | 0 | 0.0% |
| PATO | 2 | 0 | 0.0% |
| OGMS | 1 | 0 | 0.0% |

## Disease families (DRAFT — owner sign-off required)

| family | root | members (root + is_a descendants) |
|---|---|---|
| alzheimer | mondo:0004975 (Alzheimer disease) | 22 |
| parkinson | mondo:0005180 (Parkinson disease) | 25 |

Track-D presence of each family (MONDO-mapped nodes inside the family):

| family | family size | present in Track-D |
|---|---|---|
| alzheimer | 22 | 6 |
| parkinson | 25 | 11 |

Config: `configs/disease_families.yaml`. Splits (M2A-2) hold out evidence
groups touching these families; AD and PD are separate holdouts (§6.2.2).
**Limitation for sign-off:** diseases without a MONDO mapping (EFO/OBA-only)
cannot be family-attributed; the owner may expand family membership at
sign-off (e.g., name-matched EFO terms) — any edit re-triggers M2A-2.
