# 0009 — Track-R disease-family holdout maps individual MONDO ids to grouped disease nodes

- **change:** Spec §6.2.2 calls for a disease-family holdout on both tracks. The
  draft `configs/disease_families.yaml` lists individual MONDO ids (e.g.
  `mondo:0004975`). In Track-R, Alzheimer disease and Parkinson disease are
  represented as `mondo_grouped:*` nodes that aggregate many individual MONDO
  terms; the exact individual ids do not exist as disease nodes in Track-R.
- **rationale:** Without this expansion the family holdout on Track-R would
  match zero disease-protein edges and would not test zero-shot recovery. The
  grouped-node form is the Track-R representation of the family phenotype.
- **requester:** agent (implementation necessity, M2A-2).
- **impact:** Track-R family eval set is defined by the pooled grouped nodes that
  contain any family member numeric id, plus any exact-matching individual
  `mondo:*` ids. Family membership still awaits owner sign-off in
  `configs/disease_families.yaml` (§0.3.9).
- **date:** 2026-07-24
- **approved:** owner — 2026-07-24 (chat: "all docs are good")
