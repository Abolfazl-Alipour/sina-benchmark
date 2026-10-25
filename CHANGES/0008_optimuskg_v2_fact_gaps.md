# 0008 — OptimusKG v2.0 fact gaps vs spec §5.2.2/§5.2.4

- **change:** Two facts recorded in spec V2/§5.2 do not hold for the actual
  downloaded release (Dataverse doi:10.7910/DVN/IYNGEV **v2.0**, 2026-07-24):
  1. **Disease ids are not MONDO-preferred.** Of 36,345 disease nodes only
     11,724 carry MONDO ids natively; the rest are EFO (11,960), OBA (10,997),
     Orphanet, DOID, OTAR, OBI, NCIT, MP, PATO. The node `properties.xrefs`
     field that could supply a MONDO crosswalk is empty in this release.
     Canonical tables therefore keep the native CURIE and add a `mondo_id`
     column populated only for MONDO-native nodes; full DOID/EFO→MONDO
     crosswalking is deferred to Phase 2 (`configs/disease_families.yaml`
     workstream) using an external MONDO ontology mapping.
  2. **No per-gene pocket/druggable flags.** The spec §5.2.4 property list
     "structure with ligand / high-quality pocket / druggable family" is NOT
     present in v2.0 gene nodes (fields verified absent). `gene_properties.parquet`
     carries what exists: tractability (per modality), safety_liabilities,
     target_class, constraint_scores, homologues, chemical_probes.
- **rationale:** Record facts as measured (spec §0.2.8 "measured, not estimated");
  canonicalization and property extraction proceed with what the release
  actually contains.
- **requester:** agent (measured during M2 ingest).
- **impact:** Phase-2 MONDO-family holdout needs an external MONDO crosswalk
  (ontology download) — small added workstream. Phase-5 docking suitability
  gate (§9.1.3) loses the OptimusKG-assisted flags and must rely on
  PDB/AF + fpocket + binding-site conservation alone, until a future OptimusKG
  release restores them. No M1 gate criterion changes.
- **date:** 2026-07-24
- **approved:** owner — 2026-07-24 (chat: "use it all, use the full knowledge graph")
