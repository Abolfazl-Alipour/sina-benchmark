# V2 — OptimusKG facts (recorded as resolved per spec §3.3-V2)

**Status: DRAFT — awaiting owner signature**

Spec source: `SINA_Implementation_Plan_v4_FINAL.md` §3.3 item V2 and §0.0-B
(adopted after verification from primary source per Appendix C).

## Graph facts (published counts)

| Metric | Value |
|---|---|
| Nodes | **190,531** across **10 node types** |
| Edges | **21,813,816** across **26–27 relations** |
| DIS-GEN edges | 9,734,774 (44.6%; sources DisGeNET + Open Targets) |
| ANA-GEN edges | 8,787,955 (Bgee expression — see Task-G curation note §7.4) |
| GEN-GEN edges | 327,924 (PrimeKG-derived) |

## Distribution

- Versioned Harvard Dataverse releases, DOI **10.7910/DVN/IYNGEV**.
- PyPI client `optimuskg` (pinned in `docs/preflight/dependency_matrix.md`).
- Stratified per-type Parquet with Polars-native structs (no JSON-string parsing).

## Carried properties

- DIS-GEN edges: `evidence_score`, `year_initial`, `year_final`,
  specificity/pleiotropy indices.
- Gene nodes: tractability (per modality), safety liabilities, constraint scores,
  homologues (species, confidence, %identity), and per-gene flags
  "structure with ligand / high-quality pocket / druggable family".

## Publication status

- arXiv April 2026, **"In review" (not yet peer-reviewed)** — recorded here as a
  standing caveat for all downstream claims (Appendix D assumption 3/5 gates apply).

## Phase-1 reconciliation targets

- Exact counts: 190,531 nodes / 21,813,816 edges (spec §5.5 — OptimusKG reconciles
  exactly, not within tolerance).
- DIS-GEN temporal-field coverage measured in Phase 1 (Appendix D assumption 5:
  if <50% of DIS-GEN edges carry years, the temporal split downgrades to secondary).

Owner signature: ______________________  Date: __________
