<!--- track_scope: benchmark -->
# Organoid Handoff Spec Stub — CellCircuit (pre-Gate-0, 2026-08-12 Directives v2 H7)

> **Do not ship to Sadegh/CellCircuit until the wet-lab readout is defined in writing** (p-tau / α-syn / DNA-damage foci per disease panel — a target without a readout is not a validation). This stub defines the handoff shape only.

## Handoff shape (only after primary gate passes, §10)

- **Retrain once on all 678 diseases** (no holdout), query the target disease (e.g., late-onset PD `mondo:0008199` or AD3 `mondo:0011913`), rank all 17,448 GEN candidates, remove known associations.
- **Handoff = top-10 novel genes per unit** (not 20–50 — wet-lab capacity is the constraint), each with:
  - `score_delta ≥ 2× degree-baseline` (model MRR vs degree MRR on that unit)
  - `context_lines ≥ 3` (train1 lines touching the disease — gate H6)
  - Monogenic-gene rank from `m2_monogenic_rank.py` (WRN/LMNA panel, descriptive)
- **GWAS/PRS intersection:** reported as a **flag, not a gate** — a candidate outside known heritability loci is marked *risk* but not rejected. Benchmark measures knowledge-graph reconstruction; external-evidence gating would evaluate a different system. Future GWAS-layer arm is candidate future evaluated arm.
- **Lines recommended:** sporadic/idiopathic iPSC lines alongside monogenic patient lines (monogenic-island success is evidence of engine transfer, NOT proof of polygenic relevance — gap stated in every handoff document).
- **Bridge language (claim ceiling):** GEN = gene symbol for gene+protein product; drugs bind proteins. Handoff may say "candidate target prioritization" only with bridge: ranked gene → its protein is candidate target → cheap post-hoc filter (druggability, direction-of-effect e.g., WRN needs *restored*, modality). Never "validated target" / "causal gene identification" / "drug candidate".
- **Stop condition:** no handoff until readout paragraph exists and is committed.

## Placeholder (fill after Wave 1)

| unit | gene | monogenic rank | deg | context | GWAS flag | readout |
|---|---|---|---|---|---|---|
| TBD | — | — | — | — | flag | — |

