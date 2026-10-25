<!-- track_scope: benchmark -->
# C1 — Candidate universe verification (v6 §8, binding) — 2026-08-31

**Rule:** every loss, ranking, floor, and metric in the completion track uses exactly the pinned 17,448-gene list (CHANGES/0020 §1 artifact). On-disk artifact: `data/processed/m2c/completion_exports/<unit>/entity_types.txt` (original, non-ANA completion exports; type token `1`).

## Asserts (all PASS)
- `len(pinned) == 17,448` — PASS; both units' pinned lists identical.
- train2 ∪ valid ∪ test tails ⊆ pinned list, per disease — PASS (AD train2 distinct tails 17,189; PD 17,446; 0 missing in either unit).
- ANA export entity vocab regenerated (D0 defect fix, commit `a37e248`): 86,157/86,157 entities typed; GEN-type nodes 60,235 (= 17,448 pinned + 42,838 ANA-only loci).

## Biotype / namespace breakdown of the 60,235 GEN nodes (C1.5 blocking assert)
Source: `data/processed/m1/gene_properties.parquet` (61,306 genes; 60,235/60,235 matched).

| biotype | count |
|---|---:|
| protein_coding | 19,999 |
| lncRNA | 15,748 |
| processed_pseudogene | 9,456 |
| misc_RNA | 2,154 |
| unprocessed_pseudogene | 1,921 |
| transcribed_unprocessed_pseudogene | 1,586 |
| (null) | 1,329 |
| transcribed_processed_pseudogene | 1,150 |
| snoRNA | 892 |
| (other small classes) | remainder |

- **Tripwire: protein-coding 19,999 ∈ [15k, 25k] → PASS** (no STOP).
- Namespace: 100% `ensembl:` in both the pinned list and the ANA export — no prefix mismatch → PASS.
- Pinned list composition: 16,073 protein_coding + 684 lncRNA + ~318 pseudogene/miRNA/other — consistent with "protein-coding + annotated subset" → PASS.

## Consequence
- GNN loss/eval restricted to the pinned 17,448 (runner).
- D0 ANA-graph baselines re-run with candidates filtered to the pinned list before scoring (C1.4) — see `completion_baselines_ana.md` (regenerated).
- Junk loci (~42.8k) stay in the message graph, never scored (C5).
