<!-- track_scope: benchmark -->
# LODO supervision shuffles — dual negative controls (Directives v2 B3; CHANGES/0022)

Permutation negative control for the Light Gate 0 gate: the DIS-GEN supervision
(train2.txt) of 3 LODO export units is rewritten under two
degree-preserving nulls, one sibling export dir per unit and mode
(`<unit>_shuffle_<mode>/` under `data/processed/m2/lodo_exports`). train1/valid/test,
entity files and strip_list.json are copied unchanged; each new dir carries a
`shuffle_provenance.json` (mode, seed, source unit, sha256 of source and
shuffled train2). Units: `mondo:0010196`, `mondo:0008310`, `mondo:0007669`. Modes: `across_disease`, `within_disease`. Seed: 1337.

- **across_disease (popularity null):** every gene keeps its exact supervision
  degree; its disease edges are reassigned to uniformly sampled diseases
  (roster: pooled). Gene popularity survives; disease pairing is
  destroyed.
- **within_disease (specificity null):** every disease keeps its exact
  supervision degree; its edges are reassigned to uniformly sampled genes of
  the same unit's train2. Disease degree survives; gene identity is destroyed.

Expected result when the harness trains on a shuffled export: filtered MRR on
the held-out key collapses to the degree baseline (~0.01). If it does not, the
pipeline leaks supervision signal and disease-associated gene recovery claims
are blocked until the leak is found (CHANGES/0022 permutation gate).

## Verification (all rows asserted at shuffle time)

| unit | mode | before | after | genes | diseases | marginal | check | dups | disjoint | sha256 |
|---|---|---|---|---|---|---|---|---|---|---|
| mondo:0010196 | across_disease | 59,454 | 59,454 | 7,707 | 644 | gene_degree | PASS | 0 | PASS | `7b42625b1aad` |
| mondo:0010196 | within_disease | 59,454 | 59,454 | 7,704 | 616 | disease_degree | PASS | 0 | PASS | `69e5e6359983` |
| mondo:0008310 | across_disease | 58,021 | 58,021 | 7,571 | 644 | gene_degree | PASS | 0 | PASS | `b6c07626bc77` |
| mondo:0008310 | within_disease | 58,021 | 58,021 | 7,570 | 624 | disease_degree | PASS | 0 | PASS | `7bc64aff2f02` |
| mondo:0007669 | across_disease | 81,624 | 81,624 | 9,907 | 644 | gene_degree | PASS | 0 | PASS | `bd785f689c94` |
| mondo:0007669 | within_disease | 81,624 | 81,624 | 9,905 | 636 | disease_degree | PASS | 0 | PASS | `3c55a2e7e1a1` |

Leakage contract: no duplicate (disease, relation, gene) tuples in any shuffled
train2; shuffled tuples are disjoint from that unit's valid.txt/test.txt
tuples (tuple-level resampling, max 1000 attempts per edge, hard
failure beyond); the preserved marginal matches exactly; re-running with seed
1337 reproduces byte-identical train2.txt files (sorted lines).

## Reproduction

`python scripts/m2_lodo_shuffles.py --config configs/m2_lodo_shuffles_keep.yaml --resume` (KEEP-graph variant; `configs/m2_lodo_shuffles.yaml` remains the DROP-graph original)
