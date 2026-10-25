<!-- track_scope: benchmark -->
# Trimmed Task-G build (CHANGES/0016 — owner directive 2026-08-01)

Trims the background regulatory graph (BRG) AND the supervision pool to the
focused cluster (neurodegenerative + diabetes + aging). Owner's rule applied:
**all background context for held-out diseases is kept; only supervision
(DIS-GEN) edges of held-out families are cut.**

## Cluster

- MONDO descendants of neurodegenerative disease + diabetes mellitus,
  intersected with Track-D: **678 diseases**
- Plus aging traits, parental longevity, 22 progeria terms
  (v5.1 roster completion, CHANGES/0019) and 3 trap
  parent terms (in-cluster, never graded — v5.1 §5.4)
- Supervision pool: **191,953 DIS-GEN edges**; supervision
  genes: 17,448; cluster-linked drugs: 1,270
- Trimmed BRG: **1,431,945 edges** (disease context + gene PPI
  backbone + gene annotations + drug legs; ANA-GEN/EXP-*/DRG-DRG excluded)

## Exports (BioPathNet format)

| family | train1 | train2 | valid | test | held-out connected | DRG-DIS on held-out |
|---|---|---|---|---|---|---|
| trackd_family_alzheimer | 1,431,945 | 178,843 | 1,344 | 11,766 | 6/22 | 315 |
| trackd_family_parkinson | 1,431,945 | 184,707 | 708 | 6,538 | 8/25 | 247 |

## Leakage disclosure (mandatory, owner-directed)

DRG-DIS edges touching held-out family diseases are RETAINED in train1
(count per family above). Via DRG-GEN they form 2-hop disease→drug→gene proxy
paths to the supervision target. This is a known, accepted, disclosed risk
per the owner's 2026-08-01 directive, against the research memo's
recommendation (`reports/m2b/research_decision_memo.md` §2). Any claim from
these runs must carry this disclosure.

## Guards (all PASS at build time)

- No held-out family edge in train2 (supervision cut enforced).
- No DIS-GEN label anywhere in train1 (BRG never contains the supervised label).
- Every held-out family disease has degree ≥ 1 in train1 (no orphaning).

## Reproduction

`python scripts/m2_build_trimmed_taskg.py --config configs/m2_trimmed_taskg.yaml --resume`
