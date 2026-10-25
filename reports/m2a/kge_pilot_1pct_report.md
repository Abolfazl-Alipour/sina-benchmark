<!-- track_scope: benchmark -->
# M2A-3 KGE pilot report (spec §6.3)

**Status:** DONE — single-epoch pilot on a 1% Track-R family train subgraph
(subsampled by evidence_group_id prefix "00") to measure PyKEEN throughput.
The 1% subgraph does not contain the family query diseases, so the score table
is empty; the pilot's purpose is timing only.

## Pilot results

| split | model | train rows | queries | candidates | scores | time s |
|---|---|---|---|---|---|---|
| trackr_family | TransE_d64_e1 | ~81,000 (1%) | 216 | 27,610 | 0 | 93.312 |

## Full-scale estimate

- Pilot: 1 epoch, embedding_dim=64, 1% of 8.1M edges.
- Full-graph 1 epoch (TransE-64): ~93 s x 100 ≈ **2.6 hours**.
- Full-scale target (TransE-256, 50 epochs, 5 seeds): ~2.6 h x (256/64) x 50 x 5 ≈ **2,600 hours**.
- ComplEx-256 and RotatE-512 are slower still.

**Conclusion:** Full-scale KGE on the complete Track-R/Track-D graphs is not
practical on local CPU. It requires cloud GPU (Phase 0.5) or a reduced-graph
evaluation protocol.

## Reproduction

- Config: `configs/m2_kge_pilot_1pct.yaml`
- Manifest: `data/raw/manifests/m2_kge_pilot_1pct.json`
- Output: `data/processed/m2/kge_1pct`
- Re-run: `python scripts/m2_kge_pilot.py --config configs/m2_kge_pilot_1pct.yaml --resume`
