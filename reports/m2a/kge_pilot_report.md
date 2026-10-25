<!-- track_scope: benchmark -->
# M2A-3 KGE pilot report (spec §6.3)

**Status:** FAILED — local RSS budget exceeded.

## Attempt

- Script: `scripts/m2_kge_pilot.py --config configs/m2_kge_pilot.yaml`
- Model: TransE-64, 1 epoch, Track-R family holdout train (8,098,244 edges).
- Result: process exceeded the 40 GB peak-RSS budget (peak ~41 GB observed)
  and was killed after 300 s. No scores were produced.

## Interpretation

Full-scale KGE (TransE-256/ComplEx-256/RotatE-512, ≥5 seeds, ≥50 epochs) on
the full Track-R or Track-D training graphs cannot run within the local 40 GB
RSS budget. It must be run on:

- **Cloud CPU** (larger RAM), or
- **Cloud GPU** (Phase 0.5, owner-pending), or
- **A reduced subgraph** (e.g., relation-filtered or degree-filtered) if the
  evaluation question permits.

The PyKEEN wrapper (`src/sina/kge.py`) is unit-tested and ready; the bottleneck
is graph size, not code correctness.

## Next step

Schedule a cloud-CPU/GPU job for full-scale KGE, or define a reduced-graph
evaluation protocol for local piloting.
