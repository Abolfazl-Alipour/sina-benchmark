# 0001 — 12 threads instead of the spec's 14

- **change:** Spec §0.3.3 assumes "14 threads AMD / 50 GB". The execution machine has
  **12 threads** / 49 GB RAM. `torch.set_num_threads(12)` will be used when torch
  enters the env (Phase 2+); all thread-count-sensitive configs use 12.
- **rationale:** Hardware fact; cannot be changed by software.
- **requester:** agent (recorded at env bring-up, 2026-07-24).
- **impact:** Wall-clock timings for local CPU workloads will be somewhat longer than
  spec estimates; spec §0.2.8 ("measured, not estimated") already requires pilot
  timing before every uncertain-cost stage, so measured timings govern. No gate or
  claim is affected. The 40 GB RSS budget (§0.3.3) is unchanged and asserted by
  `sina.manifests.Run`.
- **date:** 2026-07-24
- **approved:** owner — 2026-07-24
