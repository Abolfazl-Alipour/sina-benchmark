# 0022 — Light Gate 0 amend (defer PrimeKG two-step, keep VRAM/permutation/resume)

- **change:** Amend v5.1 §7 Gate 0 for Light Gate 0: PrimeKG two-step reproduction (torch 2.0.1/py3.10 then 2.13.0/py3.12.13 ±10% per metric) is **deferred to supplement/post-Arms** and removed from the critical path. The **mandatory pre-spend Light Gate 0** is: (a) single-batch VRAM/throughput ladder probe on the correct 1,431,945-edge BRG (1,420,095 train1) across 300→280→250→200→150→100 (pre-registered, inside Gate 0 budget ~0.5 GPU-h), (b) degree-preserving permutation negative control on 2–3 LODO units (collapse to null), (c) preemption-resume drill via `sina.manifests.Run`. If ladder pass-time >2×70min assumption, Arms re-priced before spend (v5.1 §6). Wired scripts `m1_reproduce_biopathnet.py` / `pod/run_trimmed_long.py` remain; PrimeKG configs stay wired for deferred run.
- **rationale:** Organoid-partner target prioritization on 13 primary LODO diseases (1 seed) is the grant-critical path; 1,853 `INDICATION` rows are insufficient for direct drug repurposing, and the invalid 295k-era VRAM (142GB at batch 300) must be re-measured on the 1.43M BRG before any Arms. Full two-step is trust-anchor only and can run in parallel without blocking Step 0b (1 vs 2 passes × capped/uncapped pilot) or Arms. User approved Light Gate 0 plan 2026-08-10 (`.agents/plans/2026-08-10-light-gate0.md`).
- **requester:** owner — 2026-08-10
- **impact:** Lifts hard gate of CHANGES/0021 for Light Gate 0 only (items a–c). Step 0b and primary-only Arms (13×1 seed) proceed after Light Gate 0 gate report `reports/m2b/gate0_light.md` is committed. PrimeKG replica remains a priced contingency inside the $900 cap, not deleted.
- **date:** 2026-08-10
- **approved:** owner — 2026-08-10

---

## Light Gate 0 runbook

- **Ladder:** `conda run -n sina python scripts/m1_reproduce_biopathnet.py --config configs/m1_biopathnet_repro.yaml` single-batch forward+backward at each rung; `nvidia-smi` peak VRAM, ms/batch, extrapolated pass-time = ms/batch × (train2 rows / batch). Record max batch that fits H200.
- **Permutation:** degree-preserving DIS-GEN shuffle on 2–3 units (Werner `mondo:0010196`, HGPS, renal cysts), harness end-to-end, filtered MRR ≈ degree baseline.
- **Resume:** kill −9 mid-batch, `--resume` via `Run` manifest, hashes match.
- **Report:** `reports/m2b/gate0_light.md` (`<!-- track_scope: benchmark -->`) with ladder table, chosen batch, re-priced Arms if needed, permutation table, resume log.
