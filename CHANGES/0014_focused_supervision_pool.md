# 0014 — Task-G supervision pool focused to neurodegenerative + aging + diabetes branches

- **change:** For Task-G BioPathNet training (spec §7.4), restrict the DIS-GEN
  supervision pool from the global ~9.73M edges (all diseases) to a focused,
  mechanically defined pool of **~154,259 edges**:
  - neurodegenerative disease branch (`MONDO:0005559` descendants, 591 usable
    Track-D diseases, 131,702 edges),
  - aging trait `efo:0022597` (2,150 edges; Open Targets aging-associated
    gene set — disclosed as a *trait*, not a disease),
  - parental longevity `efo:0007796` (68 edges),
  - progeria syndromes (515 edges),
  - diabetes mellitus branch `MONDO:0005015` (26 diseases, 19,824 edges).
  Held-out AD/PD family diseases remain excluded everywhere. Batch subsampling
  (CHANGES/0012) is retained until the big-card full-pass runs; with the
  focused pool, per-epoch coverage rises from 0.04% (global) to ~3%.
- **rationale:** Owner directive (2026-08-01): training signal should be
  relevant to the target diseases (CNS/aging/diabetes) rather than diluted
  across all diseases. The exam, leakage guards, and evaluation protocol are
  untouched — this changes only which (disease, gene) pairs produce gradient
  updates. Per-batch cost is unchanged (full-graph message passing); compute
  cost is unchanged (~119 GPU-h at batch_per_epoch 1000).
- **requester:** owner 2026-08-01.
- **impact:** Task definition becomes "disease-associated gene recovery with
  neurodegeneration/aging-focused training" — disclosed in every benchmark
  report. Cluster membership is mechanical (MONDO descendants + named EFO
  traits), pre-registered here, and reproducible. Near-sibling diseases of the
  held-out families (e.g., other tauopathies) remain in the pool by design
  (same generalization principle as the published TxGNN protocol); a stricter
  exclusion is a separate owner call. Eval honesty analysis in
  `reports/m2b/training_pool_and_hardware_plan.md` §2.
- **date:** 2026-08-01
- **approved:** owner — **PENDING** (drafted for sign-off after ad_fix1 verdict)
