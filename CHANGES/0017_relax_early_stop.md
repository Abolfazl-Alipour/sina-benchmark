# 0017 — Relax early-stop for the trimmed long run: train all 20 epochs regardless of validation trajectory

- **change:** Remove the early-stop brake from the CHANGES/0015 long-run
  protocol. The run no longer stops after 3 consecutive epochs without
  validation MRR improvement, and the "flat at null level by ~epoch 5 → stop"
  kill criterion is suspended. The run trains the **full 20 epochs**
  (full passes, batch 300) regardless of the validation trajectory. The
  20-epoch ceiling, per-epoch validation scoring, per-epoch checkpointing,
  best-checkpoint selection, and the final held-out test gate are UNCHANGED.
  Implemented as `early_stop_patience: 20` in `configs/m2b_trimmed_ad_long.yaml`
  (driver `scripts/pod/run_trimmed_long.py`).
- **rationale:** owner directive 2026-08-06 ("let it learn even with 5
  non-improving epochs"). Owner's judgment: the model should be given the
  full training budget before any stop decision; the validation trajectory is
  still recorded every epoch for the report, it just doesn't gate training.
  (For the record: each epoch is a full pass over all 170,574 supervision
  pairs, so data exposure is not a limiting factor — but the owner call
  stands on giving the model the full budget.)
- **requester:** owner 2026-08-06.
- **impact:** Worst-case cost is now the committed outcome if validation never
  improves: ~23.5 GPU-h (~$74–85, within the $300 cap, spec §12). The
  trajectory table in `reports/m2b/trimmed_ad_learning_run.md` still records
  every epoch, so a flat trajectory is reported honestly at the end rather
  than acted on mid-run. If the old process (patience 3) early-stops before
  this change takes effect, training resumes from the latest checkpoint with
  the brake off and the segments are combined in the report.
- **date:** 2026-08-06
- **approved:** owner — 2026-08-06 (chat directive)
