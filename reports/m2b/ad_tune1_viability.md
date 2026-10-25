<!-- track_scope: benchmark -->
# tune1 record — terminated early (verdict moot after structural diagnosis)

**Date:** 2026-08-01
**Run:** AD seed 1338, lr 1.0e-3, batch_per_epoch 2000, batch 4,
gradient_interval 4 — the anti-overfit tuning experiment after the AD-1337
gate FAIL (`reports/m2b/ad_seed1_viability.md`).

**Outcome:** terminated after epoch 2 of 4 (best-valid checkpoint epoch 1)
when the failure diagnosis (`reports/m2b/ad_failure_diagnosis.md`) showed the
gate could not pass for structural reasons.

**Trajectory:** valid MRR 0.00949 (epoch 0) → 0.01005 (epoch 1) → 0.00667
(epoch 2). lr 1e-3 prevented the immediate collapse seen at lr 5e-3
(0.0179→0.0072) but the level stayed ≈ degree-matched null, far below the
degree baseline (0.0695). Consistent with the diagnosis: with all 6 held-out
diseases orphaned in the BRG, no optimizer setting can recover graph signal
that is not there.

**Status of artifacts:** checkpoints for epochs 1–2 exist under
`third_party/biopathnet/experiments/.../2026-07-31-14-19*` (persistent
volume); not uploaded to the release (superseded by ad_fix1).
**Verdict:** FAIL by trajectory; formally moot — superseded by the CHANGES/0013
export fix and the ad_fix1 rerun.
