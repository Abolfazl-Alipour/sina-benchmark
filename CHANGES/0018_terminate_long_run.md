# 0018 — Terminate the trimmed long run after epoch 4 (experiment called; drawing-board reset)

- **change:** End the CHANGES/0015/0017 long run after epoch 4. Epochs 5–20
  are not trained. The multi-seed benchmark program is NOT started (verdict:
  NO-GO, `reports/m2b/trimmed_ad_learning_run.md`). Next step is an owner-led
  redesign ("drawing board"), not more training on the current protocol.
- **rationale:** owner directive 2026-08-06 ("if it is not improving after the
  4th pass, we are fucked and we have to go back to the drawing board;
  document everything, commit, save everything to the synced folder, update
  the handoff"). The validation trajectory justified the call: valid MRR
  0.0162 → 0.0157 → 0.0119 → 0.0109 — monotonic decline crossing below the
  degree-matched null (0.0140) while training loss fell (0.288 → 0.239):
  memorization without transfer, the same signature as both earlier failed
  runs. Final gate on the best checkpoint (epoch 1): test MRR 0.0160 /
  Hits@10 0.0137 vs the pre-registered degree-baseline gate 0.0695 / 0.1447 —
  FAIL; Hits@10 is below even the degree-matched null (0.0297).
- **requester:** owner 2026-08-06.
- **impact:** Total spend on the H200 session ≈ 8 GPU-h (~$21 at $2.6/h).
  All artifacts preserved: logs + epoch 1–4 checkpoints on the owner-synced
  `/workspace` (Google Drive), code + report + this entry on GitHub. This
  supersedes CHANGES/0017's "train all 20 regardless" for this run (0017 was
  itself superseded within hours by the same owner's later directive).
- **date:** 2026-08-06
- **approved:** owner — 2026-08-06 (chat directive)
