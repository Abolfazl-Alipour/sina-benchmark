# 0015 — Single-seed long learning run before the multi-seed program

- **change:** Replace the next planned runs (multi-seed, 4 epochs,
  `batch_per_epoch` subsampling per CHANGES/0012) with a **single-seed
  learning-curve experiment on the AD holdout**: 1 seed (1339), the
  CHANGES/0014 focused supervision pool (~154k edges), full passes (no
  subsampling), **up to 20 epochs** on a large-memory GPU (H200 planned),
  per-epoch validation, early-stop (halt after 3 consecutive epochs without
  validation improvement), full held-out AD test at the end. The multi-seed
  benchmark program (5 seeds, mean±SD) is deferred until after this run shows
  the model learns.
- **rationale:** owner directive 2026-08-01 ("1 seed, train the fuck out of
  the model and see if it actually learns"). The subsampled 4-epoch regime was
  too thin to tell learning from noise (0.04% pool coverage/epoch globally;
  ~3% focused). A single long full-pass run is the cheapest decisive test of
  learnability before committing multi-seed budget.
- **requester:** owner 2026-08-01.
- **impact:** Run manifests record epochs, coverage (batches × batch size ≥
  pool size per epoch), validation trajectory, and early-stop events.
  Reported strictly as a learning-curve experiment — no mean±SD claims; the
  multi-seed statistics decision comes after. Estimated cost $270–650
  (H200, est. 4–9 h/epoch) — within the $300 frontier cap (§12).
- **date:** 2026-08-01
- **approved:** owner — 2026-08-01 (chat)
