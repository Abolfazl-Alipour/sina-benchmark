# 0007 — SBIR fall due date is September 5, 2026, not September 8

- **change:** Spec §3.3-V5 and the §14.2 checkpoint table use "September 8, 2026"
  as the fall SBIR submission date. Verification against NIH's official standard
  due dates (2026-07-24) shows the standard receipt date is **September 5, 2026**
  (recurring Cycle I date). Checkpoint dates shift three days earlier:
  "Task-D reproduced Aug 25" and "numbers frozen Aug 30" stand; the finalization
  window becomes **Sep 1–5**; the fallback is **January 5, 2027**.
- **rationale:** Fact correction from primary-source verification; submitting
  against a wrong (later) date would miss the real deadline.
- **requester:** agent (found during Phase-0 verification, recorded in
  `docs/preflight/V5_nih_dates.md`).
- **impact:** The September track is strictly tighter than the spec assumed; the
  §14.2 checkpoint-rule logic ("missing a checkpoint is not failure; January is
  the fallback") is unchanged. Grant-writing start (Aug 1) unchanged.
- **date:** 2026-07-24
- **approved:** owner — 2026-07-24
