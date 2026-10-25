# 0013 — Task-G export: reconnect held-out diseases in the BRG (low-leakage context edges only)

- **change:** In `scripts/m2_export_taskg_biopathnet.py`, stop stripping ALL
  BRG edges that touch held-out family diseases. Retain low-leakage context
  edges for held-out diseases: **DIS-PHE, DIS-DIS, and all hierarchy
  relations**; continue excluding **DRG-DIS and EXP-DIS** for held-out
  diseases (leakage-conservative: drug–disease indication is a 2-hop proxy of
  the DIS-GEN supervision target via DRG-GEN). DIS-GEN supervision edges for
  held-out diseases remain excluded as before (enforced by the split;
  train1 never contains the supervised label).
- **rationale:** Measured root cause of the AD gate FAIL
  (`reports/m2b/ad_failure_diagnosis.md`): all 6 held-out AD diseases have
  degree 0 in the BRG, so NBFNet cannot propagate from the query node and
  scores at degree-matched-null level regardless of training. BioPathNet's
  published PrimeKG protocol keeps held-out diseases connected via exactly
  these context relations.
- **review record (owner-required 3×):** (1) leakage: DRG-DIS/EXP-DIS kept
  excluded as close proxies; DIS-DIS hierarchy paths to non-family parents
  with train supervision are legitimate context (mirrors the published
  protocol) and are disclosed in reports; (2) protocol fidelity: retention
  set ⊂ BioPathNet's own PrimeKG context set, no unfair advantage;
  (3) experiment design: one AD rerun decides winnability before any
  switch to the KGPFLN arm.
- **requester:** owner 2026-08-01 ("i approve the new plan, only if you
  review it 3 times").
- **impact:** Export rebuild only (splits, supervision, leakage guards
  unchanged); ~30 min rebuild, then one AD seed rerun (~13 h) decides whether
  the zero-shot family exam is winnable for NBFNet on this graph.
- **date:** 2026-08-01
- **approved:** owner — 2026-08-01 (chat), conditional on the recorded 3×
  review; review outcome: APPROVED WITH REFINEMENT (DRG-DIS/EXP-DIS stay
  excluded).
