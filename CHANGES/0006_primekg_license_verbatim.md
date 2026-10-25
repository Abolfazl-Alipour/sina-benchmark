# 0006 — PrimeKG license recorded verbatim is CC0 1.0, not CC BY-NC-SA 4.0

- **change:** Spec §3.3-V1 and Appendix A assert PrimeKG's license is
  "CC BY-NC-SA 4.0". Verification against the primary source (Harvard Dataverse
  record, doi:10.7910/DVN/IXA7BM, 2026-07-24) shows the dataset license field is
  **CC0 1.0** with no custom terms. `docs/preflight/V1_licenses.md` records the
  verbatim primary-source value, per the spec's own "record verbatim" rule.
- **rationale:** Fact correction from primary-source verification; the spec
  mandates recording licenses verbatim, and the verbatim value contradicts the
  spec's assertion.
- **requester:** agent (found during Phase-0 verification).
- **impact:** Track-R's redistribution posture is *less* restricted at the
  Dataverse-record level than the spec assumed. However, PrimeKG integrates
  sources with their own restrictions (e.g. DrugBank-derived content), so the
  conservative working assumption — research/benchmark-only for Tracks R/D,
  permissive rebuild (Track-C) for commercial use — is **retained unchanged**
  pending owner sign-off on V1. No architecture, gate, or budget changes.
- **date:** 2026-07-24
- **approved:** owner — 2026-07-24
