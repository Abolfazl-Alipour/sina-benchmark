# 0003 — Heavy phase-specific dependencies deferred

- **change:** Spec §3 item 1 lists pinned binaries/deps as Phase-0 scope
  (torch, PyKEEN, RDKit, AnyBURL, smina/Vina 1.2, OpenBabel, fpocket,
  `dimorphite_dl`). These are **not installed** in the conda env yet; each enters
  `environment.yml` (pinned) when its phase starts. The dependency-compatibility
  matrix (`docs/preflight/dependency_matrix.md`) records their pinned versions,
  URLs, and SHA-256 where the publisher provides them, as spec §3 item 1 requires.
- **rationale:** Keeps the Phase 0/1 env small and reproducible; nothing in
  Phases 0–1 imports these packages. Installing unused multi-GB wheels (torch) at
  scaffold time wastes disk and lock churn.
- **requester:** agent (2026-07-24).
- **impact:** The matrix is the single source of truth for the deferred pins; when a
  phase starts, `environment.yml` is updated to the recorded pin and the lock file
  is regenerated. No gate criterion changes: Gate 0 checks are V-notes signed,
  container builds (CI), pytest green.
- **date:** 2026-07-24
- **approved:** owner — 2026-07-24
