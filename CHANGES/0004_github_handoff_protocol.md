# 0004 — GitHub + HANDOFF continuous-commit protocol

- **change:** The spec does not prescribe a VCS workflow. Adopted: a private GitHub
  repo (`Abolfazl-Alipour/sina-platform`, branch `main`) with continuous commits +
  pushes authorized by the owner as a standing workflow; `HANDOFF.md` as a living
  takeover document that must be updated in the same commit as any `src/` or
  `scripts/` change (enforced by the CI `handoff-freshness` job); CI secret-scan
  blocks credentials from entering git.
- **rationale:** Owner request (2026-07-24) for durable, auditable progress and
  agent-independent takeover capability.
- **requester:** owner.
- **impact:** Adds two CI jobs (secret-scan, handoff-freshness) and the HANDOFF
  maintenance discipline; token hygiene rule added (fine-grained PAT scoped to this
  repo only; see HANDOFF §6.1). No spec content changes.
- **date:** 2026-07-24
- **approved:** owner — 2026-07-24
