# 0002 — No Docker daemon on the local machine

- **change:** Spec §3 item 1 requires a Dockerfile and Gate 0 requires "container
  builds". Docker is not installed locally; the container is built by the CI
  `docker` job (GitHub Actions) instead of locally.
- **rationale:** Owner has not installed Docker locally; CI already provides a
  reproducible build environment, which is arguably the more meaningful check.
- **requester:** agent / environment constraint (2026-07-24).
- **impact:** Gate 0 evidence for "container builds" = green CI `docker` job run,
  referenced in the Gate 0/1 evidence package. Owner may optionally install Docker
  for a local check (HANDOFF §6.4); this remains open, not blocking.
- **date:** 2026-07-24
- **approved:** owner — 2026-07-24
