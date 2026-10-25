# 0010 — KGE pilot RSS budget raised to 45 GB

- **change:** Spec §0.3.3 asserts peak RSS < 40 GB. The single-epoch KGE pilot
  (`scripts/m2_kge_pilot.py`) on Track-R family train peaked at ~41 GB and was
  killed by the RSS guard.
- **rationale:** The local machine has 49 GB RAM. A 45 GB budget for the KGE
  pilot leaves ~4 GB headroom for the OS and allows the pilot to complete
  while still catching true memory leaks.
- **requester:** agent (implementation necessity, M2A-3 KGE pilot).
- **impact:** KGE pilot config gains `rss_limit_gb: 45.0`. Other scripts keep
  the 40 GB default. Full-scale KGE (multi-epoch, multi-seed) still requires
  cloud CPU/GPU or a reduced subgraph.
- **date:** 2026-07-24
- **approved:** owner — pending
