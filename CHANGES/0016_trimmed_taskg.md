# 0016 — Trimmed Task-G: focused BRG + focused supervision, held-out diseases keep all background context

- **change:** New pipeline (`configs/m2_trimmed_taskg.yaml`,
  `scripts/m2_build_trimmed_taskg.py`) builds the Task-G training data from a
  **trimmed graph**: both the background regulatory graph (BRG) and the
  supervision pool are restricted to the focused cluster — neurodegenerative
  disease (`MONDO:0005559`), parkinsonian disorder (`MONDO:0021095`, required
  because MONDO does not classify PD under neurodegenerative), diabetes
  mellitus (`MONDO:0005015`), aging trait (`efo:0022597`), parental longevity
  (`efo:0007796`), and progeria syndromes. Result: BRG 291,692 edges (from
  1.94M), supervision 183,684 edges (from 9.73M), 653 cluster diseases.
  **Owner's rule:** held-out family diseases keep ALL their background
  context edges in the BRG (including DRG-DIS — see disclosure); only their
  supervision edges are cut, per run (AD run cuts AD family supervision;
  PD run cuts PD family supervision; the other family remains as training
  signal). Supersedes the supervision-pool part of CHANGES/0014 (which
  proposed trimming supervision only).
- **rationale:** owner directive 2026-08-01: "trim down the background
  regulatory network and the supervision set... we are not cutting anything
  off from Alzheimer's... we keep everything [in the background network];
  on the supervision set only we cut the connections to the marked diseases."
  Full-graph training is compute-dead (measured 222 days/epoch on A40,
  `reports/m2b/gpu_smoke_estimates.md`); the trimmed graph makes full-pass
  training practical on a single large-memory card. Keeping held-out diseases
  fully connected in the BRG eliminates the orphan failure mode that killed
  the first AD run (`reports/m2b/ad_failure_diagnosis.md`).
- **requester:** owner 2026-08-01.
- **impact:** Training set is 2 orders of magnitude smaller; per-epoch time
  drops correspondingly. **DISCLOSED LEAKAGE RISK:** DRG-DIS edges touching
  held-out diseases are retained in the BRG (AD: 315, PD: 247); via DRG-GEN
  they form 2-hop disease→drug→gene proxy paths to the supervision target.
  The research memo (`reports/m2b/research_decision_memo.md` §2) recommended
  against this; the owner overrode. Every report from these runs must carry
  the disclosure. Guards at build time: no held-out supervision edges in
  train2; no DIS-GEN in the BRG; every held-out disease in Track-D is
  connected in the BRG (6/6 AD, 8/8 PD — no orphans).
- **date:** 2026-08-01
- **approved:** owner — 2026-08-01 (chat directive)
