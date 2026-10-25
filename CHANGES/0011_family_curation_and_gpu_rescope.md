# 0011 — PD family curation (parkinsonism-plus excluded) + GPU plan rescoped to OptimusKG-only

- **change 1:** `configs/disease_families.yaml` PD family reduced from 28 to 25
  members. Owner excluded three parkinsonism-*plus* terms — Kufor-Rakeb syndrome
  (`mondo:0011706`), parkinsonism due to ATP13A2 deficiency (`mondo:0017809`),
  parkinsonian-pyramidal syndrome (`mondo:0009830`) — as clinically distinct
  from classic Parkinson disease. Family config now marked owner-APPROVED.
- **change 2:** Cloud GPU plan (spec §7, §12) rescoped per owner directive
  2026-07-25: **OptimusKG Task-G only**, models = **BioPathNet + KGPFLN**
  (the in-context-learning graph foundation model, arXiv 2605.14907, evaluated
  as a benchmarked arm). The PrimeKG Task-D reproduction anchor (spec §7.2)
  and the cloud KGE baseline suite are removed from the GPU plan.
- **rationale:** (1) parkinsonism-plus syndromes dilute the zero-shot PD exam
  and leak less useful signal; (2) owner wants cloud spend concentrated on the
  novel OptimusKG result rather than the reproduction anchor.
- **requester:** owner (chat, 2026-07-25: "let's go with your recommendation"
  + "update the plan for GPU so that we only test optimuskg using our
  biopathnet and the foundational model").
- **impact:** Family splits, baselines, Task-G export, and metrics rebuilt on
  the 25-member PD family (3 dropped diseases return to training; their
  DIS-GEN edges are no longer held out). Consequence of change 2: the spec's
  Task-D reproduction trust anchor is dropped — Task-G numbers on OptimusKG
  will have no published-protocol delta to compare against; the benchmark
  stands on its own splits/nulls/baselines instead. KGPFLN needs a
  contamination audit (its 57 pretraining benchmarks may contain held-out
  diseases) before its run counts as zero-shot.
- **date:** 2026-07-25
- **approved:** owner — 2026-07-25
