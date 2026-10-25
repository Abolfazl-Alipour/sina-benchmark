# 0027 — ANA-free graph + ranking-loss epoch (owner directive 2026-09-04)

- **date:** 2026-09-04
- **status:** EXECUTED per owner directive (verbatim: "remove the ANA, do the full epoch of the new ranking loss"). Records a deviation from the frozen v6/0026 spec; formal owner signature still collected at Phase-5 sign-off.
- **supersedes:** nothing — the halted ANA wave (`checkpoints_completion_ana/`, 4 PD epochs at B=50) is preserved untouched. This entry opens a **parallel** track; it does not modify it.

## Context (measured, on origin/main)

- Probe 2 (`reports/m2c/objective_probe.md`, `4b4b62c`): masked sampled ranking loss (1 pos + 255 masked negatives, gather-down log-softmax) beats frozen multi-hot BCE +35% on train-query filtered MRR (0.5105 vs 0.3790, beyond ±20% band), +15.5% val MRR (inside band), −11% step time, zero VRAM change. Verdict: promising, not decisive at 200 batches.
- CPU nulls (`reports/m2c/cpu_nulls_and_coverage.md`, owner-quoted af53881/0993e08): ANA feature is near-constant (99.1% of pinned genes carry it, median 41/48 edges); floors ±ANA identical (degree 0.00817 vs 0.00831, CIs overlap); Arm C (ANA ablation) cancelled — ANA is inert as signal.
- ANA is 1,397,674 of the 2,829,619 train1_ana edges → removing it halves the message graph and the ~20h/epoch wall-clock is expected to drop roughly proportionally.

## Decisions executed

1. **ANA removed from the message graph.** No rebuild: the pre-D0 KEEP-only completion exports (`data/processed/m2c/completion_exports/<unit>/`) are reused as-is — train1 = 1,431,945 KEEP edges; train2/valid/test splits byte-identical to the ANA units (PD 190,647/967/339, AD 189,450/1,882/621). Runner flag `--exports-dir` (default unchanged = ANA exports).
2. **One full epoch of the rank-loss variant on the KEEP graph**, PD `mondo__0005180` first, seed 1337. Rank loss as probed: per row 1 held-out positive + K=255 seeded negatives without replacement, all known positives (train2∪valid∪test) of the row's disease masked (count-asserted), gather-down from the [B,17,448] scores, log-softmax temp 1.0, NLL. Aux 1:1 GEN→BPO/PWY multi-hot weight 1.0 unchanged; C3 filtered val unchanged; checkpoints every epoch; `--ckpt-root reports/m2c/checkpoints_completion_noana` (ANA wave dir untouched).
3. **Everything else frozen:** BioPathNet `b051a90` 32×6 distmult/PNA dependent short_cut, Adam lr 1e-3, pinned 17,448 universe (C1), batch from a fresh VRAM sanity check on the KEEP graph before launch.

## Explicitly not decided here (owner call, see HANDOFF §4)

- Whether this becomes the new wave criterion (would need this entry signed + E2 relaunch semantics).
- AD run, multi-epoch continuation (early-stop patience 5 / 3% / 20), width probe (32 vs 64), junk-loci pruning (C5 kept them; ANA removal does not prune them).

## Artifacts

- Runner changes: `git show` for `--exports-dir` / `--ckpt-root` flags (default path bit-identical; dry-run config dumps verified).
- Run artifacts: `reports/m2c/checkpoints_completion_noana/mondo__0005180/` (epoch checkpoints + per-epoch JSONL) and the post-run report `reports/m2c/noana_rank_epoch.md`.

## Amendment 2026-09-05 — owner decisions (verbatim instruction parsed)

- **ANA wave formally BURIED** (owner: "fuck ana", decision 3): never resume the halted
  ANA wave, no AD run on the ANA stack. `reports/m2c/checkpoints_completion_ana/`
  remains preserved read-only as the audit trail. The v6 Arm-E/ANA-candidate plans are
  closed; ANA edges are excluded from the criterion stack on three independent
  inertness proofs plus equal test scores (0.01887 no-ANA vs 0.01883/0.01863 ANA).
- **Criterion signature timing** (decision 1): owner agreed with agent recommendation —
  sign this entry to make rank+noANA the criterion **after the AD epoch lands and
  confirms invariance on disease #2** (AD `mondo__0004975` K=all rank epoch launched
  2026-09-05 18:46Z, ETA ~04:50Z).
- **Post-AD direction** (decision 2): **(b) close the discovery track** — Phase-5
  writeup (`reports/m2c/completion_wave.md` + paper §3.7, certified negative +
  reranker + gates as the deliverable) — **and (c) attack candidate generation on
  CPU/API now** (the proven bottleneck: degree panels cover 2/339 hidden, OT 87/339,
  reranker converts 3.3× null given candidates).
- **100GB volume** (decision 4): deferred by owner; GitHub + release asset +
  root-EBS mirror remain the persistence story (`PERSIST_RESUME.md`).
