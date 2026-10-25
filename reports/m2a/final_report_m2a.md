<!-- track_scope: benchmark -->
# SINA Platform — Phase 2 (M2A) Final Report

**Date:** 2026-07-24
**Scope:** M2A benchmark harness (spec §6) — splits, baselines, nulls, metrics,
statistics, Gate 2 evidence.
**Status:** COMPLETE (locally feasible items). Owner-pending items remain.

---

## 1. What was done

### Phase 0/1 recap (already complete)

| Milestone | Result |
|---|---|
| M0 — Phase 0 | Preflight docs, V1–V7, CHANGES 0001–0007; owner approved V1–V6 |
| M1 — Track-R (PrimeKG) | 129,312 nodes / 8,100,498 edges / 4,050,064 evidence groups |
| M2 — Track-D (OptimusKG) | EXACT 190,531 nodes / 21,813,816 edges; gene_properties; temporal split DOWNGRADED (0.25% coverage) |
| M3 — Track-C skeleton | Catalog + CI license gate; owner approved |
| M4 — Gate 1 package | `reports/m1/gate1_evidence.md`; CI green |

### M2A harness (this phase)

| Milestone | Commit | Result |
|---|---|---|
| M2A-0 — env + approvals | `6f895dd` | torch 2.13.0+cpu, pykeen 1.11.1, scipy 1.18.0; owner approvals recorded |
| M2A-1 — MONDO crosswalk + families | `2042d63` | 42.3% MONDO coverage on Track-D; AD family 22 / PD 28 drafted (awaiting owner sign-off) |
| M2A-2 — splitters + leakage guards | `621321a` | 6 split artifacts on `evidence_group_id`; `assert_split_disjoint` enforced in CI; Track-R family grouped-node expansion (CHANGES/0009) |
| M2A-3 — baselines + nulls | `81182e0` … `30147bb` | degree/PPR/null pilots; PyKEEN KGE wrapper; AnyBURL wrapper (zero-shot family unusable — unseen diseases) |
| M2A-4 — metrics + statistics | `4aeaafe` | filtered/raw MRR/Hits@k, degree-stratified metrics, mean±SD/CI, paired Wilcoxon |
| M2A-5 — Gate 2 evidence | `2a1439b` | `reports/m2a/gate2_evidence.md` |
| M2A extras | `8d80cdf`, `51325bb` | multi-seed stats runner; external-eval skeleton; KGE 1%-subgraph pilot (93 s/epoch → ~2,600 h full-scale on CPU — impractical) |

### Splits produced (all on `evidence_group_id`, leakage-guarded)

| Split | Track | Protocol | Test rows |
|---|---|---|---|
| trackr_random | R | random 90/5/5 | 405,988 |
| trackd_random | D | random 90/5/5 | 1,093,156 |
| trackr_family | R | AD/PD family holdout | 914 |
| trackd_family | D | AD/PD family holdout | 18,625 |
| trackr_txgnn_neurodegenerative | R | TxGNN area | 514 |
| trackd_temporal | D | temporal (secondary) | 252 |

### Headline metrics (filtered MRR)

| Split | Degree | Degree-matched null | PPR |
|---|---|---|---|
| trackr_family | 0.0026 | 0.0005 | 0.0000 |
| trackd_family | 0.0620 | 0.0111 | 0.0000 |

Degree beats the degree-matched null on both tracks. PPR returns no scores for
zero-shot family diseases (unseen in train).

### Key code artifacts

- `src/sina/splits.py` — random/family/TxGNN/temporal split primitives
- `src/sina/baselines.py` — degree, PPR, degree-matched null, relation-shuffled null
- `src/sina/kge.py` — PyKEEN wrapper (TransE/ComplEx/RotatE-compatible)
- `src/sina/anyburl.py` — AnyBURL wrapper (triple export, config, rule filtering, parsing)
- `src/sina/metrics.py` — filtered/raw MRR/Hits@k, degree-stratified metrics
- `src/sina/stats.py` — multi-seed aggregation, paired Wilcoxon
- `scripts/` — resumable, manifest-checkpointed pipelines for all of the above

---

## 2. What remains

### Owner-pending (cannot be done by agent)

| Item | Why pending | Action needed |
|---|---|---|
| Disease-family sign-off | `configs/disease_families.yaml` is a draft (§0.3.9 biological-semantics gate) | Owner review of AD/PD family membership |
| Phase 0.5 cloud-GPU smoke test | Needs owner cloud accounts + spend approval (≤$50 cap) | Owner provision cloud GPU |

### Locally infeasible (needs cloud or reduced protocol)

| Item | Why infeasible | Path forward |
|---|---|---|
| Full-scale KGE (TransE-256/ComplEx-256/RotatE-512, ≥5 seeds, ≥50 epochs) | ~2,600 hours on local CPU (108 days) | Cloud GPU (Phase 0.5) or reduced-graph protocol |
| AnyBURL on zero-shot family | Family diseases unseen in train; rule-based methods cannot score them | Evaluate AnyBURL on warm-start or non-family splits |
| External eval set | Data source not chosen (PENDING_DATA_SOURCE) | Owner selects source (e.g., DisGeNET 2026+, Open Targets post-cutoff) |

### Next phase (M2B frontier)

| Item | Gate | Notes |
|---|---|---|
| BioPathNet Task-D reproduction (PrimeKG) | Phase 0.5 cloud-GPU smoke test | Cloud GPU required |
| BioPathNet Task-G (first benchmark on OptimusKG) | M2A harness complete | Uses M2A splits + metrics |

---

## 3. Verification state

- **Tests:** 100 unit + integration tests green (`pytest tests/unit tests/integration -q`)
- **Lint:** ruff clean; report lint clean (10 files, 0 violations)
- **CI:** pushed to `main`, all jobs green
- **Reproducibility:** all scripts `--resume`/`--dry-run`; manifests in `data/raw/manifests/`; data artifacts hash-verified in `data/raw/REGISTRY.json`

---

## 4. Bottom line

Phase 2 (M2A) is **complete for all locally feasible work**. The benchmark harness
— splits, baselines, nulls, metrics, statistics — is built, tested, and committed.
The remaining items are either **owner sign-offs** (disease families, cloud GPU)
or **cloud-dependent** (full-scale KGE, BioPathNet). The project is ready for
Gate 2 review and for the owner to unlock Phase 0.5.
