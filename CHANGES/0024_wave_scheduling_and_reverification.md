# 0024 — Light Gate 0 + Step 0b wave + 0023 supersession (DRAFT skeleton, 2026-08-12 directives v2)

- **change:** Record supersession of CHANGES/0023 lean-10 (13→10 cut) per Directives v2 §0: frozen roster **29 =13 primary+5+11** restored and binding; cost control via **wave scheduling** (Phase D), not re-scoping. Implement Light Gate 0 (B1-B6) + Step 0b Part 1 DIS-PHE keep-vs-drop (C1, model-side) with rebuild to `train1=1,431,945` if KEEP wins (re-run baselines + dual permutation + triviality + attrition) vs confirm DROP, then Part 2 weighting×{1 vs 2 passes} (C2), then Wave 1 (13×1 seed, directional only). Core ≈$487, cap $900.
- **rationale:** Five-review round R-H…R-L (2026-08-12) found lean-10 destroyed cluster inference (min exact Wilcoxon p=0.125), dropped AD coverage, violated freeze (S1); DIS-PHE DROP was PPR-baseline-only and 29/29 islands post-DROP (S3, biggest find); PrimeKG smoke vs full claim (S4); need dual shuffles (S5); plus H1-H4 housekeeping (11.67M/train1, θ truth, causal→monogenic, +20.2% citation). Directives v2 is governing plan.
- **requester:** Directives v2 FINAL 2026-08-12 (supersedes 2026-08-10) + Task-G v5.1 as amended through 2026-08-12
- **impact:** (a) 0023 superseded — no lean-10 training, Wave 1 restores 13 primaries; (b) DIS-PHE re-verdict: **KEEP wins** (Step 0b, +73.2% mean test MRR, see §(b) below); (c) pass-time re-pricing: **not triggered** (97.8 min/pass ≤ 140 at batch 32, L40S); (d) batch ladder outcome: **stable batch 32** (29.28 GB peak); (e) KEEP triviality re-verdict: **all 13 primaries trip the frozen A3 ≤2-hop criterion on the KEEP graph (41.7%–94.2%; 19/29 roster-wide)** while no baseline reaches MRR≥0.99 — criterion-vs-empirics conflict recorded in §(b), owner decision required; (f) KEEP exports rebuilt for all 29 units (`configs/m2_lodo_exports_keep.yaml`, `reports/m2b/lodo_exports_keep.md`) — filled 2026-08-16.
- **date:** 2026-08-12
- **approved:** pending (owner sign-off after B2/C1)

---

## Placeholders (fill as results land)

### (a) 0023 supersession
CHANGES/0023 lean-10 cut superseded 2026-08-12 per Directives v2 §0/S1: saving ~$10 vs cost of destroying 4-cluster inference and signed freeze violation. Roster remains 29 (13 primary) per CHANGES/0020. Stage: done. Evidence: `CHANGES/0023` marked superseded, `HANDOFF.md` M2B-6 updated.

### (b) DIS-PHE re-verdict (C1, model-side)
- Hypothesis: KEEP restores `disease→phenotype→gene` bridge to 768k PHE-GEN web for 29 island diseases (0+1-8 context).
- Pilot: 3 units × {keep (train1=1,431,945) vs drop (1,420,095)} × pass-1 × 1 seed 1337 (Fontaine-adjacent PD + HGPS + one secondary-evidence).
- Metric: model-side Δ filtered MRR/AUPRC beyond pre-registered noise band.
- Result: **KEEP wins (2026-08-16, Step 0b on GPU, seed 1337).** keep_pass1 vs drop_pass1 filtered test MRR: Werner `0.00653 vs 0.00346` (+88.7%), HGPS `0.01015 vs 0.00555` (+82.9%), RCAD `0.00930 vs 0.00629` (+47.9%); mean +73.2%. Evidence: `reports/m2b/step0b_ablation.md` + `step0b_ablation_results.json` + `data/raw/manifests` (step0b). Keep exports rebuilt for the 13 primaries (`*_keep_disphe`, train1 = 1,431,945 incl. 11,850 PHENOTYPE_PRESENT; train2/valid/test byte-identical to frozen DROP exam; tuple-disjointness verified).
- Deviation on record: the executed Step 0b matrix was partial — `keep_pass2` and `keep_capped` arms were not run; the 1-pass freeze (C2) rests on the DROP-graph pass-1-vs-pass-2 comparisons (Werner 0.00346 vs 0.00347; RCAD 0.00629 → 0.00425 regression) plus pass-1 keep. Closing the keep_pass2 hole is optional and requires a new CHANGES entry.
- Re-runs on the KEEP graph (mandated by this entry), executed 2026-08-16 in parallel with Wave 1:
  - Baselines: `m2_run_lodo_baselines.py` (`configs/m2_lodo_baselines_keep.yaml` → `reports/m2b/lodo_baselines_keep.md`): degree MRR 0.0006–0.0128, PPR MRR 0.0009–0.0170 across the 13 primaries — **no baseline MRR≥0.99 anywhere on the KEEP graph** (DROP-graph PPR-1.0 unit `mondo:0012853` falls to PPR 0.0017).
  - Dual shuffles: `configs/m2_lodo_shuffles_keep.yaml` → `reports/m2b/lodo_shuffles_keep.md`: 3 probe units × {across_disease, within_disease}, marginals exact, disjointness PASS, seed-1337 byte-identical reproducibility. GPU collapse-control training runs still pending (post-Wave-1).
  - Attrition census (KEEP train1 context edges incident to held-out disease): all 13 primaries at **10–98 context edges** (DROP graph: 1–2, 29/29 below the ≥10 floor) — `mondo:0007669` 38 (37 PHENOTYPE_PRESENT), `mondo:0012853` 98 (96 PHE), `mondo:0010196` 37 (34 PHE), `mondo:0100165` 10 (8 PHE, exactly at floor), etc. The C1 KEEP verdict restores the context floor for every primary unit.
  - Triviality re-run (`configs/m2_triviality_keep.yaml` → `reports/m2b/triviality_audit_keep.md`): **all 13 primaries trip the frozen A3 ≤2-hop criterion on the KEEP graph** — 41.7%–94.2% of test answers reachable in ≤2 hops via `disease→PHENOTYPE_PRESENT→gene` bridges (vs DROP graph where only baseline-MRR degenerates flagged). No unit trips the baseline-MRR≥0.99 prong. Extended to all 29 units (2026-08-16, secondaries' keep exports built): **19/29 degenerate** — 13 primaries + 4 secondary-evidence (orphanet:75496 87.8%, mondo:0013523 87.5%, mondo:0014527 76.9%, mondo:0011150 75.0%) + 2 low-context (mondo:0007088 70.4%, mondo:0030089 42.1%); the remaining 10 (9 low-context + efo:0005422 skin aging at 1.5%) stay below the 10% prong. The prong flags exactly the well-connected units — consistent with phenotype bridges driving it, not baseline recoverability (max KEEP baseline: PPR 0.0799, mondo:0018911).
  - **Conflict on record (owner decision required):** literal application of the frozen A3 criterion would mark all 13 primaries degenerate (descriptive-only, excluded from the primary gate) — which would nullify the benchmark the KEEP verdict was meant to save. The empirical evidence contradicts the criterion's premise: reachability ≠ recoverability (e.g. `mondo:0010196`: 83.8% of answers ≤2 hops yet Step 0b keep_pass1 test MRR = 0.00653, at the pre-registered island band). Resolving this — e.g. amending A3 to require baseline-lift evidence rather than reachability alone, or re-scoping the degenerate flag to the MRR≥0.99 prong — requires owner sign-off / new CHANGES entry. Until then Wave 1 completes as launched and all 13 units are reported; the primary-gate stats treatment is flagged pending in HANDOFF.md.
- Gate: CHANGES/0024 filled before any Arm wave.

### (c) Pass-time re-pricing (B2)
- Assumption: 70 min/pass on buggy 295k (invalid).
- Ladder probe on 1.43M: batch 300→280→250→200→150→100; below 100 → gradient accumulation via CHANGES.
- Measured (2026-08-15, AWS g6e.2xlarge L40S 48GB, Werner export, 59,454 train2 triples): rungs 300/280/250/200/150/100/64 all OOM (>42 GB); **batch 32 fits at 29.28 GB peak, 3158 ms/batch → 97.8 min/pass**; batch 16 fits (15.82 GB, 106.4 min/pass); batch 8/4 fit but slower per pass (124/183 min). Evidence: `reports/m2b/gate0_ladder_results.json` (sha256 `6c0ceace…f71a417d5d`), `reports/m2b/gate0_light.md`.
- Rule: if >140 min/pass (2×), PAUSE, re-price: floor config primaries-only×3 seeds ≈$400, full waves recomputed, owner sign-off.
- Result: **97.8 min/pass ≤ 140 → NO pause, no re-price.** Note hardware changed vs plan: L40S 48GB (not H200) → ladder fell below 100, so per this entry the ladder rule's fallback applies: batch 32 with gradient accumulation available (effective 128/256 at ×4/×8) instead of plain batch ≥100. Owner sign-off on this entry still pending (see header).

### (d) Batch ladder outcome (B1+B6)
- Stable batch: **32** (peak 29.28 GB VRAM, L40S 48GB); conservative fallback 16 (15.82 GB). All rungs 64–300 OOM.
- Island-expectation statement verified in `reports/m2b/gate0_light.md`: post-DROP 29 islands (0 phenotype +1-8) → expected MRR 0.01-0.03, beating random IS the result. (Note: C1 KEEP verdict post-dates this statement; on the KEEP graph units regain 1-8 phenotype bridges — expected MRR range still applies as the pre-registered band.)
- Manifest: `reports/m2b/gate0_ladder_results.manifest.json` (+ `reports/m2b/gate0_resume_drill_results.manifest.json` for B4).

---

## Reference

- Directives v2 §0-§7, Hard rules 1-15, Phases B/C/D, H1-H7 housekeeping.
- Pre-registered gates: degenerate-baseline (`MRR≥0.99` OR >10% ≤2 hops → excluded, Fontaine descriptive-only), valid_baseline flag, 13 primaries exact Wilcoxon + hierarchical BCa≥1000 + cluster permutation.
