<!-- track_scope: benchmark -->
# Generator Panels Wave 2 — candidate-generation attack (CPU, 2026-09-05)

Owner directive 2026-09-05 (decision 2c, CHANGES/0027 amendment): attack candidate
generation — the proven bottleneck — on CPU/API. All panels: pinned 17,448 universe,
train2+valid positives excluded (excl-knowns, realistic "new candidates" panels),
tiebreak gene token. Script `scripts/pod/m2c_generator_wave2.py`; panels + manifest +
shas committed under `reports/m2c/panels_wave2/`.

## Coverage@K (hidden-test positives in panel; PD n=339, AD n=621)

| generator | PD K=100 | PD K=200 | PD K=500 | AD K=100 | AD K=200 | AD K=500 |
|---|---:|---:|---:|---:|---:|---:|
| ot_overall (reference) | **0.257** (87) | **0.313** (106) | 0.416 (141) | 0.084 (52) | 0.179 (111) | 0.378 (235) |
| ot_literature (europepmc+uniprot_lit harmonic) | 0.192 (65) | 0.283 (96) | 0.451 (153) | **0.151** (94) | **0.287** (178) | **0.570** (354) |
| ot_functional (impc/crispr/...) | 0.133 (45) | 0.174 (59) | 0.263 (89) | 0.077 (48) | 0.159 (99) | 0.319 (198) |
| ot_genetics (gwas/eva/burden/...) | 0.077 (26) | 0.077 (26) | 0.077 (26) | 0.042 (26) | 0.076 (47) | 0.097 (60) |
| pheno_shared (HPO, graph-internal) | 0.127 (43) | 0.186 (63) | 0.289 (98) | 0.122 (76) | 0.205 (127) | 0.235 (146) |
| drug_indication (INDICATION+TARGET) | 0.035 (12) | 0.035 (12) | 0.035 (12) | 0.029 (18) | 0.029 (18) | 0.029 (18) |
| **union_best (min-rank fusion)** | 0.195 (66) | 0.307 (104) | **0.493 (167)** | 0.106 (66) | 0.196 (122) | 0.391 (243) |

Reference points from earlier waves: degree-top-K excl-knowns PD 0.0885@100;
PPR worse than degree (shortlist_generator_probe.md); GWAS windows binned (§6).

## Read-out

1. **ot_literature is the strongest single generator** (AD 57.0%@500, PD 45.1%@500;
   beats ot_overall on AD at every K) — but it is the *circular* channel: its scores
   are built from the same literature corpus as the DIS-GEN supervision (circularity
   audit, `reports/m2c/ot_panels/circularity_audit.json`). Treat as an upper bound,
   not a deployable claim.
2. **pheno_shared is the strongest structurally-independent generator**: HPO
   annotations, not disease-gene literature. PD 12.7%@100 / AD 12.2%@100 — 1.4× the
   degree panel, from graph-internal edges only (zero external API).
3. **union_best at K=500 sets the coverage record: PD 49.3%** (167/339) vs OT 41.6%.
   But min-rank fusion *dilutes* at K=100 (66 vs OT's 87) — weak generators crowd out
   strong ones' top ranks. Fusion helps only at deep K.
4. drug_indication pools are tiny (73 PD / 36 AD genes scored) — saturated
   immediately; high precision per slot, no depth.
5. ot_genetics pools are thin too (72 PD / 256 AD scored) — genetics evidence covers
   few pinned genes at all.

## Two-stage with the reranker (no-ANA epoch_001 score vector, PD)

Every wave-2 panel converts at **two-stage recall@10 = 0.0295** — identical on all
7 generators × 3 K values, and identical to every prior panel (OT, degree, GWAS).
The reranker's hit set is panel-independent across **10+ panel constructions** now.

**The ceiling is now pinned at both ends: coverage is fixable (union gets ~50% at
K=500), rerank conversion is not (0.0295 everywhere on this supervision). The product
story that survives: broad union panel + honest ranking; the model's contribution is
mid-rank MRR lift, not top-10 discovery.**

Claim ceiling: disease-associated gene recovery / candidate prioritization only.
