<!-- track_scope: benchmark -->
# Shortlist-generator probe PD — can a better generator fix the coverage bottleneck? (zero training)

Date 2026-09-03. Unit `mondo__0005180` (PD). Zero training; read-only checkpoint scoring.
Ran entirely on CPU (`sina` env for generation, `bpn120` env with
`CUDA_VISIBLE_DEVICES=""` for model scoring); the live GPU training run was never touched.
Script: `scripts/pod/m2c_shortlist_generator_probe.py` (subcommands `generate`, `score`).
Follow-up to `reports/m2c/rerank_probe_pd.md` (Probe 1).

**Question:** Probe 1 showed the rerank bottleneck is candidate GENERATION — a
degree-top-100 panel contains held-out PD genes only 0.59% of the time (8.85% when known
positives are excluded from membership), while the model ranks covered positives at the
top. Does a better generator (disease-personalized PPR) produce panels that contain
held-out genes, and does the model still rank them top-10 when present?

## Construction

All generators: pinned 17,448-gene universe; panel membership EXCLUDES the 6,245 known
train2 (PD-head) + valid positives (the realistic not-yet-associated panel). Top-K by
generator score, tiebreak gene token (deterministic). K ∈ {100, 200, 500}.

- `degree_train1` — undirected train1 degree (Probe 1's `excl_known` panel; sha matches).
- `ppr_train1` — PPR personalized on `mondo:0005180`, train1 graph, alpha 0.5,
  max_iter 100, tol 1e-6 (`sina.baselines.ppr_baseline`, C1.4 parameters).
- `ppr_train1train2` — same on train1+train2 (matches the C1.4 baseline graph).
- DWPC-proxy generator: SKIPPED — no per-candidate score artifact exists from the C1.4
  run (`data/processed/m2c/completion_baselines_ana/` does not exist;
  `reports/m2c/completion_baselines_ana.md` carries aggregate metrics only).

Valid-split coverage note: under the primary exclusion rule valid positives are excluded
by construction (coverage structurally 0), so the table reports valid coverage under a
train2-only exclusion (valid treated as held-out) — a supplementary, clearly-labeled view.

## Input hashes (sha256)

| input | sha256 |
|---|---|
| train1.txt | `dff240c6087b304283838a8e9276fc1374c32fb7e342ee929d2de76d6cdab39f` |
| train2.txt | `18ce464837ee1e5c43a4a84c3d01f155c341bb5e5d3e36ee81a50d99d12bfeb3` |
| valid.txt | `6151e3e9a7200b5a4532aa889078fdc7c0dd9005d49d9f768c72c622794b1ba1` |
| test.txt | `18bbbcff754bf60db85620cf5b9d72cc81a3ce1a5e0cb4ac4985190590b8040e` |
| pinned_genes_17448.txt | `9dee160847afd92753eb4dcc1d0ec36320fd5d99f6114bdc74c7e7d62d2e3f95` |

Shortlist hashes (100/200/500-gene panels, token-joined): degree_train1
`0194c1c5…` / `ae50faed…` / `39aa5b15…`; ppr_train1 `eccd4b43…` / `2cf48c6d…` /
`d0aac71f…`; ppr_train1train2 `cac86a7e…` / `93886905…` / `51977bff…` (full values in
`/tmp/shortlist_gen.json`).

## Coverage of held-out positives (fraction of queries whose true gene is in the panel)

| generator | K=100 test | K=200 test | K=500 test | K=100 valid* | K=200 valid* | K=500 valid* |
|---|---:|---:|---:|---:|---:|---:|
| degree_train1 | **0.0885** (30/339) | **0.1504** | **0.2714** | 0.0496 | 0.0838 | 0.1934 |
| ppr_train1 | 0.0472 | 0.1180 | 0.2566 | 0.0341 | 0.0631 | 0.1748 |
| ppr_train1train2 | 0.0619 | 0.1180 | 0.2448 | 0.0383 | 0.0714 | 0.1675 |

\* valid under train2-only exclusion (see construction note).

**No generator beats raw train1 degree at any K.** Disease-personalized PPR concentrates
on the disease neighborhood — which is exactly where the excluded known positives live —
and after their removal PPR's remaining top genes cover held-out positives *worse* than
global degree. Even at K=500 the best panel contains only 27% of test positives.

## Model rerank on the panels (K=100, test n=339, bootstrap 1000/seed 1337)

Scoring: one NBFNet forward per checkpoint over the pinned 17,448 tails (all queries
share head `mondo:0005180` + DIS-GEN relation; Probe 1 machinery, exact). Panels contain
no known positives, so in-panel ranking needs no further masking. Null = the generator's
own score ordering reranking the SAME panel. Two-stage recovery = P(test gene in panel
AND model rank ≤ 10). `best.pth` = epoch 3 weights (Probe 1 verified).

| checkpoint | panel | coverage | two-stage | 95% CI | recall given coverage | null two-stage | null 95% CI |
|---|---|---:|---:|---|---:|---:|---|
| best.pth (ep3) | degree_train1 | 0.0885 | **0.0295** | [0.0147, 0.0472] | 0.333 (10/30) | 0.0088 | [0.0000, 0.0206] |
| best.pth (ep3) | ppr_train1train2 | 0.0619 | **0.0295** | [0.0147, 0.0501] | 0.476 (10/21) | 0.0088 | [0.0000, 0.0206] |
| epoch_002.pth | degree_train1 | 0.0885 | 0.0295 | [0.0118, 0.0472] | 0.333 (10/30) | 0.0088 | [0.0000, 0.0206] |
| epoch_002.pth | ppr_train1train2 | 0.0619 | 0.0295 | [0.0147, 0.0501] | 0.476 (10/21) | 0.0088 | [0.0000, 0.0206] |

Cross-validation: degree_train1 two-stage reproduces Probe 1's `excl_known` result
exactly (0.0295, same CI, same panel sha) — independent pipeline, identical numbers.

Hit detail (best.pth): the two panels share 7 of their 10 hits; each finds 3 uniques
(13 distinct genes total across both panels). Covered-hit model ranks sit at 1–10 while
the same genes' generator-score ranks are 1–98 — e.g. `ensembl:ENSG00000186868` model
rank 1 vs PPR rank 35; `ensembl:ENSG00000110944` model rank 9 vs degree rank 98. The
model reorders panels far better than the generator's own score: recall given coverage
0.33–0.48 vs ~0.10 for the null.

## Bottom line

**No.** No generator + current checkpoint reaches two-stage top-10 recovery ≥ 0.10.
Best observed: 0.0295 (10/339) on both the degree and PPR panels at K=100, statistically
above each panel's generator-score null at the point-estimate level (0.0088) — model CI
lower bounds 0.0147 exceed the null point estimate, though the null's bootstrap CI upper
bound (0.0206) overlaps the model CI, so CI separation is partial, not complete. The
limiter remains generation: the best K=100 panel contains 8.85% of held-out positives,
and even a K=500 degree panel caps coverage at 27.1%. Generation quality, not rerank
quality, bounds two-stage recovery for PD at this training stage.

## Honest interpretation

Candidate-generation coverage is the dominant constraint, and the obvious graph-proximity
upgrade (disease-personalized PPR) makes it *worse*, not better, once known positives are
excluded: PPR mass concentrates on the known-positive neighborhood. When a held-out gene
IS in the panel, the current checkpoints rank it top-10 a third (degree panel) to half
(PPR panel) of the time, vs about a tenth for the generator's own ordering — the reranker
already adds real value over its generator. But with 91–94% of held-out positives absent
from any 100-gene panel, two-stage recovery cannot exceed coverage (≤ 0.088 at K=100),
and the observed 0.0295 sits far below a practically useful 0.10+ bar. These numbers
support disease-associated gene recovery / candidate prioritization framing only; they
do not support stronger claims. A realistic path to higher two-stage recovery is better
generation (larger K, model-guided generation, or biological priors beyond graph
proximity) — all out of scope for this zero-training probe.

## Reproducibility

- Generation (fast: < 5 s total; PPR converges in ≤ 100 iterations on both graphs):
  ```bash
  PYTHONPATH=src ~/miniconda3/envs/sina/bin/python \
      scripts/pod/m2c_shortlist_generator_probe.py generate --out /tmp/shortlist_gen.json
  ```
- Scoring (~25 s; two forwards after cached JIT extension):
  ```bash
  P=~/miniconda3/envs/bpn120; OV=~/persist/cuda-12.8-overlay
  PATH="$OV/hostbin:$P/bin:$PATH" CUDA_HOME=$P CXX="$P/bin/x86_64-conda-linux-gnu-g++" \
  TORCH_EXTENSIONS_DIR=/tmp/torch_ext_probe CUDA_VISIBLE_DEVICES="" \
  $P/bin/python -u scripts/pod/m2c_shortlist_generator_probe.py score \
      --shortlists /tmp/shortlist_gen.json --out /tmp/shortlist_score.json
  ```
- Result blobs: `/tmp/shortlist_gen.json`, `/tmp/shortlist_score.json` (ephemeral;
  aggregates reproduced in this file).
- Zero GPU memory consumed; live PD training run (pid 31686) unaffected throughout.
