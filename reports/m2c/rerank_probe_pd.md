<!-- track_scope: benchmark -->
# Rerank probe PD — does the current Arm E PD model carry a rerank-ready signal? (zero training)

Date 2026-09-03. Unit `mondo__0005180` (PD). Zero training: read-only scoring of the live
run's checkpoints. Ran entirely on CPU (`bpn120` env, `CUDA_VISIBLE_DEVICES=""`, 6 threads);
the GPU training run was never touched. Script: `scripts/pod/m2c_rerank_probe.py`.

**Question (peer-review directed):** genome-wide filtered val MRR is 0.016–0.019 with
H@10 2–4% vs the 0.00818 degree floor. Genetics teams bring 50–100 candidates and want a
top-10. Does the current model already recover held-out PD-associated genes inside a
realistic 100-candidate shortlist at **recall@10-in-100 ≥ 0.40**, with bootstrap CI clear
of the degree-baseline null on the same shortlists?

## Shortlist construction (deterministic, fixed for all queries)

All 967 valid and 339 test queries share one head (`mondo:0005180`) and one relation
(DIS-GEN = `ASSOCIATED_WITH`), so one fixed 100-gene shortlist applies to every query.
NBFNet scoring with `all_loss=None` depends only on (h, r, graph, weights) — the tail
index is a gather — so **one forward per checkpoint over the pinned 17,448 tails** yields
every query's genome-wide and in-shortlist scores exactly.

- **primary** (pre-registered construction): top-100 pinned genes by undirected degree in
  `train1.txt` (count head+tail occurrences; tiebreak: gene token). Known train2+valid
  positives of the disease stay in the 100 but are masked out of the ranking exactly like
  the runner's filtered evaluator (`excluded = train2 + valid positives − {answer}`):
  filtered non-hits, never countable, occupying no rank position.
- **excl_known** (post-hoc sensitivity, NOT the pre-registered construction): top-100
  pinned genes by train1 degree among genes **not** already train2/valid positives — the
  realistic "100 not-yet-associated candidates" panel. Included because the primary
  construction turned out to be 94% saturated with known positives (below).

A query positive outside its shortlist is an automatic miss (intended realism).
Null: degree ranking of the same 100 through the same filtered evaluator.
Bootstrap: 1000 resamples over queries, seed 1337, percentile 95% CI.

## Input hashes (sha256)

| input | sha256 |
|---|---|
| train1.txt | `dff240c6087b304283838a8e9276fc1374c32fb7e342ee929d2de76d6cdab39f` |
| train2.txt | `18ce464837ee1e5c43a4a84c3d01f155c341bb5e5d3e36ee81a50d99d12bfeb3` |
| valid.txt | `6151e3e9a7200b5a4532aa889078fdc7c0dd9005d49d9f768c72c622794b1ba1` |
| test.txt | `18bbbcff754bf60db85620cf5b9d72cc81a3ce1a5e0cb4ac4985190590b8040e` |
| pinned_genes_17448.txt | `9dee160847afd92753eb4dcc1d0ec36320fd5d99f6114bdc74c7e7d62d2e3f95` |
| shortlist primary (100 tokens) | `d066f30c543c9fd9a8bf76def2215ec554b3f45875c15b59ace7ced81cfc764c` |
| shortlist excl_known (100 tokens) | `0194c1c5103cc666ece82dc6e9e1066cd0314b8013605b339a72060c160fbb1a` |
| best.pth (epoch 3 weights) | `1c245bb0fc31a1e9abf6cb1b867918963aa78ce6481e1e14cd18b03f175cfd73` |
| epoch_002.pth (epoch 2) | `745cb3a84b8fca6d41159d2043cb9eeac61073eadcc4403c7bdcdf2b42fc6e32` |
| epoch_003.pth (epoch 3) | `67bd4fe5025aefb602a41bbd81b9a54bf5d7a5a1477ddedfb696263f30d8377d` |
| last.pth (epoch 3) | `67bd4fe5025aefb602a41bbd81b9a54bf5d7a5a1477ddedfb696263f30d8377d` |

Checkpoint identity note: `best.pth` currently holds **epoch 3** weights (val loss
0.000658 beat epoch 1 at epoch 3, so the runner overwrote it); its model tensors are
verified elementwise-identical to `epoch_003.pth`/`last.pth` (file bytes differ only via
a separate `torch.save` of the same state). The mission header's "best.pth = epoch 1" is
stale. Scored checkpoints are therefore epochs {2, 3}, with epoch 3 × 3 aliases.

Machinery cross-check: the probe's genome-wide filtered **valid** metrics reproduce the
runner's logged values exactly (epoch 2: MRR 0.0158 / H@10 0.0372; epoch 3: MRR 0.0168 /
H@10 0.0259), confirming scoring and filtered-ranking equivalence with the trainer.
Genome-wide filtered **test**: MRR 0.0186, H@10 0.0295 (both epochs; floor 0.00818).

## Coverage (the binding constraint)

| shortlist | scorable candidates (of 100) | test coverage | valid coverage |
|---|---:|---:|---:|
| primary | 6 | **2/339 = 0.59%** | 13/967 = 1.34% |
| excl_known | 100 | 30/339 = 8.85% | 0/967 (structural: valid positives are known by construction) |

The degree-top-100 shortlist for PD contains **94 already-known train2/valid positives**;
only 6 scorable candidates remain, and only 2 of the 339 held-out test positives are in
the list at all. PD test positives are low-degree, not-yet-heavily-studied genes;
train1 degree is effectively a proxy for "already documented for PD".

## Results — primary (pre-registered) shortlist, test split (n = 339)

| checkpoint | epoch | recall@10-in-100 | 95% CI | MRR-in-100 | degree null (CI) | calib slope (hits) | med gw rank | med in-100 | med shift |
|---|---:|---:|---|---:|---|---|---:|---:|---:|
| best.pth | 3 | 0.0059 | [0.0000, 0.0147] | 0.0044 | 0.0059 ([0.0000, 0.0147]) | 1.11 (2) | 174 | 1.5 | 3 |
| epoch_002.pth | 2 | 0.0059 | [0.0000, 0.0147] | 0.0044 | 0.0059 ([0.0000, 0.0147]) | 1.33 (2) | 175 | 1.5 | 13 |
| epoch_003.pth | 3 | 0.0059 | [0.0000, 0.0147] | 0.0044 | 0.0059 ([0.0000, 0.0147]) | 1.11 (2) | 174 | 1.5 | 3 |
| last.pth | 3 | 0.0059 | [0.0000, 0.0147] | 0.0044 | 0.0059 ([0.0000, 0.0147]) | 1.11 (2) | 174 | 1.5 | 3 |

Both covered test positives (`ensembl:ENSG00000177606` genome-wide rank 3 → in-100 rank 1;
`ensembl:ENSG00000188906` gw rank 6 → in-100 rank 2) are ranked at the very top — but the
degree null ranks them top as well (they are the highest-degree scorable genes), so model
recall exactly equals the null. Valid split mirrors this: recall 0.0134 = valid coverage =
null (13/967), MRR-in-100 0.0062 (ep3).

## Results — excl_known sensitivity shortlist, test split (n = 339)

| checkpoint | epoch | recall@10-in-100 | 95% CI | MRR-in-100 | degree null (CI) | calib slope (hits) | med gw rank | med in-100 | med shift |
|---|---:|---:|---|---:|---|---|---:|---:|---:|
| best.pth | 3 | 0.0295 | [0.0147, 0.0472] | 0.0118 | 0.0088 ([0.0000, 0.0206]) | 0.62 (30) | 174 | 15.5 | 69 |
| epoch_002.pth | 2 | 0.0295 | [0.0118, 0.0472] | 0.0118 | 0.0088 ([0.0000, 0.0206]) | 0.69 (30) | 175 | 15.5 | 75.5 |
| epoch_003.pth | 3 | 0.0295 | [0.0147, 0.0472] | 0.0118 | 0.0088 ([0.0000, 0.0206]) | 0.62 (30) | 174 | 15.5 | 69 |
| last.pth | 3 | 0.0295 | [0.0147, 0.0472] | 0.0118 | 0.0088 ([0.0000, 0.0206]) | 0.62 (30) | 174 | 15.5 | 69 |

Of the 30 covered test positives, the model puts 10 in the top-10 (genome-wide ranks of
the hits: 1, 3, 5, 6, 7, 11, 13, 14, 15, 27) vs 3 for the degree ordering — 3.3× the null,
and the model CI lower bound (0.0147) sits above the null point estimate (0.0088).

Calibration slope: pooled logistic fit of hit indicator on raw model score over
test queries × filtered shortlist candidates (Newton–Raphson). Positive slopes throughout
(scores carry hit information) but hit counts are small (2 primary / 30 sensitivity), so
slopes are indicative, not precise.

## Verdict

**FAILS the pre-set bar.** On the pre-registered degree-top-100 construction,
recall@10-in-100 = 0.0059 for every checkpoint vs the required ≥ 0.40, and the CI
[0.0000, 0.0147] is not clear of the degree null (0.0059 — identical by construction).
On the not-yet-associated sensitivity panel the model does beat the degree null
(0.0295 vs 0.0088, CI lower bound above the null point) but remains an order of
magnitude under the bar.

## Honest interpretation

The "genome-wide 2–4% H@10 routinely becomes 30–50% within-100" expectation does not
materialize for PD here, and the dominant reason is **coverage, not ranking quality**.
A degree-built top-100 for this disease is 94% occupied by genes already documented as
PD-associated in train2/valid; the held-out test positives are mostly low-degree genes
that no degree-based shortlist contains (2/339 covered). Where a test positive is in the
shortlist, the model ranks it at or near the top (ranks 1–2 of 6 scorable), but so does
the degree ordering — the two signals are confounded at the top of the degree
distribution. On the sensitivity panel of not-yet-associated genes the model shows a
real but small prioritization signal: 3.3× the degree ordering at top-10 (10/339 vs
3/339), bootstrap CI above the null point estimate, positive calibration slope, and a
median genome-wide→in-100 rank shift of ~69–76 positions for covered positives. That is
far short of a usable 40% top-10 recovery. Claim ceiling for these numbers:
disease-associated gene recovery / candidate prioritization only; nothing here supports
stronger claims, and the shortlist-construction sensitivity (coverage 0.6% vs 8.9%
depending on design) is itself a key finding for how a rerank product would have to
source candidates.

## Reproducibility

- Script: `scripts/pod/m2c_rerank_probe.py` (ruff-clean). Full results blob:
  `/tmp/rerank_probe_pd.json` (ephemeral, not committed; all aggregates are in this file).
- Run line (CPU, ~55 s total; ~12 s per checkpoint forward after one-time JIT compile):
  ```bash
  P=~/miniconda3/envs/bpn120; OV=~/persist/cuda-12.8-overlay
  PATH="$OV/hostbin:$P/bin:$PATH" CUDA_HOME=$P CXX="$P/bin/x86_64-conda-linux-gnu-g++" \
  TORCH_EXTENSIONS_DIR=/tmp/torch_ext_probe CUDA_VISIBLE_DEVICES="" \
  $P/bin/python -u scripts/pod/m2c_rerank_probe.py --out /tmp/rerank_probe_pd.json
  ```
  `TORCH_EXTENSIONS_DIR` is a probe-private JIT cache: with CUDA hidden, torchdrug's
  spmm/rspmm extensions rebuild as CPU-only under a different version hash, so ninja (in
  `$P/bin`) must be on PATH; the training run's extension cache and the GPU are untouched.
- Smoke verification before the full run: `--smoke 3` processed 3 valid + 3 test queries
  end-to-end per checkpoint (all four checkpoints), then the full run scored all
  339 test + 967 valid queries per checkpoint.
- Zero GPU memory consumed; the live PD training run (pid 31686) was unaffected.
