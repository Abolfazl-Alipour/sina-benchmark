<!-- track_scope: benchmark -->
# E2 PD — hidden test set scored once per checkpoint (loss-best vs rank-best)

Date 2026-09-04. Unit `mondo__0005180` (PD). The wave is halted by owner directive; this
is a one-shot read-only scoring of existing checkpoints — no training, no relaunch.
Ran on CPU (`bpn120`, `CUDA_VISIBLE_DEVICES=""`, 6 threads, ~12 s per checkpoint
forward); the CPU path is the Probe-1-validated one (it reproduces the trainer's logged
val MRR/H@10 exactly). Script: `scripts/pod/m2c_pd_test_once.py`.

Protocol: existing filtered evaluator over the pinned 17,448-gene universe (candidates =
pinned minus the 6,245 train2+valid positives of the PD head → 11,203 scorable
candidates), optimistic-tie bisect ranking — identical semantics to the trainer's
validation loop. All 339 hidden test queries share head `mondo:0005180` and relation
DIS-GEN, so one forward over the pinned tails per checkpoint scores every query exactly.
Test edges were never used in training, message passing, or filtering. Bootstrap CIs:
1000 resamples over the 339 queries, seed 1337.

## Checkpoint identity check (post-epoch-4 state)

`best.pth` as it stands after epoch 4 carries **epoch 3 weights**: epoch field 3, and all
model tensors are elementwise-identical to `epoch_003.pth` (shared tensor-hash
`c246451c30c5eb05c3766d141b809dd5…`; `epoch_004.pth` meta confirms `best_epoch: 3`).
Epoch 4 (val loss 0.0006516, −1.0% < 3% min-delta) did not improve on epoch 3's best
0.000658, so `best.pth` was not overwritten at epoch 4.

## Results (PD test, n = 339, filtered, pinned 17,448)

| checkpoint | role | epoch | test MRR [95% CI] | test Hits@10 [95% CI] | median rank | sha256 |
|---|---|---:|---:|---:|---:|---|
| epoch_001.pth | rank-best (val MRR 0.0195) | 1 | 0.01883 [0.01300, 0.02695] | 0.02950 [0.01180, 0.05015] | 171 | `a64d231d919589bf791464adc971f8d90b0aec5c7c141b86ca0c3cbd4b446ed3` |
| best.pth | loss-best (val loss 0.000658; = epoch_003.pth weights, verified above) | 3 | 0.01863 [0.01298, 0.02641] | 0.02950 [0.01475, 0.04727] | 174 | `1c245bb0fc31a1e9abf6cb1b867918963aa78ce6481e1e14cd18b03f175cfd73` |
| supervision-only co-count null (context) | no-graph null, train2-only counting | — | 0.01484 [0.00884, 0.02289] | 0.02950 [0.01180, 0.05015] | — | — |

Null context line quoted exactly from `reports/m2c/cpu_nulls_and_coverage.md` §4:
Score(PD,g) = Σ over d′∈D(g) of |G_pd ∩ genes(d′)|, train2 only, 1.2 s CPU.

## Reading

- The loss-best vs rank-best pair is statistically indistinguishable on the hidden test
  set: MRR 0.0186–0.0188 (CIs overlap almost completely), identical H@10 0.02950
  (10/339). The epoch-1 → epoch-3 val-loss improvement (0.000686 → 0.000658) does not
  translate into any test-ranking difference.
- Both checkpoints sit ~0.004 MRR above the supervision-only co-count null point
  estimate (0.01484) and ~2.3× the degree floor (0.00818), with all CIs overlapping;
  H@10 is exactly equal to the null (0.02950). The earlier caveat in
  `cpu_nulls_and_coverage.md` §4 ("GNN test unscored") is now resolved: on test the GNN
  edges ahead of the counting null on MRR point estimate only, and ties it on H@10.
- Cross-validation: the best.pth row matches Probe 1's genome-wide test derivation
  (MRR 0.0186, H@10 0.0295) from an independent code path, and matches the epoch_003.pth
  scoring bit-for-bit (same weights).

Claim ceiling: disease-associated gene recovery / candidate prioritization only.

## Reproducibility

```bash
P=~/miniconda3/envs/bpn120; OV=~/persist/cuda-12.8-overlay
PATH="$OV/hostbin:$P/bin:$PATH" CUDA_HOME=$P CXX="$P/bin/x86_64-conda-linux-gnu-g++" \
TORCH_EXTENSIONS_DIR=/tmp/torch_ext_probe CUDA_VISIBLE_DEVICES="" \
$P/bin/python -u scripts/pod/m2c_pd_test_once.py --out /tmp/pd_test_once.json
```

Result blob: `/tmp/pd_test_once.json` (ephemeral; aggregates reproduced above).
Zero GPU memory consumed.
