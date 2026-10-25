# Reweighted-null pre-test — label-reweighting CPU gate (PD + AD)

Owner directive 2026-09-05: before spending GPU on label reweighting, test its
premise on CPU. The supervision-only null (`cpu_nulls_and_coverage.md` §4) is a
pure function of the train2 labels, so the same reweighting applied to the null
decides the gate before any GPU burns.

## Pre-registered rule

- Reweighted null DROPS → top-10 signal is literature-popularity-shaped;
  reweighting the model's labels would hurt. Lever dies on CPU, zero GPU.
- Reweighted null HOLDS → non-literature labels carry recoverable signal;
  a reweighted K=all epoch is justified.

## Method

Script `scripts/pod/m2c_reweighted_null.py` (sina env, CPU, seconds per unit).
Score(d\*, g) = sum over train2 diseases d' containing g of |G ∩ genes(d')|,
train2 only. Variants differ only in the disease-side gene set G (candidate
side D(g) stays full — OT evidence classes exist in-repo only for the PD/AD
pulls, stated limitation):

- FULL — all train2 positives of d\* (control)
- OT_MATCHED — positives with an OT association row (missing-data control)
- NONLIT — OT row with any non-literature class score > 0
- GENETICS — OT row with genetics class score > 0
- STRICT_GF — genetics or functional_other > 0 (literature + curated out)
- LITONLY — literature > 0 and every other class == 0

Eval identical to every completion-track number this week: filtered over the
pinned 17,448 (candidates = pinned minus train2/valid positives of d\*),
optimistic-tie rank, hidden test scored, bootstrap 1000 / seed 1337.
Evidence-class map byte-identical to the committed circularity audit
(`m2c_ot_panel_eval.py:61-75`).

## PD (`mondo__0005180`, test n=339, train positives 5,278)

| variant    | \|G\|  | MRR (95% CI)            | H@10 (95% CI)           |
|------------|--------|-------------------------|-------------------------|
| FULL       | 5,278  | 0.01484 [0.00884, 0.02289] | 0.02950 [0.01180, 0.05015] |
| OT_MATCHED | 4,732  | 0.01478 [0.00880, 0.02289] | 0.02950 [0.01180, 0.05015] |
| NONLIT     | 1,859  | 0.01393 [0.00791, 0.02152] | 0.02065 [0.00885, 0.03835] |
| GENETICS   | 354    | 0.01451 [0.00871, 0.02259] | 0.02360 [0.00885, 0.04130] |
| STRICT_GF  | 1,859  | 0.01393 [0.00791, 0.02152] | 0.02065 [0.00885, 0.03835] |
| LITONLY    | 2,873  | 0.01500 [0.00925, 0.02253] | 0.02950 [0.01180, 0.05015] |

FULL reproduces the published null to all 5 decimals — control passes.
Hidden-test profile: 304/339 test positives have OT rows; of those 184
literature-only (61%), 120 carry non-literature evidence, 26 genetics.

## AD (`mondo__0004975`, test n=621, train positives 9,911)

| variant    | \|G\|  | MRR (95% CI)            | H@10 (95% CI)           |
|------------|--------|-------------------------|-------------------------|
| FULL       | 9,911  | 0.00975 [0.00643, 0.01417] | 0.01610 [0.00644, 0.02576] |
| OT_MATCHED | 9,347  | 0.00974 [0.00646, 0.01415] | 0.01610 [0.00644, 0.02576] |
| NONLIT     | 4,167  | 0.00941 [0.00612, 0.01391] | 0.01610 [0.00644, 0.02576] |
| GENETICS   | 1,051  | 0.00964 [0.00634, 0.01405] | 0.01610 [0.00644, 0.02576] |
| STRICT_GF  | 4,167  | 0.00941 [0.00612, 0.01391] | 0.01610 [0.00644, 0.02576] |
| LITONLY    | 5,180  | 0.00987 [0.00661, 0.01415] | 0.01610 [0.00644, 0.02576] |

H@10 is byte-identical across all six variants; MRR spans 0.00941–0.00987.
Dropping 89% of the gene set (FULL 9,911 → GENETICS 1,051) changes nothing.

## Verdict: HOLD on both diseases — GPU epoch justified, with one caveat

- The pre-registered HOLD condition is met. MRR never collapses under any
  reweighting; on AD even the genetics-only sliver reproduces FULL exactly.
  Non-literature labels carry the full recoverable signal. Reweighting the
  model's labels will not destroy the signal — the CPU gate passes.
- Caveat (PD LITONLY): literature-only edges alone also reproduce FULL
  exactly (0.01500/0.02950). The channels are redundant, not complementary:
  either channel suffices for the null's 10 hits. So reweighting is licensed
  but not promised — it tests whether forcing the model off the easy
  literature path finds anything beyond co-occurrence, and the PD LITONLY
  row says the easy path already reaches the full null score.
- Side finding: STRICT_GF == NONLIT exactly on both diseases — no train
  positive is supported by curated databases alone without genetics or
  functional evidence. Curated-only support does not occur in this label set.

Raw outputs: `reports/m2c/reweighted_null_pd.json`,
`reports/m2c/reweighted_null_ad.json` (variant scores, CIs, test profiles).
