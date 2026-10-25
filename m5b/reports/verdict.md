# M5B verdict: LOSS — nothing built here beats Open Targets L2G

Deciding metric: covered PR AUC diff (ours minus OT), locus-blocked
bootstrap CI. OT covered PR AUC: **0.4767** (n=491 loci).

## Arms (chr folds; gene folds worse everywhere)

| Arm | Our PR AUC | Diff vs OT + 95% CI | Top-1 diff |
|---|---|---|---|
| v4 base (fresh-safe, no graph) | 0.253 | -0.222 [-0.271, -0.176] | -0.134 |
| + fine-mapping overlaps | 0.253 | -0.222 (no gain) | -0.165 |
| + ClinVar | 0.276 | -0.198 (+0.024) | -0.143 |
| + regulatory elements | 0.250 | -0.225 (no gain) | -0.153 |
| + eQTL overlaps | 0.253 | -0.221 (no gain) | -0.157 |
| all five tracks | 0.300 | -0.174 [-0.225, -0.126] | -0.128 |
| HPO/model-selected | 0.255 | -0.218 (worse than fixed) | -0.139 |

Every CI excludes zero on the negative side. LOSS on all arms.

## What each track earned (PR AUC deltas on chr folds)

- Fine-mapping posteriors: ~0.000. Distance already says what they say.
- ClinVar: ~+0.024. Only track with a measurable pulse.
- Regulatory elements: ~0.000.
- eQTL overlaps: ~0.000 at 1.4% row coverage (exact-match sparsity;
  3,159 genes, 70k links, 66 gene errors recorded).
- Combined: +0.047 over base. The gap is 0.22. Not close.

## HPO post-mortem (kept: fixed LightGBM)

Inner-fold AUPRC during selection was 0.54-0.60; pooled outer PR AUC came
out 0.255, below the untuned fixed model (0.300). Per the pre-registered
rule (outer trails inner by > 5pts -> keep the simpler model), the HPO
scores are discarded for the verdict. Likely causes, stated as hypotheses
not findings: per-fold winner-take-all mixed XGBoost and LightGBM score
scales in one pooled ranking, and 50-trial selection overfit inner folds.
The fixed LightGBM v5-all stands as our best: 0.300 vs 0.477.

## Leakage statement

No current-graph features. No OT scores as inputs. OT fine-mapping variant
lists used as public evidence (variant-level, never gene scores). HPO tuned
inside folds only. Colocalization engine never built (disease-stats gate
failed at 33.7%; correctly demoted, not faked).

## Limitations

- Single-run, one seed family (1337/1339/5027) unless noted.
- eQTL overlap is exact-variant-match only (no LD expansion, no ABF).
- Prospective 4,104-locus set curated but not yet evaluated (Track P).
- SINA GNN shows no unseen-gene signal (0.16/0.30); untouched by this work.

## Bottom line

Five public datasets, honest splits, tuned and untuned: our best is 0.300
against OT's 0.477. OT's functional-evidence advantage survives everything
we could legally add on CPU. Beating them requires either their missing
coverage (unprocessed spots), future discoveries (prospective set), or a
modeling breakthrough — not more public features.
