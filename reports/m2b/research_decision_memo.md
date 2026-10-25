<!-- track_scope: benchmark -->
# Research decision memo — zero-shot family holdout winnability on OptimusKG (owner-requested, research-only)

**Date:** 2026-08-01 · **Scope:** research only; no training run, no src/scripts/configs changes.
**Inputs:** HANDOFF §10/§10.5, `reports/m2b/ad_failure_diagnosis.md`, CHANGES/0012+0013,
`reports/m2b/ad_seed1_viability.md`, `reports/m2b/gpu_smoke_estimates.md`.
**Decision requested for:** ad_fix1 in flight (verdict gate ~2026-08-01 19:15Z).

---

## 1. Winnability — is the exam winnable at 3–23 context edges per held-out disease?

**Answer: yes, with caveats — and our connectivity is the SAME regime as the published protocol, not thinner.**

The key fact correction from primary sources: TxGNN's held-out diseases are not richly connected either.
From the published SI ([Nature Medicine 2024, Fig. S6, p. S15](https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-024-03233-x/MediaObjects/41591_2024_3233_MOESM1_ESM.pdf)):
held-out diseases average **18.8 first-hop connections** (cell-proliferation split) and **10.1** (anemia split),
vs ~60 for training diseases. Our CHANGES/0013 retention (3–23 edges/disease, mean ~10 across the 6 AD
diseases: DRG-DIS 315 + DIS-PHE 47 + DIS-DIS 13 + EXP-DIS 10 stripped, only DIS-PHE+DIS-DIS+hierarchies
retained) sits **inside the published 10–19 band at the top of our range and below it at the bottom**.

What the published protocol actually keeps for held-out diseases (TxGNN repo,
[`datasplit.py get_edge_group`](https://github.com/mims-harvard/TxGNN)): removes *only*
indication/contraindication/off-label edges of held-out diseases + a random sample of the 2-hop
neighborhood (5% of the *global* edge count) — **no per-relation-type removal**. Disease–protein,
disease–phenotype, disease–disease, exposure edges all survive. Note: their supervised task is
drug–disease, so disease–protein is *context* for them; for us disease–gene (DIS-GEN) is the
supervised target and must be removed. What remains as legitimate context for our held-out diseases
is exactly what CHANGES/0013 retains — **our retention set matches the protocol's intent; it is
thinner only because OptimusKG diseases carry fewer context edge types than PrimeKG diseases**
(we have no disease–protein context to keep — it is our label).

Supporting evidence that this regime is winnable for NBFNet-class models:

- Semi-inductive k-shot benchmark ([Kochsiek et al., Findings of EMNLP 2023, Table 7](https://aclanthology.org/2023.findings-emnlp.713/)):
  structure-only aggregation at **3–10 observed edges reaches H@1 0.29–0.40**, at 10 edges
  matching a transductive ComplEx average. Our 3–23 range covers that whole regime.
- BioPathNet's own zero-shot drug-indication results ([Nature BME 2026](https://www.nature.com/articles/s41551-025-01598-z)):
  NBFNet-class with sparse query connectivity + auxiliary BRG message passing, peer-reviewed, works.
- BioPathNet's PNA aggregator is degree-scale-robust by design.
- Caveat: NBFNet scores are path sums — edge *count* matters less than whether the retained edges
  connect the disease to well-connected neighbors forming multi-hop paths to genes. DIS-PHE/DIS-DIS
  are exactly the bridge types the published protocol relies on.

**Verdict:** the exam is winnable in principle at 3–23 edges; at the 3–5 edge end it is doubtful.
ad_fix1 is the correct decisive experiment.

## 2. DRG-DIS (drug–disease indication edges for held-out diseases in the BRG) — yes or no?

**Answer: NO for the primary protocol — it is a reviewer kill-shot. Acceptable only as an explicitly disclosed secondary ablation.**

- **No published protocol includes them.** TxGNN removes all indication/contraindication/off-label
  edges of held-out diseases ([repo code + README](https://github.com/mims-harvard/TxGNN));
  BioPathNet same ([Methods](https://www.nature.com/articles/s41551-025-01598-z)); COMIC uses
  TxGNN's splits verbatim with "complete exclusion of approved drugs from the training data"
  ([bioRxiv 2025.03.26.645428, §4.1](https://www.biorxiv.org/content/10.1101/2025.03.26.645428v1.full));
  KGML-xDTD avoids disease holdout entirely by design ([arXiv 2212.01384](https://arxiv.org/pdf/2212.01384v2)).
- **The leakage mechanism is concrete, not hypothetical:** OptimusKG has DRG-GEN (drug→gene
  target/enzyme, 20,694 edges). DRG-DIS + DRG-GEN gives a 2-hop *disease→drug→gene* path that
  directly proxies the DIS-GEN supervised answer — the model would rank genes by "genes targeted
  by the disease's known drugs," which is close to reading the answer key.
- Withholding is also cheap: DRG-DIS is 315 of the 385 stripped edges (82%), but the published
  regime shows the task is winnable *without* it (§1).

## 3. Alternatives

### 3a. KGPFN (arXiv 2605.14907 — note: the model's actual name is "KGPFN")

**Not a rescue for orphan diseases; viable only as an audited arm. Three hard gates before any number counts.**

- **Same connectivity floor.** KGPFN's per-candidate signal still runs NBFNet propagation from the
  query head ([§4.2–4.3](https://arxiv.org/html/2605.14907v1)); an isolated head collapses to
  uniform scores (its own §5.3 ablation shows monotonic degradation as neighborhood radius → 0).
  Its relation-level global context (retrieves 20 positive + 60–80 negative instances of the query
  relation from the target graph) is degree-independent, which degrades more gracefully than bare
  NBFNet — but no orphan-node evaluation exists in the paper.
- **Contamination audit (mandatory, from the research):**
  1. Claimed pretraining = FB15K-237 + WN18RR + CoDEx-Medium (general domain, no MONDO) — clean on paper.
  2. BUT the released checkpoint's actual training mixture is **undocumented** (repo ships
     `train_all.yaml` → all 57 graphs incl. Hetionet; HF checkpoint has no config card). If `all`
     was used, **Hetionet Disease–associates–Gene supervision is inside the checkpoint** — must be
     resolved (inspect checkpoint config / authors).
  3. The frozen ULTRA `ultra_50g` backbone's 50-graph mixture is **not public** — Hetionet is
     excluded from ULTRA's documented 10-graph transductive set, but the full mixture needs
     verification (Wikidata-QID overlap vs OptimusKG).
  4. Hetionet↔OptimusKG: map DOID→MONDO (DOID:10652 ↔ mondo:0004975; DOID:14330 ↔ mondo:0005180)
     and quantify triple-level overlap with our held-out families if any exposure is confirmed.
  5. **Protocol trap:** KGPFN retrieves in-context positives/negatives from the *target graph* at
     inference — our harness must exclude held-out family edges from that context pool or the model
     reads test supervision. Single most important protocol detail.
- **License blocker:** the [repo](https://github.com/HKUST-KnowComp/KGPFN) has **no LICENSE file**
  (all rights reserved by default) — cannot be used in the platform until clarified.
- **Biomedical track record is its worst tier:** Hetionet zero-shot MRR **0.195** (its weakest
  benchmark class, vs 0.432 average) and hard failure on structurally foreign graphs (HM splits,
  MRR 0.04–0.07) ([Appendix B.2, Tables 5–6](https://arxiv.org/html/2605.14907v1)).

### 3b. Other approaches for isolated query nodes

- Structure-only (ULTRA, InGram, NodePiece, TRIX, MOTIF, KG-ICL): **all share the ≥1-edge floor**;
  none has an isolated-node mechanism ([ULTRA arXiv 2310.04562](https://arxiv.org/abs/2310.04562);
  [InGram arXiv 2305.19987](https://arxiv.org/abs/2305.19987); [NodePiece arXiv 2106.12144](https://arxiv.org/abs/2106.12144)).
- The SI-LP benchmark is explicit: at 0-shot, structure-only methods **cannot run**; only text-anchored
  models (SimKGC, KGT5) work, "when textual information is provided" ([§4](https://aclanthology.org/2023.findings-emnlp.713/)).
  We do not have reliable text descriptions for all OptimusKG nodes — text-anchored fallback is a
  data project of its own, not a drop-in.
- New-node-prediction work targeting truly isolated nodes exists but needs node features
  ([arXiv 2401.05468](https://arxiv.org/abs/2401.05468)).

## 4. The 18.7% orphan answer genes — how do benchmarks handle unreachable candidates?

**Default: they never face the question.** Mainstream inductive benchmarks (GraIL, ILPC, InGram,
ULTRA evals) are constructed so every evaluated entity is connected in the inference graph
([InGram §5.4 requires a spanning tree](https://ar5iv.labs.arxiv.org/html/2305.19987)).
The accepted practices around it:

- **Silent exclusion is condemned:** LibKGE documents that removing unseen-entity test triples
  inflates results (~+3pp MRR) and calls it inconsistent ([README](https://github.com/uma-pi1/kge)).
- **Transparent stratification is published good practice:** KRACL reports MRR/Hits@10 across
  **7 in-degree bins** ([KDD 2022, Table 4](https://arxiv.org/abs/2208.07622)); Pujara et al.
  (EMNLP 2017) is the canonical sparsity-degradation analysis ([D17-1184](https://aclanthology.org/D17-1184.pdf));
  Hamaguchi et al.'s out-of-KB protocol discards unsupported entities **at construction time,
  pre-registered**.
- **Defensible form for us (pre-registered now):** (1) the headline metric is always the FULL test
  set (all 11,766 triples, orphans included, rank-floor for unreachable); (2) the stratum is
  structural — "degree 0 in the declared BRG," fixed before results; (3) stratified numbers are
  reported **alongside**, never instead of; (4) the unreachable stratum is reported as a coverage
  ceiling (18.7% structural cap on any purely path-based method). This mirrors KRACL + LibKGE +
  Hamaguchi precedent and is novel as a benchmark contribution — no published benchmark reports
  this exact split.

## 5. Recommendation (pick ONE of HANDOFF §10.5 options)

**Pick option (b) — richer BRG context, executed as a staged ladder with ad_fix1 as stage 1. Do NOT pick (a) DRG-DIS as primary.**

Execution plan:

1. **Stage 1 (running):** ad_fix1 verdict (CHANGES/0013, DIS-PHE+DIS-DIS+hierarchies, ~19:15Z).
   - PASS (beats degree baseline) → protocol settled; proceed to PD with the same export policy.
   - NEAR-MISS (beats null, within 2× of baseline) → add reachability-stratified metrics (§4) and
     report the connectivity gradient; then Stage 2.
   - FAIL (null level) → Stage 2.
2. **Stage 2 (fallback, disclosed ablation):** DRG-DIS inclusion, run ONCE on AD, explicitly labeled
   "leakage-disclosed secondary analysis" in every report — its numbers never headline; they exist
   to measure how much of the gap is context vs proxy. If it dramatically outperforms, that itself
   is reported as evidence of proxy leakage, not as a win.
3. **Stage 3 (parallel, no GPU):** implement reachability-stratified metrics in
   `scripts/m2_compute_metrics.py` per §4's pre-registered form (data work, ~1 day).
4. **Stage 4 (only if Stage 1–2 both fail):** KGPFN arm — gated on the audit in §3a (checkpoint
   mixture, ULTRA mixture, license) and the context-pool exclusion protocol. Not before.
5. **Stage 5 (always prepared):** if the honest outcome is negative, publish it: "NBFNet-class
   models on sparsely-connected disease holdouts" with the diagnosis, the connectivity-gradient
   stratification, and the published-regime comparison (TxGNN SI Fig. S6). That is a citable
   benchmark finding, and the claim ceiling language ("disease-associated gene recovery
   evaluation") already fits it.

**Honest expected outcome:** the published regime (TxGNN 10–19 edges works) says ad_fix1 has a real
chance; the thin end (3–5 edges on some of our 6 diseases) says a uniform PASS is unlikely. Most
probable outcome: **partial win** — beats null clearly, beats or approaches degree baseline on
better-connected diseases, with a visible connectivity gradient in the stratified metrics. That
outcome is publishable either way and leaves the DRG-DIS disclosure and KGPFN arm as documented
fallbacks rather than failures.

---

## Citations index (primary sources)

- TxGNN: [Nature Medicine 2024](https://www.nature.com/articles/s41591-024-03233-x) + [SI PDF](https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-024-03233-x/MediaObjects/41591_2024_3233_MOESM1_ESM.pdf) + [repo](https://github.com/mims-harvard/TxGNN)
- BioPathNet: [Nature BME 2026](https://www.nature.com/articles/s41551-025-01598-z)
- COMIC: [bioRxiv 2025.03.26.645428](https://www.biorxiv.org/content/10.1101/2025.03.26.645428v1.full)
- KGPFN: [arXiv 2605.14907](https://arxiv.org/abs/2605.14907) + [HTML](https://arxiv.org/html/2605.14907v1) + [repo](https://github.com/HKUST-KnowComp/KGPFN) + [HF checkpoints](https://huggingface.co/Eason-nuo/KGPFN)
- ULTRA: [arXiv 2310.04562](https://arxiv.org/abs/2310.04562) + [repo](https://github.com/DeepGraphLearning/ULTRA)
- InGram: [arXiv 2305.19987](https://arxiv.org/abs/2305.19987) · NodePiece: [arXiv 2106.12144](https://arxiv.org/abs/2106.12144)
- NBFNet: [arXiv 2106.06935](https://arxiv.org/abs/2106.06935)
- SI-LP k-shot benchmark: [Findings of EMNLP 2023](https://aclanthology.org/2023.findings-emnlp.713/)
- KRACL degree bins: [KDD 2022](https://arxiv.org/abs/2208.07622) · Sparsity: [EMNLP 2017](https://aclanthology.org/D17-1184.pdf) · LibKGE exclusion warning: [README](https://github.com/uma-pi1/kge)
- KGML-xDTD: [arXiv 2212.01384](https://arxiv.org/pdf/2212.01384v2)
