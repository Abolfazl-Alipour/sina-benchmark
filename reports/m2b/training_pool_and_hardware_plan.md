<!-- track_scope: benchmark -->
# Training-pool focus + hardware plan — research findings (owner-requested, 2026-08-01)

**Scope:** findings only; **no protocol changes are approved by this document**.
The pool change itself is drafted as `CHANGES/0014` (pending owner sign-off).
Claim ceiling applies: "disease-associated gene recovery" language only.

## 1. Why full epochs are impossible (the measured mechanics)

The per-batch cost of BioPathNet is dominated by **full-graph message passing**
(6 NBFNet layers over the `train1 = 1,420,095`-triple LODO graph — 11.67M is the source OptimusKG count, never the per-unit GPU graph), not by the supervision pool
size. Measured on the A40 pod (`reports/m2b/gpu_smoke_estimates.md`): 7.9 s/batch
at batch 4 (memory-capped), backward scatter kernels = 75% of GPU time. Pool
size changes how many batches an epoch contains, never the per-batch cost:

| supervision pool | edges | full-epoch time (batch 4, 7.9 s) |
|---|---|---|
| global (current) | 9,734,774 | 222 days |
| nervous system disorder branch | 754,996 | 17 days |
| CNS disorder branch | 326,093 | 7.5 days |
| **neurodegenerative branch (proposed focus)** | **131,702** | **3.0 days** |
| diabetes mellitus branch | 19,824 | 10.9 h |

Conclusion: **no biologically meaningful pool makes full epochs viable on the
A40**; subsampling (CHANGES/0012) or a bigger-memory GPU (§3) is required.

## 2. Focused-pool candidates (computed from canonical Track-D, held-out AD/PD families excluded)

| Cluster (MONDO branch, mechanical definition) | diseases usable | DIS-GEN edges | per-epoch coverage at 1000 batches |
|---|---|---|---|
| global | ~10,000 | 9,734,774 | 0.04% |
| nervous system disorder (`MONDO:0005071`) | 2,995 | 754,996 | 0.53% |
| CNS disorder (`MONDO:0002602`) | 1,301 | 326,093 | 1.2% |
| **neurodegenerative disease (`MONDO:0005559`)** | **591** | **131,702** | **3.0%** |
| diabetes mellitus (`MONDO:0005015`) | 26 | 19,824 | — |

### Aging anchors (owner question: "is there an aging cluster?")

MONDO has **no aging disease branch**. What exists in the graph:

| anchor | edges | note |
|---|---|---|
| `efo:0022597` "aging" (Open Targets trait) | **2,150 DIS-GEN edges / 2,150 genes** | the real aging-associated gene set |
| `efo:0007796` "parental longevity" | 68 | longevity genetics |
| progeria syndromes (5 MONDO terms) | 515 | premature-aging diseases |
| GO cellular-senescence processes | 0 | no gene edges in this graph — unusable |
| HPO "Mortality/Aging" (`hp:0040006`) | 2 | unusable |

**Proposed focused pool (CHANGES/0014 draft, needs owner sign-off):**
neurodegenerative branch (131,702) + aging trait (2,150) + parental longevity
(68) + progeria (515) + diabetes branch (19,824) ≈ **154,259 supervision
edges** ≈ 3% per-epoch coverage (~75× the global relevance per batch, same cost).
Disclosure: the EFO "aging" node is a *trait*, not a disease — recorded as such.
Honesty: exam untouched (held-out families excluded); cluster is mechanical,
pre-registered, disclosed.

## 3. Hardware plan (measured-bottleneck-driven)

Bottleneck = memory bandwidth (sparse scatter ops) + capacity, not compute.

| card | capacity | bandwidth | $/hr | batch | full-pool epoch (154k edges) | program cost (AD+PD, 4 epochs × 5 seeds) |
|---|---|---|---|---|---|---|
| A40 (current) | 46 GB | ~0.7 TB/s | 0.40–0.70 | 4 | 3.0 days | dead (≈1.5 yr) |
| **A100-80GB** | 80 GB | ~2 TB/s | 1.1–2.0 | 16 | **~7 h** | **~$200–560** |
| **H200** | 141 GB | ~4.8 TB/s | 3.13–3.60 | 32 | **~1.5–2 h** | **~$200–400** (similar $, ~4–7 days wall-clock) |
| B200 | 180 GB | ~8 TB/s | ~6 | — | — | **rejected**: needs CUDA 12.4+/torch ≥2.4 (our verified stack is torch 2.0.1+cu118); capacity we cannot use; worst $/answer |

Notes: batch 32 deviates from the authors' effective-batch-16 recipe (minor,
documented). Multi-A40 data-parallel (~4× throughput) is viable but adds
engineering/failure modes for no accuracy gain. H200 is the speed pick at
equal total dollars; A100-80GB is the floor-price pick with zero surprises.

## 4. A*Net verdict (owner question: did BioPathNet test it?)

- **The A*Net inventor is on the BioPathNet team**: Zhaocheng Zhu (first
  author of NBFNet 2021 and A*Net 2022) is a BioPathNet co-author
  ([author list](https://www.nature.com/articles/s41551-025-01598-z)).
  No public head-to-head or "we tried it" statement exists.
- **A*Net's own Limitations section** ([arXiv 2206.04798](https://arxiv.org/pdf/2206.04798v1.pdf)):
  (1) in the inductive setting it "fails to match the performance of the
  original NBFNet" — an accuracy tax in exactly our zero-shot setting;
  (2) it is "algorithm design rather than system design," so real-world
  speedup is less than the theoretical message-count improvement.
- **Verdict:** A*Net is demoted from "obvious upgrade" to **fallback with a
  documented accuracy tax**, to be self-validated before any number counts.
  The literature-aligned path is full NBFNet on a bigger card.

## 5. Longevity rubric (clarification, not training)

The "longevity rubric" (spec §8.2.B) is a Phase-4 candidate-triage dimension:
a hand-curated 20–40 entry longevity list scored after the benchmark
(blinded curation, 25% double-curation, exploratory weight ≤0.15), cross-checked
against OptimusKG's embedded homologue field. It is **not** a training signal
and is separate from the aging-trait supervision edges in §2.

## 6. Standing order of operations

1. ad_fix1 verdict (CHANGES/0013 rerun) — is the exam winnable at all.
2. If winnable/near: owner signs `CHANGES/0014` (focused pool) + pick card
   (H200 for speed, A100-80GB for price) → full-pass runs (no subsampling)
   on the 154k-edge pool.
3. A*Net stays documented as fallback with its accuracy tax; BRG-shrink stays
   documented as an ablation dial.
