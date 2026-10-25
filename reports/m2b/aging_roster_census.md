<!-- track_scope: benchmark -->
# Aging roster census + exam-design direction (owner discussion, 2026-08-07)

**Purpose.** Record the census of aging-related disease/trait nodes in
Track-D (OptimusKG) and the exam-design direction discussed with the owner
on 2026-08-07. Owner decision pending; nothing here is implemented yet.
CPU-only analysis of `data/processed/m1/trackd/{nodes,edges}.parquet`
(36,345 disease nodes; 9,734,774 DIS-GEN edges) and the current trimmed
supervision pool.

## 1. Aging-related nodes in OptimusKG (the census)

Name search over all disease nodes (`aging|senesc|longev|proger|frailt|
werner|hutchinson|lifespan|elderly|centenar`, false positives like
"im**aging**" removed by hand) → **38 nodes, 14,831 DIS-GEN edges**, in
three buckets:

**Bucket 1 — aging traits (normal-aging population genetics):**

| node | DIS-GEN edges |
|---|---|
| efo:0022597 aging | 2,150 |
| efo:0005422 skin aging | 368 |
| oba:0005494 aging rate | 150 |
| efo:0007796 parental longevity | 68 |
| efo:0009885 frailty measurement | 11 |
| oba:VT0005372 life span trait | 6 |

**Bucket 2 — progeroid syndromes (monogenic premature-aging diseases,
genetically specific):**

| node | DIS-GEN edges |
|---|---|
| mondo:0010196 Werner syndrome | 857 |
| mondo:0008310 Hutchinson-Gilford progeria syndrome | 468 |
| mondo:0012853 Fontaine progeroid syndrome | 230 |
| orphanet:75496 Ehlers-Danlos syndrome, progeroid type | 45 |
| mondo:0018203 LMNA-related cardiocutaneous progeria | 41 |
| mondo:0013523 Nestor-Guillermo progeria syndrome | 35 |
| mondo:0014157 mandibular hypoplasia-deafness-progeroid | 60 |
| mondo:0019321 atypical Werner syndrome | 27 |
| mondo:0957266 RECON progeroid syndrome | 17 |
| mondo:0030880 mandibuloacral dysplasia progeroid | 16 |
| ~13 rarer terms | 1–29 each |

**Bucket 3 — TRAP terms (parent categories, never grade on these):**
`mondo:0019303 premature aging syndrome` (4,049 edges), `mondo:0015333
progeroid syndrome` (1,589), `mondo:0020732 progeria` (706) are MONDO
**parent classes** (verified PARENT hierarchy edges to the specific
syndromes). Their large gene lists are inherited aggregates — the same
mega-key pathology documented for generic Alzheimer disease in
`trimmed_regrade_autopsy.md` §3. Keep them in the graph as context;
exclude them from evaluation.

## 2. Gap vs the current trimmed graph (CHANGES/0016)

The trimmed supervision pool currently contains only `aging`,
`parental longevity`, Hutchinson-Gilford progeria, pseudoprogeria, and
progeria-short stature-pigmented nevi. **Missing:** Werner syndrome,
progeria, premature aging syndrome, progeroid syndrome, Fontaine, the
LMNA/Nestor-Guillermo/mandibular syndromes, aging rate, frailty, skin
aging, life span trait, hippocampal sclerosis of aging, and the remaining
rare progerias. Completing the roster = extend `cluster_roots`/
`aging_traits`/`progeria` in `configs/m2_trimmed_taskg.yaml` + rebuild
(~30 min, CPU). Bucket-3 terms belong in the background graph only.

## 3. Exam-design direction (discussed 2026-08-07)

Owner proposal: evaluate across many diseases instead of the single AD
family. Two standard designs:

- **Warm random-edge holdout** (hide ~10% of every disease's links):
  large and representative, but the test diseases are already known to the
  model — measures gap-filling, not zero-shot recovery. Most published
  numbers use this design. Cheap to add as a secondary metric.
- **Cold many-disease holdout** (hide ALL links of many diseases, one at a
  time, macro-averaged per disease): keeps zero-shot honesty, prevents any
  single noisy key from dominating (the `trimmed_regrade_autopsy.md`
  finding), and makes the aging traits + specific progerias first-class
  test cases.

**Direction recommended:** cold many-disease holdout inside the focused
cluster (neurodegenerative + parkinsonian + diabetes + aging — the
owner-chosen scope), aging roster completed per §2, disease-macro grading,
warm holdout as a secondary reported metric. Expanding to all 36,345
diseases = back to the full 9.7M-edge training problem (weeks of GPU);
the focused cluster keeps a training pass at ~70 min on one H200.

**Planned sequence (pending owner go-ahead):** complete aging roster →
rebuild trimmed graph → cold many-disease exam → short training run
(1–2 passes, the empirically observed peak) → ranked aging candidate
prioritization with mechanistic paths. Any grading/protocol change gets a
CHANGES entry first.
