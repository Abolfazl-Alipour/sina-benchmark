# 0020 — Task-G v5.1 freeze package (pre-registration numbers, roster, thresholds)

- **change:** Freezes the Task-G exam per `SINA_TaskG_Redesign_v5_1.md` Step 4a
  with actual measured numbers from Steps 1–3 (roster completion + rebuild,
  freeze-prep audits, F@k + provenance — all re-run on the corrected
  1,431,945-edge BRG after the annotation-orientation bug fix, commit
  `34c1d32`). This entry freezes: the roster, co-holdout clusters, key-size
  ceiling θ, floors, seeds, resolved training config, BRG composition,
  baselines, statistics + power statement, arm structure, and license
  posture. Two decision points are explicitly reserved for the owner (§12).
  No training on any non-frozen exam; this entry IS the freeze.
- **rationale:** v5.1 §6.1 ordering: roster rebuild → re-measured counts →
  overlap audit → ceiling classification → freeze → power/MDE → Gate 0. All
  inputs are now final and published (`reports/m2b/lodo_freeze_audits.md`,
  `reports/m2b/fk_drug_paths_and_provenance.md`,
  `data/processed/m2/lodo_audits/`).
- **requester:** v5.1 protocol Step 4a (owner-adopted via CHANGES/0019).
- **impact:** Everything downstream (LODO export builder, baselines, Gate 0,
  Step 0b, Arms A/B, statistics) is bound to the frozen values below. Any
  change requires a new CHANGES entry. GPU budget cap $900 unchanged.
- **date:** 2026-08-08
- **approved:** owner — 2026-08-08 (chat directive: "D1: do your
  recommendation. D2. do your recommendation. approved"). Both §12
  recommendations adopted: D1 = three-panel roster structure (primary
  confirmatory n=13; secondary-evidence n=5; secondary-low-context n=9+2;
  aging-rate/parental-longevity/frailty/life-span supervision-only);
  D2 = parkinsonism-plus sibling-strip (Kufor-Rakeb + PSP-parkinsonism
  mondo:0009839 DIS-GEN edges stripped when grading PD units, per-unit strip
  lists published).

---

## 1. Data and graph pins

| Pin | Value |
|---|---|
| Base graph | OptimusKG, Dataverse doi:10.7910/DVN/IYNGEV **version 2.0** (`configs/m1_trackd_ingest.yaml`), client `optimuskg==1.0.0`, 190,531 nodes / 21,813,816 edges / **27 relation labels** (ingest-verified 2026-08-08; DRG-BPO present) |
| Cluster | 678 diseases (neuro/aging/metabolic), roster completed per CHANGES/0019 Step 1 |
| Supervision | 191,953 DIS-GEN edges; 17,448 supervision genes; 1,270 cluster drugs |
| BRG (corrected) | **1,431,945 edges**: PHE-GEN 768,844 · GEN-GEN 275,220 · BPO-GEN 147,477 · CCO-GEN 94,150 · MFN-GEN 82,362 · PWY-GEN 43,754 · DIS-PHE 11,850 · DRG-GEN 5,310 · DRG-DIS 1,853 · DIS-DIS 1,125 |
| Candidate universe | 17,448 genes (all key-size percentages computed against this) |
| EFO/OBA pinning | all 6 aging traits resolve as direct Track-D CURIEs (`efo:0022597`, `efo:0007796`, `oba:0005494`, `efo:0005422`, `efo:0009885`, `oba:VT0005372`); none has a MONDO crosswalk xref; census method per `reports/m2b/aging_roster_census.md` (name search + manual false-positive removal; near-synonym misses unbounded) |
| MONDO | release 2026-07-06, hash-verified; ancestors taken from Track-D node `properties.ancestors` |

## 2. Frozen roster (47 scored units; mondo:0011149 dropped — 0 DIS-GEN edges by source)

Tiers (v5.1 §6.1 + §4 floors; **the connectivity-floor structure is owner
decision D1, §12** — the tiers below are the recommended three-panel
structure):

### Primary panel (key ≥50, context ≥10, key ≤5%) — 13 units, confirmatory

| unit | name | key | key % | ctx | F@3 |
|---|---|---|---|---|---|
| mondo:0010196 | Werner syndrome | 857 | 4.91 | 37 | 0.000 |
| mondo:0008310 | Hutchinson-Gilford progeria | 468 | 2.68 | 23 | 1.000 |
| mondo:0012853 | Fontaine progeroid syndrome | 230 | 1.32 | 98 | 0.000 |
| mondo:0014157 | mandibular hypoplasia-deafness-progeroid | 60 | 0.34 | 29 | 0.000 |
| mondo:0011913 | Alzheimer disease 3 | 75 | 0.43 | 23 | 0.000 |
| mondo:0008199 | late-onset Parkinson disease | 66 | 0.38 | 32 | 0.000 |
| mondo:0957576 | Parkinson disease 25, autosomal recessive | 151 | 0.87 | 13 | 0.000 |
| mondo:0007669 | renal cysts and diabetes syndrome | 307 | 1.76 | 38 | 0.000 |
| mondo:0010785 | maternally-inherited diabetes and deafness | 110 | 0.63 | 17 | 0.000 |
| mondo:0012192 | permanent neonatal diabetes mellitus | 335 | 1.92 | 33 | 0.000 |
| mondo:0012522 | diabetes mellitus, transient neonatal | 302 | 1.73 | 12 | 0.000 |
| mondo:0100164 | permanent neonatal diabetes mellitus | 532 | 3.05 | 15 | 0.000 |
| mondo:0100165 | permanent neonatal diabetes mellitus 1 | 381 | 2.18 | 10 | 0.000 |

### Secondary panel — evidence band (key 20–49, context ≥10) — 5 units, descriptive

mondo:0007088 AD1 (30/10) · mondo:0011150 acroosteolysis-keloid-like (29/65) ·
mondo:0013523 Nestor-Guillermo progeria (35/52) · mondo:0014527 progeroid
features-HCC (29/26) · orphanet:75496 EDS-progeroid (45/88).
All F@3 = 0.000.

### Secondary panel — low-context (key ≥50, context 2–9) — 9 units, descriptive

mondo:0000828 juvenile-onset PD (207/3) · mondo:0017279 young-onset PD
(226/2) · mondo:0100087 familial AD (405/4) · mondo:0015140 early-onset AD
(167/4) · mondo:0018911 MODY (389/9) · mondo:0020525 transient neonatal DM
(601/3) · mondo:0012480 transient neonatal DM (214/4) · mondo:0030089
permanent neonatal DM 4 (467/9) · efo:0005422 skin aging (368/8, F@3 0.976).
Plus double-secondary (key 20–49 AND low context): mondo:0018203
LMNA-cardiocutaneous progeria (41/2), mondo:0019321 atypical Werner (27/1).

### Context-only (key >5% ceiling, §4) — 6 units, NEVER graded

efo:0022597 aging (2,150; 12.3%) · mondo:0004975 Alzheimer disease (12,414;
71.1%) · mondo:0005180 Parkinson disease (6,584; 37.7%) · mondo:0005147 T1D
(5,664; 32.5%) · mondo:0005148 T2D (9,197; 52.7%) · mondo:0016391 neonatal
diabetes mellitus (1,161; 6.7%). Plus the 3 named trap terms (CHANGES/0019:
progeria / premature aging syndrome / progeroid syndrome).

### Supervision-only (never graded) — evidence floor or ~zero context

Below evidence floor (<20 keys): efo:0009885 frailty (11), oba:VT0005372 life
span (6), mondo:0014265 AD18 (19), mondo:0008311 (12), mondo:0008694 (6),
mondo:0012590 (8), mondo:0030880 (16), mondo:0957266 (8), mondo:0957953 (17),
orphanet:2963 (10). Effectively ungradeable on context (≤1 BRG edge):
efo:0007796 parental longevity (68/1), oba:0005494 aging rate (150/1).

## 3. Co-holdout clusters (v5.1 §5.3, measured + qualifier)

**Frozen rule:** pairwise overlap **coefficient ≥ 0.3 AND absolute
intersection ≥ 20 genes** vs any graded unit or supervision disease →
co-hold-out (the ≥20-intersection qualifier stops tiny keys from
chain-merging the roster; below it, 0.15–0.3 = disclosed in the per-unit
appendix; sensitivity at 0.2/0.5 reported). The full matrix is published at
`data/processed/m2/lodo_audits/overlap_matrix.parquet` (218,791 pairs).
MONDO ancestors of held-out members are stripped per the ancestor-strip rule
(v5.1 §5.3) regardless of overlap. Measured high-overlap **non-ancestor**
supervision diseases (e.g. Kufor-Rakeb 0.952 vs juvenile-onset PD;
PSP-parkinsonism mondo:0009839, 1.000 vs late-onset PD) are handled by a
per-unit **sibling-strip list** in the export builder — same mechanism as
ancestor-strip, published per unit (subject to owner decision D2, §12).

**Frozen clusters (each = ONE statistical unit):**

1. **Progeria nucleus** — Werner + atypical Werner + HGPS +
   LMNA-cardiocutaneous + Fontaine + skin aging. (Werner↔HGPS 0.346/int 162;
   Werner↔Fontaine 0.326/75; skin aging↔Werner 0.329/121 — the audit's
   non-obvious catch: skin aging nests into the progeria cluster, not the
   trait stratum.)
2. **AD familial** — familial AD + early-onset AD + AD3 + AD1 (+AD18
   supervision-only). All pairwise 1.000 containment into familial AD.
3. **PD graded** — juvenile-onset + young-onset + late-onset + PD25
   (juvenile≡young-onset 1.000/207; late-onset↔young-onset 0.682/45;
   PD25↔young-onset 0.570/86). One unit.
4. **MODY** — MODY + renal cysts and diabetes (1.000/307). One unit.
5. **NDM** — all 7 neonatal/transient/permanent diabetes subtypes
   (0012192, 0012522, 0020525, 0012480, 0030089, 0100164, 0100165; pairwise
   0.71–1.000). One unit.
6. Singletons: mandib-hyp-deaf-progeroid · MIDD · Nestor-Guillermo (hand-curated
   separate from LMNA; verified: max int to any graded unit 17 < 20) ·
   acroosteolysis-keloid-like · EDS-progeroid · progeroid features-HCC.

**Gradeable statistical units: 11** (6 progeria-area + AD + PD + MODY + MIDD +
NDM). Primary-panel-only: 13 diseases → **7 independent clusters**.

## 4. Key-size ceiling: θ = 5% (frozen, judgment call — disclosed)

The v5.1 working range was 2–5% "from the BFS-separation diagnostic." After
the BRG fix, that diagnostic **no longer shows a collapse**: sep@3 band means
+14.1 (below 2%) / +11.0 (2–5%) / +7.2 (above 5%) — monotonic, no threshold.
θ=5% is therefore frozen on the surviving evidence: (a) the mega-key
pathology itself (generic AD = 71.1% of the gene universe with sep@3 = −2.5,
the only clearly negative unit; the pilot autopsy that motivated the rule);
(b) face validity — every unit above 5% is an umbrella term, every unit below
is a specific disease; (c) robustness — the only boundary case is Werner
(4.91%): θ=4% would drop the single most important positive-control disease,
which we read as the line being arbitrary ±1%, disclosed, not as a reason to
move it. Sensitivity at θ=4%/6% reported with the results.

## 5. BRG composition (settled by measurement, v5.1 §4)

- **DIS-PHE: DROP.** The pre-registered PPR ± DIS-PHE ablation on the
  corrected BRG: ΔMRR −0.0152, ΔHits@10 −0.0310 over 42 informative units —
  with PHE-GEN exits restored, phenotype context dilutes gene signal.
  (Verdict flipped from KEEP after the annotation-orientation bug fix;
  both measurements published in `reports/m2b/lodo_freeze_audits.md` §7.)
  PHE-GEN gene-annotation edges STAY; only disease→phenotype context leaves.
- **DRG-BPO: EXCLUDE** from the BRG for this program. It is absent from the
  trimmed build's `message_passing_labels`; the F@k drug-circularity
  measurement was computed without it, so including it would invalidate the
  frozen F@k table and could only increase drug-mediated reachability.
  Recorded as a limitation (v5.1 §2.1 note); reconsider at v6.
- Included (measured composition, §1): GEN-GEN, BPO/MFN/CCO/PWY/PHE-GEN,
  DIS-DIS, DRG-DIS, DRG-GEN. ANA-GEN and EXP-* remain excluded (v5.1 §4).

## 6. Training protocol (frozen)

| Setting | Frozen value |
|---|---|
| Model | BioPathNet NBFNet, repo `b051a90`, committed patches (strict-negative guard, streaming eval) |
| Passes | **1 (pass-1 only), fixed, no checkpoint selection**; valid scored post-pass, reporting-only; intra-pass valid probes every ~150 batches, diagnostic-only |
| Config | input_dim 32; hidden_dims [32×6]; message distmult; aggregate pna; short_cut yes; layer_norm yes; dependent yes; criterion bce; num_negative 64; strict_negative yes; adversarial_temperature 1; heterogeneous_negative yes; heterogeneous_evaluation yes; full_batch_eval yes; Adam lr **1e-3** (frozen un-retuned; Step 0b re-confirms on the rebuilt graph); metric mrr |
| Batch | 300, fallbacks 280 → 250 pre-registered. **VRAM must be re-measured in Gate 0/Step 0b: the 142/144 GB measurement was taken on the 295k-edge BRG; the corrected BRG is 1.43M edges. If 250 OOMs, the fallback ladder extends by CHANGES entry — never ad hoc.** |
| Seeds | {1337, 1339, 5027}; 5-seed extension if primary-gate Wilcoxon p ∈ [0.05, 0.15] |
| Splits | per unit: 15% valid / 85% test on evidence-group IDs (SHA-1 `relation\|a\|b`); reciprocal twins co-assigned; ancestor-strip + sibling-strip lists applied per unit; unit context edges retained in BRG |
| Guards (Phase F, CI) | ancestor-strip verification; zero-vector GEN check (no isolated GEN in fact graph); `remove_pos: no` on warm-exam predict; no DIS-GEN in BRG; connectivity floor |

## 7. Baselines (frozen, same-graph rule)

Per LODO export, seeing exactly the model's training graph: (1) degree
baseline, pessimistic ties; (2) degree-matched null ≥100 resamples; (3) PPR
from the query disease over BRG-minus-DIS-PHE + training DIS-GEN, restart
**0.5** (sensitivity 0.15/0.85), symmetric normalization; (4) RotatE via
pinned PyKEEN on all export edge types, GEN-restricted eval, identical
filter/tie policy, d=256, 2 seeds (scheduled Step 4b, isolated CPU);
(5) engine-protocol popularity, both directions, pessimistic ties. AnyBURL:
warm-exam only. Orphan answers: included in primary ranking + excluded in a
secondary reachable-subset row, coverage statistic alongside.

## 8. Statistics and power (frozen, honest arithmetic)

- Co-primary: filtered MRR + filtered AUPRC, disease-macro, per-stratum;
  Hits@{1,3,10} reported; positive-control panel = exact causal-gene ranks on
  all monogenic units, gate = median rank ≤ 50 (chance ≈ 50/17,448 = 0.29%).
- Inference: hierarchical BCa bootstrap ≥1000 (strata → units); paired
  Wilcoxon per baseline comparison, two-sided α=0.05, BH across comparisons;
  effect sizes beside p-values; cluster = one unit.
- **Minimum detectable effect (normal approximation, paired, 80% power,
  α=0.05 two-sided, Wilcoxon efficiency penalty ~5%):**
  primary panel 13 units → dz ≈ 0.8; at the 7-cluster level → dz ≈ 1.1;
  all 29 graded diseases → dz ≈ 0.55; at the 11-cluster level → dz ≈ 0.9.
  **Interpretation (pre-registered): a NO-GO means "no very large effect at
  the cluster level," never "no effect"; CI excluding 0 with BH-corrected
  p > 0.05 = underpowered → triggers the 5-seed extension before any NO-GO.**
  Unit non-independence is disclosed; the hierarchical bootstrap is the
  mitigation; per-stratum co-primary reporting stands.

## 9. Drug-edge arms and exam integrity (frozen from measurement)

- **Exam-integrity verdict: PASS** — roster-median F@3 = 0.000 vs the 0.50
  pharmacological-reconstruction trigger; the full-roster Arm B contingency
  stays LOCKED (unlock condition was ≥0.10).
- F-covariate bands across graded units: F@3 ≥10% = {HGPS 1.000, skin aging
  0.976}; all other graded units <10% (all 0.000). F@2 ≤ 0.06 everywhere;
  exclusivity 0.997–1.0 (drug paths redundant with the backbone).
- **Arm B always runs on: the monogenic panel + skin aging.** Reconciliation
  recorded: v5.1 §6.3 says "monogenic panel + top-5 drug-rich units," but 4
  of the measured top-5 (generic PD, T1D, T2D + skin aging the exception) are
  above the §5.4 ceiling and therefore never graded; grading them in Arm B
  would violate the ceiling rule. The ceiling rule wins; the F-covariate
  regression uses the two gradeable ≥10% units vs the <10% rest, disclosed as
  thin (2-vs-many) in Limitations. HGPS remains always dual-arm.
- A/B per-answer rank-shift analysis on all dual-arm units, reporting-only.

## 10. License posture (frozen, auditable basis)

Supervision pool 191,953 edges: **99.36% Open Targets only** (190,734),
**0.64% DisGeNET+Open Targets** (1,219), 0 DisGeNET-only. Highest per-unit
any-DisGeNET: familial-AD member 7.65%, juvenile-onset PD 5.31%, MODY 3.60%.
Drug/clinical-token provenance (CTD/CGI via sources_indirect): pool 0.49%.
Full table: `reports/m2b/fk_drug_paths_and_provenance.md` §5 +
`license_provenance.parquet`. Research-use only; released harness inherits
NC-SA; commercial regeneration on Track-C; Track-C transferability gap
acknowledged (v5.1 §11.8) with conditional Step 8.

## 11. Execution schedule after sign-off (unchanged, v5.1 §13)

Phase F (local CPU): LODO export builder + strip-list guards + baseline
battery (~1 week wall-clock). Then owner-triggered GPU: Gate 0 (two-step
PrimeKG reproduction + permutation negative control + preemption-resume
drill + **VRAM re-measurement on the 1.43M-edge BRG**) → Step 0b pilot →
Arm A → Arm B subset → statistics/release. Cap $900.

## 12. OWNER DECISIONS — DECIDED 2026-08-08 (both recommendations adopted)

- **D1 — connectivity floor / roster tiers.** Recommendation (adopted in the
  tier structure above): three panels — primary confirmatory (context ≥10,
  n=13), secondary-evidence (key 20–49, n=5), secondary-low-context (context
  2–9, n=9+2) — secondary panels descriptive only, never pooled; aging-rate /
  parental-longevity / frailty / life-span become supervision-only (0–1
  context edges). Alternative: floor ≥2 for the primary exam (n≈22 primary)
  with per-unit disclosure — more units but grades diseases that are nearly
  isolated in the message-passing graph; we advise against for the
  confirmatory exam.
- **D2 — parkinsonism-plus sibling-strip.** Kufor-Rakeb (0.952 vs
  juvenile-onset PD) and PSP-parkinsonism mondo:0009839 (1.000 vs late-onset
  PD) are in the cluster via MONDO closure despite the CHANGES/0011
  parkinsonism-plus exclusion, and are NOT MONDO ancestors (ancestor-strip
  does not reach them). Recommendation: strip their DIS-GEN edges via the
  per-unit sibling-strip list when grading PD units (consistent with the
  measured co-holdout rule). Alternative: keep as context and disclose.

Sign-off on this entry with no amendments = both recommendations adopted.
**Owner signed off 2026-08-08 with no amendments — both adopted.**
