# Phase B Gate 0 Light on GPU Report

<!-- track_scope: benchmark -->

**Run Date**: 2026-08-15  
**Hardware**: AWS `g6e.2xlarge` (NVIDIA L40S 48GB GPU, 30 GB RAM, 4 vCPUs)  
**Governing Directives**: SINA Task G Redesign v5.1 + Coding Agent Directives v2 (2026-08-12) + CHANGES/0022  
**Status**: **ALL GATES PASS**

---

## 1. Summary of Gate 0 Light Outcomes

| Gate | Description | Target / Requirement | Outcome | Status |
|---|---|---|---|---|
| **B1** | VRAM / Throughput Ladder Probe | Rungs on 1.43M BRG export | Micro-batch 32 fits 29.28 GB VRAM, 3158 ms/batch | **PASS** |
| **B2** | Pass Time Bound | $\le 140$ min/pass ceiling | 97.8 min/pass at batch 32 | **PASS** |
| **B3** | Dual Shuffles Negative-Control | Across + within disease permutations | 3 units $\times$ 2 modes generated with marginals preserved | **PASS** |
| **B4** | Resume Drill Mid-Epoch | Save / reload / weight bitwise fidelity | Max parameter diff: 0.00e+00, loss step verified | **PASS** |
| **B5** | PrimeKG Smoke | Single-area smoke | Deferred to supplement per CHANGES/0022 | **DEFERRED** |
| **B6** | Gate 0 Light Report | Island expectation statement | Included verbatim below | **PASS** |

---

## 2. Benchmark Expectation & Context

> 29 islands (0 phenotype + 1-8) -> expected MRR 0.01-0.03, beating random IS the result.

The evaluation protocol tests zero-shot link prediction generalization across 29 frozen disease units. Due to severe topological isolation (0 to 8 phenotype connections), expected MRRs fall naturally in the 0.01–0.03 range. Beating random permutations and degree-preserving shuffles confirms signal extraction.

---

## 3. Detailed Results

### 3.1 B1/B2 VRAM & Throughput Ladder (NVIDIA L40S 48GB)

Measured on Werner syndrome LODO export (`mondo__0010196`, 59,454 training triples, 1,420,095 BRG edges):

| Batch Size | Status | Peak VRAM (GB) | Mean Latency (ms/batch) | Pass Time (min) | Note |
|:---:|:---:|:---:|:---:|:---:|:---|
| 300 | OOM | 45.51 | - | - | Initial ladder rung |
| 280 | OOM | 42.52 | - | - | Ladder rung |
| 250 | OOM | 45.95 | - | - | Ladder rung |
| 200 | OOM | 43.77 | - | - | Ladder rung |
| 150 | OOM | 42.69 | - | - | Ladder rung |
| 100 | OOM | 41.99 | - | - | Ladder rung |
| 64  | OOM | 44.45 | - | - | Extension rung |
| **32**  | **OK**  | **29.28** | **3158.1** | **97.8** | **Optimal training configuration** |
| 16  | OK  | 15.82 | 1718.6 | 106.4 | Conservative fallback rung |
| 8   | OK  | 9.12  | 1001.6 | 124.1 | Low memory footprint |
| 4   | OK  | 5.77  | 737.4  | 182.7 | Below throughput floor |

**Key Finding**: Micro-batch size 32 runs in 97.8 minutes per pass, well beneath the 140 min/pass ceiling. With gradient accumulation factor 4 or 8, effective batch sizes of 128 to 256 are achieved without NeighborLoader or memory overflow.

### 3.2 B3 Dual Shuffles Verification

Negative-control permutations generated on 3 probe units:
- `mondo:0010196` (Werner syndrome, 59,454 edges)
- `mondo:0008310` (Hutchinson-Gilford progeria syndrome, 58,021 edges)
- `mondo:0007669` (Renal cysts and diabetes syndrome, 81,624 edges)

Modes generated:
1. **Across-Disease (Popularity Null)**: Exact per-gene marginal degrees preserved.
2. **Within-Disease (Specificity Null)**: Exact per-disease marginal degrees preserved.

All generated splits passed strict tuple-disjointness and marginal-preservation guards.

### 3.3 B4 Resume Drill Verification

- **Checkpoint Size**: 88.28 MB
- **Parameter Fidelity**: Bitwise identical across solver reload (`max_param_diff = 0.00e+00`).
- **Optimization Continuation**: Pre-save loss zsh.7153 \rightarrow$ Post-resume step loss zsh.6316$.
- **Status**: **PASS**.

---

## 4. Manifests & Artifacts

- Ladder Results: `reports/m2b/gate0_ladder_results.json` (SHA256: `6c0ceace1d86d1f1e6d06eda7c7691186b2b79b48ab27e1a97af62f71a417d5d`)
- Ladder Manifest: `reports/m2b/gate0_ladder_results.manifest.json`
- Resume Drill Results: `reports/m2b/gate0_resume_drill_results.json`
- Resume Manifest: `reports/m2b/gate0_resume_drill_results.manifest.json`
- LODO Shuffles Report: `reports/m2b/lodo_shuffles.md`
- LODO Shuffles Manifest: `data/raw/manifests/m2_lodo_shuffles.json`
- Freeze Audits Report: `reports/m2b/lodo_freeze_audits.md`
- LODO Exports Report: `reports/m2b/lodo_exports.md`
- LODO Exports Manifest: `data/raw/manifests/m2_build_lodo_exports.json`
