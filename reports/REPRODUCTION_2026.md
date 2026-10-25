<!-- track_scope: benchmark -->
# REPRODUCTION_2026 — clean verification log (T2 of the 2026-10-06 audit handoff)

**Date:** 2026-10-06/07. **Machine:** beast (owner's local box, 12 threads / 44GB RAM, CPU-only).
**Repo state:** `main` at `d5caa5c` (paper correction pass). **Env:** existing `sina` conda env (python 3.12.13, torch 2.13.0+cpu, polars 1.43.0, duckdb 1.5.5). This is the environment the numbers were reproduced *in*; a from-scratch `environment.lock.yml` rebuild was not performed (env rebuild is a heavier exercise; recorded as a gap, not silently skipped).
**Rule for this log:** every claim below states *how* it was verified. Anything not re-executed is marked VERIFIED-BY-ARTIFACT (traced to committed files) or NOT-RUN-HERE (with the reason). No number was patched to match.

## Verdicts

| # | Item | How verified | Result |
|---|---|---|---|
| 1 | ruff | `ruff check src tests scripts` | PASS (0) |
| 2 | Unit tests | `pytest tests/unit` | **211 passed**, 0 failed |
| 3 | Integration tests | `pytest tests/integration` | **47 passed**, 0 skipped, 0 failed (235s) |
| 4 | Claim lint | `python -m sina.lint_reports reports/ --config configs/report_lint.yaml` | **89 files, 0 violations** |
| 5 | Per-locus gates (PD 0.904 / AD 0.954) | `m2c_loci_headtohead.py gate` on committed score vectors + checkpoints | PASS both, exact (0.904/0.954) |
| 6 | Head-to-head, all 4 panels | full `run` re-execution, both units × {train2, test-knowns} | **EXACT MATCH** to committed JSONs (n, chance draws, all strata + bootstrap CIs; seeded bootstrap 1000/1337 is deterministic) |
| 7 | M4B v3 report | `tabular_baseline.py report --featureset v3` from committed `train_results_v3.json` | **byte-identical** regeneration (after fixing one hardcoded banned-term string in the script — see findings) |
| 8 | M5B verdict numbers | parsed `m5b/reports/v5_results.json` | covered PR AUCs match `verdict.md` (0.253/0.253/0.276/0.250/0.253/0.300 vs OT 0.4767; all-five diff −0.1741 [−0.2249,−0.1257], n=491) |
| 9 | Prospective probes | re-ran `m2c_prospective_null.py` + `m2c_prospective_model.py` (CPU) | PD model MRR 0.00373 H@10 0/29; AD 0.00300, 2/239; nulls PRE 0/82 and 0/1086 — all match §3.7 of the paper |
| 10 | Release asset `data-m2c-2026-09-05` | downloaded via GitHub API (58MB) | **SHA-256 `a0bb706c…806a42` OK** vs PERSIST_RESUME.md |
| 11 | Full export counts | local `data/processed/m2c/full_exports/` line counts | 9,709,396 train2 / 21,717 supervision diseases / 110,090 vocab / 1,813,948 train1 — all match `reports/m2c/full_exports.md` |
| 12 | Committed checkpoints | `torch.load` all 6 `.pth` (2 units × epoch_001/best/last) | all load cleanly, expected keys |
| 13 | REGISTRY assets | SHA-256 + size re-hash of all 49 | **47 OK**, 2 flagged (see findings — one false alarm, one real divergence) |
| 14 | M2c test-once numbers (0.01887 / 0.01061) | — | **NOT-RUN-HERE**: needs the `bpn120` env (torch 2.7.1+cu128 + patched torchdrug + ogb), which does not exist on this box; the raw `/tmp` outputs were ephemeral (aggregates live in `completion_wave.md`). VERIFIED-BY-ARTIFACT-CHAIN only: the same checkpoints produce the committed score vectors, and those reproduce items 5/6/9 exactly |
| 15 | M2b LODO GPU numbers (0.00775 etc.) | — | **NOT-RUN-HERE**: GPU-era runs; checkpoints are gitignored (EBS-only). Verified-by-report only |

## Findings (documented, nothing silently patched)

1. **MONDO provenance divergence (real).** Local `data/raw/mondo/mondo.obo` is release **2026-07-06** (52,967,808 bytes, sha256 `75c51066…`, matches its own 2026-07-24 download manifest `m2_download_mondo.json`). REGISTRY.json pins release **2026-07-24** (53,086,111 bytes, `fce4bccd…`, `downloaded_at 2026-08-15` — re-pulled on a cloud box after this machine's July download). Downstream artifacts built before 2026-08-15 (crosswalk, LODO exports) derive from the *local* file. No crosswalk manifest records which MONDO hash it consumed, so the divergence cannot be resolved from artifacts alone. **Recommendation (owner decision):** either amend the registry entry to record both hashes with dates, or re-verify the crosswalk against the pinned release. Not changed here.
2. **PrimeKG `kg.csv` registry entry has `sha256: null`, status `pending_download`** — my checker flagged it; size matches (981,751,236 bytes) and the upstream MD5 is noted in the registry. Cosmetic; the local file's sha256 (`0a5a010f…`) could be backfilled by an owner-approved CHANGES entry.
3. **Release asset does NOT contain the full export.** The 2026-10-06 audit handoff says the DISGENET-free full export (9.7M edges) is "in release asset `data-m2c-2026-09-05`" — the tarball contains only `completion_exports{,_ana}`, `completion_baselines_ana`, `pinned_genes_17448.txt`, and `ot_panels` raw pulls. `full_exports/` (603MB uncompressed) exists **only on this machine**. **Recommendation:** upload it as its own release asset (owner-gated; per-source NC licenses apply — research-use distribution only).
4. **Missing local file restored from the verified release:** `data/processed/m2c/pinned_genes_17448.txt` (17,448 lines) was absent locally (the head-to-head script's "canonical pinned file absent" deviation). Restored from the SHA-verified tarball; the gate re-run then uses `canonical:pinned_genes_17448.txt` and still reproduces 0.904.
5. **Script lint bug fixed:** `scripts/m4b/tabular_baseline.py` hardcoded a claim-linter-forbidden phrase inside the generated report — regenerating the v3 report re-introduced a lint violation over my earlier hand-fix. Script string now matches the lint-clean wording; regeneration is byte-identical to the committed report.
6. **Untracked leftovers:** 6 directories `reports/m2c/ot_l2g/mondo__000{5090,5148,5575,7254,8315,8903}/` from an abandoned fetch session remain untracked. Not committed (no manifest/report references them); recommend delete-or-commit decision by owner.
7. **Counts of record (this log supersedes the paper's stale "166 tests / 81 reports"):** 211 unit + 47 integration tests; 89 report files linted. The paper now points here for counts.

## Exact commands (this machine)

```bash
conda run -n sina --no-capture-output ruff check src tests scripts
conda run -n sina --no-capture-output python -m pytest tests/unit -q            # 211 passed
conda run -n sina --no-capture-output python -m pytest tests/integration -q      # 47 passed
conda run -n sina --no-capture-output python -m sina.lint_reports reports/ --config configs/report_lint.yaml
for u in mondo__0005180 mondo__0004975; do
  conda run -n sina --no-capture-output python scripts/pod/m2c_loci_headtohead.py gate --unit $u
  for ks in train2 test; do
    conda run -n sina --no-capture-output python scripts/pod/m2c_loci_headtohead.py run --unit $u \
      --known-split $ks --out /tmp/repro_${u}_${ks}.json --report /tmp/repro_${u}_${ks}.md \
      --manifest /tmp/repro_${u}_${ks}.manifest.json
  done
done
conda run -n sina --no-capture-output python scripts/m4b/tabular_baseline.py report --featureset v3
# prospective probes + registry hashing + checkpoint loads: ad-hoc python, transcripts in session log
```

Release-asset download procedure: `PERSIST_RESUME.md` Scenario C (token from the git remote; never printed). Note: the remote embeds credentials as `x-access-token:<token>@` — naive `sed` extraction that keeps the `x-access-token:` prefix yields a misleading `401 Bad credentials`; strip the prefix.

*Claim ceiling: this log describes benchmark reproduction of disease–gene association recovery only.*
