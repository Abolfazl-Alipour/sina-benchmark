# 0023 — Lean target-prioritization scope (10 primaries, 1 seed, Light Gate only)

- **change:** Scope Arms to **target prioritization only** (disease→gene, 1 seed 1337, 1 pass) on **10 primary diseases** (lean subset of the 13-primary freeze; 3 lowest-evidence primaries descriptive-only) + **Light Gate 0** mandatory. De-scope all secondary panels (16 diseases: 5 evidence + 11 low-context) from GPU (folders stay on disk, no training). De-scope drug repurposing direct `drug→disease` training (1,853 `INDICATION` rows insufficient; 5,310 `DRG-GEN` not a platform) — drug mapping stays as post-hoc `gene→drug` lookup on the 5k `TARGET` edges (descriptive, no GPU). Lock to 6-layer NBFNet (6 hops, `training_pool_and_hardware_plan.md`), 1 pass (overfit 0.0162→0.0109 observed), H200 batch ladder 300→100 probed on 1.43M. Bear-minimum GPU = Light Gate ~0.5h (~$1.50) + 10×~70min (~11.7h) = **~12.2h ≈ $34** at $2.80/H200 (worst 2× assumption: ~24h ≈ $68). Inside $900 cap. Step 0b 1-vs-2 pilot skipped by owner directive (lock to 1); re-enable via new CHANGES if Light Gate shows >2× pass-time.
- **rationale:** 3D organoid partner validates gene targets via CRISPR in organoids (not HTS drug screen); 192k `DIS-GEN` supervision is dense vs 1,853 drug indications. Grant narrative needs zero-shot gene finding for PD/AD/familial forms, not a sparse repurposing claim. Lean scope saves ~$135 vs 29×3 and ~$10 vs 13 primaries; 6-hop not 3-hop (correct per 11.67M-edge ×6 layers). Secondaries and TxGNN/foundation sidecar remain priced contingencies inside cap, not deleted.
- **requester:** owner — 2026-08-10 (conversation: "only the primary confirmatory diseases and the light gate", "how much bear minimum")
- **impact:** Amends CHANGES/0020 roster interpretation (13→10 graded, 3→descriptive) and v5.1 §6 pass-count (lock to 1) without touching frozen 29 export folders. Hands off: GPU clone needs `git pull` + `scripts/pod/cold_start.sh` (1-2h) then Light Gate ladder. Gate report `reports/m2b/gate0_light.md` still mandatory before Arms.
- **date:** 2026-08-10
- **approved:** owner — 2026-08-10

---

## Primary lean roster (10, graded)

Picked as the 10 highest-evidence primaries from the 13 (drop the 3 smallest `test` splits to save GPU while keeping family coverage):

- c1_progeria: `mondo:0010196` Werner (730 test), `mondo:0008310` HGPS (394), `mondo:0012853` Fontaine (201) — 3
- c3_pd: `mondo:0008199` late-onset PD (57), `mondo:0957576` PD25 (127) — 2 (young/juvenile PD were secondary, already out)
- c4_mody: `mondo:0007669` renal cysts+diabetes (258) — 1 (MODY 319 is secondary, out)
- c5_ndm: `mondo:0012192` PNDM panc+cerebellar (281), `mondo:0012522` TNDM3 (251), `mondo:0100164` PNDM (459), `mondo:0100165` PNDM1 (329) — 4

Dropped to descriptive: `mondo:0014157` mandibular hypo (52), `mondo:0010785` maternal diabetes+deafness (96) — smallest primaries, still built but not trained in lean scope. Restores to 13 via new CHANGES if budget allows.
