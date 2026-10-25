<!-- track_scope: benchmark -->
# Full-graph product export (single-train, no holdouts)

Generated 2026-09-12T01:32:54Z.
Supervision: 9709396 DIS-GEN / 21717 diseases. Background: 1813948 (GEN-GEN, PHE-GEN, BPO-GEN, CCO-GEN, MFN-GEN, PWY-GEN, DRG-GEN, DRG-DIS, DIS-PHE, DIS-DIS).
Vocab: 110090 entities (32074 genes, 36121 diseases).
DISGENET-tagged edges removed from every kept label (counts in manifest; DISGENET-only-unique-loss is zero by construction — all tagged edges co-sourced).
valid.txt/test.txt empty by design: one model trains on everything; evaluation is recovery disclosure, not holdout.
