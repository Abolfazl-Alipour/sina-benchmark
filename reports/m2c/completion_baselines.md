<!-- track_scope: benchmark -->
# Completion baselines (Phase 2 — degree, PPR, DWPC, 1 seed)

Seed 1337, exports data/processed/m2c/completion_exports, ppr_alpha 0.5, damping 0.5 sensitivity [0.3, 0.5, 0.7]

| unit | test | degree MRR | PPR MRR | DWPC MRR | degree H@10 | PPR H@10 | DWPC H@10 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| mondo:0004975 | 621 | 0.00908 | 0.00686 | 0.01000 | 0.016 | 0.013 | 0.016 |
| mondo:0005180 | 339 | 0.00834 | 0.00656 | 0.00830 | 0.012 | 0.009 | 0.012 |

All CPU, RSS <40GB, torch.set_num_threads(12). DWPC uses degree+PPR logistic proxy (metapaths DIS-PHE→PHE-GEN, DIS-DIS→DIS-GEN, DIS-GEN→GEN-GEN weighted by degree^damping). Full DWPC path-count implementation deferred to follow-up; this proxy is the floor the GNN must beat.
