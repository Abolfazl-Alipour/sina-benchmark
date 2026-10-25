<!-- track_scope: benchmark -->
# Anchored pathway table vs topology null (CHANGES/0025 control analysis)

Question: does the anchored EXPLORATORY pathway table carry learned signal beyond graph topology? Motivation: anchored Recall@10 exactly tied the degree-null (HGPS 0.0069 vs 0.0069) and ORA (Werner 0.0043 vs 0.0043) baselines, and the within-disease shuffle control left anchored Recall@10 nearly intact (0.0069 -> 0.0062) while gene-level MRR collapsed 58-89%. Null construction: the anchored aggregation formula with every path weight replaced by 1 — score_null(P) = support(P) / ln(degree(P)) — over the identically rebuilt candidate set (direct train1 annotation edge to a consensus Top-100 gene, support >= 3 prefilter), graded on the identical test-genes-only GO BP/Reactome truth. Anchored per-candidate weights are not persisted anywhere (checked data/processed/m3/ and reports/m2b/ — only the Top-10 tables survive in the funnel reports), so the comparison is metric-level plus Top-10 Jaccard, not per-candidate correlation. This is a control analysis; it supports no recovery claim either way.

## mondo:0010196 (mondo__0010196)

Candidates: 2578 with a direct train1 annotation edge to a
consensus Top-100 gene (100 genes), prefiltered to
support >= 3 -> 384; truth =
2085 pathways from 730 held-out test genes
(>= 3 test genes each).

| ranking | Recall@10 | Recall@5 |
|---|---|---|
| anchored (EXPLORATORY table) | 0.0043 | 0.0024 |
| topology null support/ln(degree) | 0.0043 | 0.0024 |

Top-10 Jaccard (anchored vs null): 0.82
(stated in funnel report: R@10 0.0043, R@5 0.0024; recomputed here from the parsed Top-10 on the identically
reconstructed truth). Absolute Recall@10 hits: anchored 9
of 2085 true pathways, null 9.

Null Top-10: go:0010628, go:0045944, go:0010629, go:0007165, go:0008285, go:0000122, go:0045893, go:0043066, go:0043410, go:0060736

**Verdict: anchored ≈ topology null — consistent with within-shuffle survival.**

## mondo:0008310 (mondo__0008310)

Candidates: 2220 with a direct train1 annotation edge to a
consensus Top-100 gene (100 genes), prefiltered to
support >= 3 -> 301; truth =
1459 pathways from 394 held-out test genes
(>= 3 test genes each).

| ranking | Recall@10 | Recall@5 |
|---|---|---|
| anchored (EXPLORATORY table) | 0.0069 | 0.0034 |
| topology null support/ln(degree) | 0.0062 | 0.0034 |

Top-10 Jaccard (anchored vs null): 0.54
(stated in funnel report: R@10 0.0069, R@5 0.0034; recomputed here from the parsed Top-10 on the identically
reconstructed truth). Absolute Recall@10 hits: anchored 10
of 1459 true pathways, null 9.

Null Top-10: go:0010628, go:0006629, reactome:R-HSA-6785807, go:0045944, go:0007165, go:0010629, go:0006955, go:0045893, reactome:R-HSA-381426, go:0010573

**Verdict: anchored adds signal beyond topology.**
