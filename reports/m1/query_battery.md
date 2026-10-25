<!-- track_scope: benchmark -->
# M1 DuckDB query battery (spec §5.2.5)

**Status: DRAFT — reviewer-auditable evidence (Gate 1)**

Identical battery over DuckDB views of each track's canonical tables.

## Track-R (PrimeKG)

- nodes: **129,312** (10 types)
- edges: **8,100,498** (30 relations)
- evidence groups: **4,050,064**
  (groups with exactly 2 rows: 4,049,879 —
  on Track-R these are reciprocal twins, because kg.csv stores both directions
  of every edge; on Track-D undirected edges are stored once (§5.2.3), so
  0 is the expected value)
- FK integrity (dangling endpoints): **0**

Top relations:

| relation | edges |
|---|---|
| anatomy_protein_present | 3036406 |
| drug_drug | 2672628 |
| protein_protein | 642150 |
| disease_phenotype_positive | 300634 |
| bioprocess_protein | 289610 |
| cellcomp_protein | 166804 |
| disease_protein | 160822 |
| molfunc_protein | 139060 |
| drug_effect | 129568 |
| bioprocess_bioprocess | 105772 |

Top-degree nodes:

| node_id | name | degree |
|---|---|---|
| uberon:2106 | spleen | 33200 |
| uberon:1043 | esophagus | 32952 |
| uberon:2108 | small intestine | 33710 |
| uberon:468 | multicellular organism | 34738 |
| uberon:3889 | fallopian tube | 33484 |
| uberon:9834 | dorsolateral prefrontal cortex | 33026 |
| uberon:473 | testis | 33486 |
| uberon:160 | intestine | 33084 |
| uberon:945 | stomach | 32946 |
| uberon:2367 | prostate gland | 33238 |

## Track-D (OptimusKG)

- nodes: **190,531** (10 types)
- edges: **21,813,816** (36 relations)
- evidence groups: **21,813,816**
  (groups with exactly 2 rows: 0 —
  on Track-R these are reciprocal twins, because kg.csv stores both directions
  of every edge; on Track-D undirected edges are stored once (§5.2.3), so
  0 is the expected value)
- FK integrity (dangling endpoints): **0**

Top relations:

| relation | edges |
|---|---|
| ASSOCIATED_WITH | 10531730 |
| EXPRESSION_PRESENT | 6616463 |
| EXPRESSION_ABSENT | 2171492 |
| SYNERGISTIC_INTERACTION | 1341086 |
| INTERACTS_WITH | 734862 |
| PHENOTYPE_PRESENT | 157144 |
| PARENT | 95711 |
| IS_A | 61720 |
| INDICATION | 58690 |
| CONTRAINDICATION | 19997 |

Top-degree nodes:

| node_id | name | degree |
|---|---|---|
| uberon:0000178 | blood | 55767 |
| uberon:0001323 | tibial nerve | 54070 |
| uberon:0001388 | gastrocnemius | 54608 |
| uberon:0001416 | skin of abdomen | 54266 |
| uberon:0001511 | skin of leg | 55211 |
| uberon:0002190 | subcutaneous adipose tissue | 53883 |
| uberon:0004533 | left testis | 55841 |
| uberon:0007610 | tibial artery | 54183 |
| uberon:0010414 | omental fat pad | 54269 |
| uberon:0035833 | lower esophagus muscularis layer | 54458 |
