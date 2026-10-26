# SINA Benchmark

A leakage-certified, leave-one-disease-out (LODO) benchmark for disease-associated
gene recovery on OptimusKG, with a graph neural network (BioPathNet/NBFNet)
evaluated against same-graph degree, personalized PageRank, and counting nulls.

**Headline result (a certified negative):** under strict LODO zero-shot conditions
the model reaches mean MRR 0.00775 across 13 primary diseases (3 seeds, Wilcoxon
p = 0.34, not significant); degree-preserving shuffle controls collapse 58 to 89
percent, certifying the evaluation leakage-free. On a seen-disease completion
track (Parkinson's, Alzheimer's) the model ties a supervision-only counting null
on Hits@10. Prospective prediction of post-2020 GWAS associations collapses for
model and null alike. On this supervision, the ceiling is the labels, not the model.

The full manuscript is `reports/SINA_PAPER.md`. Per-locus head-to-head results
against Open Targets L2G are in `reports/m2c/` (train-knowns and leakage-controlled
test-knowns panels, with committed score vectors in `reports/m2c/ot_panels/`).

## Reproduce

```bash
conda env create -f environment.lock.yml
python -m pytest tests/unit -q
python -m sina.lint_reports reports/ --config configs/report_lint.yaml
```

## Data and licenses

Code: MIT (see LICENSE). Data is NOT redistributed here: OptimusKG
(doi:10.7910/DVN/IYNGEV) and PrimeKG (doi:10.7910/DVN/IXA7BM, CC BY-NC-SA 4.0)
carry per-source licenses, some non-commercial. Download them from their DOIs;
`data/raw/REGISTRY.json` pins the exact assets by SHA-256.

## Frozen artifacts (rev4)

Release `data-rev4-2026-10-25`
([tarball](https://github.com/Abolfazl-Alipour/sina-benchmark/releases/download/data-rev4-2026-10-25/sina-benchmark-data-rev4-2026-10-25.tar.gz),
[MANIFEST.json](https://github.com/Abolfazl-Alipour/sina-benchmark/releases/download/data-rev4-2026-10-25/MANIFEST.json),
60,742,782 bytes): KEEP background graph (1,431,945 edges), LODO roster plus
per-unit splits for all 29 units, and completion-track splits for
Parkinson disease and Alzheimer disease. Derived splits and edge list only;
no raw source content. Every file is SHA-256 pinned in the manifest.

## Citation

Alipour, A. SINA: A Leakage-Certified, Pathway-Aware Benchmark for
Disease-Associated Gene Recovery on OptimusKG. Preprint, 2026.
