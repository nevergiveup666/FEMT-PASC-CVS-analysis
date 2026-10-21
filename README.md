# FEMT–PASC-CVS Network Pharmacology Analysis

This repository contains the R code used to generate all computational
results reported in the manuscript:

> **Five-Element Music Therapy targets overlap with severe acute COVID-19
> transcriptomic signatures: A hypothesis-generating network-pharmacology
> framework for candidate gene discovery relevant to PASC-CVS**

The study is purely computational and hypothesis-generating. No causal
relationships or clinical efficacy of Five-Element Music Therapy (FEMT)
for post-acute sequelae of SARS-CoV-2 infection with cardiovascular
manifestations (PASC-CVS) can be inferred from these results. All
predicted candidate genes require prospective validation in independent
PASC-CVS cohorts and wet-laboratory experiments.

---

## 1. Repository contents

| File | Description |
|------|-------------|
| `FEMT_PASC_CVS_analysis.R` | Main analysis script (all steps, Figures 1–6, Tables 1–2, Supplementary Tables S1–S3) |
| `GSE171110_series_matrix.txt` | GEO series matrix for GSE171110, used only to derive sample groups |
| `README.md` | This file |

The two expression datasets are **not** redistributed here. They are
publicly available from GEO (see Section 3) and should be downloaded by
the user.

---

## 2. Runtime environment

- R version: **4.2.2**
- Key packages (versions used in the manuscript):

| Package | Version |
|---------|---------|
| limma | 3.50.3 |
| clusterProfiler | 4.4.4 |
| org.Hs.eg.db | 3.16 |
| pROC | 1.18.0 |
| caret | 6.0-94 |
| igraph | 1.3.5 |
| tidyverse | ≥ 1.3.2 |
| ggvenn | ≥ 0.1.10 |
| glmnet, randomForest, patchwork, ggrepel, RColorBrewer, gridExtra | any recent version |

Install missing packages from CRAN / Bioconductor, e.g.:

```r
install.packages(c("tidyverse", "pROC", "caret", "igraph", "ggvenn",
                   "glmnet", "randomForest", "patchwork", "ggrepel",
                   "RColorBrewer", "gridExtra"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("limma", "clusterProfiler", "AnnotationDbi",
                       "org.Hs.eg.db"))
