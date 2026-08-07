# perturbMatch

[![Test R-universe](https://github.com/soulj/perturbMatch/actions/workflows/r-universe.yml/badge.svg)](https://github.com/soulj/perturbMatch/actions/workflows/r-universe.yml)

**Match a transcriptomics signature to single-gene perturbation reference databases.**

`perturbMatch` matches a query gene signature (e.g. a differential-expression
result) against a reference database of perturbation signatures and scores
which perturbed regulators best explain your data, using cosine, GSEA and
signed-Jaccard similarity.

## Installation

```r
# Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Development version from GitHub
BiocManager::install("soulj/perturbMatch")
```
The vignettes additionally need packages listed under `Suggests`, which are
not installed automatically:

```r
BiocManager::install(c(
  # main vignette
  "airway", "DESeq2", "GO.db", "tibble",
  # single-cell vignette
  "muscData", "SingleCellExperiment", "scuttle", "scater", "singscore",
  "scales"
))
```

## Getting started

```r
library(perturbMatch)

# Load the reference data. Downloaded and cached on first use (0.7-3.1 GB)
# add demo = TRUE for the small bundled demonstration database instead
ref <- getReferenceDatabase()

# Turn a DE result (DESeq2 / limma / topTable) into a ranked query
query <- prepareQuery(dds, reference = ref)

# Match the query against the perturbation reference
result <- calcEnrichment(ref, expressionMatrix = query, method = "cosine")

# Visualise the top mimicking / opposing regulators
plotSimilarity(result, topN = 20)

# GO enrichment of the mimicking / reversing regulators
plotRegulatorGSEA(regulatorGSEA(result))
```
## Learn more

- Full walkthrough: the package vignette (`browseVignettes("perturbMatch")`).
- Function help: `?calcEnrichment`, `?prepareQuery`, `?plotSimilarity`.
