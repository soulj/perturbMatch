# perturbMatch

[![Test R-universe](https://github.com/soulj/perturbMatch/actions/workflows/r-universe.yml/badge.svg)](https://github.com/soulj/perturbMatch/actions/workflows/r-universe.yml)

**Match a transcriptomics signature to single-gene perturbation reference databases.**

`perturbMatch` matches a query gene signature (e.g. a differential-expression
result) against a reference database of perturbation signatures and scores
which perturbed regulators best explain your data, using cosine, GSEA and
signed-Jaccard similarity.

<img src="man/figures/README-similarity.png" alt="Lollipop plot ranking the top mimicking and reversing gene perturbations by similarity Z score" width="600" />

*Regulators whose perturbation signature most resembles (orange) or opposes
(blue) the query, with marker shape showing gain- or loss-of-function.*

## Installation

These packages are required for building the vignettes, that are not installed automatically:

```r
# Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  # main vignette
  "airway", "DESeq2", "GO.db", "tibble",
  # single-cell vignette
  "muscData", "SingleCellExperiment", "scuttle", "scater", "singscore",
  "scales"
))
```

```r
# Development version from GitHub
BiocManager::install("soulj/perturbMatch", build_vignettes = TRUE)
```

## Getting started

```r
library(perturbMatch)
library(airway)
library(DESeq2)

# A query signature: dexamethasone-treated versus untreated airway smooth
# muscle cells
data(airway)
de <- DESeqDataSet(airway, design = ~ cell + dex) |>
    DESeq() |>
    results(contrast = c("dex", "trt", "untrt"))

# The 21-signature demonstration reference bundled with the package.
# getReferenceDatabase() with no arguments fetches a full reference instead,
# downloaded and cached on first use (0.7-3.1 GB)
ref <- getReferenceDatabase(demo = TRUE)

# Map the query onto the reference, then score every signature against it
query <- prepareQuery(de, reference = ref)
result <- calcEnrichment(ref, expressionMatrix = query, method = "cosine")

# The best-matching regulators, and a lollipop plot of them
cols <- c("pert_symbol", "pert_type", "score", "zscore")
topRegulators(result, n = 5)[, cols]
plotSimilarity(result, topN = 20)
```

The demonstration reference is too small to nominate real regulators.
Swap in `getReferenceDatabase()` for that. `prepareQuery()` also takes a
`limma` `topTable()` or any data frame with gene identifiers and a
log2 fold-change column, and `calcEnrichment()` accepts plain up- and
down-regulated gene lists instead of a ranked query. GO enrichment of the
matched regulators is available through `regulatorGSEA()` and
`plotRegulatorGSEA()`.

## Learn more

- Full walkthrough:
  [bulk RNA-seq](https://biocstaging.r-universe.dev/articles/perturbMatch/perturbMatch.html),
  [single-cell](https://biocstaging.r-universe.dev/articles/perturbMatch/perturbMatch_singlecell_vignette.html)
- Function help: `?calcEnrichment`, `?prepareQuery`, `?plotSimilarity`.

## Citation

If you use perturbMatch, please cite:

> Soul J, Young DA (2026). Automated generation of a gene perturbation
> transcriptomic atlas using large language models. *bioRxiv*.
> [doi:10.64898/2026.08.08.743502](https://doi.org/10.64898/2026.08.08.743502)
