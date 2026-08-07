library(testthat)
library(ggplot2)
library(dplyr)

#' Minimal 3-column data frame for testing
make_data <- function(n = 100, seed = 1) {
    set.seed(seed)
    data.frame(
        gene = paste0("Gene", seq_len(n)),
        logFC = rnorm(n, sd = 2),
        padj = rbeta(n, 0.5, 5),
        stringsAsFactors = FALSE
    )
}

test_that("volcanoPlot errors on fewer than 3 columns", {
    bad <- data.frame(gene = letters, logFC = seq_along(letters))
    expect_error(volcanoPlot(bad), "ID, foldchange and FDR columns")
})

test_that("volcanoPlot errors when fold-change column is non-numeric", {
    bad <- data.frame(
        gene = letters[1:5],
        logFC = letters[1:5], # character, not numeric
        padj = rep(0.01, 5),
        stringsAsFactors = FALSE
    )
    expect_error(volcanoPlot(bad), "numeric fold changes")
})

test_that("volcanoPlot errors when FDR column is non-numeric", {
    bad <- data.frame(
        gene = letters[1:5],
        logFC = seq_len(5),
        padj = letters[1:5], # character, not numeric
        stringsAsFactors = FALSE
    )
    expect_error(volcanoPlot(bad), "numeric FDRs")
})

test_that("volcanoPlot accepts a data.frame-like input (tibble)", {
    df <- tibble::tibble(
        gene  = paste0("G", 1:50),
        logFC = rnorm(50),
        padj  = runif(50)
    )
    expect_s3_class(volcanoPlot(df), "ggplot")
})

test_that("volcanoPlot returns a ggplot object", {
    p <- volcanoPlot(make_data())
    expect_s3_class(p, "ggplot")
})

test_that("genes above thresholds are classed Up-regulated", {
    df <- data.frame(
        gene = "GeneA",
        logFC = 2, # > log2(1.5)
        padj = 0.001, # < 0.05
        stringsAsFactors = FALSE
    )
    p <- volcanoPlot(df, numberPoints = 1)
    pdat <- ggplot2::layer_data(p, 1) # geom_point layer
    # Up-regulated → firebrick3
    expect_true(any(pdat$colour == "firebrick3"))
})

test_that("genes below negative threshold are classed Down-regulated", {
    df <- data.frame(
        gene = "GeneB",
        logFC = -2,
        padj = 0.001,
        stringsAsFactors = FALSE
    )
    p <- volcanoPlot(df)
    pdat <- ggplot2::layer_data(p, 1)
    expect_true(any(pdat$colour == "dodgerblue3"))
})

test_that("non-significant genes are coloured black (Unchanged)", {
    df <- data.frame(
        gene  = "GeneC",
        logFC = 2,
        padj  = 0.9
    )
    p <- volcanoPlot(df)
    pdat <- ggplot2::layer_data(p, 1)
    expect_true(any(pdat$colour == "black"))
})


test_that("rows with NA padj are silently dropped", {
    df <- make_data(20)
    df$padj[c(1, 5, 10)] <- NA
    expect_no_error(volcanoPlot(df))
    p <- volcanoPlot(df)
    pdat <- ggplot2::layer_data(p, 1)
    expect_equal(nrow(pdat), 17L)
})

test_that("p-values of zero (Inf -log10) do not cause errors", {
    df <- make_data(10)
    df$padj[1] <- 0 # will produce Inf after -log10
    expect_no_error(volcanoPlot(df))
    p <- volcanoPlot(df)
    pdat <- ggplot2::layer_data(p, 1)
    expect_false(any(is.infinite(pdat$y)))
})

test_that("numberPoints = 0 produces no text-repel layer", {
    p <- volcanoPlot(make_data(), numberPoints = 0)
    layers <- sapply(p$layers, function(l) class(l$geom)[1])
    expect_false("GeomTextRepel" %in% layers)
})

test_that("numberPoints > 0 adds a text-repel layer", {
    p <- volcanoPlot(make_data(), numberPoints = 3)
    layers <- sapply(p$layers, function(l) class(l$geom)[1])
    expect_true("GeomTextRepel" %in% layers)
})

test_that("numberPoints = NULL suppresses automatic labels without erroring", {
    p <- volcanoPlot(make_data(), numberPoints = NULL)
    layers <- sapply(p$layers, function(l) class(l$geom)[1])
    expect_false("GeomTextRepel" %in% layers)
})

test_that("numberPoints = NULL still allows selectedPoints to be labelled", {
    p <- volcanoPlot(make_data(), numberPoints = NULL, selectedPoints = c("Gene1", "Gene2"))
    layers <- sapply(p$layers, function(l) class(l$geom)[1])
    expect_true("GeomTextRepel" %in% layers)
})


test_that("stricter logFCThreshold reduces up/down-regulated count", {
    set.seed(7)
    df <- make_data(200)

    p_loose <- volcanoPlot(df, logFCThreshold = 0)
    p_strict <- volcanoPlot(df, logFCThreshold = 2)

    count_regulated <- function(p) {
        pdat <- ggplot2::layer_data(p, 1)
        sum(pdat$colour %in% c("firebrick3", "dodgerblue3"))
    }

    expect_gt(count_regulated(p_loose), count_regulated(p_strict))
})

test_that("stricter fdrThreshold reduces up/down-regulated count", {
    set.seed(8)
    df <- make_data(200)

    p_loose <- volcanoPlot(df, fdrThreshold = 0.9)
    p_strict <- volcanoPlot(df, fdrThreshold = 1e-5)

    count_regulated <- function(p) {
        pdat <- ggplot2::layer_data(p, 1)
        sum(pdat$colour %in% c("firebrick3", "dodgerblue3"))
    }

    expect_gt(count_regulated(p_loose), count_regulated(p_strict))
})


#' Minimal SummarizedExperiment for testing
make_se <- function(n_genes = 50, n_samples = 3, seed = 1) {
    set.seed(seed)
    logFC_mat <- matrix(
        rnorm(n_genes * n_samples, sd = 2),
        nrow = n_genes,
        dimnames = list(paste0("Gene", seq_len(n_genes)), paste0("S", seq_len(n_samples)))
    )
    padj_mat <- matrix(
        rbeta(n_genes * n_samples, 0.5, 5),
        nrow     = n_genes,
        dimnames = dimnames(logFC_mat)
    )
    SummarizedExperiment::SummarizedExperiment(
        assays = list(logFC = logFC_mat, padj = padj_mat)
    )
}

test_that("plotVolcano returns a ggplot", {
    se <- make_se()
    p <- plotVolcano(se, signature = "S1")
    expect_s3_class(p, "ggplot")
})

test_that("plotVolcano errors on missing assay", {
    se <- make_se()
    expect_error(
        plotVolcano(se, signature = "S1", pvalueAssay = "nonexistent"),
        "requires assays"
    )
})

test_that("plotVolcano errors on missing signature", {
    se <- make_se()
    expect_error(
        plotVolcano(se, signature = "S99"),
        "Signature not found"
    )
})

test_that("plotVolcano respects custom assay names", {
    set.seed(2)
    n <- 30
    mat <- matrix(rnorm(n * 2), nrow = n,
        dimnames = list(paste0("G", seq_len(n)), c("A", "B")))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(effect = mat, pval = abs(mat) / max(abs(mat)))
    )
    p <- plotVolcano(
        se,
        signature   = "A",
        effectAssay = "effect",
        pvalueAssay = "pval"
    )
    expect_s3_class(p, "ggplot")
})
