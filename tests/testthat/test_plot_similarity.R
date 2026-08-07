library(testthat)
library(SummarizedExperiment)
library(ggplot2)

# Mock data helper
# plotSimilarity() dispatches on SummarizedExperiment and calls
# perturbMatch::result(), which is only defined for "PerturbMatch" objects, so the mock
# must be a PerturbMatch object with an "enrichment" assay (first column becomes
# "score", remaining columns such as "zscore" are kept as-is).

make_plot_object <- function() {
    mat <- cbind(
        SignedJaccard = c(0.9, -0.7, -0.4, 0.5, -0.2),
        zscore        = c(1.8, -1.4, -0.8, 1.0, -0.4)
    )
    rownames(mat) <- paste0("sig", 1:5)

    se <- SummarizedExperiment(
        assays = list(enrichment = mat),
        rowData = DataFrame(
            pert_symbol = c("GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E"),
            pert_type   = c("OE", "KD", "KO", "OE", "KD")
        )
    )

    new("PerturbMatch", se)
}

# .inferActivityScore

test_that(".inferActivityScore flips sign for KD/KO but not OE", {
    df <- data.frame(
        pert_symbol = c("A", "B", "C", "D"),
        pert_type   = c("KD", "OE", "KO", "OE"),
        zscore      = c(-1, 2, -3, 0.5)
    )

    flipped <- .inferActivityScore(df)

    expect_equal(flipped$zscore, c(1, 2, 3, 0.5))
})

test_that(".inferActivityScore is case-insensitive for pert_type", {
    df <- data.frame(pert_symbol = "A", pert_type = "kd", zscore = -2)

    flipped <- .inferActivityScore(df)

    expect_equal(flipped$zscore, 2)
})

# .preparePlotData

test_that(".preparePlotData sorts by zscore and labels direction", {
    df <- data.frame(
        pert_symbol = c("A", "B", "C", "D"),
        pert_type   = c("OE", "OE", "OE", "OE"),
        zscore      = c(-1, 2, -3, 0.5)
    )

    res <- .preparePlotData(df, topN = 10, inferActivity = FALSE)

    expect_equal(res$pert_symbol, c("B", "D", "A", "C"))
    expect_equal(res$Direction, c(
        "Mimic (Positive)", "Mimic (Positive)",
        "Reverse (Negative)", "Reverse (Negative)"
    ))
})

test_that(".preparePlotData restricts to top and bottom n when there are many rows", {
    df <- data.frame(
        pert_symbol = paste0("G", 1:10),
        pert_type   = "OE",
        zscore      = c(5, 4, 3, 2, 1, -1, -2, -3, -4, -5)
    )

    res <- .preparePlotData(df, topN = 2, inferActivity = FALSE)

    expect_equal(nrow(res), 4)
    expect_equal(res$pert_symbol, c("G1", "G2", "G9", "G10"))
})

test_that(".preparePlotData keeps all rows when there are fewer than 2 * topN", {
    df <- data.frame(
        pert_symbol = paste0("G", 1:4),
        pert_type   = "OE",
        zscore      = c(2, 1, -1, -2)
    )

    res <- .preparePlotData(df, topN = 5, inferActivity = FALSE)

    expect_equal(nrow(res), 4)
})

test_that(".preparePlotData applies activity inference before sorting", {
    df <- data.frame(
        pert_symbol = c("A", "B"),
        pert_type   = c("KD", "OE"),
        zscore      = c(-5, 1)
    )

    res <- .preparePlotData(df, topN = 10, inferActivity = TRUE)

    # A's score flips from -5 to 5, so it now ranks above B
    expect_equal(res$pert_symbol, c("A", "B"))
    expect_equal(res$zscore, c(5, 1))
})

# plotSimilarity

test_that("plotSimilarity returns a ggplot object", {
    pm <- make_plot_object()

    p <- plotSimilarity(pm, topN = 3)

    expect_s3_class(p, "ggplot")
})

test_that("plotSimilarity labels the y-axis for transcriptional similarity by default", {
    pm <- make_plot_object()

    p <- plotSimilarity(pm, topN = 3)

    expect_equal(p$labels$y, "Similarity Z Score")
})

test_that("plotSimilarity labels the y-axis for inferred activity when requested", {
    pm <- make_plot_object()

    p <- plotSimilarity(pm, topN = 3, inferActivity = TRUE)

    expect_equal(p$labels$y, "Inferred Activity Z Score")
})

test_that("plotSimilarity respects topN", {
    pm <- make_plot_object()

    p <- plotSimilarity(pm, topN = 1)
    pdat <- ggplot2::layer_data(p, 3) # geom_point layer

    expect_equal(nrow(pdat), 2) # top 1 positive + top 1 negative
})
