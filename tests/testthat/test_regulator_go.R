library(testthat)
library(SummarizedExperiment)
library(ggplot2)

# Mock data helper
# regulatorGSEA() works from a PerturbMatch result: it reads pert_symbol / pert_type
# from rowData and the "zscore" column of the enrichment assay. Real gene
# symbols are used so they map to Entrez IDs via org.Hs.eg.db. TP53 appears
# twice (a KD and an OE) to exercise per-gene aggregation.

make_regulator_object <- function() {
    mat <- cbind(
        SignedJaccard = c(-0.8, 0.4, 0.6, -0.3, 0.5, -0.2),
        zscore        = c(-2.0, 1.0, 1.5, -0.6, 1.2, -0.4)
    )
    rownames(mat) <- paste0("sig", 1:6)

    se <- SummarizedExperiment(
        assays  = list(enrichment = mat),
        rowData = DataFrame(
            pert_symbol = c("TP53", "TP53", "MYC", "EGFR", "STAT3", "JUN"),
            pert_type   = c("KD", "OE", "OE", "KO", "OE", "KD")
        )
    )

    new("PerturbMatch", se)
}

# .regulatorGeneStats

test_that(".regulatorGeneStats flips KD/KO and aggregates by mean", {
    pm <- make_regulator_object()

    stats <- .regulatorGeneStats(
        pm, aggregate = "mean", inferActivity = TRUE, org.Hs.eg.db::org.Hs.eg.db
    )

    # TP53 (Entrez 7157): KD zscore -2 flips to +2, OE zscore +1 stays;
    # mean = 1.5.
    expect_equal(unname(stats["7157"]), 1.5)
    expect_true(all(!is.na(names(stats))))
})

test_that(".regulatorGeneStats 'extreme' keeps the largest-magnitude score", {
    pm <- make_regulator_object()

    stats <- .regulatorGeneStats(
        pm, aggregate = "extreme", inferActivity = TRUE,
        org.Hs.eg.db::org.Hs.eg.db
    )

    # After the KD flip TP53 has scores {+2, +1}; the extreme is +2.
    expect_equal(unname(stats["7157"]), 2)
})

test_that(".regulatorGeneStats without inferActivity keeps raw signs", {
    pm <- make_regulator_object()

    stats <- .regulatorGeneStats(
        pm, aggregate = "mean", inferActivity = FALSE,
        org.Hs.eg.db::org.Hs.eg.db
    )

    # TP53 raw scores {-2, +1}; mean = -0.5.
    expect_equal(unname(stats["7157"]), -0.5)
})

test_that(".regulatorGeneStats errors without pert_symbol/pert_type", {
    mat <- cbind(x = 1:3, zscore = c(0.1, 0.2, 0.3))
    rownames(mat) <- paste0("sig", 1:3)
    se <- SummarizedExperiment(assays = list(enrichment = mat))
    pm <- new("PerturbMatch", se)

    expect_error(
        .regulatorGeneStats(
            pm, aggregate = "mean", inferActivity = TRUE,
            org.Hs.eg.db::org.Hs.eg.db
        ),
        "pert_symbol"
    )
})

# regulatorGSEA

test_that("regulatorGSEA returns a tidy GO table", {
    pm <- make_regulator_object()

    go <- suppressWarnings(regulatorGSEA(pm, ontology = "BP", minSize = 1))

    expect_s3_class(go, "data.frame")
    expect_setequal(
        colnames(go),
        c("ID", "term", "NES", "pval", "padj", "size", "leadingEdge")
    )
    expect_true(all(grepl("^GO:", go$ID)))
    # padj is sorted ascending
    expect_false(is.unsorted(go$padj))
})

# plotRegulatorGSEA

test_that("plotRegulatorGSEA returns a ggplot", {
    pm <- make_regulator_object()
    go <- suppressWarnings(regulatorGSEA(pm, ontology = "BP", minSize = 1))

    p <- plotRegulatorGSEA(go, pvalueCutoff = 1)

    expect_s3_class(p, "ggplot")
})

test_that("plotRegulatorGSEA errors when nothing passes the cutoff", {
    df <- data.frame(
        ID = "GO:0000001", term = "t", NES = 1, pval = 0.9,
        padj = 0.9, size = 3, leadingEdge = "TP53"
    )

    expect_error(plotRegulatorGSEA(df, pvalueCutoff = 0.05), "pass")
})
