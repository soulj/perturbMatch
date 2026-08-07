library(testthat)
library(SummarizedExperiment)
library(DelayedArray)
library(BiocParallel)

# Mock data helpers
make_test_se <- function() {
    mat <- matrix(
        c(
            10, 2,
            9,  3,
            8,  4,
            1,  7,
            2,  8,
            3,  9
        ),
        nrow = 6,
        byrow = TRUE
    )

    rownames(mat) <- paste0("gene", 1:6)
    colnames(mat) <- c("sig1", "sig2")

    se <- SummarizedExperiment(
        assays = list(chrdir = mat),
        colData = DataFrame(
            pert_id     = c("KO1", "OE2"),
            pert_iname  = c("KO1", "OE1"),
            pert_type   = c("KO", "OE"),
            pert_symbol = c("A", "B"),
            sig_id      = c("s1", "s2"),
            technology  = c("RNASeq", "RNASeq"),
            PGR         = c("0.92", "0.6")
        )
    )

    S4Vectors::metadata(se) <- list(ref_metric = "chrdir", default_metric = "chrdir")

    se
}

make_expression_matrix <- function() {
    data.frame(
        gene = paste0("gene", 1:6),
        logFC = c(5, 4, 3, -3, -4, -5),
        stringsAsFactors = FALSE
    )
}

# Reference with two selectable metrics (mirrors a limma-style database:
# ref_metric lists everything available, default_metric picks the fallback)
make_limma_se <- function() {
    log2fc <- matrix(
        c(5, 4, 3, -3, -4, -5, 1, 1, 1, -1, -1, -1),
        nrow = 6,
        dimnames = list(paste0("gene", 1:6), c("sig1", "sig2"))
    )
    padj <- matrix(0.01, nrow = 6, ncol = 2, dimnames = dimnames(log2fc))
    pi <- log2fc * -log10(padj)

    se <- SummarizedExperiment(
        assays  = list(log2fc = log2fc, pi = pi),
        colData = DataFrame(pert_id = c("KO1", "OE2"))
    )
    S4Vectors::metadata(se) <- list(ref_metric = c("pi", "log2fc"), default_metric = "pi")

    se
}

# Validation tests

test_that("calcEnrichment validates gene set inputs", {
    se <- make_test_se()

    expect_error(
        calcEnrichment(
            se,
            upRegulated = 1:3,
            downRegulated = c("gene4", "gene5"),
            background = paste0("gene", 1:6)
        ),
        "upRegulated must be a character vector"
    )

    expect_error(
        calcEnrichment(
            se,
            upRegulated = c("gene1", "gene2"),
            downRegulated = 1:3,
            background = paste0("gene", 1:6)
        ),
        "downRegulated must be a character vector"
    )

    expect_error(
        calcEnrichment(
            se,
            upRegulated = c("gene1", "gene2"),
            downRegulated = c("gene4", "gene5"),
            background = 1:6
        ),
        "background must be a character vector"
    )
})


test_that("expressionMatrix is required for cosine methods", {
    se <- make_test_se()

    expect_error(
        calcEnrichment(
            se,
            upRegulated = c("gene1", "gene2"),
            downRegulated = c("gene4", "gene5"),
            background = paste0("gene", 1:6),
            method = "cosine"
        ),
        "An expressionMatrix with score values must be provided"
    )
})


test_that("calcEnrichment errors when gene sets do not overlap the reference", {
    se <- make_test_se()

    expect_error(
        calcEnrichment(
            se,
            upRegulated = c("notagene1", "notagene2"),
            downRegulated = c("gene5", "gene6"),
            background = paste0("gene", 1:6),
            method = "SignedJaccard"
        ),
        "None of the 'upRegulated' genes were found in the reference object"
    )

    expect_error(
        calcEnrichment(
            se,
            upRegulated = c("gene1", "gene2"),
            downRegulated = c("notagene5", "notagene6"),
            background = paste0("gene", 1:6),
            method = "SignedJaccard"
        ),
        "None of the 'downRegulated' genes were found in the reference object"
    )

    expect_error(
        calcEnrichment(
            se,
            upRegulated = c("gene1", "gene2"),
            downRegulated = c("gene5", "gene6"),
            background = c("notagene1", "notagene2"),
            method = "SignedJaccard"
        ),
        "None of the 'background' genes were found in the reference object"
    )
})


test_that("calcEnrichment validates topN, topSig, and chunkSize are positive whole numbers", {
    se <- make_test_se()

    expect_error(
        calcEnrichment(se, upRegulated = c("gene1"), downRegulated = c("gene5"),
            background = paste0("gene", 1:6), method = "SignedJaccard", topN = -5),
        "'topN' must be a single positive whole number"
    )

    expect_error(
        calcEnrichment(se, expressionMatrix = make_expression_matrix(),
            method = "cosine", topSig = 0),
        "'topSig' must be a single positive whole number"
    )

    expect_error(
        calcEnrichment(se, upRegulated = c("gene1"), downRegulated = c("gene5"),
            background = paste0("gene", 1:6), method = "SignedJaccard", chunkSize = 0),
        "'chunkSize' must be a single positive whole number"
    )
})


test_that("calcEnrichment validates threshold is non-negative", {
    se <- make_test_se()

    expect_error(
        calcEnrichment(se, upRegulated = c("gene1"), downRegulated = c("gene5"),
            background = paste0("gene", 1:6), method = "SignedJaccard", threshold = -1),
        "'threshold' must be NULL or a single non-negative number"
    )
})

# Functional tests

test_that("calcEnrichment returns a PerturbMatch object for SignedJaccard", {
    se <- make_test_se()

    result <- calcEnrichment(
        se,
        upRegulated = c("gene1", "gene2"),
        downRegulated = c("gene5", "gene6"),
        background = paste0("gene", 1:6),
        method = "SignedJaccard",
        topN = 2,
        BPPARAM = SerialParam(progressbar = FALSE),
        chunkSize = 1
    )

    expect_s4_class(result, "PerturbMatch")
    expect_equal(nrow(assay(result)), 2)
    expect_true("SignedJaccard" %in% colnames(assay(result)))
    expect_true("overlap" %in% colnames(assay(result)))
})


test_that("calcEnrichment records the parameter that defined the gene sets", {
    se <- make_test_se()

    score <- function(...) {
        calcEnrichment(
            se,
            upRegulated = c("gene1", "gene2"),
            downRegulated = c("gene5", "gene6"),
            background = paste0("gene", 1:6),
            method = "SignedJaccard",
            BPPARAM = SerialParam(progressbar = FALSE),
            chunkSize = 1,
            ...
        )
    }

    byThreshold <- score(threshold = 5)
    expect_equal(metadata(byThreshold)$threshold, 5)
    expect_null(metadata(byThreshold)$topN)

    byTopN <- score(topN = 2)
    expect_equal(metadata(byTopN)$topN, 2)
    expect_null(metadata(byTopN)$threshold)
})


test_that("calcEnrichment returns expected columns for GSEA", {
    se <- make_test_se()

    result <- calcEnrichment(
        se,
        upRegulated = c("gene1", "gene2"),
        downRegulated = c("gene5", "gene6"),
        background = paste0("gene", 1:6),
        method = "GSEA",
        BPPARAM = SerialParam(progressbar = FALSE),
        chunkSize = 1
    )

    expect_s4_class(result, "PerturbMatch")

    expect_true(all(
        c("GSEAup", "GSEAdown", "GSEAScore") %in%
            colnames(assay(result))
    ))
})


test_that("GSEA ranks on the combined GSEAScore, not up-enrichment alone", {
    # Two signatures with identical up-enrichment but opposite down-enrichment:
    # only the combined GSEAScore (up - down) separates them. The primary
    # 'score' must reflect that, not GSEAup.
    genes <- paste0("g", 1:12)
    mat <- cbind(
        # query-down genes (g11, g12) are UP here -> discordant, low score
        sigA = c(3, 2.5, rep(0, 8), 2.4, 2.2),
        # query-down genes are DOWN here -> concordant, high score
        sigB = c(3, 2.5, rep(0, 8), -2.4, -2.2)
    )
    rownames(mat) <- genes

    se <- SummarizedExperiment(
        assays = list(chrdir = mat),
        colData = DataFrame(sig_id = c("sigA", "sigB"))
    )
    S4Vectors::metadata(se) <-
        list(ref_metric = "chrdir", default_metric = "chrdir")

    res <- calcEnrichment(
        se,
        upRegulated = c("g1", "g2"),
        downRegulated = c("g11", "g12"),
        background = genes,
        method = "GSEA",
        BPPARAM = SerialParam(progressbar = FALSE),
        chunkSize = 1
    )

    # primary score column is GSEAScore, and it is not GSEAup
    expect_equal(unname(regulatorTable(res)$score), unname(assay(res)[, "GSEAScore"]))
    expect_false(
        isTRUE(all.equal(assay(res)[, "GSEAScore"], assay(res)[, "GSEAup"]))
    )

    # the fully concordant signature (sigB) ranks first
    expect_equal(topRegulators(res, n = 1)$sig_id, "sigB")
})


test_that("calcEnrichment derives signatures from expressionMatrix", {
    se <- make_test_se()
    expr <- make_expression_matrix()

    result <- calcEnrichment(
        se,
        expressionMatrix = expr,
        method = "cosine",
        topSig = 2,
        BPPARAM = SerialParam(progressbar = FALSE),
        chunkSize = 1
    )

    expect_s4_class(result, "PerturbMatch")
    expect_equal(upRegulated(result), c("gene1", "gene2"))
    expect_equal(downRegulated(result), c("gene6", "gene5"))
})


test_that("calcEnrichment returns expected columns for cosineExtreme", {
    se <- make_test_se()
    expr <- make_expression_matrix()

    result <- calcEnrichment(
        se,
        expressionMatrix = expr,
        method = "cosineExtreme",
        BPPARAM = SerialParam(progressbar = FALSE),
        chunkSize = 1
    )

    expect_s4_class(result, "PerturbMatch")
    expect_true("cosineExtreme" %in% colnames(assay(result)))
})


# refMetric resolution

test_that("calcEnrichment uses default_metric when refMetric is not supplied", {
    se <- make_limma_se()

    result <- calcEnrichment(
        se,
        upRegulated = c("gene1", "gene2"),
        downRegulated = c("gene5", "gene6"),
        background = paste0("gene", 1:6),
        method = "SignedJaccard",
        topN = 2,
        BPPARAM = SerialParam(progressbar = FALSE),
        chunkSize = 1
    )

    expect_equal(metadata(result)$ref_metric, "pi")
})


test_that("calcEnrichment uses the explicitly requested refMetric", {
    se <- make_limma_se()

    result <- calcEnrichment(
        se,
        upRegulated = c("gene1", "gene2"),
        downRegulated = c("gene5", "gene6"),
        background = paste0("gene", 1:6),
        method = "SignedJaccard",
        refMetric = "log2fc",
        topN = 2,
        BPPARAM = SerialParam(progressbar = FALSE),
        chunkSize = 1
    )

    expect_equal(metadata(result)$ref_metric, "log2fc")
})


test_that("calcEnrichment errors when refMetric is not one of the available metrics", {
    se <- make_limma_se()

    expect_error(
        calcEnrichment(
            se,
            upRegulated = c("gene1", "gene2"),
            downRegulated = c("gene5", "gene6"),
            background = paste0("gene", 1:6),
            method = "SignedJaccard",
            refMetric = "not_a_metric",
            BPPARAM = SerialParam(progressbar = FALSE),
            chunkSize = 1
        ),
        "'refMetric' must be one of"
    )
})


test_that("calcEnrichment errors when the reference object has no ref_metric set", {
    se <- make_test_se()
    metadata(se) <- list()

    expect_error(
        calcEnrichment(
            se,
            upRegulated = c("gene1", "gene2"),
            downRegulated = c("gene5", "gene6"),
            background = paste0("gene", 1:6),
            method = "SignedJaccard",
            BPPARAM = SerialParam(progressbar = FALSE),
            chunkSize = 1
        ),
        "ref_metric.*is not set"
    )
})


test_that("calcEnrichment errors when default_metric is not set", {
    se <- make_test_se()
    metadata(se) <- list(ref_metric = "chrdir") # default_metric missing

    expect_error(
        calcEnrichment(
            se,
            upRegulated = c("gene1", "gene2"),
            downRegulated = c("gene5", "gene6"),
            background = paste0("gene", 1:6),
            method = "SignedJaccard",
            BPPARAM = SerialParam(progressbar = FALSE),
            chunkSize = 1
        ),
        "default_metric.*is not set"
    )
})


test_that("calcEnrichment records whether the reference was the demo database", {
    se <- make_test_se()

    score <- function(reference) {
        calcEnrichment(
            reference,
            upRegulated = c("gene1", "gene2"),
            downRegulated = c("gene5", "gene6"),
            background = paste0("gene", 1:6),
            method = "SignedJaccard",
            topN = 2,
            BPPARAM = SerialParam(progressbar = FALSE),
            chunkSize = 1
        )
    }

    expect_false(metadata(score(se))$demo)

    S4Vectors::metadata(se)$demo <- TRUE
    expect_true(metadata(score(se))$demo)
})


test_that(".runScoringMethod rejects an unknown method", {
    # match.arg() guards the public entry point, so this defensive branch is
    # only reachable by calling the block-level dispatcher directly.
    block <- matrix(1, nrow = 1, ncol = 1, dimnames = list("g1", "s1"))

    expect_error(
        perturbMatch:::.runScoringMethod(
            block, "NotAMethod", "g1", "g1", NULL, NULL, 1
        ),
        "Unknown method"
    )
})
