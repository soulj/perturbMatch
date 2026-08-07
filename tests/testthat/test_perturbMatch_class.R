library(testthat)
library(SummarizedExperiment)

# Mock data helper

make_object <- function() {
    se <- SummarizedExperiment(
        assays = list(
            enrichment = matrix(
                c(0.9, 0.2, -0.5),
                ncol = 1,
                dimnames = list(c("sig1", "sig2", "sig3"), "enrichment")
            )
        ),
        rowData = DataFrame(
            pert_id    = c("p1", "p2", "p3"),
            pert_iname = c("gene1_KO", "gene2_KO", "gene3_KO")
        )
    )

    S4Vectors::metadata(se) <- list(
        up_regulated   = c("gene1", "gene2"),
        down_regulated = c("gene5", "gene6"),
        background     = paste0("gene", 1:10),
        threshold      = 0.5,
        topN           = 250L,
        topSig         = 100L
    )

    new("PerturbMatch", se)
}


# Class tests

test_that("PerturbMatch object can be created", {
    obj <- make_object()

    expect_s4_class(obj, "PerturbMatch")
    expect_s4_class(obj, "SummarizedExperiment")
})

test_that("an object without an enrichment assay is rejected", {
    se <- SummarizedExperiment(
        metadata = list(up_regulated = c("gene1", "gene2"))
    )

    expect_error(new("PerturbMatch", se), "enrichment")
})

test_that("gene sets must be character vectors", {
    se <- SummarizedExperiment(
        assays = list(
            enrichment = matrix(0.5, dimnames = list("sig1", "score"))
        ),
        metadata = list(up_regulated = 1:3)
    )

    expect_error(new("PerturbMatch", se), "up_regulated")
})

test_that("score_column must name a column of the enrichment assay", {
    se <- SummarizedExperiment(
        assays = list(
            enrichment = matrix(0.5, dimnames = list("sig1", "score"))
        ),
        metadata = list(score_column = "SignedJaccard")
    )

    expect_error(new("PerturbMatch", se), "score_column")
})

test_that("scoring parameters must be single numbers", {
    se <- SummarizedExperiment(
        assays = list(
            enrichment = matrix(0.5, dimnames = list("sig1", "score"))
        ),
        metadata = list(topN = c(100L, 250L))
    )

    expect_error(new("PerturbMatch", se), "topN")
})

test_that("provenance accessors return values set at construction", {
    obj <- make_object()

    expect_equal(upRegulated(obj), c("gene1", "gene2"))
    expect_equal(downRegulated(obj), c("gene5", "gene6"))
    expect_equal(background(obj), paste0("gene", 1:10))
})

test_that("scoring parameters are readable from metadata()", {
    # threshold/topN/topSig have no accessor of their own: they are the
    # arguments the caller passed, echoed back through metadata() and show()
    obj <- make_object()
    meta <- S4Vectors::metadata(obj)

    expect_equal(meta$threshold, 0.5)
    expect_equal(meta$topN, 250L)
    expect_equal(meta$topSig, 100L)
})

test_that("provenance accessors have no replacement methods", {
    expect_false(existsMethod("upRegulated<-", "PerturbMatch"))
    expect_false(existsMethod("downRegulated<-", "PerturbMatch"))
    expect_false(existsMethod("background<-", "PerturbMatch"))
    expect_false(existsMethod("threshold<-", "PerturbMatch"))
    expect_false(existsMethod("topN<-", "PerturbMatch"))
    expect_false(existsMethod("topSig<-", "PerturbMatch"))
})

test_that("enrichment accessor returns enrichment assay", {
    obj <- make_object()

    enrich <- enrichment(obj)

    expect_true(is.matrix(enrich))
    expect_equal(nrow(enrich), 3)
    expect_equal(colnames(enrich), "enrichment")
})

test_that("result combines rowData and enrichment scores", {
    obj <- make_object()

    res <- regulatorTable(obj)

    expect_s3_class(res, "data.frame")
    expect_true(all(c("pert_id", "pert_iname", "score") %in% colnames(res)))
    expect_equal(res$score, c(0.9, 0.2, -0.5))
})

test_that("topRegulators returns highest scoring entries ranked by absolute score", {
    obj <- make_object()

    tt <- topRegulators(obj, n = 2)

    expect_equal(nrow(tt), 2)
    expect_true(abs(tt$score[1]) >= abs(tt$score[2]))
    # sig1 (0.9) and sig3 (-0.5) have the largest magnitude scores, ahead of
    # sig2 (0.2), even though sig3's score is negative
    expect_equal(tt$pert_id, c("p1", "p3"))
})

test_that("topRegulators handles n larger than available rows", {
    obj <- make_object()

    tt <- topRegulators(obj, n = 10)

    expect_equal(nrow(tt), 3)
})

test_that("show method prints gene sets and parameters", {
    obj <- make_object()

    expect_output(show(obj), "up_regulated")
    expect_output(show(obj), "down_regulated")
    expect_output(show(obj), "background")
    expect_output(show(obj), "parameters")
    expect_output(show(obj), "threshold")
})

test_that("show method flags results built on the demo database", {
    obj <- make_object()

    expect_false(any(grepl("DEMONSTRATION", capture.output(show(obj)))))

    S4Vectors::metadata(obj)$demo <- TRUE
    expect_output(show(obj), "DEMONSTRATION")
})

test_that("provenance accessors return NULL when not set at construction", {
    obj <- perturbMatch:::.newPerturbMatch(up = c("GENE1", "GENE2"))

    expect_null(downRegulated(obj))
    expect_null(background(obj))

    meta <- S4Vectors::metadata(obj)
    expect_null(meta$threshold)
    expect_null(meta$topN)
    expect_null(meta$topSig)
})

test_that("show method prints NULL for unset gene sets", {
    obj <- perturbMatch:::.newPerturbMatch(up = c("GENE1", "GENE2"))

    expect_output(show(obj), "down_regulated: NULL")
    expect_output(show(obj), "background: NULL")
})
