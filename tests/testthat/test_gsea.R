test_that(".gseaSingleColumnBoth returns 0 for both sets when no query genes overlap", {
    stats <- c(1, 2, 3, 4, 5)
    inUp <- rep(FALSE, 5) # none overlap
    inDown <- rep(FALSE, 5) # none overlap

    result <- .gseaSingleColumnBoth(stats, inUp, inDown)

    expect_equal(unname(result["up"]), 0)
    expect_equal(unname(result["down"]), 0)
})

test_that(".gseaSingleColumnBoth returns a named up/down numeric vector", {
    stats <- c(3, 1, 4, 1, 5)
    inUp <- c(TRUE, FALSE, TRUE, FALSE, FALSE)
    inDown <- c(FALSE, TRUE, FALSE, TRUE, FALSE)

    result <- .gseaSingleColumnBoth(stats, inUp, inDown)

    expect_length(result, 2)
    expect_named(result, c("up", "down"))
    expect_type(result, "double")
})

test_that(".gseaSingleColumnBoth silently drops NA values before scoring", {
    # position 2 (geneB) is NA and gets dropped; remaining genes are in the
    # same order in both cases
    stats_with_na <- c(3, NA, 4, 1, 5)
    stats_without_na <- c(3, 4, 1, 5)
    inUp_with_na <- c(TRUE, FALSE, TRUE, FALSE, FALSE) # geneA, geneC
    inUp_without_na <- c(TRUE, TRUE, FALSE, FALSE)
    inDown_with_na <- rep(FALSE, 5)
    inDown_without_na <- rep(FALSE, 4)

    result_na <- .gseaSingleColumnBoth(stats_with_na, inUp_with_na, inDown_with_na)
    result_no_na <- .gseaSingleColumnBoth(stats_without_na, inUp_without_na, inDown_without_na)

    expect_equal(result_na, result_no_na)
})

test_that(".gseaSingleColumnBoth up-score is higher when up-genes have higher expression", {
    # Function sorts stats descending internally before scoring.
    # Favourable: the up-genes have the highest expression -> top of ranked list
    # Unfavourable: the up-genes have the lowest expression -> bottom of ranked list
    stats_favourable <- c(4, 3, 2, 1, 0)
    stats_unfavourable <- c(0, 1, 2, 3, 4)
    inUp <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
    inDown <- rep(FALSE, 5)

    score_high <- .gseaSingleColumnBoth(stats_favourable, inUp, inDown)["up"]
    score_low <- .gseaSingleColumnBoth(stats_unfavourable, inUp, inDown)["up"]

    expect_gt(score_high, score_low)
})

test_that(".gseaSingleColumnBoth down-score is more negative when down-genes cluster at the bottom of the ranking", {
    # calcGseaStat's running-sum ES is negative when a gene set is enriched at
    # the low end of a descending-sorted ranking, so a good match for
    # down-regulated genes (low expression) produces a negative "down" score.
    stats_favourable <- c(0, 1, 2, 3, 4) # down-genes at the bottom of the ranking
    stats_unfavourable <- c(4, 3, 2, 1, 0) # down-genes at the top of the ranking
    inUp <- rep(FALSE, 5)
    inDown <- c(TRUE, TRUE, FALSE, FALSE, FALSE)

    score_matched <- .gseaSingleColumnBoth(stats_favourable, inUp, inDown)["down"]
    score_mismatched <- .gseaSingleColumnBoth(stats_unfavourable, inUp, inDown)["down"]

    expect_lt(score_matched, score_mismatched)
})

# .getGSEAScores

test_that(".getGSEAScores returns a named list of two numeric vectors", {
    set.seed(1)
    expr_mat <- matrix(
        rnorm(20), nrow = 5, ncol = 4,
        dimnames = list(
            c("geneA", "geneB", "geneC", "geneD", "geneE"),
            c("s1", "s2", "s3", "s4")
        )
    )
    up_genes <- c("geneA", "geneB")
    down_genes <- c("geneD", "geneE")

    result <- .getGSEAScores(expr_mat, up_genes, down_genes)

    expect_type(result, "list")
    expect_named(result, c("up", "down"))
    expect_length(result$up, ncol(expr_mat))
    expect_length(result$down, ncol(expr_mat))
    expect_type(result$up, "double")
    expect_type(result$down, "double")
})

test_that(".getGSEAScores up and down scores differ when gene sets are distinct", {
    set.seed(7)
    expr_mat <- matrix(
        rnorm(30), nrow = 6, ncol = 5,
        dimnames = list(
            c("g1", "g2", "g3", "g4", "g5", "g6"),
            paste0("s", 1:5)
        )
    )
    up_genes <- c("g1", "g2")
    down_genes <- c("g5", "g6")

    result <- .getGSEAScores(expr_mat, up_genes, down_genes)

    expect_false(identical(result$up, result$down))
})

test_that(".getGSEAScores returns zeros for both sets when no genes overlap", {
    expr_mat <- matrix(
        rnorm(20), nrow = 4, ncol = 5,
        dimnames = list(
            c("geneA", "geneB", "geneC", "geneD"),
            paste0("s", 1:5)
        )
    )
    up_genes <- c("geneX", "geneY")
    down_genes <- c("geneP", "geneQ")

    result <- .getGSEAScores(expr_mat, up_genes, down_genes)

    expect_equal(result$up, c(s1 = 0, s2 = 0, s3 = 0, s4 = 0, s5 = 0))
    expect_equal(result$down, c(s1 = 0, s2 = 0, s3 = 0, s4 = 0, s5 = 0))
})
