library(testthat)
library(DelayedArray)

# Mock Data
# Reference matrix: 3 genes (A, B, C) and 2 gene expression signatures (S1, S2)
# S1: [1, 2, 3]
# S2: [-1, -2, -3] (Perfectly anti-correlated to S1)
ref_mat <- matrix(c(1, 2, 3, -1, -2, -3),
    nrow = 3,
    dimnames = list(c("A", "B", "C"), c("S1", "S2")))
reference <- DelayedArray::DelayedArray(ref_mat)

# Query signature: Gene B is up, Gene C is down
# Gene D isn't in the reference to test filtering
query_df <- data.frame(
    gene = c("A", "B", "C", "D"),
    logFC = c(0, 5, -5, 100) # D should be ignored
)

up_genes <- "B"
down_genes <- "C"

test_that("perturbMatch:::.getCosineScores handles standard cosine (topOnly = FALSE)", {
    # Logic: All genes in reference (A, B, C) are used.
    # S1 vector: [1, 2, 3]
    # Query vector aligned: [0, 5, -5]
    # Dot product: (1*0) + (2*5) + (3*-5) = 10 - 15 = -5

    results <- perturbMatch:::.getCosineScores(
        reference = reference,
        expressionMatrix = query_df,
        upRegulated = up_genes,
        downRegulated = down_genes,
        topOnly = FALSE
    )

    expect_type(results, "double")
    expect_length(results, 2)
    expect_named(results, c("S1", "S2"))

    # S1 and S2 should be exact opposites
    expect_equal(results[["S1"]], -results[["S2"]])

    # Manual calculation check for S1:
    # norm_ref = sqrt(1^2 + 2^2 + 3^2) = sqrt(14)
    # norm_query = sqrt(0^2 + 5^2 + (-5)^2) = sqrt(50)
    # cosine = -5 / (sqrt(14) * sqrt(50))
    expected_s1 <- -5 / (sqrt(14) * sqrt(50))
    expect_equal(results[["S1"]], expected_s1)
})

test_that("perturbMatch:::.getCosineScores handles extreme cosine (topOnly = TRUE)", {
    # Logic: Only genes B and C are used. Gene A is dropped.
    # S1 vector: [2, 3]
    # Query vector: [5, -5]
    # Dot product: (2*5) + (3*-5) = 10 - 15 = -5

    results <- perturbMatch:::.getCosineScores(
        reference = reference,
        expressionMatrix = query_df,
        upRegulated = up_genes,
        downRegulated = down_genes,
        topOnly = TRUE
    )

    # Manual calculation check for S1:
    # norm_ref = sqrt(2^2 + 3^2) = sqrt(13)
    # norm_query = sqrt(5^2 + (-5)^2) = sqrt(50)
    expected_s1_extreme <- -5 / (sqrt(13) * sqrt(50))

    expect_equal(results[["S1"]], expected_s1_extreme)
})

test_that("perturbMatch:::.getCosineScores handles missing genes gracefully", {
    # Query with no overlapping genes
    query_empty <- data.frame(gene = c("X", "Y"), logFC = c(1, 2))

    # Depending on how .cosine handles empty colSums, this usually returns NaN
    # but we want to ensure it doesn't crash.
    results <- perturbMatch:::.getCosineScores(reference, query_empty, "X", "Y", topOnly = FALSE)

    expect_true(all(is.na(results)))
})

test_that(".alignExpressionMatrix correctly orders and filters", {
    # Purpose: Ensure query genes are in the SAME order as reference rows
    # even if the input data frame is scrambled.
    scrambled_query <- query_df[c(3, 1, 2, 4), ] # C, A, B, D

    aligned <- perturbMatch:::.alignExpressionMatrix(
        scrambled_query, reference, up_genes, down_genes, topOnly = FALSE
    )

    # Should match reference rownames: A, B, C
    expect_equal(aligned[, 1], rownames(reference))
})
