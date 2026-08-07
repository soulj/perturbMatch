library(testthat)
library(DelayedArray)

# Setup mock data.
# rankObject holds ascending ranks (ties.method = "min"): rank 1 is the
# smallest reference value, so the largest ranks mark the reference "up" set.
# Genes A-D carry the largest ranks (reference up) and G-J the smallest
# (reference down).
rankObject <- DelayedArray::DelayedArray(seed = data.frame(rank = 10:1))
rownames(rankObject) <- LETTERS[1:10]

# `reference` carries the raw metric values and is only consulted by the
# threshold branch (values 1..10, so >= 8 selects H, I, J).
reference <- DelayedArray::DelayedArray(seed = data.frame(rank = 1:10))
rownames(reference) <- LETTERS[1:10]
topN <- 4

# Mock query sets
query_up <- LETTERS[1:4]
query_down <- LETTERS[7:10]
background <- LETTERS[1:10]

test_that(".getSignedJaccards returns expected values for perfect match", {
    # query_up matches reference_up and query_down matches reference_down
    # signed = (1 + 1 - 0 - 0) / 2 = 1
    # overlap = 4 + 4 + 0 + 0 = 8

    result <- perturbMatch:::.getSignedJaccards(
        reference = reference,
        rankObject = rankObject,
        upRegulated = query_up,
        downRegulated = query_down,
        threshold = NULL,
        topN = 4
    )

    expect_equal(result[[1]], c(rank = 1))
    expect_equal(result[[2]], c(rank = 8))
})

test_that(".getSignedJaccards returns expected values for perfect opposite", {
    # When query_up matches reference_down and query_down matches reference_up
    # signed = (0 + 0 - 1 - 1) / 2 = -1

    result <- perturbMatch:::.getSignedJaccards(
        reference = reference,
        rankObject = rankObject,
        upRegulated = query_down, # Swapped
        downRegulated = query_up, # Swapped
        threshold = NULL,
        topN = 4
    )

    expect_equal(result[[1]], c(rank = -1))
    expect_equal(result[[2]], c(rank = 8))
})

test_that(".getSignedJaccards handles the threshold parameter correctly", {
    # If threshold is 8:
    # reference_up (>= 8) is H, I, J (3 genes)
    # reference_down (<= -8) is empty (0 genes)

    # If we query with H, I, J as upRegulated
    result <- perturbMatch:::.getSignedJaccards(
        reference = reference,
        rankObject = rankObject,
        upRegulated = c("H", "I", "J"),
        downRegulated = c("A"),
        threshold = 8,
        topN = NULL
    )

    # upup should be 1 (match), others 0
    expect_equal(result[[1]], c(rank = 0.5)) # (1 + 0 - 0 - 0) / 2
    expect_equal(result[[2]], c(rank = 3)) # Only the 3 up genes match
})

test_that(".getSignedJaccards returns zero for non-overlapping sets", {

    result <- perturbMatch:::.getSignedJaccards(
        reference = reference,
        rankObject = rankObject,
        upRegulated = c("E", "F"),
        downRegulated = c("E", "F"),
        threshold = NULL,
        topN = 4
    )

    expect_equal(unname(result[[1]]), 0)
    expect_equal(unname(result[[2]]), 0)
})

test_that(".getSignedJaccards errors when topN is too large for the gene count", {
    # rankObject has 10 genes; topN = 6 makes the top-6 and bottom-6 sets
    # overlap (2 * 6 > 10).
    expect_error(
        perturbMatch:::.getSignedJaccards(
            reference = reference,
            rankObject = rankObject,
            upRegulated = query_up,
            downRegulated = query_down,
            threshold = NULL,
            topN = 6
        ),
        "too large"
    )
})

test_that(".getSignedJaccards thresholds each column by its own gene count", {
    # Regression test: with NAs the number of measured genes varies per
    # signature, so the top-topN cut (nGenes - topN) differs per column and
    # must be broadcast down the columns. A naive `rank > (nGenes - topN)`
    # recycles the length-ncol vector element-wise and mis-selects the up set.
    m <- matrix(
        NA_integer_, nrow = 10, ncol = 2,
        dimnames = list(LETTERS[1:10], c("c1", "c2"))
    )
    m[, 1] <- 1:10       # c1: all 10 genes measured (A smallest ... J largest)
    m[1:6, 2] <- 1:6     # c2: only A-F measured, G-J are NA
    rankObj <- DelayedArray::DelayedArray(m)
    ref <- rankObj

    result <- perturbMatch:::.getSignedJaccards(
        reference = ref,
        rankObject = rankObj,
        upRegulated = c("D", "E", "F"),   # c2 up-set (ranks 4-6); not c1's
        downRegulated = c("A", "B", "C"), # down-set for both columns
        threshold = NULL,
        topN = 3
    )

    # c1 up-set is H,I,J so query up D,E,F only hits the down side -> 0.5;
    # c2 up-set is D,E,F so both up and down match perfectly -> 1.
    expect_equal(result[[1]], c(c1 = 0.5, c2 = 1))
})

test_that(".getSignedJaccards allows topN at exactly half the gene count", {
    # 2 * 5 == 10: the top-5 and bottom-5 sets partition the genes exactly,
    # so this boundary must remain valid.
    expect_no_error(
        perturbMatch:::.getSignedJaccards(
            reference = reference,
            rankObject = rankObject,
            upRegulated = query_up,
            downRegulated = query_down,
            threshold = NULL,
            topN = 5
        )
    )
})
