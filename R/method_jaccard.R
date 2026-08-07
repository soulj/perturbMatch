#' Calculate signed Jaccard similarity between a reference matrix and query
#' gene sets
#'
#' @param reference A (Delayed)Matrix of expression values
#' @param rankObject Column-rank matrix of \code{reference}, ranked ascending
#'   so rank 1 is the smallest value
#' @param upRegulated Character vector of upregulated query genes
#' @param downRegulated Character vector of downregulated query genes
#' @param threshold Optional numeric threshold. If \code{NULL} \code{topN}
#'   is used instead
#' @param topN Number of top/bottom genes to define reference gene sets
#'
#' @return A list with two numeric vectors: signed Jaccard scores and overlap
#'   counts
#' @noRd

.getSignedJaccards <- function(
    reference,
    rankObject,
    upRegulated,
    downRegulated,
    threshold = NULL,
    topN = 250
) {
    if (is.null(threshold)) {
        nGenes <- DelayedArray::colSums(!is.na(rankObject))
        # Keep the top-topN and bottom-topN sets separate: once 2 * topN
        # exceeds the measured genes, the middle ranks fall in both.
        if (any(2 * topN > nGenes)) {
            stop(
                "'topN' (",
                topN,
                ") is too large for the number of ",
                "reference genes: the top and bottom gene sets would ",
                "overlap. Use topN <= genes / 2, or set a 'threshold'.",
                call. = FALSE
            )
        }
        referenceDown <- rankObject <= topN
        referenceUp <- .thresholdByColumn(rankObject, nGenes - topN)
    } else {
        referenceUp <- reference >= threshold
        referenceDown <- reference <= -threshold
    }

    dimnames(referenceUp) <- dimnames(reference)
    dimnames(referenceDown) <- dimnames(reference)

    upup <- .jaccardDelayedMatrix(referenceUp, upRegulated)
    downdown <- .jaccardDelayedMatrix(referenceDown, downRegulated)
    updown <- .jaccardDelayedMatrix(referenceUp, downRegulated)
    downup <- .jaccardDelayedMatrix(referenceDown, upRegulated)

    signedJaccard <- (upup[[1]] + downdown[[1]] - updown[[1]] - downup[[1]]) / 2
    overlap <- upup[[2]] + downdown[[2]] + updown[[2]] + downup[[2]]

    list(signedJaccard, overlap)
}


#' Test each matrix element against a per-column threshold, \code{x[i, j] >
#' cut[j]}
#' @noRd

.thresholdByColumn <- function(x, cut) {
    if (is.matrix(x)) {
        x > cut[col(x)]
    } else {
        t(t(x) > cut)
    }
}


#' Compute Jaccard similarity between a binary matrix and a query gene set
#'
#' @param reference A logical (Delayed)Matrix where \code{TRUE} indicates
#'   membership in a gene set
#' @param query Character vector of query genes
#'
#' @return A list: (1) Jaccard scores per column, (2) intersection sizes
#' @noRd

.jaccardDelayedMatrix <- function(reference, query) {
    overlap <- DelayedArray::colSums(
        reference[query, , drop = FALSE],
        na.rm = TRUE
    )
    setSize <- DelayedArray::colSums(reference, na.rm = TRUE)
    union <- setSize + length(query) - overlap
    list(overlap / union, overlap)
}
