#' Calculate cosine similarity between reference gene expression profiles
#' and a query signature
#'
#' @param reference A (Delayed)Matrix of expression values
#' @param expressionMatrix Data frame with gene IDs in column 1 and
#'   fold-change values in column 2
#' @param upRegulated Character vector of upregulated query genes
#' @param downRegulated Character vector of downregulated query genes
#' @param topOnly Logical. If \code{TRUE} restrict to up/down genes only
#'   (\code{"cosineExtreme"} variant)
#'
#' @return Numeric vector of cosine similarity scores, one per column of
#'   \code{reference}
#' @noRd

.getCosineScores <- function(
    reference,
    expressionMatrix,
    upRegulated,
    downRegulated,
    topOnly = FALSE
    ) {
    expressionMatrix <- .alignExpressionMatrix(
        expressionMatrix, reference, upRegulated, downRegulated, topOnly
    )

    if (topOnly) {
        keep <- rownames(reference) %in% expressionMatrix[, 1]
        reference <- reference[keep, , drop = FALSE]
    }

    unlist(.cosine(reference, expressionMatrix[, 2]))
}


#' Vectorised cosine similarity between matrix columns and a single vector
#'
#' Each reference column is compared to \code{vec} over the genes measured in
#' both. Missing values (\code{NA}) are excluded from the dot product and from
#' both norms alike, so a signature's norm never counts genes that cannot
#' contribute to the numerator.
#'
#' @param mat Numeric (Delayed)Matrix of reference values.
#' @param vec Numeric vector aligned to the rows of \code{mat}. \code{NA} marks
#'   genes the query did not score.
#'
#' @return Numeric vector of cosine similarities, one per column of \code{mat}.
#' @noRd

.cosine <- function(mat, vec) {
    # Drop genes the query never scored so they don't inflate the norms.
    present <- !is.na(vec)
    mat <- mat[present, , drop = FALSE]
    vec <- vec[present]

    numerator <- MatrixGenerics::colSums(mat * vec, na.rm = TRUE)
    # vecNorm2 zeroes vec^2 wherever mat is NA, keeping each column's two norms
    # and its dot product on the same set of measured genes.
    matNorm2 <- MatrixGenerics::colSums(mat^2, na.rm = TRUE)
    vecNorm2 <- MatrixGenerics::colSums((!is.na(mat)) * vec^2, na.rm = TRUE)

    numerator / sqrt(matNorm2 * vecNorm2)
}


#' Align and optionally subset an expression matrix to reference row names
#'
#' @param expressionMatrix Data frame (gene IDs column 1, values column 2)
#' @param reference A matrix whose rownames define the gene universe
#' @param upRegulated Character vector of upregulated genes
#' @param downRegulated Character vector of downregulated genes
#' @param topOnly Logical. If \code{TRUE} subset to extreme genes
#'
#' @return Filtered and row-aligned expression matrix
#' @noRd

.alignExpressionMatrix <- function(
    expressionMatrix,
    reference,
    upRegulated,
    downRegulated,
    topOnly
    ) {
    expressionMatrix <- expressionMatrix[
        expressionMatrix[, 1] %in% rownames(reference),
    ]
    expressionMatrix <- expressionMatrix[
        match(rownames(reference), expressionMatrix[, 1]),
    ]

    if (topOnly) {
        expressionMatrix <- expressionMatrix[
            expressionMatrix[, 1] %in% c(upRegulated, downRegulated),
        ]
    }

    expressionMatrix
}
