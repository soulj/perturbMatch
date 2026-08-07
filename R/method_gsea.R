#' Score the query gene sets against every reference signature by GSEA
#'
#' @param rankObject A (Delayed)Matrix of reference values, one column per
#'   signature, used directly as the ranking statistic
#' @param upRegulated Character vector of upregulated query genes
#' @param downRegulated Character vector of downregulated query genes
#'
#' @return A list of two numeric vectors, the up and down enrichment statistics
#'   of each column
#' @noRd

.getGSEAScores <- function(rankObject, upRegulated, downRegulated) {
    # Resolve membership once and reuse it across every column.
    genes <- rownames(rankObject)
    inUp <- genes %in% upRegulated
    inDown <- genes %in% downRegulated

    scores <- apply(rankObject, 2, .gseaSingleColumnBoth, inUp, inDown)

    list(up = scores["up", ], down = scores["down", ])
}

#' Enrichment statistics of both query gene sets in one reference signature
#'
#' @param stats Numeric vector of reference values for a single signature
#' @param inUp,inDown Logical vectors marking which elements of \code{stats}
#'   belong to the up and down query gene sets
#'
#' @return A named numeric vector of length two, \code{up} and \code{down}
#' @noRd

.gseaSingleColumnBoth <- function(stats, inUp, inDown) {
    keep <- !is.na(stats)
    stats <- stats[keep]
    inUp <- inUp[keep]
    inDown <- inDown[keep]

    ord <- order(stats, decreasing = TRUE)
    statsSorted <- stats[ord]

    selUp <- which(inUp[ord])
    selDown <- which(inDown[ord])

    c(
        up = if (length(selUp) == 0) {
            0
        } else {
            fgsea::calcGseaStat(statsSorted, selUp)
        },
        down = if (length(selDown) == 0) {
            0
        } else {
            fgsea::calcGseaStat(statsSorted, selDown)
        }
    )
}
