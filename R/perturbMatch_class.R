#' PerturbMatch enrichment results
#'
#' @description
#' An S4 class extending
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} that holds the
#' enrichment scores produced by \code{\link{calcEnrichment}}. Scores are
#' available through \code{enrichment()}, or as a tidy table through
#' \code{\link{regulatorTable}()}. The usual \code{SummarizedExperiment}
#' accessors also work.
#'
#' Details of analysis are stored in \code{metadata()} when
#' the object is built. The query gene sets have accessors of their own
#' (\code{\link{upRegulated}}, \code{\link{downRegulated}},
#' \code{\link{background}}). The scoring parameters are read straight from
#' \code{metadata()} and are printed by \code{show()}. The reference gene sets
#' come from either \code{threshold} or \code{topN}.
#'
#' The class contains an \code{"enrichment"}
#' assay of numeric scores, gene sets held as character vectors, and a
#' \code{score_column} entry. Objects
#' failing any of these are rejected by \code{\link[methods]{validObject}}.
#'
#' @param object A \code{PerturbMatch} object returned by
#'   \code{calcEnrichment}.
#' @return
#' \itemize{
#'   \item \code{enrichment}: matrix of enrichment scores.
#'   \item \code{regulatorTable}: tidy \code{data.frame} of scores plus row
#'     metadata.
#' }
#'
#' @name PerturbMatch-class
#' @exportClass PerturbMatch
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @import SummarizedExperiment
#' @importFrom methods callNextMethod is new validObject
#' @importFrom S4Vectors metadata
setClass(
    "PerturbMatch",
    contains = "SummarizedExperiment"
)


.isSingleString <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x)
}


.validGeneSets <- function(meta) {
    fields <- c("up_regulated", "down_regulated", "background")
    invalid <- vapply(
        fields,
        function(field) {
            genes <- meta[[field]]
            !is.null(genes) && (!is.character(genes) || anyNA(genes))
        },
        logical(1)
    )

    sprintf(
        "metadata()$%s must be a character vector with no NAs",
        fields[invalid]
    )
}


# regulatorTable() renames this column to "score". A name that matches
# nothing would otherwise surface as a missing-column error far from here.
.validScoreColumn <- function(meta, enrichmentAssay) {
    scoreColumn <- meta$score_column
    if (is.null(scoreColumn)) {
        return(character())
    }
    if (!.isSingleString(scoreColumn)) {
        return("metadata()$score_column must be a single string")
    }
    if (
        is.null(enrichmentAssay) ||
            scoreColumn %in% colnames(enrichmentAssay)
    ) {
        return(character())
    }

    paste0(
        "metadata()$score_column '",
        scoreColumn,
        "' names no column of the 'enrichment' assay"
    )
}


.validScalarParams <- function(meta) {
    fields <- c("threshold", "topN", "topSig")
    invalid <- vapply(
        fields,
        function(field) {
            value <- meta[[field]]
            !is.null(value) &&
                (!is.numeric(value) || length(value) != 1L || anyNA(value))
        },
        logical(1)
    )

    sprintf("metadata()$%s must be a single number", fields[invalid])
}


.validPerturbMatch <- function(object) {
    errors <- character()

    enrichmentAssay <- NULL
    if (!"enrichment" %in% assayNames(object)) {
        errors <- c(errors, "an assay named 'enrichment' is required")
    } else {
        enrichmentAssay <- assay(object, "enrichment")
        if (!is.numeric(enrichmentAssay)) {
            errors <- c(errors, "the 'enrichment' assay must be numeric")
        }
    }

    meta <- S4Vectors::metadata(object)
    errors <- c(
        errors,
        .validGeneSets(meta),
        .validScoreColumn(meta, enrichmentAssay),
        .validScalarParams(meta)
    )

    if (length(errors) > 0) errors else TRUE
}

setValidity("PerturbMatch", .validPerturbMatch)


# Internal constructor used by the tests.
.newPerturbMatch <- function(
    up = character(),
    down = character(),
    bg = character(),
    threshold = NULL,
    topN = NULL,
    topSig = NULL
) {
    params <- Filter(
        Negate(is.null),
        list(
            up_regulated = if (length(up) > 0) up else NULL,
            down_regulated = if (length(down) > 0) down else NULL,
            background = if (length(bg) > 0) bg else NULL,
            threshold = threshold,
            topN = topN,
            topSig = topSig
        )
    )

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(enrichment = matrix(numeric(0), nrow = 0, ncol = 0))
    )
    S4Vectors::metadata(se) <- params

    new("PerturbMatch", se)
}


#' Show a PerturbMatch object
#'
#' Prints the inherited \code{SummarizedExperiment} summary followed by the
#' gene sets and parameters recorded when the object was built.
#'
#' @param object A \code{PerturbMatch} object.
#' @return Invisibly returns \code{NULL}. Called for its printing side effect.
#' @examples
#' pm <- methods::new(
#'     "PerturbMatch",
#'     SummarizedExperiment::SummarizedExperiment(
#'         assays = list(
#'             enrichment = matrix(0.5, dimnames = list("SIG1", "score"))
#'         ),
#'         metadata = list(
#'             up_regulated   = c("GENE1", "GENE2"),
#'             down_regulated = "GENE3",
#'             background     = c("GENE1", "GENE2", "GENE3"),
#'             threshold      = 0.5,
#'             topN           = 250L,
#'             topSig         = 100L
#'         )
#'     )
#' )
#' show(pm)
#' @rdname show-PerturbMatch
#' @export
setMethod("show", "PerturbMatch", function(object) {
    callNextMethod()

    meta <- S4Vectors::metadata(object)

    if (isTRUE(meta$demo)) {
        cat(
            "NOTE: scored against the bundled DEMONSTRATION database",
            "(21 signatures), not a full reference\n"
        )
    }

    .printGenes <- function(label, genes) {
        n <- length(genes)
        if (n > 0) {
            preview <- if (n > 5) {
                paste0(paste(head(genes, 5), collapse = ", "), " ...")
            } else {
                paste(genes, collapse = ", ")
            }
            cat(sprintf("%s(%d): %s\n", label, n, preview))
        } else {
            cat(sprintf("%s: NULL\n", label))
        }
    }

    .printGenes("up_regulated", meta$up_regulated)
    .printGenes("down_regulated", meta$down_regulated)
    .printGenes("background", meta$background)

    params <- meta[c("threshold", "topN", "topSig")]
    params <- Filter(Negate(is.null), params)
    if (length(params) > 0) {
        cat(
            "parameters:",
            paste(
                paste0(names(params), " = ", unlist(params)),
                collapse = ", "
            ),
            "\n"
        )
    }
})


#' Up-regulated genes from a PerturbMatch object
#'
#' Returns the up-regulated gene set recorded in \code{metadata()} when the
#' object was built by \code{calcEnrichment}.
#'
#' @param object A \code{PerturbMatch} object.
#' @return A character vector of gene identifiers, or \code{NULL} if unset.
#' @examples
#' pm <- methods::new(
#'     "PerturbMatch",
#'     SummarizedExperiment::SummarizedExperiment(
#'         assays = list(
#'             enrichment = matrix(0.5, dimnames = list("SIG1", "score"))
#'         ),
#'         metadata = list(up_regulated = c("GENE1", "GENE2"),
#'             background = c("GENE1", "GENE2"))))
#' upRegulated(pm)
#' @rdname upRegulated
#' @aliases upRegulated,PerturbMatch-method
#' @export
#' @exportMethod upRegulated
setGeneric("upRegulated", function(object) standardGeneric("upRegulated"))

setMethod("upRegulated", "PerturbMatch", function(object) {
    S4Vectors::metadata(object)$up_regulated
})


#' Down-regulated genes from a PerturbMatch object
#'
#' Returns the down-regulated gene set recorded in \code{metadata()} when the
#' object was built by \code{calcEnrichment}.
#'
#' @param object A \code{PerturbMatch} object.
#' @return A character vector of gene identifiers, or \code{NULL} if unset.
#' @examples
#' pm <- methods::new(
#'     "PerturbMatch",
#'     SummarizedExperiment::SummarizedExperiment(
#'         assays = list(
#'             enrichment = matrix(0.5, dimnames = list("SIG1", "score"))
#'         ),
#'         metadata = list(down_regulated = c("GENE3", "GENE4"),
#'             background = c("GENE3", "GENE4"))))
#' downRegulated(pm)
#' @rdname downRegulated
#' @aliases downRegulated,PerturbMatch-method
#' @export
#' @exportMethod downRegulated
setGeneric("downRegulated", function(object) standardGeneric("downRegulated"))

setMethod("downRegulated", "PerturbMatch", function(object) {
    S4Vectors::metadata(object)$down_regulated
})


#' Background genes from a PerturbMatch object
#'
#' Returns the background gene set recorded in \code{metadata()} when the
#' object was built by \code{calcEnrichment}. This is the universe of genes
#' against which enrichment was scored.
#'
#' @param object A \code{PerturbMatch} object.
#' @return A character vector of gene identifiers, or \code{NULL} if unset.
#' @examples
#' pm <- methods::new(
#'     "PerturbMatch",
#'     SummarizedExperiment::SummarizedExperiment(
#'         assays = list(
#'             enrichment = matrix(0.5, dimnames = list("SIG1", "score"))
#'         ),
#'         metadata = list(background = c("GENE1", "GENE2", "GENE3"))))
#' background(pm)
#' @rdname background
#' @aliases background,PerturbMatch-method
#' @export
#' @exportMethod background
setGeneric("background", function(object) standardGeneric("background"))

setMethod("background", "PerturbMatch", function(object) {
    S4Vectors::metadata(object)$background
})


#' Enrichment scores from a PerturbMatch object
#'
#' Returns the matrix stored in the \code{"enrichment"} assay.
#'
#' @param object A \code{PerturbMatch} object.
#' @return A matrix of enrichment scores (rows = features, columns =
#'   signatures).
#' @examples
#' mat <- matrix(runif(3), nrow = 3, dimnames = list(
#'     c("GENE1", "GENE2", "GENE3"), "score"
#' ))
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(enrichment = mat)
#' )
#' S4Vectors::metadata(se) <- list(
#'     up_regulated   = "GENE1",
#'     down_regulated = "GENE3",
#'     background     = c("GENE1", "GENE2", "GENE3")
#' )
#' pm <- new("PerturbMatch", se)
#' enrichment(pm)
#' @rdname enrichment
#' @aliases enrichment,PerturbMatch-method
#' @export
#' @exportMethod enrichment
setGeneric("enrichment", function(object) standardGeneric("enrichment"))

setMethod("enrichment", "PerturbMatch", function(object) {
    assay(object, "enrichment")
})


#' Tidy table of regulator scores from a PerturbMatch object
#'
#' Combines \code{rowData} with the enrichment scores into a single
#' \code{data.frame}, one row per reference signature.
#'
#' The name deliberately avoids \code{result}/\code{results}, which would sit
#' one letter from \code{\link[DESeq2]{results}} in a script that uses both.
#'
#' @param object A \code{PerturbMatch} object.
#' @return A \code{data.frame} with one row per feature holding every
#'   \code{rowData} column plus a \code{score} column.
#' @examples
#' mat <- matrix(runif(3), nrow = 3, dimnames = list(
#'     c("GENE1", "GENE2", "GENE3"), "score"
#' ))
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(enrichment = mat)
#' )
#' S4Vectors::metadata(se) <- list(
#'     up_regulated   = "GENE1",
#'     down_regulated = "GENE3",
#'     background     = c("GENE1", "GENE2", "GENE3")
#' )
#' pm <- new("PerturbMatch", se)
#' regulatorTable(pm)
#' @rdname regulatorTable
#' @aliases regulatorTable,PerturbMatch-method
#' @export
#' @exportMethod regulatorTable
setGeneric("regulatorTable", function(object) standardGeneric("regulatorTable"))

setMethod("regulatorTable", "PerturbMatch", function(object) {
    df <- as.data.frame(rowData(object))
    scoringResults <- enrichment(object)

    # calcEnrichment() records which column holds the primary score, under
    # whatever name the scoring method gave it. Objects assembled by hand carry
    # no such record, so their first column is taken instead.
    scoreColumn <- S4Vectors::metadata(object)$score_column
    if (is.null(scoreColumn)) {
        scoreColumn <- colnames(scoringResults)[1L]
    }
    colnames(scoringResults)[colnames(scoringResults) == scoreColumn] <- "score"

    cbind(df, scoringResults)
})


#' Top-scoring regulator signatures
#'
#' Returns the \code{n} highest-scoring rows of \code{\link{regulatorTable}},
#' ranked by
#' absolute enrichment score.
#'
#' @param object A \code{PerturbMatch} object.
#' @param ...    Additional arguments passed to methods.
#' @param n      Number of signatures to return, capped at \code{nrow(object)}.
#'   Default \code{5}.
#' @return A \code{data.frame} with at most \code{n} rows, ordered by descending
#'   absolute \code{score}.
#' @examples
#' mat <- matrix(runif(5), nrow = 5, dimnames = list(
#'     paste0("GENE", seq_len(5)), "score"
#' ))
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(enrichment = mat)
#' )
#' S4Vectors::metadata(se) <- list(
#'     up_regulated   = "GENE1",
#'     down_regulated = "GENE5",
#'     background     = paste0("GENE", seq_len(5))
#' )
#' pm <- new("PerturbMatch", se)
#' topRegulators(pm, n = 3)
#' @rdname topRegulators
#' @aliases topRegulators,PerturbMatch-method
#' @export
#' @exportMethod topRegulators
setGeneric(
    "topRegulators",
    function(object, n = 5, ...) standardGeneric("topRegulators")
)

setMethod("topRegulators", "PerturbMatch", function(object, n = 5, ...) {
    tbl <- regulatorTable(object)
    n <- min(n, nrow(tbl))
    tbl <- tbl[order(abs(tbl[, "score"]), decreasing = TRUE), , drop = FALSE]
    tbl[seq_len(n), , drop = FALSE]
})
