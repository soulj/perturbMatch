#' @title Calculate enrichment
#' @description Calculates the enrichment of a regulator expression signature
#'   in the query expression signature. This measure is used to prioritise
#'   which upstream regulators are driving the observed gene expression changes.
#'
#' @section Input signature:
#'   Two input formats are supported:
#'   \itemize{
#'     \item \strong{Ranked list}: supply \code{expressionMatrix} (a data frame
#'       of gene IDs and score values). \code{upRegulated},
#'       \code{downRegulated}, and \code{background} are ignored.
#'     \item \strong{Gene sets}: supply \code{upRegulated},
#'       \code{downRegulated}, and \code{background} as character vectors.
#'       \code{expressionMatrix} should be \code{NULL}.
#'   }
#'   If both are provided, \code{expressionMatrix} takes precedence.
#'
#' @section Reference metric:
#'   Two metadata slots record which metrics a reference supports:
#'   \itemize{
#'     \item \code{metadata(object)$ref_metric}: the available metrics, each
#'       naming a stored assay. limma references hold \code{c("pi", "log2fc")}
#'       and chrdir references hold \code{"chrdir"}.
#'     \item \code{metadata(object)$default_metric}: the metric used when
#'       \code{refMetric = NULL}.
#'   }
#'
#' @param object A \code{SummarizedExperiment} reference database with
#'   \code{metadata(object)$ref_metric} set to one of \code{"chrdir"},
#'   \code{"log2fc"}, or \code{"pi"}.
#' @param expressionMatrix Optional data frame with gene IDs in column 1 and
#'   score values in column 2. Required for \code{"cosine"} and
#'   \code{"cosineExtreme"} methods. If supplied, \code{upRegulated},
#'   \code{downRegulated}, and \code{background} are not used.
#' @param upRegulated Character vector of upregulated genes in the experiment.
#'   Only used when \code{expressionMatrix} is \code{NULL}.
#' @param downRegulated Character vector of downregulated genes in the
#'   experiment. Only used when \code{expressionMatrix} is \code{NULL}.
#' @param background Character vector of all genes expressed/measured in the
#'   experiment. Only used when \code{expressionMatrix} is \code{NULL}.
#' @param method Enrichment method to use. One of \code{"SignedJaccard"},
#'   \code{"GSEA"}, \code{"cosine"}, \code{"cosineExtreme"}.
#'   Defaults to \code{"SignedJaccard"}.
#' @param refMetric Reference metric to score against, one of
#'   \code{metadata(object)$ref_metric}. When \code{NULL} (default),
#'   \code{metadata(object)$default_metric} is used.
#' @param threshold Optional numeric threshold applied to reference metric
#'   values to define reference gene sets (\code{|metric| > threshold}).
#'   When \code{NULL}, \code{topN} is used instead. Default: \code{NULL}.
#' @param topN Number of top and bottom genes per reference signature to use,
#'   ranked by the reference metric. Used only when \code{threshold} is
#'   \code{NULL}. Default: \code{250}.
#' @param topSig Number of top/bottom genes from the query signature to use
#'   when deriving gene sets from \code{expressionMatrix}. Default: \code{250}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object
#'   controlling parallelisation. Defaults to
#'   \code{BiocParallel::SerialParam(progressbar = TRUE)}.
#' @param chunkSize Number of reference signatures (columns) per processing
#'   block. When \code{NULL} (default) a value giving roughly one block per
#'   worker in \code{BPPARAM} is chosen. Larger values use more memory, and
#'   \code{chunkSize >= ncol(object)} leaves a single block that cannot be
#'   parallelised.
#'
#' @return A \link[=PerturbMatch-class]{PerturbMatch} object containing the
#'   enrichment scores.
#'
#' @examples
#' # Load the small bundled demo reference database
#' se <- getReferenceDatabase(demo = TRUE)
#'
#' # Score a query gene set against every reference signature. Here we just
#' # use a handful of reference genes to stand in for a real up/down
#' # signature.
#' result <- calcEnrichment(
#'     se,
#'     upRegulated   = rownames(se)[seq_len(3)],
#'     downRegulated = rownames(se)[seq(4, 5)],
#'     background    = rownames(se),
#'     method        = "SignedJaccard",
#'     BPPARAM       = BiocParallel::SerialParam(progressbar = FALSE)
#' )
#'
#' cols <- c("pert_symbol", "pert_type", "score", "zscore")
#' topRegulators(result, n = 3)[, cols]
#'
#' @rdname calcEnrichment
#' @aliases calcEnrichment,SummarizedExperiment-method
#' @export
#' @exportMethod calcEnrichment
#' @importFrom fgsea calcGseaStat
#' @importFrom stats na.omit reorder setNames
#' @importFrom utils head tail
setGeneric(
    "calcEnrichment",
    function(
        object,
        expressionMatrix = NULL,
        upRegulated = NULL,
        downRegulated = NULL,
        background = NULL,
        method = c("SignedJaccard", "GSEA", "cosine", "cosineExtreme"),
        refMetric = NULL,
        threshold = NULL,
        topN = 250L,
        topSig = 250L,
        BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
        chunkSize = NULL
    ) {
        standardGeneric("calcEnrichment")
    }
)

setMethod(
    "calcEnrichment",
    "SummarizedExperiment",
    function(
        object,
        expressionMatrix = NULL,
        upRegulated = NULL,
        downRegulated = NULL,
        background = NULL,
        method = c("SignedJaccard", "GSEA", "cosine", "cosineExtreme"),
        refMetric = NULL,
        threshold = NULL,
        topN = 250L,
        topSig = 250L,
        BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
        chunkSize = NULL
    ) {
        .calcEnrichmentInternal(
            object = object,
            expressionMatrix = expressionMatrix,
            upRegulated = upRegulated,
            downRegulated = downRegulated,
            background = background,
            method = method,
            refMetric = refMetric,
            threshold = threshold,
            topN = topN,
            topSig = topSig,
            BPPARAM = BPPARAM,
            chunkSize = chunkSize
        )
    }
)


# Internal helpers: reference scoring assay
#' Resolve the refMetric to use for scoring.
#'
#' If refMetric is NULL, falls back to metadata(object)$default_metric.
#' Validates the resolved value against metadata(object)$ref_metric, the
#' vector of metrics the reference object supports.
#' @noRd
.resolveRefMetric <- function(object, refMetric) {
    available <- metadata(object)$ref_metric
    default <- metadata(object)$default_metric

    if (is.null(available)) {
        stop(
            "metadata(object)$ref_metric is not set. ",
            "Please use getReferenceDatabase() to load a valid ",
            "reference object.",
            call. = FALSE
        )
    }

    if (is.null(default)) {
        stop(
            "metadata(object)$default_metric is not set. ",
            "Please use getReferenceDatabase() to load a valid ",
            "reference object.",
            call. = FALSE
        )
    }

    resolved <- if (is.null(refMetric)) default else refMetric

    if (!resolved %in% available) {
        stop(
            "'refMetric' must be one of: ",
            paste(available, collapse = ", "),
            ". Got: '",
            resolved,
            "'.",
            call. = FALSE
        )
    }

    resolved
}


# Internal helpers: query gene-set preparation

.prepareGeneSets <- function(object, expressionMatrix, topSig) {
    expressionMatrix <- expressionMatrix[
        expressionMatrix[, 1] %in% rownames(object),
        ,
        drop = FALSE
    ]
    list(
        background = expressionMatrix[, 1],
        upRegulated = expressionMatrix[
            order(-expressionMatrix[, 2])[seq_len(topSig)],
            1
        ],
        downRegulated = expressionMatrix[
            order(expressionMatrix[, 2])[seq_len(topSig)],
            1
        ]
    )
}


# Internal helper: blockApply scoring

.scoreChunks <- function(
    scoringAssay,
    method,
    upRegulated,
    downRegulated,
    background,
    expressionMatrix,
    threshold,
    topN,
    chunkSize,
    BPPARAM
) {
    # Restrict to background genes before chunking, so the grid reflects the
    # rows actually scored.
    keepIdx <- rownames(scoringAssay) %in% background
    scoringAssay <- scoringAssay[keepIdx, , drop = FALSE]

    grid <- DelayedArray::colAutoGrid(
        scoringAssay,
        ncol = min(ncol(scoringAssay), chunkSize)
    )

    results <- DelayedArray::blockApply(
        scoringAssay,
        .runScoringMethod,
        method,
        upRegulated,
        downRegulated,
        expressionMatrix,
        threshold,
        topN,
        grid = grid,
        BPPARAM = BPPARAM
    )

    dplyr::bind_rows(results)
}


# Query gene-set resolution

#' Resolve the up/down/background query gene sets and restrict them to genes
#' present in the reference object.
#'
#' If \code{expressionMatrix} is supplied, the gene sets are derived from it
#' (overriding any \code{upRegulated}/\code{downRegulated}/\code{background}
#' arguments). Otherwise the supplied gene sets are used directly.
#' @noRd
.resolveQueryGeneSets <- function(
    object,
    expressionMatrix,
    upRegulated,
    downRegulated,
    background,
    topSig
) {
    if (!is.null(expressionMatrix)) {
        sets <- .prepareGeneSets(object, expressionMatrix, topSig)
        background <- sets$background
        upRegulated <- sets$upRegulated
        downRegulated <- sets$downRegulated
    }

    .validateGeneSets(upRegulated, downRegulated, background)

    resolved <- list(
        upRegulated = na.omit(unique(upRegulated[
            upRegulated %in% rownames(object)
        ])),
        downRegulated = na.omit(unique(downRegulated[
            downRegulated %in% rownames(object)
        ])),
        background = background
    )

    .validateGeneSetOverlap(resolved, object)

    resolved
}


# Core internal function

.calcEnrichmentInternal <- function(
    object,
    expressionMatrix = NULL,
    upRegulated = NULL,
    downRegulated = NULL,
    background = NULL,
    method = c("SignedJaccard", "GSEA", "cosine", "cosineExtreme"),
    refMetric = NULL,
    threshold = NULL,
    topN = 250L,
    topSig = 250L,
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
    chunkSize = NULL
) {
    method <- match.arg(method)
    refMetric <- .resolveRefMetric(object, refMetric)
    .validateMethodInputs(method, expressionMatrix)
    if (is.null(chunkSize)) {
        chunkSize <- .defaultChunkSize(object, BPPARAM)
    }
    .validateScoringParams(topN, topSig, chunkSize, threshold)

    # Retrieve scoring assay once, outside blockApply.
    scoringAssay <- assay(object, refMetric)

    sets <- .resolveQueryGeneSets(
        object, expressionMatrix, upRegulated, downRegulated, background,
        topSig
    )

    results <- .scoreChunks(
        scoringAssay, method, sets$upRegulated, sets$downRegulated,
        sets$background, expressionMatrix, threshold, topN, chunkSize,
        BPPARAM
    )

    .buildPerturbMatch(
        results, object, sets, threshold, topN, topSig, refMetric
    )
}


# Validation helpers

.validateMethodInputs <- function(method, expressionMatrix) {
    matrixRequired <- c("cosine", "cosineExtreme")
    if (method %in% matrixRequired && is.null(expressionMatrix)) {
        stop(
            "An expressionMatrix with score values must be provided ",
            "for the '",
            method,
            "' method.",
            call. = FALSE
        )
    }
}

.validateGeneSets <- function(upRegulated, downRegulated, background) {
    if (!is(upRegulated, "character")) {
        stop("upRegulated must be a character vector.", call. = FALSE)
    }
    if (!is(downRegulated, "character")) {
        stop("downRegulated must be a character vector.", call. = FALSE)
    }
    if (!is(background, "character")) {
        stop("background must be a character vector.", call. = FALSE)
    }
}

# Fail early if the query gene sets don't overlap the reference genes
.validateGeneSetOverlap <- function(sets, object) {
    .checkOverlap <- function(genes, name) {
        if (length(genes) == 0) {
            stop(
                "None of the '",
                name,
                "' genes were found in the reference ",
                "object (rownames(object)). Check that gene ID types match ",
                "between your query and the reference (e.g. Entrez vs Ensembl ",
                "vs gene symbol).",
                call. = FALSE
            )
        }
    }

    inReference <- sets$background %in% rownames(object)
    backgroundInReference <- sets$background[inReference]

    .checkOverlap(sets$upRegulated, "upRegulated")
    .checkOverlap(sets$downRegulated, "downRegulated")
    .checkOverlap(backgroundInReference, "background")
}

#' Choose a default column-chunk size for block processing.
#'
#' Targets roughly one block per worker
#' @noRd
.defaultChunkSize <- function(object, BPPARAM) {
    nworkers <- BiocParallel::bpnworkers(BPPARAM)
    ncols <- ncol(object)
    max(500L, as.integer(ceiling(ncols / nworkers)))
}

.validateScoringParams <- function(topN, topSig, chunkSize, threshold) {
    .validatePositiveInteger(topN, "topN")
    .validatePositiveInteger(topSig, "topSig")
    .validatePositiveInteger(chunkSize, "chunkSize")
    .validateThreshold(threshold)
}

.validatePositiveInteger <- function(value, name) {
    ok <- is.numeric(value) &&
        length(value) == 1 &&
        !is.na(value) &&
        value == as.integer(value) &&
        value >= 1
    if (!ok) {
        stop(
            "'",
            name,
            "' must be a single positive whole number.",
            call. = FALSE
        )
    }
}

.validateThreshold <- function(threshold) {
    if (is.null(threshold)) {
        return(invisible(NULL))
    }
    ok <- is.numeric(threshold) &&
        length(threshold) == 1 &&
        !is.na(threshold) &&
        threshold >= 0
    if (!ok) {
        stop(
            "'threshold' must be NULL or a single non-negative number.",
            call. = FALSE
        )
    }
}


# Scoring method dispatcher (block-level)

.scoreCosine <- function(block, expressionMatrix, upRegulated, downRegulated) {
    scores <- .getCosineScores(
        block,
        expressionMatrix,
        upRegulated,
        downRegulated,
        topOnly = FALSE
    )
    data.frame(cosine = scores)
}

.scoreCosineExtreme <- function(
    block,
    expressionMatrix,
    upRegulated,
    downRegulated
) {
    scores <- .getCosineScores(
        block,
        expressionMatrix,
        upRegulated,
        downRegulated,
        topOnly = TRUE
    )
    data.frame(cosineExtreme = scores)
}

.scoreSignedJaccard <- function(
    block,
    rankBlock,
    upRegulated,
    downRegulated,
    threshold,
    topN
) {
    res <- .getSignedJaccards(
        block,
        rankBlock,
        upRegulated,
        downRegulated,
        threshold,
        topN
    )
    data.frame(SignedJaccard = res[[1]], overlap = res[[2]])
}

.scoreGSEA <- function(block, upRegulated, downRegulated) {
    res <- .getGSEAScores(block, upRegulated, downRegulated)
    data.frame(
        GSEAScore = res[[1]] - res[[2]],
        GSEAup = res[[1]],
        GSEAdown = res[[2]]
    )
}

.runScoringMethod <- function(
    block,
    method,
    upRegulated,
    downRegulated,
    expressionMatrix,
    threshold,
    topN
) {
    rankBlock <- NULL
    if (method == "SignedJaccard") {
        rankBlock <- DelayedMatrixStats::colRanks(
            block,
            ties.method = "min",
            preserveShape = TRUE
        )
        dimnames(rankBlock) <- dimnames(block)
    }

    switch(
        method,
        cosine = .scoreCosine(
            block,
            expressionMatrix,
            upRegulated,
            downRegulated
        ),
        cosineExtreme = .scoreCosineExtreme(
            block,
            expressionMatrix,
            upRegulated,
            downRegulated
        ),
        SignedJaccard = .scoreSignedJaccard(
            block,
            rankBlock,
            upRegulated,
            downRegulated,
            threshold,
            topN
        ),
        GSEA = .scoreGSEA(block, upRegulated, downRegulated),
        stop("Unknown method: '", method, "'", call. = FALSE)
    )
}


# Result builder

.buildPerturbMatch <- function(
    results,
    object,
    sets,
    threshold,
    topN,
    topSig,
    refMetric
) {
    upRegulated <- sets$upRegulated
    downRegulated <- sets$downRegulated
    background <- sets$background

    # Each scoring method puts its primary score first and any secondary
    # columns after it.
    scoreColumn <- colnames(results)[1L]
    results$zscore <- as.vector(scale(results[[scoreColumn]]))

    scoringResults <- SummarizedExperiment::SummarizedExperiment(
        assays = list(enrichment = as.matrix(results))
    )

    rownames(scoringResults) <- SummarizedExperiment::colData(object)[, 1]

    SummarizedExperiment::rowData(scoringResults) <-
        SummarizedExperiment::colData(object)

    scoringResults <- new("PerturbMatch", scoringResults)

    metadata(scoringResults)$up_regulated <- upRegulated
    metadata(scoringResults)$down_regulated <- downRegulated
    metadata(scoringResults)$background <- background

    # Recording both would leave the unused one looking like it shaped the
    # scores. Assigning NULL drops the entry.
    metadata(scoringResults)$threshold <- threshold
    metadata(scoringResults)$topN <- if (is.null(threshold)) topN else NULL

    metadata(scoringResults)$topSig <- topSig
    metadata(scoringResults)$ref_metric <- refMetric
    metadata(scoringResults)$score_column <- scoreColumn
    metadata(scoringResults)$demo <- isTRUE(metadata(object)$demo)
    metadata(scoringResults)$source <- metadata(object)$source

    # metadata<- bypasses validity, so check once the object is fully assembled.
    validObject(scoringResults)

    scoringResults
}
