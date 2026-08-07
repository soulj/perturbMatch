#' Volcano Plot from a SummarizedExperiment Object
#'
#' Extracts effect size and p-value assays for a given signature from a
#' \code{SummarizedExperiment} and delegates to \code{\link{volcanoPlot}}.
#'
#' @param object A \code{SummarizedExperiment} object. Must contain assays
#'   named by \code{effectAssay} and \code{pvalueAssay}, and \code{rownames}
#'   corresponding to gene identifiers.
#' @param signature A single character string matching one of
#'   \code{colnames(object)}, selecting the column to plot.
#' @param effectAssay Character string giving the name of the assay that holds
#'   log2 fold-change values. Default: \code{"logFC"}.
#' @param pvalueAssay Character string giving the name of the assay that holds
#'   adjusted p-values (FDR). Default: \code{"padj"}.
#' @param numberPoints Integer. Number of top up- and down-regulated genes to
#'   label automatically (ranked by p-value, then |logFC|). Set to \code{NULL}
#'   to suppress automatic labelling. Default: \code{5}.
#' @param selectedPoints Character vector of gene IDs to label regardless of
#'   significance. Combined with \code{numberPoints} labels when both are
#'   supplied. Default: \code{NULL}.
#' @param logFCThreshold Numeric. Absolute log2 fold-change cut-off used to
#'   classify genes as regulated. Default: \code{log2(1.5)}.
#' @param fdrThreshold Numeric. Adjusted p-value cut-off. Default: \code{0.05}.
#' @param pointSize Numeric. Size of plotted points passed to
#'   \code{\link[ggplot2]{geom_point}}. Default: \code{2}.
#' @param labSize Numeric. Text size for gene labels passed to
#'   \code{\link[ggrepel]{geom_text_repel}}. Default: \code{4}.
#' @param truncateInfPvalue Logical. If \code{TRUE}, genes with an adjusted
#'   p-value of exactly zero (an infinite \eqn{-\log_{10}} FDR) are drawn at
#'   the largest finite \eqn{-\log_{10}} FDR in the signature, so they appear
#'   as a flat row along the top of the plot. If \code{FALSE}, these genes are
#'   left at their true infinite value and dropped by the plotting layer.
#'   Default: \code{TRUE}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{volcanoPlot}}
#'
#' @examples
#' set.seed(1)
#' n <- 30
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(
#'         logFC = matrix(rnorm(n * 2), nrow = n,
#'             dimnames = list(paste0("Gene", seq_len(n)), c("A", "B"))),
#'         padj = matrix(runif(n * 2), nrow = n,
#'             dimnames = list(paste0("Gene", seq_len(n)), c("A", "B")))
#'     )
#' )
#'
#' plotVolcano(se, signature = "A")
#' plotVolcano(se, signature = "A", selectedPoints = c("Gene1", "Gene2"))
#'
#' @export
plotVolcano <- function(
    object,
    signature,
    effectAssay = "logFC",
    pvalueAssay = "padj",
    numberPoints = 5,
    selectedPoints = NULL,
    logFCThreshold = log2(1.5),
    fdrThreshold = 0.05,
    pointSize = 2,
    labSize = 4,
    truncateInfPvalue = TRUE
    ) {

    idx <- .validateVolcanoSE(object, signature, effectAssay, pvalueAssay)

    data <- data.frame(
        gene = rownames(object),
        logFC = assay(object, effectAssay)[, idx],
        padj = assay(object, pvalueAssay)[, idx],
        stringsAsFactors = FALSE
    )

    volcanoPlot(
        data,
        numberPoints   = numberPoints,
        selectedPoints = selectedPoints,
        logFCThreshold = logFCThreshold,
        fdrThreshold   = fdrThreshold,
        pointSize      = pointSize,
        labSize        = labSize,
        truncateInfPvalue = truncateInfPvalue
    )
}

# Validate that a SummarizedExperiment has the assays and signature
# plotVolcano needs, returning the signature's column index.
.validateVolcanoSE <- function(object, signature, effectAssay, pvalueAssay) {
    required <- c(effectAssay, pvalueAssay)
    missingAssays <- setdiff(required, assayNames(object))

    if (length(missingAssays)) {
        stop(
            "Volcano plotting requires assays: ",
            paste(required, collapse = ", "),
            call. = FALSE
        )
    }

    if (!signature %in% colnames(object)) {
        stop("Signature not found in object.", call. = FALSE)
    }

    match(signature, colnames(object))
}


#' Volcano Plot of Differential Expression Results
#'
#' Renders a volcano plot of log2 fold-change against \eqn{-\log_{10}}(FDR),
#' colouring genes as up-regulated, down-regulated, or unchanged relative to
#' supplied thresholds. The most significant genes at each tail are labelled
#' automatically, and arbitrary genes can be highlighted via
#' \code{selectedPoints}.
#'
#' Infinite \eqn{-\log_{10}} p-values (arising from \eqn{p = 0}) are, by
#' default, truncated to the largest finite value to keep the plot axis
#' bounded. See \code{truncateInfPvalue} to switch this off.
#'
#' @param data A data frame with at least three columns, in this order:
#'   \enumerate{
#'     \item Gene / feature identifier (character).
#'     \item Log2 fold-change (numeric).
#'     \item Adjusted p-value / FDR (numeric, in the interval \eqn{[0, 1]}).
#'   }
#'   Additional columns are ignored.
#' @param numberPoints Integer. Number of top up- and down-regulated genes to
#'   label (ranked first by p-value, then by \code{|logFC|}). Pass \code{NULL}
#'   to disable automatic labelling. Default: \code{5}.
#' @param selectedPoints Character vector of gene identifiers to label
#'   regardless of their statistics. Merged with automatic labels when
#'   \code{numberPoints} is also set. Default: \code{NULL}.
#' @param logFCThreshold Numeric. Minimum absolute log2 fold-change required
#'   for a gene to be called regulated. Default: \code{log2(1.5)} (\eqn{\approx
#'   0.585}).
#' @param fdrThreshold Numeric. Maximum adjusted p-value for a gene to be
#'   called regulated. Shown as a horizontal dashed line. Default: \code{0.05}.
#' @param pointSize Numeric. Passed to \code{\link[ggplot2]{geom_point}} as
#'   \code{size}. Default: \code{2}.
#' @param labSize Numeric. Passed to \code{\link[ggrepel]{geom_text_repel}} as
#'   \code{size}. Default: \code{4}.
#' @param truncateInfPvalue Logical. If \code{TRUE}, genes with \eqn{p = 0}
#'   (an infinite \eqn{-\log_{10}} FDR) are placed at the largest finite
#'   \eqn{-\log_{10}} FDR, forming a flat row along the top of the plot. If
#'   \code{FALSE}, they keep their infinite value and are dropped by the
#'   plotting layer. Default: \code{TRUE}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object. Colour scale:
#'   \describe{
#'     \item{firebrick3}{Up-regulated (\code{logFC >= logFCThreshold} and
#'       \code{padj <= fdrThreshold}).}
#'     \item{dodgerblue3}{Down-regulated (\code{logFC <= -logFCThreshold} and
#'       \code{padj <= fdrThreshold}).}
#'     \item{black}{Unchanged (all other genes).}
#'   }
#'
#' @importFrom dplyr %>% mutate filter arrange bind_rows desc
#' @importFrom rlang .data sym
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual geom_hline
#'   geom_vline theme xlab ylab
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot theme_cowplot
#'
#' @examples
#' set.seed(42)
#' n <- 500
#' df <- data.frame(
#'     gene  = paste0("Gene", seq_len(n)),
#'     logFC = rnorm(n, sd = 2),
#'     padj  = rbeta(n, 0.5, 5)
#' )
#' volcanoPlot(df)
#'
#' # Highlight specific genes
#' volcanoPlot(df, selectedPoints = c("Gene1", "Gene42"))
#'
#' # Stricter thresholds, no automatic labels
#' volcanoPlot(df, numberPoints = NULL,
#'     logFCThreshold = log2(2), fdrThreshold = 0.01)
#'
#' @export
volcanoPlot <- function(
    data,
    numberPoints = 5,
    selectedPoints = NULL,
    logFCThreshold = log2(1.5),
    fdrThreshold = 0.05,
    pointSize = 2,
    labSize = 4,
    truncateInfPvalue = TRUE
    ) {

    data <- as.data.frame(data)
    cols <- .validateVolcanoColumns(data)

    data <- .classifyVolcanoExpression(
        data, cols$fcCol, cols$pvalCol, logFCThreshold, fdrThreshold,
        truncateInfPvalue
    )

    p <- .volcanoBasePlot(
        data, cols$fcCol, cols$pvalCol, cols$idCol,
        logFCThreshold, fdrThreshold, pointSize
    )

    topGenes <- .volcanoTopGenes(
        data, cols$fcCol, cols$pvalCol, cols$idCol, numberPoints, selectedPoints
    )

    .addVolcanoLabels(p, topGenes, cols$idCol, labSize)
}

# Validate and identify the ID/fold-change/FDR columns of volcano input data
.validateVolcanoColumns <- function(data) {
    if (ncol(data) < 3) {
        stop("The dataset needs ID, foldchange and FDR columns", call. = FALSE)
    }

    idCol <- colnames(data)[1]
    fcCol <- colnames(data)[2]
    pvalCol <- colnames(data)[3]

    if (!all(is.numeric(data[, fcCol]))) {
        stop("The 2nd column should be numeric fold changes", call. = FALSE)
    }
    if (!all(is.numeric(data[, pvalCol]))) {
        stop("The 3rd column should be numeric FDRs", call. = FALSE)
    }

    list(idCol = idCol, fcCol = fcCol, pvalCol = pvalCol)
}

# Add -log10(FDR), drop NA FDRs, and classify each gene's expression
# direction relative to the fold-change/FDR thresholds
.classifyVolcanoExpression <- function(
    data, fcCol, pvalCol, logFCThreshold, fdrThreshold,
    truncateInfPvalue = TRUE
) {
    data <- data %>%
        dplyr::mutate(logpval = -log10(!!sym(pvalCol))) %>%
        dplyr::filter(!is.na(!!sym(pvalCol))) %>%
        dplyr::mutate(
            Expression = dplyr::case_when(
                !!sym(fcCol) >= logFCThreshold &
                    !!sym(pvalCol) <= fdrThreshold ~ "Up-regulated",
                !!sym(fcCol) <= -logFCThreshold &
                    !!sym(pvalCol) <= fdrThreshold ~ "Down-regulated",
                TRUE ~ "Unchanged"
            )
        )

    data$Expression <- factor(
        data$Expression,
        levels = c("Down-regulated", "Unchanged", "Up-regulated")
    )

    # Optionally truncate Inf -log10(p) values (p == 0) to the largest finite
    # value so they remain visible along the top of the plot rather than being
    # dropped for being non-finite.
    if (truncateInfPvalue && any(is.infinite(data$logpval))) {
        maxNonInf <- max(data[is.finite(data$logpval), "logpval"], na.rm = TRUE)
        data[is.infinite(data$logpval), "logpval"] <- maxNonInf
    }

    data
}

# Build the base scatter plot (points, threshold lines, colour scale, axes)
.volcanoBasePlot <- function(
    data, fcCol, pvalCol, idCol, logFCThreshold, fdrThreshold, pointSize
) {
    myColors <- c("dodgerblue3", "black", "firebrick3")
    names(myColors) <- levels(data$Expression)

    ggplot2::ggplot(
        data,
        ggplot2::aes(x = !!sym(fcCol), y = .data$logpval)
    ) +
        ggplot2::geom_point(
            ggplot2::aes(color = .data$Expression), size = pointSize
        ) +
        cowplot::theme_cowplot(font_size = 16) +
        ggplot2::scale_color_manual(name = "Expression", values = myColors) +
        ggplot2::geom_hline(
            yintercept = -log10(fdrThreshold), linetype = "dashed"
        ) +
        ggplot2::geom_vline(
            xintercept = c(-logFCThreshold, logFCThreshold),
            linetype = "dashed"
        ) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::xlab(expression("log"[2] * "FC")) +
        ggplot2::ylab(expression("-log"[10] * "FDR"))
}

# Build the label set: automatic top genes plus any user-selected genes
.volcanoTopGenes <- function(
    data, fcCol, pvalCol, idCol, numberPoints, selectedPoints
) {
    # NULL means "no automatic labels". head() itself only accepts a number
    n <- if (is.null(numberPoints)) 0 else numberPoints

    topGenes <- dplyr::bind_rows(
        data %>%
            dplyr::filter(.data$Expression == "Up-regulated") %>%
            dplyr::arrange(!!sym(pvalCol), dplyr::desc(abs(!!sym(fcCol)))) %>%
            head(n),
        data %>%
            dplyr::filter(.data$Expression == "Down-regulated") %>%
            dplyr::arrange(!!sym(pvalCol), dplyr::desc(abs(!!sym(fcCol)))) %>%
            head(n)
    )

    if (!is.null(selectedPoints)) {
        selectedRows <- data %>% dplyr::filter(!!sym(idCol) %in% selectedPoints)
        topGenes <- dplyr::bind_rows(selectedRows, topGenes)
    }

    topGenes
}

# Add repelled gene-label text to a volcano plot, if there are any to add
.addVolcanoLabels <- function(p, topGenes, idCol, labSize) {
    if (nrow(topGenes) == 0) {
        return(p)
    }

    p +
        ggrepel::geom_text_repel(
            data = topGenes,
            ggplot2::aes(label = !!sym(idCol)),
            size = labSize,
            max.overlaps = Inf,
            min.segment.length = 0,
            seed = 42,
            force = 5
        )
}
