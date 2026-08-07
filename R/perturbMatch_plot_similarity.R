#' Plot signature similarity scores
#'
#' Draws the top positive and top negative signatures as a lollipop plot.
#'
#' @param object A \code{SummarizedExperiment} containing scores in its colData.
#' @param topN Integer. Number of top positive and top negative hits to show.
#' @param inferActivity Logical. If \code{TRUE}, flip the sign of the score for
#'   KD/KO perturbations so it reflects predicted gene activity rather than
#'   transcriptional similarity. If \code{FALSE} (default), the score keeps its
#'   similarity meaning and the perturbation direction is shown by point shape:
#'   up-pointing triangles for gain-of-function (OE), down-pointing triangles
#'   for loss-of-function (KD/KO).
#' @return A \code{ggplot} object.
#' @examples
#' # A PerturbMatch object as returned by calcEnrichment(): one row per
#' # reference signature, scores in the 'enrichment' assay and the signature
#' # annotation in rowData().
#' enr <- cbind(
#'     score  = c(0.9, -0.7, -0.4, 0.5, -0.2),
#'     zscore = c(1.6, -1.1, -0.6, 0.8, -0.3)
#' )
#' rownames(enr) <- paste0("SIG", seq_len(5))
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays  = list(enrichment = enr),
#'     rowData = S4Vectors::DataFrame(
#'         pert_symbol = c("GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E"),
#'         pert_type   = c("OE", "KD", "KO", "OE", "KD")
#'     )
#' )
#' pm <- methods::new("PerturbMatch", se)
#' plotSimilarity(pm, topN = 3)
#' plotSimilarity(pm, topN = 3, inferActivity = TRUE)
#' @aliases plotSimilarity,PerturbMatch-method
#' @import ggplot2
#' @export
#' @exportMethod plotSimilarity
setGeneric(
    "plotSimilarity",
    function(object, topN = 10, inferActivity = FALSE) {
        standardGeneric("plotSimilarity")
    }
)

setMethod(
    "plotSimilarity",
    "PerturbMatch",
    function(object, topN = 10, inferActivity = FALSE) {
        plotData <- .preparePlotData(
            regulatorTable(object), topN, inferActivity
        )
        yLabel <- if (inferActivity) {
            "Inferred Activity Z Score"
        } else {
            "Similarity Z Score"
        }

        p <- ggplot(
            plotData,
            aes(x = reorder(.data$pert_symbol, .data$zscore), y = .data$zscore)
        ) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
            geom_segment(
                aes(xend = .data$pert_symbol, yend = 0),
                color = "grey70",
                linewidth = 0.8
            )

        # With activity inference the direction is already in the score, so use
        # a single circle. Otherwise map it to the point shape.
        p <- if (inferActivity) {
            p +
                geom_point(
                    aes(fill = .data$Direction),
                    shape = 21,
                    size = 4.5,
                    color = "black",
                    stroke = 0.6
                )
        } else {
            p +
                geom_point(
                    aes(fill = .data$Direction, shape = .data$Perturbation),
                    size = 4.5,
                    color = "black",
                    stroke = 0.6
                )
        }

        .lollipopScales(p, yLabel, inferActivity)
    }
)

# Flip scores for KD/KO so score reflects inferred gene activity
.inferActivityScore <- function(data) {
    pert <- toupper(data$pert_type)
    isKnockdown <- pert %in% c("KD", "KO")
    data$zscore <- ifelse(isKnockdown, data$zscore * -1L, data$zscore)
    data
}

# Sort by zscore, restrict to the top/bottom topN hits, and label direction
.preparePlotData <- function(data, topN, inferActivity) {
    if (inferActivity) {
        data <- .inferActivityScore(data)
    }

    data <- data[order(data$zscore, decreasing = TRUE), ]
    if (nrow(data) > 2 * topN) {
        data <- rbind(head(data, topN), tail(data, topN))
    }

    data$Direction <- ifelse(
        data$zscore >= 0,
        "Mimic (Positive)",
        "Reverse (Negative)"
    )

    # KD/KO are loss-of-function, everything else is gain-of-function.
    pert <- toupper(data$pert_type)
    data$Perturbation <- factor(
        ifelse(
            pert %in% c("KD", "KO"),
            "Loss-of-function (KD/KO)",
            "Gain-of-function (OE)"
        ),
        levels = c("Gain-of-function (OE)", "Loss-of-function (KD/KO)")
    )
    data
}

# Add scales, guides, labels, and theme to a lollipop base plot
.lollipopScales <- function(
    p,
    yLabel = "Similarity Score",
    inferActivity = FALSE
) {
    okabeIto <- c(
        "Mimic (Positive)" = "#E69F00",
        "Reverse (Negative)" = "#0072B2"
    )

    p <- p + scale_fill_manual(values = okabeIto)

    # Filled triangles for the perturbation shape (24 = up = OE,
    # 25 = down = KD/KO). The fill legend keeps a circle to show only sign.
    if (inferActivity) {
        p <- p + guides(fill = guide_legend(override.aes = list(shape = 21)))
    } else {
        pertShapes <- c(
            "Gain-of-function (OE)" = 24,
            "Loss-of-function (KD/KO)" = 25
        )
        p <- p +
            scale_shape_manual(values = pertShapes, drop = FALSE) +
            guides(
                fill = guide_legend(override.aes = list(shape = 21)),
                shape = guide_legend(override.aes = list(fill = "grey70"))
            ) +
            labs(shape = "Perturbation")
    }

    p +
        coord_flip() +
        labs(
            title = "Top positive and negative signatures",
            x = "Gene Perturbation",
            y = yLabel,
            fill = "Direction Profile"
        ) +
        theme_minimal(base_size = 14) +
        theme(
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(face = "bold"),
            legend.position = "bottom",
            legend.box = "vertical",
            legend.box.just = "left"
        )
}
