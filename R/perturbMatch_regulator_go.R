#' GO enrichment of perturbMatch regulator hits by GSEA
#'
#' @description
#' Runs gene set enrichment analysis (\pkg{fgsea}) over Gene Ontology terms
#' using the reference regulators ranked by their perturbMatch similarity
#' score. This asks which biological processes are enriched among the
#' regulators whose perturbation signatures most closely \emph{mimic} the
#' query (positive enrichment) or most strongly \emph{reverse} it (negative
#' enrichment).
#'
#' @details
#' A reference database usually holds several signatures per regulator gene
#' (different cell lines, and knockdown versus overexpression). These are
#' collapsed to a single ranking statistic per gene before scoring:
#' \enumerate{
#'   \item scores are put on a consistent inferred-activity scale by flipping
#'     the sign of the z-score for KD/KO perturbations (as in
#'     \code{plotSimilarity(inferActivity = TRUE)}), controlled by
#'     \code{inferActivity}.
#'   \item the per-signature scores of each gene are aggregated with
#'     \code{aggregate} (the single \code{"extreme"} value with the largest
#'     magnitude, or the \code{"mean"} of the activity z-scores).
#' }
#' The resulting one-value-per-gene ranking is scored against GO gene sets, so
#' the positive end of the ranking (mimicking regulators) and the negative end
#' (reversing regulators) are tested together and separated by the sign of each
#' term's normalised enrichment score (\code{NES}).
#'
#' Gene sets are built from \code{orgDb} (using the ancestor-propagated
#' \code{GOALL} mapping restricted to \code{ontology}), which also supplies
#' the regulator symbol-to-Entrez mapping. Human-readable term names are added
#' from \pkg{GO.db} when it is installed. Otherwise \code{term} falls back to
#' the GO accession.
#'
#' @param object A \link[=PerturbMatch-class]{PerturbMatch} object from
#'   \code{\link{calcEnrichment}}.
#' @param ontology GO ontology to test: one of \code{"BP"} (default),
#'   \code{"MF"} or \code{"CC"}.
#' @param aggregate How to collapse the multiple signatures of a regulator into
#'   one ranking statistic. \code{"extreme"} (default) keeps the single
#'   largest-magnitude z-score. \code{"mean"} averages the activity z-scores.
#' @param inferActivity Logical. When \code{TRUE} (default) the z-score sign is
#'   flipped for KD/KO perturbations before aggregating, so the ranking
#'   reflects inferred gene activity rather than raw transcriptional
#'   similarity.
#' @param orgDb An \pkg{AnnotationDbi} \code{OrgDb} object, or the name of an
#'   installed \code{OrgDb} package, supplying the GO annotation and the
#'   symbol/Entrez mapping. Default \code{"org.Hs.eg.db"}, matching the species
#'   of the packaged reference databases.
#' @param minSize,maxSize Minimum and maximum GO gene-set sizes to test, passed
#'   to \code{\link[fgsea]{fgsea}}.
#' @param ... Further arguments passed to \code{\link[fgsea]{fgsea}}.
#'
#' @return A \code{data.frame}, one row per GO term ordered by adjusted
#'   p-value, with columns \code{ID} (GO accession), \code{term}, \code{NES}
#'   (positive = enriched among mimicking regulators, negative = among
#'   reversing regulators), \code{pval}, \code{padj}, \code{size} and
#'   \code{leadingEdge} (comma-separated regulator gene symbols).
#'
#' @examples
#' se <- getReferenceDatabase(demo = TRUE)
#' res <- calcEnrichment(
#'     se,
#'     upRegulated   = rownames(se)[seq_len(5)],
#'     downRegulated = rownames(se)[seq(6, 10)],
#'     background    = rownames(se),
#'     BPPARAM       = BiocParallel::SerialParam(progressbar = FALSE)
#' )
#' # minSize is lowered because the demo database is tiny. Keep the default
#' # on a full reference.
#' go <- regulatorGSEA(res, minSize = 1)
#' head(go)
#'
#' @seealso \code{\link{plotRegulatorGSEA}}, \code{\link{calcEnrichment}}
#' @rdname regulatorGSEA
#' @export
#' @importFrom fgsea fgsea
#' @importFrom AnnotationDbi mapIds select
regulatorGSEA <- function(
    object,
    ontology = c("BP", "MF", "CC"),
    aggregate = c("extreme", "mean"),
    inferActivity = TRUE,
    orgDb = "org.Hs.eg.db",
    minSize = 15L,
    maxSize = 500L,
    ...
) {
    ontology <- match.arg(ontology)
    aggregate <- match.arg(aggregate)
    db <- .resolveOrgDb(orgDb)

    stats <- .regulatorGeneStats(object, aggregate, inferActivity, db)
    pathways <- .buildGOSets(names(stats), ontology, db)

    res <- fgsea::fgsea(
        pathways = pathways,
        stats = stats,
        minSize = minSize,
        maxSize = maxSize,
        ...
    )
    .tidyGSEAResult(res, stats, db)
}

# Collapse the per-signature z-scores of a perturbMatch result to one
# activity-scaled statistic per regulator gene, keyed by Entrez ID (the
# namespace of the GO gene sets). Direction harmonisation reuses the KD/KO
# sign-flip logic used by plotSimilarity(inferActivity = TRUE).
.regulatorGeneStats <- function(object, aggregate, inferActivity, db) {
    tbl <- regulatorTable(object)
    if (!all(c("pert_symbol", "pert_type") %in% colnames(tbl))) {
        stop(
            "'object' must have 'pert_symbol' and 'pert_type' columns in ",
            "rowData(); is it a perturbMatch result from calcEnrichment()?",
            call. = FALSE
        )
    }

    z <- tbl$zscore
    if (inferActivity) {
        isKnockdown <- toupper(tbl$pert_type) %in% c("KD", "KO")
        z <- ifelse(isKnockdown, -z, z)
    }

    perGene <- if (aggregate == "mean") {
        vapply(split(z, tbl$pert_symbol), mean, numeric(1))
    } else {
        vapply(
            split(z, tbl$pert_symbol),
            function(v) v[which.max(abs(v))],
            numeric(1)
        )
    }

    entrez <- AnnotationDbi::mapIds(
        db,
        keys = names(perGene),
        column = "ENTREZID",
        keytype = "SYMBOL",
        multiVals = "first"
    )
    perGene <- perGene[!is.na(entrez)]
    entrez <- entrez[!is.na(entrez)]

    # If two symbols map to the same Entrez ID, keep the larger-magnitude score.
    ord <- order(abs(perGene), decreasing = TRUE)
    perGene <- perGene[ord]
    entrez <- entrez[ord]
    keep <- !duplicated(entrez)

    stats::setNames(perGene[keep], entrez[keep])
}

# GO gene sets for the supplied Entrez genes, using the GOALL mapping
.buildGOSets <- function(entrez, ontology, db) {
    # A gene belongs to many GO terms and a term to many genes, so select()
    # always reports a many:many mapping. It has no quiet switch, and the
    # notice would fire once per regulatorGSEA() call, so it is suppressed.
    map <- suppressMessages(AnnotationDbi::select(
        db,
        keys = entrez,
        columns = c("GOALL", "ONTOLOGYALL"),
        keytype = "ENTREZID"
    ))
    map <- map[!is.na(map$GOALL) & map$ONTOLOGYALL == ontology, ]
    map <- unique(map[, c("ENTREZID", "GOALL")])
    if (nrow(map) == 0) {
        stop(
            "No ",
            ontology,
            " GO annotations found for the regulator genes.",
            call. = FALSE
        )
    }
    split(map$ENTREZID, map$GOALL)
}

# Turn the fgsea data.table into a tidy data.frame: readable term names (from
# GO.db when available), leading-edge Entrez IDs mapped back to symbols, ordered
# by adjusted p-value.
.tidyGSEAResult <- function(res, stats, db) {
    res <- as.data.frame(res)
    if (nrow(res) == 0) {
        return(data.frame(
            ID = character(),
            term = character(),
            NES = numeric(),
            pval = numeric(),
            padj = numeric(),
            size = integer(),
            leadingEdge = character(),
            stringsAsFactors = FALSE
        ))
    }

    res <- res[order(res$padj), , drop = FALSE]

    entrez2symbol <- AnnotationDbi::mapIds(
        db,
        keys = names(stats),
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first"
    )
    leadingEdge <- vapply(
        res$leadingEdge,
        function(ids) paste(entrez2symbol[ids], collapse = ", "),
        character(1)
    )

    data.frame(
        ID = res$pathway,
        term = .goTerms(res$pathway),
        NES = res$NES,
        pval = res$pval,
        padj = res$padj,
        size = res$size,
        leadingEdge = leadingEdge,
        row.names = NULL,
        stringsAsFactors = FALSE
    )
}

# Human-readable GO term names, from GO.db when installed. Falls back to the GO
# accession so the function still works without the optional dependency.
.goTerms <- function(ids) {
    if (!requireNamespace("GO.db", quietly = TRUE)) {
        return(ids)
    }
    labels <- AnnotationDbi::Term(GO.db::GOTERM)[ids]
    labels[is.na(labels)] <- ids[is.na(labels)]
    unname(labels)
}


#' Dot plot of regulator GO enrichment
#'
#' Draws a dot plot of the GO terms returned by
#' \code{\link{regulatorGSEA}}. Terms are placed by normalised enrichment
#' score, so those enriched among mimicking regulators (positive \code{NES})
#' sit to the right of zero and those enriched among reversing regulators
#' (negative \code{NES}) to the left. Point size shows the GO set size and
#' colour the adjusted p-value.
#'
#' @param x A \code{data.frame} returned by \code{\link{regulatorGSEA}}.
#' @param n Maximum number of terms to show. Taken by largest \code{|NES|} so
#'   both directions are represented. Default \code{15}.
#' @param pvalueCutoff Only terms with \code{padj <= pvalueCutoff} are shown.
#'   Default \code{0.1}.
#' @param wrapWidth Wrap long GO term labels to this many characters per
#'   line. \code{NULL} disables wrapping. Default \code{40}.
#' @param fontSize Base font size passed to
#'   \code{\link[cowplot]{theme_cowplot}}. Default \code{12}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' se <- getReferenceDatabase(demo = TRUE)
#' res <- calcEnrichment(
#'     se,
#'     upRegulated   = rownames(se)[seq_len(5)],
#'     downRegulated = rownames(se)[seq(6, 10)],
#'     background    = rownames(se),
#'     BPPARAM       = BiocParallel::SerialParam(progressbar = FALSE)
#' )
#' go <- regulatorGSEA(res, minSize = 1)
#' if (any(go$padj <= 0.5)) plotRegulatorGSEA(go, pvalueCutoff = 0.5)
#'
#' @seealso \code{\link{regulatorGSEA}}
#' @rdname plotRegulatorGSEA
#' @export
#' @import ggplot2
#' @importFrom cowplot theme_cowplot
#' @importFrom rlang .data
plotRegulatorGSEA <- function(
    x,
    n = 15L,
    pvalueCutoff = 0.1,
    wrapWidth = 40,
    fontSize = 12
) {
    stopifnot(is.data.frame(x))
    df <- x[!is.na(x$padj) & x$padj <= pvalueCutoff, , drop = FALSE]
    if (nrow(df) == 0) {
        stop(
            "No GO terms pass padj <= ",
            pvalueCutoff,
            " to plot.",
            call. = FALSE
        )
    }

    # Keep the top terms by absolute NES so mimicking and reversing regulators
    # are both represented, then order the axis by signed NES.
    df <- df[order(abs(df$NES), decreasing = TRUE), , drop = FALSE]
    df <- utils::head(df, n)
    df$term <- .wrapLabels(df$term, wrapWidth)
    df$term <- factor(df$term, levels = df$term[order(df$NES)])

    ggplot(df, aes(x = .data$NES, y = .data$term)) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
        geom_point(aes(size = .data$size, colour = .data$padj)) +
        scale_colour_viridis_c(option = "plasma", name = "Adjusted\np-value") +
        scale_size(range = c(3, 8), name = "Set size") +
        guides(
            size = guide_legend(order = 1),
            colour = guide_colorbar(order = 2)
        ) +
        labs(
            x = "Normalised enrichment score (NES)",
            y = NULL,
            title = "GO enrichment of regulator hits"
        ) +
        theme_cowplot(font_size = fontSize)
}

# Wrap long GO term labels onto multiple lines
.wrapLabels <- function(labels, width) {
    if (is.null(width)) {
        return(labels)
    }
    vapply(
        labels,
        function(s) paste(strwrap(s, width = width), collapse = "\n"),
        character(1),
        USE.NAMES = FALSE
    )
}
