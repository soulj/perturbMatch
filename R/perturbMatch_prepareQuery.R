#' Prepare a differential expression result for calcEnrichment
#'
#' Extracts gene identifiers and log2 fold changes from a differential
#' expression result, maps the identifiers to the Entrez gene IDs used by
#' the perturbMatch reference databases, and reports the fraction
#' of the query that was successfully mapped. The returned data frame can be
#' passed directly to \code{\link{calcEnrichment}} as \code{expressionMatrix}.
#'
#' @param object A differential expression result. One of:
#'   \itemize{
#'     \item a \pkg{DESeq2} \code{DESeqResults} object.
#'     \item a \code{data.frame} (e.g. \pkg{limma} \code{topTable()}
#'       output, or a plain table read from a CSV/TSV) with gene identifiers
#'       in its row names or in a column, and a fold-change column named
#'       \code{"log2FoldChange"}, \code{"logFC"} or \code{"log2FC"}. An
#'       adjusted p-value column named \code{"padj"} or \code{"adj.P.Val"},
#'       if present, is carried through.
#'   }
#'
#'   Columns are matched by name, not position, and any further columns are
#'   ignored. Only \emph{adjusted} p-values are recognised, so a
#'   \code{"pvalue"} or \code{"P.Value"} column does not drive
#'   \code{dropLowExpressed}.
#'
#'   Gene identifiers are taken from the row names when those are real names,
#'   as in \code{DESeqResults} and \code{topTable()} output. Positional row
#'   names, as produced by \code{data.frame()} and \code{read.csv()}, are
#'   ignored, as are all-digit row names when a column names itself as an
#'   identifier. The identifiers are then read from a column named like an
#'   identifier (\code{"gene"}, \code{"gene_id"}, \code{"symbol"},
#'   \code{"entrezid"}, \code{"id"}, \ldots) if present, otherwise from the
#'   first column that is neither the fold-change nor the adjusted p-value
#'   column. Which of the two was used is reported in a message.
#' @param orgDb An \pkg{AnnotationDbi} \code{OrgDb} object, or the name of an
#'   installed \code{OrgDb} package, used for ID mapping. Default
#'   \code{"org.Hs.eg.db"}.
#' @param keytype The type of the query gene identifiers (e.g.
#'   \code{"ENSEMBL"}, \code{"SYMBOL"}), passed to
#'   \code{\link[AnnotationDbi]{mapIds}}. When \code{NULL} (default) it is
#'   guessed from the identifiers. Identifiers that are already Entrez IDs
#'   are used as-is without mapping.
#' @param reference Optional reference database (a
#'   \code{SummarizedExperiment}). When supplied, the fraction of the mapped
#'   query present in \code{rownames(reference)} is also reported.
#' @param minMappedFraction Numeric in \eqn{[0, 1]}. A warning is emitted if
#'   the fraction of query genes mapped to Entrez falls below this value.
#'   Default \code{0.5}.
#' @param dropLowExpressed Logical. When the input carries an adjusted p-value
#'   column, tools such as \pkg{DESeq2} set it to \code{NA} for genes that were
#'   filtered out before testing, typically because they are too lowly
#'   expressed. When \code{TRUE} (default) these genes are dropped and the
#'   number removed is reported. Has no effect when the input has no adjusted
#'   p-value column.
#'
#' @return A \code{data.frame} with a \code{gene} column (Entrez IDs) and a
#'   \code{log2FC} column, plus a \code{padj} column when one was available
#'   in \code{object}. Rows with a missing fold change or an unmapped
#'   identifier are dropped, and identifiers mapping to the same Entrez ID
#'   are collapsed to the one with the largest absolute fold change.
#'
#' @examples
#' # A limma topTable-style data frame (gene IDs in the row names)
#' tt <- data.frame(
#'     logFC = c(2.5, -1.8, 0.3),
#'     adj.P.Val = c(0.001, 0.02, 0.4),
#'     row.names = c("100", "1000", "10000")
#' )
#' # IDs are already Entrez, so no OrgDb is needed
#' prepareQuery(tt)
#'
#' @seealso \code{\link{calcEnrichment}}
#' @importFrom AnnotationDbi mapIds keytypes
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @export
prepareQuery <- function(
    object,
    orgDb = "org.Hs.eg.db",
    keytype = NULL,
    reference = NULL,
    minMappedFraction = 0.5,
    dropLowExpressed = TRUE
) {
    de <- .extractDEResult(object)
    de <- de[!is.na(de$log2FC), , drop = FALSE]

    if (dropLowExpressed) {
        de <- .dropLowExpressed(de)
    }

    de$entrez <- .mapToEntrez(de$gene, orgDb, keytype)
    .reportMapping(de$entrez, minMappedFraction)

    de <- de[!is.na(de$entrez), , drop = FALSE]
    de <- .collapseDuplicateEntrez(de)

    if (!is.null(reference)) {
        inReference <- mean(de$entrez %in% rownames(reference))
        message(sprintf(
            "%.1f%% of mapped genes are present in the reference database.",
            100 * inReference
        ))
    }

    out <- data.frame(
        gene = de$entrez,
        log2FC = de$log2FC,
        stringsAsFactors = FALSE
    )
    if ("padj" %in% colnames(de)) {
        out$padj <- de$padj
    }
    out
}

# A low mapping rate almost always means the wrong keytype was guessed or
# supplied, so the rate is reported either way and warned about below the
# caller's threshold.
.reportMapping <- function(entrez, minMappedFraction) {
    mapped <- !is.na(entrez)
    mappedFraction <- mean(mapped)
    message(sprintf(
        "Mapped %.1f%% of query genes to Entrez (%d / %d).",
        100 * mappedFraction,
        sum(mapped),
        length(entrez)
    ))
    if (mappedFraction < minMappedFraction) {
        warning(
            "Only ",
            round(100 * mappedFraction),
            "% of query genes mapped ",
            "to Entrez. Check that 'keytype' matches your gene IDs.",
            call. = FALSE
        )
    }
    invisible(mappedFraction)
}

# When the input carries adjusted p-values, tools such as DESeq2 set them to NA
# for genes filtered out before testing (typically too lowly expressed). Those
# carry no reliable signal, so drop them and report how many. A no-op when the
# input has no adjusted p-value column.
.dropLowExpressed <- function(de) {
    if (!"padj" %in% colnames(de)) {
        return(de)
    }
    nDropped <- sum(is.na(de$padj))
    if (nDropped == 0) {
        return(de)
    }
    de <- de[!is.na(de$padj), , drop = FALSE]
    msg <- sprintf(
        paste(
            "Dropped %d genes with a missing adjusted p-value",
            "(not tested, typically lowly expressed). Set",
            "dropLowExpressed = FALSE to keep them."
        ),
        nDropped
    )
    message(msg)
    de
}

# Pull (gene, log2FC[, padj]) out of a DE result. as.data.frame() covers
# DESeqResults (a DataFrame subclass) and limma topTable output alike. Gene
# IDs come from the row names, falling back to the first column. Columns are
# matched by name, so their order is irrelevant and any extra columns are
# dropped.
.extractDEResult <- function(object) {
    df <- as.data.frame(object)

    fcCol <- intersect(c("log2FoldChange", "logFC", "log2FC"), colnames(df))
    if (length(fcCol) == 0) {
        stop(
            "Could not find a fold-change column ('log2FoldChange', ",
            "'logFC' or 'log2FC') in the input.",
            call. = FALSE
        )
    }
    pCol <- intersect(c("padj", "adj.P.Val"), colnames(df))

    out <- data.frame(
        gene = .extractGeneIds(df, c(fcCol, pCol)),
        log2FC = df[[fcCol[1]]],
        stringsAsFactors = FALSE
    )
    if (length(pCol) > 0) {
        out$padj <- df[[pCol[1]]]
    }
    out
}

# Gene IDs live either in the row names or in a column. A data frame built by
# data.frame() or read.csv() carries automatic row names, so a column is used
# instead: one with an ID-like name if there is one, otherwise the first column
# not already claimed above.
.extractGeneIds <- function(df, valueCols) {
    candidates <- setdiff(colnames(df), valueCols)
    named <- candidates[
        tolower(candidates) %in%
            c(
                "gene",
                "genes",
                "gene_id",
                "geneid",
                "gene_name",
                "genename",
                "id",
                "symbol",
                "gene_symbol",
                "entrez",
                "entrezid",
                "ensembl"
            )
    ]

    rn <- attr(df, "row.names")
    if (.rowNamesAreIds(rn, length(named) > 0)) {
        message("Taking gene identifiers from the row names.")
        return(rn)
    }

    idCol <- if (length(named) > 0) named[1] else candidates[1]
    if (is.na(idCol)) {
        stop(
            "Could not find gene identifiers in the input: the row names are ",
            "positional (1, 2, 3, ...) and every column holds values rather ",
            "than identifiers. Put the gene IDs in the row names or in a ",
            "column of their own.",
            call. = FALSE
        )
    }
    message(sprintf("Taking gene identifiers from column '%s'.", idCol))
    as.character(df[[idCol]])
}

# Automatic row names are positions rather than identifiers. They are usually
# an integer vector, but a write.csv() / read.csv(row.names = 1) can
# store them as the strings "1", "2", ..., and a row subset leaves gaps in
# those ("1", "3", ...), both of which read as Entrez IDs. All-digit row names
# therefore lose to a column that names itself as an identifier.
.rowNamesAreIds <- function(rn, hasNamedIdCol) {
    if (!is.character(rn) || identical(rn, as.character(seq_along(rn)))) {
        return(FALSE)
    }
    !(hasNamedIdCol && all(grepl("^[0-9]+$", rn)))
}

# Map query IDs to Entrez. Already-Entrez IDs are returned unchanged.
# Ensembl version suffixes (".12") are stripped before mapping.
.mapToEntrez <- function(ids, orgDb, keytype) {
    ids <- as.character(ids)
    if (is.null(keytype)) {
        keytype <- .guessKeytype(ids)
    }
    if (identical(keytype, "ENTREZID")) {
        return(ids)
    }
    if (identical(keytype, "ENSEMBL")) {
        ids <- sub("\\..*$", "", ids)
    }
    db <- .resolveOrgDb(orgDb)
    if (!keytype %in% AnnotationDbi::keytypes(db)) {
        stop(
            "'",
            keytype,
            "' is not a valid keytype for the supplied ",
            "OrgDb. See AnnotationDbi::keytypes().",
            call. = FALSE
        )
    }
    # mapIds() errors (rather than returning all NA) when none of the keys
    # are valid. Degrade that to all-unmapped so the low-mapping warning in
    # prepareQuery() reports it clearly instead.
    mapped <- tryCatch(
        AnnotationDbi::mapIds(
            db,
            keys = ids,
            column = "ENTREZID",
            keytype = keytype,
            multiVals = "first"
        ),
        error = function(e) rep(NA_character_, length(ids))
    )
    unname(mapped)
}

.resolveOrgDb <- function(orgDb) {
    if (!is.character(orgDb)) {
        return(orgDb)
    }
    if (!requireNamespace(orgDb, quietly = TRUE)) {
        stop(
            "Package '",
            orgDb,
            "' is required for gene ID mapping but is ",
            "not installed.\n",
            "  Install it with: BiocManager::install(\"",
            orgDb,
            "\")\n",
            "  Or pass a different OrgDb via 'orgDb' (e.g. \"org.Mm.eg.db\" ",
            "for mouse). If your gene IDs are already Entrez IDs, set ",
            "keytype = \"ENTREZID\" to skip mapping entirely.",
            call. = FALSE
        )
    }
    getExportedValue(orgDb, orgDb)
}

.guessKeytype <- function(ids) {
    ids <- ids[!is.na(ids)]
    if (length(ids) == 0) {
        return("SYMBOL")
    }
    if (all(grepl("^[0-9]+$", ids))) {
        return("ENTREZID")
    }
    if (all(grepl("^ENS[A-Z]*G[0-9]+", ids))) {
        return("ENSEMBL")
    }
    "SYMBOL"
}

# When several query IDs map to the same Entrez ID, keep the one with the
# largest absolute fold change.
.collapseDuplicateEntrez <- function(de) {
    if (anyDuplicated(de$entrez) == 0) {
        return(de)
    }
    ord <- order(abs(de$log2FC), decreasing = TRUE)
    de <- de[ord, , drop = FALSE]
    de[!duplicated(de$entrez), , drop = FALSE]
}
