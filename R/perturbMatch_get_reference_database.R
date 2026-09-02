#' Load the perturbMatch reference database
#'
#' Retrieve a perturbMatch reference expression database as a
#' [SummarizedExperiment::SummarizedExperiment-class] object.
#'
#' @param demo Logical scalar indicating whether to load the small bundled
#'   demonstration dataset instead of a full reference. Defaults to
#'   \code{FALSE}. The demo ignores \code{type} and \code{batchCorrection}. See
#'   Details for what it is and is not good for.
#' @param type Character scalar specifying which reference database to load.
#'   One of \code{"chrdir"} (characteristic direction) or \code{"limma"}
#'   (log2 fold-change and adjusted p-values). Defaults to \code{"chrdir"}.
#' @param batchCorrection Logical scalar selecting the batch-corrected variant
#'   of a reference. Available for both \code{type = "chrdir"} and
#'   \code{type = "limma"}: \code{FALSE} (default) loads the uncorrected
#'   signatures and \code{TRUE} loads signatures computed with study/batch
#'   effects removed (SVA).
#' @param source Character scalar selecting where the full reference is
#'   retrieved from. \code{"experimenthub"} takes it from Bioconductor's
#'   \pkg{ExperimentHub} and errors if the records are not registered there.
#'   \code{"zenodo"} downloads it from the Zenodo deposit instead.
#'   \code{"auto"} (the default) uses ExperimentHub, falling back to Zenodo
#'   only when the records are absent from the hub, and announces which of the
#'   two it used. The records are registered for Bioconductor 3.24 and later,
#'   so an older release reaches the data through the Zenodo fallback.
#'   Ignored when \code{demo = TRUE}.
#' @param cacheDir Directory holding the \pkg{BiocFileCache} cache of
#'   downloaded reference files. Defaults to a per-user cache directory
#'   (\code{tools::R_user_dir("perturbMatch", "cache")}). Only used by the
#'   Zenodo download path, as ExperimentHub manages its own cache.
#' @param ask Logical scalar controlling whether \pkg{BiocFileCache} asks
#'   before creating \code{cacheDir}. Defaults to \code{interactive()}, so a
#'   session with a user present is prompted and a scripted one creates the
#'   directory. Declining the prompt, or passing \code{TRUE} outside an
#'   interactive session, puts the cache under \code{tempdir()} and so
#'   re-downloads the reference in every session.
#'
#' @return A [SummarizedExperiment::SummarizedExperiment-class] object with
#'   assay(s) stored as HDF5-backed delayed arrays. The object carries
#'   \code{metadata(object)$ref_metric}, which tells
#'   \code{\link{calcEnrichment}} how to compute the scoring matrix from
#'   the assay(s):
#'   \describe{
#'     \item{\code{"chrdir"}}{Single assay named \code{"chrdir"}, used
#'       directly.}
#'     \item{\code{"log2fc"}}{Assay named \code{"log2fc"}, used directly.}
#'     \item{\code{"pi"}}{Combined \eqn{\pi = \log_2FC \times -\log_{10}(padj)}
#'       score stored alongside the \code{"log2fc"} and \code{"padj"} assays.}
#'   }
#'   \code{ref_metric} is fixed when the reference database is built, travels
#'   with the data, and does not need to be set or modified by the user.
#'   Alongside it are \code{method}, \code{batch_correction},
#'   \code{default_metric} and the \code{build_version} and \code{source_doi}
#'   of the deposit the reference came from.
#'   \code{metadata(object)$source} records how it was retrieved
#'   (\code{"experimenthub"}, \code{"zenodo"} or \code{"bundled"}), and
#'   \code{metadata(object)$demo} whether it is the bundled demo database.
#'   All are carried through to the
#'   \link[=PerturbMatch-class]{PerturbMatch} object returned by
#'   \code{\link{calcEnrichment}}.
#'
#' @details
#' Each reference is hosted as three files: an HDF5 file of the assay matrices
#' (0.7 GB for \code{chrdir_sva} up to 3.1 GB for \code{limma}), a serialised
#' \code{colData} of the signature annotation, and one gene annotation shared
#' by all four references. They are fetched once, cached, and reassembled into
#' a \code{SummarizedExperiment} whose assays are read from disk on demand.
#'
#' The resources are registered on ExperimentHub, and documented by
#' the \pkg{perturbMatchData} package, which does not need to be installed to
#' use them.
#' The signatures are derived from public
#' NCBI GEO RNA-seq series and released under CC-BY-4.0,
#' \doi{10.5281/zenodo.21792179}.
#'
#' The bundled demo mirrors the structure of the full databases so the package
#' can be used offline. It is a 21-signature subset tailored to the
#' vignette examples. The rankings it produces are
#' illustrative only. Any real analysis needs a full reference
#' (\code{demo = FALSE}, the default).
#'
#' @examples
#' # Load the bundled demonstration database (no download required)
#' se <- getReferenceDatabase(demo = TRUE)
#' se
#' S4Vectors::metadata(se)$ref_metric
#'
#' # A full reference is a multi-gigabyte one-time download, so it is fetched
#' # only in an interactive session
#' if (interactive()) {
#'     seChrdir <- getReferenceDatabase(type = "chrdir")
#'
#'     # chrdir reference, with batch-effect correction
#'     seChrdirBc <- getReferenceDatabase(
#'         type = "chrdir",
#'         batchCorrection = TRUE
#'     )
#'
#'     # limma reference, without batch-effect correction
#'     seLimma <- getReferenceDatabase(type = "limma", batchCorrection = FALSE)
#' }
#'
#' @importFrom HDF5Array loadHDF5SummarizedExperiment HDF5Array
#' @importFrom BiocFileCache BiocFileCache bfcinfo bfcrpath
#' @importFrom ExperimentHub ExperimentHub package
#' @importFrom tools R_user_dir
#' @export
getReferenceDatabase <- function(
    demo = FALSE,
    type = c("chrdir", "limma"),
    batchCorrection = FALSE,
    source = c("auto", "experimenthub", "zenodo"),
    cacheDir = tools::R_user_dir("perturbMatch", "cache"),
    ask = interactive()
) {
    typeGiven <- !missing(type)
    type <- match.arg(type)
    source <- match.arg(source)

    if (demo) {
        # The demo is always the uncorrected limma ref
        otherBuildRequested <- (typeGiven && type != "limma") ||
            isTRUE(batchCorrection)
        return(.loadDemo(otherBuildRequested))
    }

    key <- .referenceKey(type, batchCorrection)

    # Under "auto" only an absent record falls through to Zenodo
    parts <- NULL
    if (source != "zenodo") {
        parts <- .fetchExperimentHub(key, mustWork = source == "experimenthub")
    }

    used <- "experimenthub"
    if (is.null(parts)) {
        if (source == "auto") {
            message(
                "The perturbMatchData records are not on the ExperimentHub ",
                "of this Bioconductor version. Falling back to the Zenodo ",
                "deposit. Pass source = \"zenodo\" to select it directly."
            )
        }
        parts <- .fetchDownload(key, cacheDir, ask)
        used <- "zenodo"
    }

    # Both routes supply the same three files, so they assemble the same way
    se <- .buildReference(key, parts)

    # How the object was retrieved, unlike everything else in metadata(), is a
    # property of this call rather than of the reference
    S4Vectors::metadata(se)$demo <- FALSE
    S4Vectors::metadata(se)$source <- used
    se
}


# The demo scores to a ranked table with plausible-looking z-scores, so it is
# flagged on load here and again by show() on any result derived from it.
.loadDemo <- function(otherBuildRequested) {
    if (otherBuildRequested) {
        warning(
            "'type' and 'batchCorrection' are ignored when demo = TRUE: the ",
            "demo database is always the uncorrected limma build.",
            call. = FALSE
        )
    }

    se <- loadHDF5SummarizedExperiment(
        system.file("extdata", "demo_database", package = "perturbMatch")
    )
    S4Vectors::metadata(se)$demo <- TRUE
    S4Vectors::metadata(se)$source <- "bundled"

    message(
        "Loaded the bundled DEMONSTRATION database (",
        nrow(se),
        " genes x ",
        ncol(se),
        " signatures). It is a small illustrative subset of the ",
        "limma reference and should not be used for real analysis. Use ",
        "getReferenceDatabase(demo = FALSE) for a full reference database."
    )

    se
}


.referenceKey <- function(type, batchCorrection) {
    if (
        !is.logical(batchCorrection) ||
            length(batchCorrection) != 1L ||
            is.na(batchCorrection)
    ) {
        stop("'batchCorrection' must be a single TRUE or FALSE.", call. = FALSE)
    }

    if (batchCorrection) paste0(type, "_sva") else type
}

# The three hosted files a reference is assembled from. The gene annotation is
# shared by all four references, so it is named and cached once.
.referenceFiles <- function(key) {
    c(
        assays = paste0(key, "_assays.h5"),
        colData = paste0(key, "_colData.rds"),
        rowData = "perturbMatch_rowData.rds"
    )
}


# ExperimentHub record titles, one per hosted file. These are fixed when
# perturbMatchData is ingested, so they are the stable handle on the resources.
.referenceTitles <- function(key) {
    c(
        assays = paste0("perturbMatch_", key, "_assays"),
        colData = paste0("perturbMatch_", key, "_colData"),
        rowData = "perturbMatch_rowData"
    )
}


# The reference databases are registered on ExperimentHub by perturbMatchData,
# which documents and cites the deposit rather than serving it: the records are
# resolved here so that nothing beyond ExperimentHub has to be installed.
# Returns NULL when the records are not on the hub so the caller can fall back.
# Every other failure is raised.
.fetchExperimentHub <- function(key, mustWork) {
    titles <- .referenceTitles(key)
    hub <- ExperimentHub::ExperimentHub()

    # Subsetting on the preparer package is exact, where matching titles across
    # the whole hub would also match another package's record of the same name
    records <- hub[ExperimentHub::package(hub) == "perturbMatchData"]

    # A reference needs all three records, so a partial ingest is treated the
    # same as no records at all.
    absent <- titles[!titles %in% records$title]
    if (length(absent)) {
        if (mustWork) {
            stop(
                "No ExperimentHub record titled ",
                paste0("'", absent, "'", collapse = ", "),
                " found. The records are registered for Bioconductor 3.24 ",
                "and later, so an older release will not see them. Use ",
                "source = \"zenodo\" to take the deposit directly.",
                call. = FALSE
            )
        }
        return(NULL)
    }

    # DispatchClass H5File returns the path of the cached file, Rds returns
    # the deserialised object.
    lapply(titles, function(title) {
        record <- records[records$title == title]
        if (length(record) != 1L) {
            stop(
                "The title '",
                title,
                "' matches ",
                length(record),
                " perturbMatchData records, so the right one cannot be ",
                "chosen. Please report this against perturbMatchData.",
                call. = FALSE
            )
        }
        record[[1L]]
    })
}


# Base URL of the Zenodo record hosting the reference files. The record id is
# tied to a specific version of the deposit, not the concept DOI.
.zenodoUrl <- "https://zenodo.org/records/21825475/files"

# Zenodo URL of each file backing a reference, named by the part it supplies.
# paste0() drops names, so they are reapplied: everything downstream addresses
# the parts by name rather than by position.
.referenceUrls <- function(key) {
    files <- .referenceFiles(key)
    setNames(paste0(.zenodoUrl, "/", files), names(files))
}

# Download the files backing a reference on first use, then read them.
# bfcrpath() adds each URL to the cache and downloads it, or returns the path
# of the already-cached copy.
.fetchDownload <- function(key, cacheDir, ask = interactive()) {
    urls <- .referenceUrls(key)

    # BiocFileCache answers its own prompt with "no" when there is nobody to
    # ask, and silently caches under tempdir() instead. That costs a multi-GB
    # download per session, so it is worth saying so before it happens.
    if (ask && !interactive() && !dir.exists(cacheDir)) {
        warning(
            "'ask = TRUE' outside an interactive session leaves the ",
            "reference cached under tempdir(), so it is downloaded again ",
            "in every session. Pass ask = FALSE to cache it at ",
            cacheDir,
            " instead.",
            call. = FALSE
        )
    }

    bfc <- BiocFileCache(cacheDir, ask = ask)

    # bfcrpath() names its return after the rname, so it is stripped
    fetch <- function(part) unname(bfcrpath(bfc, urls[[part]]))

    # The annotation is a few MB and carries the size of the assay file, so it
    # is read first and the real figure reported before the large download
    colData <- readRDS(fetch("colData"))
    rowData <- readRDS(fetch("rowData"))

    if (!urls[["assays"]] %in% bfcinfo(bfc)$rname) {
        size <- .referenceMetadata(colData, key)$assay_bytes
        message(
            "Downloading the '",
            key,
            "' reference database (",
            format(structure(size, class = "object_size"), units = "auto"),
            "). This is a one-time download; it ",
            "will be cached at ",
            cacheDir,
            " and reused on subsequent calls."
        )
    }

    list(
        assays = fetch("assays"),
        colData = colData,
        rowData = rowData
    )
}


# What distinguishes the four references is written into the hosted colData
# when they are built, so it is read back rather than restated here. A resource
# cached before the metadata was added deserialises fine and then fails
# obscurely when the assays are named from it, hence the check on read.
.referenceMetadata <- function(colData, key) {
    md <- S4Vectors::metadata(colData)
    required <- c("assay_names", "ref_metric", "default_metric")

    if (!all(required %in% names(md))) {
        stop(
            "The '",
            key,
            "' signature annotation is missing the reference metadata (",
            paste(setdiff(required, names(md)), collapse = ", "),
            "). A cached copy may be left over from an earlier release; ",
            "clear the cache directory and retry.",
            call. = FALSE
        )
    }

    md
}


# Assemble the hosted files into a SummarizedExperiment. The HDF5 datasets
# carry no dimnames, they come from the row and column annotation.
.buildReference <- function(key, parts) {
    md <- .referenceMetadata(parts$colData, key)

    assays <- lapply(md$assay_names, function(name) {
        HDF5Array::HDF5Array(parts$assays, name)
    })
    names(assays) <- md$assay_names

    dims <- dim(assays[[1L]])
    if (
        dims[[1L]] != nrow(parts$rowData) ||
            dims[[2L]] != nrow(parts$colData)
    ) {
        stop(
            "The '",
            key,
            "' assay matrix (",
            dims[[1L]],
            " x ",
            dims[[2L]],
            ") does not match its annotation (",
            nrow(parts$rowData),
            " genes, ",
            nrow(parts$colData),
            " signatures). A cached file may ",
            "be incomplete or left over from an earlier release.",
            call. = FALSE
        )
    }

    # Moved onto the assembled object so the description is not left repeated
    # on the colData inside it
    colData <- parts$colData
    S4Vectors::metadata(colData) <- list()

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = assays,
        rowData = parts$rowData,
        colData = colData
    )

    S4Vectors::metadata(se) <- md
    se
}
