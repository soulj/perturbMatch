# The four reference builds. Nothing beyond their names is held in this
# package any more, so the keys are stated here rather than read back from it.
referenceKeys <- c("chrdir", "chrdir_sva", "limma", "limma_sva")

# What distinguishes a build travels with the hosted colData, so a stand-in for
# a downloaded annotation has to carry it too.
stubColData <- function(key,
                        sigs = paste0("sig", seq_len(4)),
                        assayNames = c("log2fc", "padj", "pi"),
                        bytes = 1234567) {
    cd <- S4Vectors::DataFrame(pert_symbol = sigs, row.names = sigs)

    S4Vectors::metadata(cd) <- list(
        dataset = key,
        method = if (startsWith(key, "limma")) "limma" else "chrdir",
        batch_correction = endsWith(key, "_sva"),
        assay_names = assayNames,
        ref_metric = if (identical(assayNames, "chrdir")) {
            "chrdir"
        } else {
            c("pi", "log2fc")
        },
        default_metric = if (identical(assayNames, "chrdir")) {
            "chrdir"
        } else {
            "pi"
        },
        assay_bytes = bytes
    )

    cd
}

stubReference <- function() {
    SummarizedExperiment::SummarizedExperiment(
        assays = list(chrdir = matrix(0, nrow = 1, ncol = 1))
    )
}


test_that("getReferenceDatabase loads demo data correctly", {

    se <- getReferenceDatabase(demo = TRUE)

    expect_s4_class(se, "SummarizedExperiment")

    expect_equal(nrow(se), 6500)
    expect_equal(ncol(se), 21)

    expect_s4_class(assays(se)[[1]], "DelayedMatrix")
    expect_true(all(c("log2fc", "padj", "pi") %in% assayNames(se)))

    expect_true("method" %in% names(metadata(se)))
    expect_true("Symbol" %in% colnames(rowData(se)))
    expect_true("cell_line" %in% colnames(colData(se)))
})

test_that("getReferenceDatabase sets the ref_metric/default_metric contract used by calcEnrichment", {
    se <- getReferenceDatabase(demo = TRUE)

    expect_true(all(c("pi", "log2fc") %in% metadata(se)$ref_metric))
    expect_equal(metadata(se)$default_metric, "pi")
    expect_true(all(metadata(se)$ref_metric %in% assayNames(se)))
})

test_that("getReferenceDatabase errors on an invalid type", {
    expect_error(getReferenceDatabase(demo = TRUE, type = "not_a_type"))
})

test_that("getReferenceDatabase does not load the demo database by default", {
    # A silent demo in a script is indistinguishable from a real reference,
    # so it has to be asked for
    expect_false(eval(formals(getReferenceDatabase)$demo))
})

test_that("the demo database announces itself and is flagged as demo", {
    expect_message(
        se <- getReferenceDatabase(demo = TRUE),
        "DEMONSTRATION"
    )

    expect_true(metadata(se)$demo)
    expect_equal(metadata(se)$source, "bundled")
})

test_that("the demo describes itself the way a full reference does", {
    # calcEnrichment() reads the same keys whichever it is handed
    md <- metadata(getReferenceDatabase(demo = TRUE))

    expect_true(all(
        c("dataset", "method", "batch_correction", "assay_names",
          "ref_metric", "default_metric") %in% names(md)
    ))
    expect_equal(md$assay_names, assayNames(getReferenceDatabase(demo = TRUE)))
})

test_that("demo = TRUE warns that type and batchCorrection are ignored", {
    expect_warning(
        getReferenceDatabase(demo = TRUE, type = "chrdir"),
        "ignored when demo = TRUE"
    )
    expect_warning(
        getReferenceDatabase(demo = TRUE, batchCorrection = TRUE),
        "ignored when demo = TRUE"
    )

    expect_warning(getReferenceDatabase(demo = TRUE), NA)
})

test_that("demo = TRUE does not warn about arguments the demo already is", {
    # The demo is the uncorrected limma build, so asking for exactly that is
    # not a mismatch worth warning about
    expect_warning(getReferenceDatabase(demo = TRUE, type = "limma"), NA)
    expect_warning(
        getReferenceDatabase(demo = TRUE, batchCorrection = FALSE), NA
    )
})

test_that(".referenceKey maps type/batchCorrection to distinct reference keys", {
    expect_equal(perturbMatch:::.referenceKey("chrdir", FALSE), "chrdir")
    expect_equal(perturbMatch:::.referenceKey("chrdir", TRUE), "chrdir_sva")
    expect_equal(perturbMatch:::.referenceKey("limma", FALSE), "limma")
    expect_equal(perturbMatch:::.referenceKey("limma", TRUE), "limma_sva")

    # Every combination has to land on a distinct build
    keys <- unlist(lapply(c("chrdir", "limma"), function(type) {
        vapply(
            c(FALSE, TRUE),
            function(bc) perturbMatch:::.referenceKey(type, bc),
            character(1L)
        )
    }))

    expect_setequal(keys, referenceKeys)
    expect_length(unique(keys), length(referenceKeys))
})

test_that(".referenceKey validates batchCorrection", {
    expect_error(perturbMatch:::.referenceKey("limma", NA))
    expect_error(perturbMatch:::.referenceKey("limma", "yes"))
})

test_that(".referenceMetadata reads the description off the annotation", {
    md <- perturbMatch:::.referenceMetadata(stubColData("limma_sva"), "limma_sva")

    expect_equal(md$default_metric, "pi")
    expect_true(all(md$ref_metric %in% md$assay_names))
})

test_that(".referenceMetadata rejects annotation from before the metadata", {
    # A cache left over from an earlier deposit deserialises fine and then
    # fails obscurely when the assays are named from it
    stale <- S4Vectors::DataFrame(pert_symbol = "sig1")

    expect_error(
        perturbMatch:::.referenceMetadata(stale, "chrdir"),
        "missing the reference metadata"
    )
})

test_that("each reference resolves to its three hosted files", {
    for (key in referenceKeys) {
        files <- perturbMatch:::.referenceFiles(key)

        expect_named(files, c("assays", "colData", "rowData"))
        expect_equal(files[["assays"]], paste0(key, "_assays.h5"))
        expect_equal(files[["rowData"]], "perturbMatch_rowData.rds")
    }
})

test_that("each reference resolves to three named download URLs", {
    for (key in referenceKeys) {
        urls <- perturbMatch:::.referenceUrls(key)

        # .fetchDownload() addresses the parts by name, so losing the names
        # here breaks the download only after the file has been fetched
        expect_named(urls, c("assays", "colData", "rowData"))
        expect_true(all(startsWith(urls, "https://zenodo.org/")))
        expect_true(endsWith(urls[["assays"]], paste0(key, "_assays.h5")))
    }
})

test_that(".fetchDownload returns the parts keyed by name", {
    cacheDir <- tempfile()
    dir.create(cacheDir)
    rds <- file.path(cacheDir, "annotation.rds")
    saveRDS(stubColData("chrdir", assayNames = "chrdir"), rds)

    # Stand in for the cache so the test neither downloads nor writes to the
    # user cache. bfcrpath() names its return after the rname it was given.
    requested <- character()
    local_mocked_bindings(
        BiocFileCache = function(...) "cache",
        bfcinfo = function(...) data.frame(rname = character()),
        bfcrpath = function(bfc, rname) {
            requested <<- c(requested, rname)
            path <- if (endsWith(rname, ".h5")) "assays.h5" else rds
            stats::setNames(path, rname)
        },
        .package = "perturbMatch"
    )

    expect_message(
        parts <- perturbMatch:::.fetchDownload("chrdir", cacheDir),
        "one-time download"
    )

    expect_named(parts, c("assays", "colData", "rowData"))
    expect_equal(parts$assays, "assays.h5")
    expect_s4_class(parts$colData, "DataFrame")
    expect_s4_class(parts$rowData, "DataFrame")

    # The annotation is read before the assays so the size can be reported
    urls <- perturbMatch:::.referenceUrls("chrdir")
    expect_equal(requested, unname(urls[c("colData", "rowData", "assays")]))
})

test_that(".fetchDownload reports the size recorded in the deposit", {
    cacheDir <- tempfile()
    dir.create(cacheDir)
    rds <- file.path(cacheDir, "annotation.rds")
    saveRDS(
        stubColData("chrdir", assayNames = "chrdir", bytes = 789743060),
        rds
    )

    local_mocked_bindings(
        BiocFileCache = function(...) "cache",
        bfcinfo = function(...) data.frame(rname = character()),
        bfcrpath = function(bfc, rname) {
            path <- if (endsWith(rname, ".h5")) "assays.h5" else rds
            stats::setNames(path, rname)
        },
        .package = "perturbMatch"
    )

    expect_message(
        perturbMatch:::.fetchDownload("chrdir", cacheDir),
        "753.2 Mb"
    )
})

test_that("each reference resolves to three ExperimentHub record titles", {
    for (key in referenceKeys) {
        titles <- perturbMatch:::.referenceTitles(key)

        expect_named(titles, c("assays", "colData", "rowData"))
        expect_equal(
            titles[["assays"]], paste0("perturbMatch_", key, "_assays")
        )
        expect_equal(titles[["rowData"]], "perturbMatch_rowData")
    }
})

# The record titles are the contract with perturbMatchData. Every other hub
# test here is mocked, so this is what would notice one being renamed.
# skip_if_offline() also skips on CRAN, which would take this out of the
# Bioconductor builds, so the connectivity check is made directly.
test_that("every reference resolves against the hosted records", {
    skip_if_not_installed("curl")
    skip_if(!curl::has_internet(), "no internet connection")

    hub <- ExperimentHub::ExperimentHub()
    hosted <- hub[ExperimentHub::package(hub) == "perturbMatchData"]$title

    # The hub only offers records up to the running Bioconductor version, so an
    # older release sees none of them. That is not a renamed title.
    skip_if(
        length(hosted) == 0L,
        "perturbMatchData is not on the hub of this Bioconductor version"
    )

    for (key in referenceKeys) {
        titles <- perturbMatch:::.referenceTitles(key)
        expect_true(all(titles %in% hosted), info = key)
    }
})

test_that("source = 'auto' falls back to Zenodo only when records are absent", {
    fellBack <- FALSE
    local_mocked_bindings(
        .fetchExperimentHub = function(key, mustWork) NULL,
        .fetchDownload = function(key, cacheDir, ask) {
            fellBack <<- TRUE
            list(assays = "a", colData = "c", rowData = "r")
        },
        .buildReference = function(key, parts) stubReference(),
        .package = "perturbMatch"
    )

    expect_message(
        se <- getReferenceDatabase(type = "chrdir"),
        "Falling back to the Zenodo deposit"
    )

    expect_true(fellBack)
    expect_equal(metadata(se)$source, "zenodo")
    expect_false(metadata(se)$demo)
})

test_that("source = 'auto' does not fall back when the hub itself fails", {
    # An unreachable hub is not the same as an absent record: falling back
    # would silently change where a result's reference came from.
    local_mocked_bindings(
        .fetchExperimentHub = function(key, mustWork) stop("hub unreachable"),
        .fetchDownload = function(key, cacheDir, ask) {
            stop("must not reach the Zenodo fallback")
        },
        .package = "perturbMatch"
    )

    expect_error(getReferenceDatabase(type = "chrdir"), "hub unreachable")
})

test_that("the reference records where it came from", {
    local_mocked_bindings(
        .fetchExperimentHub = function(key, mustWork) {
            list(assays = "a", colData = "c", rowData = "r")
        },
        .buildReference = function(key, parts) stubReference(),
        .package = "perturbMatch"
    )

    expect_equal(
        metadata(getReferenceDatabase(type = "chrdir"))$source,
        "experimenthub"
    )
    expect_equal(
        metadata(getReferenceDatabase(demo = TRUE))$source,
        "bundled"
    )
})

test_that("source = 'zenodo' skips ExperimentHub entirely", {
    local_mocked_bindings(
        .fetchExperimentHub = function(key, mustWork) {
            stop("must not consult ExperimentHub")
        },
        .fetchDownload = function(key, cacheDir, ask) {
            list(assays = "a", colData = "c", rowData = "r")
        },
        .buildReference = function(key, parts) stubReference(),
        .package = "perturbMatch"
    )

    expect_silent(se <- getReferenceDatabase(source = "zenodo"))
    expect_equal(metadata(se)$source, "zenodo")
})

test_that(".buildReference assembles the hosted parts into a reference", {
    genes <- paste0("gene", seq_len(20))
    sigs <- paste0("sig", seq_len(4))

    h5 <- tempfile(fileext = ".h5")
    mat <- matrix(seq_len(80) / 10, nrow = 20, ncol = 4)
    for (name in c("log2fc", "padj", "pi")) {
        HDF5Array::writeHDF5Array(mat, h5, name, with.dimnames = FALSE)
    }

    parts <- list(
        assays = h5,
        rowData = S4Vectors::DataFrame(
            Symbol = genes,
            row.names = genes
        ),
        colData = stubColData("limma_sva", sigs = sigs)
    )

    se <- perturbMatch:::.buildReference("limma_sva", parts)

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(assayNames(se), c("log2fc", "padj", "pi"))
    expect_s4_class(assay(se, "pi"), "DelayedMatrix")

    # The HDF5 datasets carry no dimnames, so they have to come from the
    # annotation
    expect_equal(rownames(se), genes)
    expect_equal(colnames(se), sigs)

    expect_equal(metadata(se)$default_metric, "pi")
    expect_true(metadata(se)$batch_correction)

    # The description belongs to the object, not to the annotation inside it
    expect_length(S4Vectors::metadata(colData(se)), 0L)
})

test_that(".buildReference names the assays the deposit says it holds", {
    h5 <- tempfile(fileext = ".h5")
    HDF5Array::writeHDF5Array(
        matrix(1, nrow = 20, ncol = 4),
        h5,
        "chrdir",
        with.dimnames = FALSE
    )

    parts <- list(
        assays = h5,
        rowData = S4Vectors::DataFrame(row.names = paste0("gene", seq_len(20))),
        colData = stubColData(
            "chrdir",
            sigs = paste0("sig", seq_len(4)),
            assayNames = "chrdir"
        )
    )

    se <- perturbMatch:::.buildReference("chrdir", parts)

    expect_equal(assayNames(se), "chrdir")
    expect_equal(metadata(se)$ref_metric, "chrdir")
})

test_that(".buildReference rejects annotation that does not match the assay", {
    h5 <- tempfile(fileext = ".h5")
    HDF5Array::writeHDF5Array(
        matrix(1, nrow = 20, ncol = 4),
        h5,
        "chrdir",
        with.dimnames = FALSE
    )

    parts <- list(
        assays = h5,
        rowData = S4Vectors::DataFrame(row.names = paste0("gene", seq_len(19))),
        colData = stubColData(
            "chrdir",
            sigs = paste0("sig", seq_len(4)),
            assayNames = "chrdir"
        )
    )

    expect_error(
        perturbMatch:::.buildReference("chrdir", parts),
        "does not match its annotation"
    )
})
