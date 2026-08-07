library(testthat)
library(SummarizedExperiment)

# .guessKeytype

test_that(".guessKeytype recognises Entrez, Ensembl, and symbols", {
    expect_equal(perturbMatch:::.guessKeytype(c("100", "2000")), "ENTREZID")
    expect_equal(
        perturbMatch:::.guessKeytype(c("ENSG00000141510.12", "ENSG00000012048")),
        "ENSEMBL"
    )
    expect_equal(perturbMatch:::.guessKeytype(c("TP53", "BRCA1")), "SYMBOL")
})

test_that(".guessKeytype falls back to SYMBOL when all IDs are NA", {
    expect_equal(
        perturbMatch:::.guessKeytype(c(NA_character_, NA_character_)),
        "SYMBOL"
    )
})

# .resolveOrgDb

test_that(".resolveOrgDb returns a non-character OrgDb object unchanged", {
    fakeDb <- structure(list(), class = "OrgDb")
    expect_identical(perturbMatch:::.resolveOrgDb(fakeDb), fakeDb)
})

test_that(".resolveOrgDb errors when the named package is not installed", {
    expect_error(
        perturbMatch:::.resolveOrgDb("no.such.orgdb.pkg.xyz"),
        "is required for gene ID mapping but is not installed"
    )
})

# .extractDEResult

test_that(".extractDEResult reads a limma topTable-style data frame", {
    tt <- data.frame(
        logFC = c(2.5, -1.8, 0.3),
        adj.P.Val = c(0.001, 0.02, 0.4),
        row.names = c("100", "1000", "10000")
    )

    de <- perturbMatch:::.extractDEResult(tt)

    expect_equal(de$gene, c("100", "1000", "10000"))
    expect_equal(de$log2FC, c(2.5, -1.8, 0.3))
    expect_equal(de$padj, c(0.001, 0.02, 0.4))
})

test_that(".extractDEResult reads a DESeq2-style data frame", {
    df <- data.frame(
        log2FoldChange = c(1, -2),
        padj = c(0.01, 0.2),
        row.names = c("g1", "g2")
    )

    de <- perturbMatch:::.extractDEResult(df)

    expect_equal(de$gene, c("g1", "g2"))
    expect_equal(de$log2FC, c(1, -2))
    expect_equal(de$padj, c(0.01, 0.2))
})

test_that(".extractDEResult takes gene IDs from the first column when row names are default", {
    df <- data.frame(
        gene = c("g1", "g2"),
        logFC = c(1, -2),
        stringsAsFactors = FALSE
    )

    de <- perturbMatch:::.extractDEResult(df)

    expect_equal(de$gene, c("g1", "g2"))
    expect_false("padj" %in% colnames(de))
})

test_that(".extractDEResult errors when no fold-change column is present", {
    expect_error(
        perturbMatch:::.extractDEResult(data.frame(x = 1:3, y = 4:6)),
        "fold-change column"
    )
})

test_that(".extractDEResult matches columns by name, not position", {
    df <- data.frame(
        padj = c(0.01, 0.2),
        baseMean = c(10, 20),
        gene_id = c("g1", "g2"),
        pvalue = c(0.001, 0.1),
        logFC = c(1, -2),
        stringsAsFactors = FALSE
    )

    de <- perturbMatch:::.extractDEResult(df)

    expect_equal(colnames(de), c("gene", "log2FC", "padj"))
    expect_equal(de$gene, c("g1", "g2"))
    expect_equal(de$log2FC, c(1, -2))
    expect_equal(de$padj, c(0.01, 0.2))
})

test_that(".extractDEResult ignores positional row names left by a row subset", {
    df <- data.frame(
        gene = c("g1", "g2", "g3"),
        logFC = c(1, -2, 3),
        padj = c(0.01, 0.5, 0.02),
        stringsAsFactors = FALSE
    )

    # Subsetting rows leaves the automatic row names with gaps ("1", "3"),
    # which must not be mistaken for gene IDs.
    de <- perturbMatch:::.extractDEResult(df[df$padj < 0.05, ])

    expect_equal(de$gene, c("g1", "g3"))
    expect_equal(de$log2FC, c(1, 3))
})

test_that(".extractDEResult keeps real row names through a row subset", {
    df <- data.frame(
        logFC = c(1, -2, 3),
        padj = c(0.01, 0.5, 0.02),
        row.names = c("g1", "g2", "g3")
    )

    de <- perturbMatch:::.extractDEResult(df[df$padj < 0.05, ])

    expect_equal(de$gene, c("g1", "g3"))
})

test_that(".extractDEResult errors when no column holds gene identifiers", {
    expect_error(
        perturbMatch:::.extractDEResult(data.frame(logFC = c(1, -2), padj = c(0.01, 0.2))),
        "Could not find gene identifiers"
    )
})

test_that(".extractDEResult prefers an ID column to all-digit row names", {
    df <- data.frame(
        ENTREZID = c(7105, 64102, 8813),
        log2FoldChange = c(-0.36, -0.07, 0.05)
    )
    # A write.csv() / read.csv(row.names = 1) round trip stores the automatic
    # row names as characters, which would otherwise read as Entrez IDs.
    rownames(df) <- as.character(seq_len(nrow(df)))

    de <- perturbMatch:::.extractDEResult(df)

    expect_equal(de$gene, c("7105", "64102", "8813"))
})

test_that(".extractDEResult ignores character positional row names", {
    df <- data.frame(gene = c("g1", "g2"), logFC = c(1, -2))
    rownames(df) <- c("1", "2")

    expect_equal(perturbMatch:::.extractDEResult(df)$gene, c("g1", "g2"))
})

test_that(".extractDEResult keeps all-digit row names when no ID column exists", {
    df <- data.frame(
        baseMean = c(10, 20),
        log2FoldChange = c(1, -2),
        row.names = c("7105", "64102")
    )

    expect_equal(perturbMatch:::.extractDEResult(df)$gene, c("7105", "64102"))
})

# .collapseDuplicateEntrez

test_that(".collapseDuplicateEntrez keeps the largest absolute fold change", {
    de <- data.frame(
        gene = c("a", "b", "c"),
        log2FC = c(1, -3, 2),
        entrez = c("10", "10", "20"),
        stringsAsFactors = FALSE
    )

    collapsed <- perturbMatch:::.collapseDuplicateEntrez(de)

    expect_equal(nrow(collapsed), 2)
    # entrez 10 keeps gene "b" (|-3| > |1|)
    expect_equal(collapsed$log2FC[collapsed$entrez == "10"], -3)
})

# prepareQuery (Entrez input: no OrgDb needed)

make_entrez_query <- function() {
    data.frame(
        logFC = c(2.5, -1.8, 0.3, 1.0),
        adj.P.Val = c(0.001, 0.02, 0.4, 0.03),
        row.names = c("100", "1000", "10000", "9999")
    )
}

test_that("prepareQuery returns gene/log2FC/padj for already-Entrez input", {
    out <- suppressMessages(prepareQuery(make_entrez_query()))

    expect_s3_class(out, "data.frame")
    expect_equal(colnames(out), c("gene", "log2FC", "padj"))
    expect_equal(out$gene, c("100", "1000", "10000", "9999"))
})

test_that("prepareQuery reports the mapped fraction as a message", {
    expect_message(
        prepareQuery(make_entrez_query()),
        "Mapped 100.0% of query genes to Entrez"
    )
})

test_that("prepareQuery drops rows with a missing fold change", {
    df <- data.frame(
        logFC = c(2, NA, 1),
        row.names = c("100", "1000", "10000")
    )

    out <- suppressMessages(prepareQuery(df))

    expect_equal(out$gene, c("100", "10000"))
})

test_that("prepareQuery drops NA-padj genes by default with a message", {
    df <- data.frame(
        logFC = c(2, -1, 0.5, 1.2),
        padj = c(0.001, NA, 0.02, NA),
        row.names = c("100", "1000", "10000", "9999")
    )

    expect_message(
        out <- prepareQuery(df),
        "Dropped 2 genes with a missing adjusted p-value"
    )
    expect_equal(out$gene, c("100", "10000"))
    expect_false(anyNA(out$padj))
})

test_that("prepareQuery keeps NA-padj genes when dropLowExpressed = FALSE", {
    df <- data.frame(
        logFC = c(2, -1, 0.5),
        padj = c(0.001, NA, 0.02),
        row.names = c("100", "1000", "10000")
    )

    out <- suppressMessages(prepareQuery(df, dropLowExpressed = FALSE))
    expect_equal(out$gene, c("100", "1000", "10000"))
    expect_true(anyNA(out$padj))
})

test_that("prepareQuery reports reference overlap when a reference is given", {
    ref <- SummarizedExperiment(
        assays = list(m = matrix(
            0, nrow = 2, ncol = 1,
            dimnames = list(c("100", "1000"), "s1")
        ))
    )

    expect_message(
        prepareQuery(make_entrez_query(), reference = ref),
        "present in the reference database"
    )
})

# prepareQuery (real ID mapping)

test_that("prepareQuery maps gene symbols to Entrez via an OrgDb", {
    skip_if_not_installed("org.Hs.eg.db")

    df <- data.frame(logFC = c(2, -1), row.names = c("TP53", "BRCA1"))

    out <- suppressMessages(prepareQuery(df, keytype = "SYMBOL"))

    expect_equal(out$gene, c("7157", "672"))
})

test_that("prepareQuery strips Ensembl version suffixes before mapping", {
    skip_if_not_installed("org.Hs.eg.db")

    # ENSG00000141510 is TP53 (Entrez 7157); the ".11" version suffix must
    # be dropped for the lookup to succeed.
    df <- data.frame(logFC = 2, row.names = "ENSG00000141510.11")

    out <- suppressMessages(prepareQuery(df, keytype = "ENSEMBL"))

    expect_equal(out$gene, "7157")
})

test_that("prepareQuery warns rather than errors when nothing maps", {
    skip_if_not_installed("org.Hs.eg.db")

    df <- data.frame(
        logFC = c(1, 2),
        row.names = c("NOTAGENE1", "NOTAGENE2")
    )

    expect_warning(
        suppressMessages(prepareQuery(df, keytype = "SYMBOL")),
        "Only 0% of query genes mapped"
    )
})

test_that("prepareQuery errors on an invalid keytype", {
    skip_if_not_installed("org.Hs.eg.db")

    df <- data.frame(logFC = c(1, 2), row.names = c("TP53", "BRCA1"))

    expect_error(
        prepareQuery(df, keytype = "NONSENSE"),
        "not a valid keytype"
    )
})
