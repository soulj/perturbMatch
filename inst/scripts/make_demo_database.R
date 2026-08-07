# make_demo_database.R
#
# Builds the small bundled demonstration reference database shipped in
# inst/extdata/demo_database.
#
# Source data
#
# The demo is a subset of the published full limma (batch-uncorrected)
# perturbMatch reference, taken from the Zenodo deposit that also backs
# getReferenceDatabase():
#
#   https://doi.org/10.5281/zenodo.21792179  (files limma_assays.h5,
#   limma_colData.rds, perturbMatch_rowData.rds)
#
#
# The reference signatures are derived from public RNA-seq series in NCBI GEO
# (the sig_id of every signature carries its GSE accession) and are released
# under CC-BY-4.0, doi:10.5281/zenodo.21792179.
#
# Subsetting
#
# Signatures are chosen so that the demo contains informative hits for BOTH
# vignette queries:
#
#   * the bulk airway / dexamethasone example (expected regulator NR3C1), and
#   * the single-cell Kang PBMC interferon example (antiviral / IFN regulators
#     such as STING1, IFI16, IRF1).
#
# Genes (rows) are reduced to keep the object within Bioconductor data-size
# limits.

suppressPackageStartupMessages({
    library(perturbMatch)
    library(SummarizedExperiment)
    library(HDF5Array)
    library(airway)
    library(DESeq2)
    library(muscData)
    library(SingleCellExperiment)
    library(scuttle)
})

set.seed(1)

OUT_DIR <- "inst/extdata/demo_database"
TARGET_ROWS <- 6500L

full <- getReferenceDatabase(
    type = "limma",
    batchCorrection = FALSE,
    source = "zenodo"
)

## signatures to keep
selected_sigs <- c(
    # bulk airway / dexamethasone hits
    "GSE200312_NR3C1_KD_1", # NR3C1 -- expected dexamethasone regulator
    "GSE197354_AFF4_KD_1",
    "GSE125167_ZBTB16_KD_1",
    "GSE179474_SMAD4_KD_1",
    "GSE197118_TAPT1_KO_1",
    "GSE227512_LINC01638_KD_1",
    "GSE215822_CCDC88A_KO_2",
    "GSE263122_VGLL4_OE_1",
    "GSE87831_NUP153_KD_2",
    "GSE194080_TBX5_KD_1",
    "GSE165963_ESRRG_KO_1",
    "GSE213710_FOSL2_OE_1",
    # single-cell / interferon hits
    "GSE129436_STING1_KO_1", # STING1 KO -- top IFN hit / singscore target
    "GSE165910_STING1_OE_1", # STING1 OE -- opposite-direction companion
    "GSE163705_IFI16_KD_1",
    "GSE114284_IRF1_KO_1",
    "GSE137155_STAT2_KD_1",
    "GSE182754_ICP0_OE_1",
    "GSE182754_E4ORF3_OE_1",
    "GSE202869_ORF9b_OE_1",
    "GSE202869_S_OE_1"
)

missing_sigs <- setdiff(selected_sigs, colnames(full))
if (length(missing_sigs)) {
    stop(
        "Signatures not found in full reference: ",
        paste(missing_sigs, collapse = ", ")
    )
}
sub <- full[, selected_sigs]

# Bulk airway query
data(airway)
dds <- DESeqDataSet(airway, design = ~ cell + dex) |>
    DESeq() |>
    results(contrast = c("dex", "trt", "untrt"))
q_airway <- prepareQuery(dds, reference = full)

# Single-cell Kang PBMC pseudobulk query
sce <- Kang18_8vs8()
sce$donor <- paste0("Donor", sce$ind)
sce <- sce[, sce$multiplets == "singlet"]
sce <- sce[, !is.na(sce$cell)]
sce <- sce[, sce$cell == "CD14+ Monocytes"]
pb <- aggregateAcrossCells(
    sce,
    ids = DataFrame(donor = sce$donor, stim = sce$stim)
)
dds_sc <- DESeqDataSet(pb, design = ~ donor + stim) |> DESeq()
res_sc <- results(dds_sc, contrast = c("stim", "stim", "ctrl"))
q_kang <- prepareQuery(res_sc, reference = full)

sig_airway <- q_airway$gene[!is.na(q_airway$padj) & q_airway$padj < 0.05]
sig_kang <- q_kang$gene[!is.na(q_kang$padj) & q_kang$padj < 0.05]

# perturbed target genes
targets <- as.character(colData(sub)$entrez_gene_id)
targets <- targets[!is.na(targets)]

must_keep <- intersect(
    unique(c(sig_airway, sig_kang, targets)),
    rownames(sub)
)

## pad with the most variable remaining genes
pool <- setdiff(rownames(sub), must_keep)
if (length(must_keep) < TARGET_ROWS && length(pool) > 0) {
    pi_mat <- as.matrix(assay(sub, "pi")[pool, , drop = FALSE])
    vars <- MatrixGenerics::rowVars(pi_mat, na.rm = TRUE)
    vars[is.na(vars)] <- 0
    n_pad <- min(TARGET_ROWS - length(must_keep), length(pool))
    pad <- pool[order(vars, decreasing = TRUE)][seq_len(n_pad)]
    keep <- c(must_keep, pad)
} else {
    keep <- must_keep
}
keep <- rownames(sub)[rownames(sub) %in% keep]
sub <- sub[keep, ]

## strip rowData to the essentials
rowData(sub) <- rowData(sub)[, c("GeneID", "Symbol")]

metadata(sub)$demo <- TRUE
metadata(sub)$source <- "bundled"

metadata(sub)$dataset <- "limma"
metadata(sub)$assay_names <- assayNames(sub)

assays(sub) <- lapply(assays(sub), as.matrix)

setHDF5DumpCompressionLevel(9L)
if (dir.exists(OUT_DIR)) {
    unlink(OUT_DIR, recursive = TRUE)
}
saveHDF5SummarizedExperiment(sub, OUT_DIR, replace = TRUE, level = 9L)

cat("Demo database written to", OUT_DIR, "\n")
cat("dim:", nrow(sub), "genes x", ncol(sub), "signatures\n")
cat("must_keep genes:", length(must_keep), " total kept:", nrow(sub), "\n")
print(as.data.frame(colData(sub)[, c(
    "pert_symbol",
    "pert_type",
    "entrez_gene_id",
    "target_gene_percentile"
)]))
