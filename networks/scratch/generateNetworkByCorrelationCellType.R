# -*- coding: utf-8 -*-
# ## -----------------------------------------------------------------------------
# title: Create Celltype Specific Networks for Single Islet Samples 
# Authorship: Adam Klie, 11/30/2021
# Description: Notebook to generate cell type specific links between CREs and genes
# Usage: Rscript --vanilla multiomeGenerateSingleSampleNetwork.R $sample
# Note: You must have the file directory structure noted in the README.txt 
# in order to successfully run this script
# ## -----------------------------------------------------------------------------

# # -----------------------------------------------------------------------------
# Set-up
# # -----------------------------------------------------------------------------

# Grab args and sample name
args <- commandArgs(trailingOnly=TRUE)
sample <- args[1]  # Sample identifier
print(sprintf("Creating CRE-gene links for celltypes in sample %s", sample))
wd <- sprintf(sample)

# Set-up reticulate for running Python in R
print("Setting up environment with proper libraries")
Sys.setenv(RETICULATE_PYTHON="/cellar/users/aklie/opt/miniconda3/envs/Renv/bin/python")
library(reticulate)
reticulate::use_python("/cellar/users/aklie/opt/miniconda3/envs/Renv/bin/python")
reticulate::use_condaenv("/cellar/users/aklie/opt/miniconda3/envs/Renv")
reticulate::py_module_available(module='leidenalg') #needs to be TRUE
reticulate::import('leidenalg') #good to make sure this doesn't error

# Load libraries
suppressMessages(library(hdf5r))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratDisk))
suppressMessages(library(SeuratData))
suppressMessages(library(Signac))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(harmony))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))
warnLevel <- getOption('warn')
options(warn = -1)
set.seed(1234)

# Set up multithreading
suppressMessages(library(future))
plan("multicore", workers = 16)  # Change based on availability, same for below
options(future.globals.maxSize = 128 * 1024 ^ 3)

# Load in genomic annotations from the interwebs. Takes about 5 mins
print("Loading hg38 genome annotations from Ensembl")
annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- 'hg38'

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# # -----------------------------------------------------------------------------
# Load Sample
# # -----------------------------------------------------------------------------

# Load in Single 10X Multiome Sample
wd <- sprintf(sample)  # Working directory
frag.file <- file.path(wd, 'atac_fragments.tsv.gz')  # Fragment path needed for chromatin object

# Print to console
print(sprintf("Loading from %s", file.path(wd, 'raw_feature_bc_matrix.h5')))
print(sprintf("Frag file in %s", frag.file))

# Load in the h5 raw Seurat object. Takes about a minute
inputdata.10x <- Read10X_h5(file.path(wd, 'raw_feature_bc_matrix.h5'))

# Grab the counts for gene expression and atac
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$`Peaks`

# # -----------------------------------------------------------------------------
# QC and Preprocess Gene Expression (GEX) and ATAC Peaks
# # -----------------------------------------------------------------------------

# Create an object for the rna data
print("Creating Seurat RNA object and filtering cells based on standard metrics")
adata <- CreateSeuratObject(counts=rna_counts, project = sample)
print(paste(length(colnames(adata[["RNA"]])), "total BCs"))

# Calculate mitochondrial count percentage and fix 0 count barcodes
adata[['percent.mt']] <- PercentageFeatureSet(adata, pattern = '^MT-')
adata$percent.mt[is.nan(adata$percent.mt)] <- 0

# Filter cells based on rna metrics
adata_sub <- subset(
  x = adata,
  subset = (nFeature_RNA >= 500) & (nFeature_RNA <= 5000) & (percent.mt <= 5)
)
print(paste(length(colnames(adata_sub[["RNA"]])), "total BCs after RNA filtering"))

# Add per barcode metric metadata to object. About 30s
print("Removing putative multiplets called by cellranger")
qc <- read.table(file.path(wd, 'per_barcode_metrics.csv'), sep=',', header=TRUE, stringsAsFactors=1)
qc <- as.data.frame(qc)
rownames(qc) <- qc$gex_barcode
qc <- qc[Cells(adata_sub), 1:length(colnames(qc))]
adata_sub <- AddMetaData(adata_sub, qc)
gc()

# Multiplet count
print(paste(table(adata_sub$excluded_reason)[2][[1]], "putative multiplets"))

# Removal of multiplets
adata_sub_multiplet <- subset(
  x = adata_sub,
  subset = excluded_reason != 1 
)
print(paste(length(colnames(adata_sub_multiplet[["RNA"]])), "total BCs after multiplet removal"))

# Initial ATAC filter
print("Performing an intial ATAC filter based on fragment numbers")
adata_sub_atac <- subset(
  x = adata_sub_multiplet,
  subset = (atac_fragments >= 1000) & (atac_fragments <= 60000) & (atac_mitochondrial_reads <= 5000)
)
print(paste(length(colnames(adata_sub_atac[["RNA"]])), "total BCs after first ATAC thresholds"))

# Clean up the ATAC matrix
print("Generating a ChromatinAssay object and adding to RNA object")
atac_counts <- atac_counts[,colnames(adata_sub_atac)]  # Grab only the left-over cells from the preprocessing
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(':', '-'))  # Fix for annotations
genome(grange.counts) <- "hg38"
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)  # Keep only peaks within chromosomes
atac_counts <- atac_counts[as.vector(grange.use), ]
print(paste(table(grange.use)[1][[1]], "ranges removed not on chroms"))

# Create a chromatin assay object, add annotations and add to adata object as 'ATAC' assay. Takes about 30s
chrom_assay <- CreateChromatinAssay(counts=atac_counts, sep=c(':', '-'), genome='hg38', fragments=frag.file, min.cells=0, min.features=0)
Annotation(chrom_assay) <- annotations
adata_sub_atac[['ATAC']] <- chrom_assay

# Second ATAC filter
print("Performing a second ATAC filter based on nucleosome signal and TSS enrichment. Takes a about 7mins")
DefaultAssay(adata_sub_atac) <- "ATAC"
adata_sub_atac <-  NucleosomeSignal(object = adata_sub_atac)
adata_sub_atac <- TSSEnrichment(object = adata_sub_atac, fast = FALSE)
adata_sub_atac <- subset(
  x = adata_sub_atac,
  subset = nucleosome_signal < 2 & TSS.enrichment > 1
)
print(paste(length(colnames(adata_sub_atac[["RNA"]])),"total BCs after nucelosome and tss thresholds"))

# # -----------------------------------------------------------------------------
# Data Transformations and Embeddings
# # -----------------------------------------------------------------------------

# Variable switch
adata <- adata_sub_atac

# RNA data transformation and PCA
print("Performing SCTransform and PCA on RNA data while regressing out coverage, mitochondrial percentage and cell cycle")
DefaultAssay(adata) <- 'RNA'
adata <- NormalizeData(adata)
adata <- CellCycleScoring(adata, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
adata <- SCTransform(adata, vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'), verbose = FALSE)
adata <- RunPCA(adata)

# # -----------------------------------------------------------------------------
# Cell Type Labeling from Reference
# # -----------------------------------------------------------------------------

print("Using the Panc8 reference transcriptome to transfer cell type labels to query data")

# Load in the reference Seurat object previously generated
reference <- LoadH5Seurat("../reference/panc8.sctransform.smartseq2-reference-integrated.h5seurat")

# Change the default assays and make sure features are the same
DefaultAssay(reference) <- "SCT"
DefaultAssay(adata) <- "SCT"
VariableFeatures(reference) <- VariableFeatures(adata)

# Transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = adata,
  normalization.method = "SCT",
  reference.reduction = "pca",
  recompute.residuals = FALSE,
  dims = 1:50
)

# Make predictions for cell type based on pca
predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype,
  weight.reduction = adata[['pca']],
  dims = 1:50
)

# Add those predictins to the metadata
adata <- AddMetaData(
  object = adata,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Idents(adata) <- "predicted.id"

# Check the cell type labeling and save the result
table(adata$predicted.id)
num_celltypes <- data.frame(unclass(table(adata$predicted.id)))
colnames(num_celltypes) <- "number_of_cells"
print(sprintf("Saving celltypes to %s", file.path("indv_sample_networks", sprintf("%s.indv.num.celltypes.csv",sample))))
write.csv(num_celltypes, file.path("indv_sample_networks", sprintf("%s.indv.num.celltypes.csv",sample)), row.names = TRUE)

# # -----------------------------------------------------------------------------
# Add MACS Peaks for each cell type
# # -----------------------------------------------------------------------------

# call peaks using MACS2
DefaultAssay(adata) <- "ATAC"
peaks <- CallPeaks(
    object = adata, 
    macs2.path = "/cellar/users/aklie/opt/miniconda3/envs/seq_tools_dev/bin/macs2", 
    group.by = "predicted.id",
    combine.peaks = FALSE,
    outdir = file.path(wd, "macs_peaks"),
    fragment.tempdir = file.path(wd, "macs_peaks"),
    cleanup = FALSE)


alpha.peaks <- NULL
beta.peaks <- NULL
delta.peaks <- NULL
for (gr in peaks) {
    print(gr$ident[1])
    if (is.null(gr$ident[1])) {
        next
    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    curr.peaks <- keepStandardChromosomes(gr, pruning.mode = "coarse")
    curr.peaks <- subsetByOverlaps(x = curr.peaks, ranges = blacklist_hg38_unified, invert = TRUE)
    if (gr$ident[1] == "alpha") {
        alpha.peaks <- curr.peaks
    } else if (gr$ident[1] == "beta") {
        beta.peaks <- curr.peaks
    } else if (gr$ident[1] == "delta") {
        delta.peaks <- curr.peaks
    }
}

# Alpha cells
alpha.ids <- colnames(adata)[adata$predicted.id == "alpha"]

# quantify counts in each peak: https://satijalab.org/signac/reference/featurematrix. Process_n increases speed with trade-off for memory
alpha.counts <- FeatureMatrix(
    fragments = Fragments(adata),
    features = alpha.peaks,
    cells = alpha.ids,
    process_n = 2000000
)

# Subset
alpha.cells <- subset(x = adata, idents = "alpha")

# create a new assay using the MACS2 peak set and add it to the Seurat object
alpha.cells[["peaks"]] <- CreateChromatinAssay(
    counts = alpha.counts,
    fragments = frag.file,
    annotation = annotations,
    ranges = alpha.peaks,
    genome = "hg38"
)

# Beta cells
beta.ids <- colnames(adata)[adata$predicted.id == "beta"]

# quantify counts in each peak: https://satijalab.org/signac/reference/featurematrix. Process_n increases speed with trade-off for memory
beta.counts <- FeatureMatrix(
    fragments = Fragments(adata),
    features = beta.peaks,
    cells = beta.ids,
    process_n = 2000000
)
# Subset
beta.cells <- subset(x = adata, idents = "beta")

# create a new assay using the MACS2 peak set and add it to the Seurat object
beta.cells[["peaks"]] <- CreateChromatinAssay(
    counts = beta.counts,
    fragments = frag.file,
    annotation = annotations,
    ranges = beta.peaks,
    genome = "hg38"
)

# Delta cells
delta.ids <- colnames(adata)[adata$predicted.id == "delta"]

# quantify counts in each peak: https://satijalab.org/signac/reference/featurematrix. Process_n increases speed with trade-off for memory
delta.counts <- FeatureMatrix(
    fragments = Fragments(adata),
    features = delta.peaks,
    cells = delta.ids,
    process_n = 2000000
)

# Subset
delta.cells <- subset(x = adata, idents = "delta")

# create a new assay using the MACS2 peak set and add it to the Seurat object
delta.cells[["peaks"]] <- CreateChromatinAssay(
    counts = delta.counts,
    fragments = frag.file,
    annotation = annotations,
    ranges = delta.peaks,
    genome = "hg38"
)

# # -----------------------------------------------------------------------------
# Link Peaks to Genes
# # -----------------------------------------------------------------------------

print("Calculating correlations between normalized expression and MACS peak counts by celltype")

# Alpha cells
DefaultAssay(alpha.cells) <- "peaks"
alpha.cells <- RegionStats(alpha.cells, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
alpha.cells <- LinkPeaks(
    object = alpha.cells,
    peak.assay = "peaks",
    expression.assay = "SCT",
    method="pearson",
    pvalue_cutoff = 2,
    score_cutoff = 0
)

# Save links
gr.links <- Links(alpha.cells)
linked.peaks <- gr.links$peak
chr.start.end <- unlist(strsplit(linked.peaks, "-"))
chr.start.end <- matrix(chr.start.end, length(chr.start.end)/3, 3, byrow = T)
seqs <- chr.start.end[, 1]
starts <- chr.start.end[, 2]
ends <- chr.start.end[, 3]
links.df <- data.frame(
    peak=gr.links$peak,
    gene=gr.links$gene,
    score=score(gr.links),
    zscore=gr.links$zscore,
    pvalue=gr.links$pvalue)
gr.peaks <- granges(alpha.cells)
gr.peaks$peak <- paste0(seqnames(gr.peaks), "-", start(gr.peaks), "-", end(gr.peaks))
annnotated.links <- merge(gr.peaks, links.df, by="peak", suffixes=c("_peaks", "_links"))
write.table(annnotated.links, 
            file=file.path("indv_sample_networks", "correlation_links_celltype", sprintf("%s.alpha.links.tsv", sample)),
            quote=F, sep="\t", row.names=F, col.names=T)

# Beta cells
DefaultAssay(beta.cells) <- "peaks"
beta.cells <- RegionStats(beta.cells, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
beta.cells <- LinkPeaks(
    object = beta.cells,
    peak.assay = "peaks",
    expression.assay = "SCT",
    method="pearson",
    pvalue_cutoff = 2,
    score_cutoff = 0
)

# Save links
gr.links <- Links(beta.cells)
linked.peaks <- gr.links$peak
chr.start.end <- unlist(strsplit(linked.peaks, "-"))
chr.start.end <- matrix(chr.start.end, length(chr.start.end)/3, 3, byrow = T)
seqs <- chr.start.end[, 1]
starts <- chr.start.end[, 2]
ends <- chr.start.end[, 3]
links.df <- data.frame(
    peak=gr.links$peak,
    gene=gr.links$gene,
    score=score(gr.links),
    zscore=gr.links$zscore,
    pvalue=gr.links$pvalue)
gr.peaks <- granges(beta.cells)
gr.peaks$peak <- paste0(seqnames(gr.peaks), "-", start(gr.peaks), "-", end(gr.peaks))
annnotated.links <- merge(gr.peaks, links.df, by="peak", suffixes=c("_peaks", "_links"))
write.table(annnotated.links, 
            file=file.path("indv_sample_networks", "correlation_links_celltype", sprintf("%s.beta.links.tsv", sample)),
            quote=F, sep="\t", row.names=F, col.names=T)

# Delta cells
DefaultAssay(delta.cells) <- "peaks"
delta.cells <- RegionStats(delta.cells, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
delta.cells <- LinkPeaks(
    object = delta.cells,
    peak.assay = "peaks",
    expression.assay = "SCT",
    method="pearson",
    pvalue_cutoff = 2,
    score_cutoff = 0
)

# Save links
gr.links <- Links(delta.cells)
linked.peaks <- gr.links$peak
chr.start.end <- unlist(strsplit(linked.peaks, "-"))
chr.start.end <- matrix(chr.start.end, length(chr.start.end)/3, 3, byrow = T)
seqs <- chr.start.end[, 1]
starts <- chr.start.end[, 2]
ends <- chr.start.end[, 3]
links.df <- data.frame(
    peak=gr.links$peak,
    gene=gr.links$gene,
    score=score(gr.links),
    zscore=gr.links$zscore,
    pvalue=gr.links$pvalue)
gr.peaks <- granges(delta.cells)
gr.peaks$peak <- paste0(seqnames(gr.peaks), "-", start(gr.peaks), "-", end(gr.peaks))
annnotated.links <- merge(gr.peaks, links.df, by="peak", suffixes=c("_peaks", "_links"))
write.table(annnotated.links, 
            file=file.path("indv_sample_networks", "correlation_links_celltype", sprintf("%s.delta.links.tsv", sample)),
            quote=F, sep="\t", row.names=F, col.names=T)

# # -----------------------------------------------------------------------------
# Save Analyzed Objects
# # -----------------------------------------------------------------------------

# All cells
saveRDS(alpha.cells, file = file.path("indv_sample_networks", "correlation_links_celltype", sprintf("%s.alpha.linked.rds", sample)))
saveRDS(beta.cells, file = file.path("indv_sample_networks", "correlation_links_celltype", sprintf("%s.beta.linked.rds", sample)))
saveRDS(delta.cells, file = file.path("indv_sample_networks", "correlation_links_celltype", sprintf("%s.delta.linked.rds", sample)))