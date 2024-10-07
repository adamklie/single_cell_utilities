# -*- coding: utf-8 -*-
# ## -----------------------------------------------------------------------------
# title: Create Cell-type Specific Networks for Single Islet Samples
# Authorship: Adam Klie, 12/11/2021
# Description: Notebook to generate cell type specific links between CREs and genes
# Usage: Rscript --vanilla generatePeakToGeneLinks.R $sample
# Note: You must have the file directory structure noted in the README.txt 
# in order to successfully run this script
# ## -----------------------------------------------------------------------------

# # -----------------------------------------------------------------------------
# Set-up
# # -----------------------------------------------------------------------------

# Grab args and sample name
args <- commandArgs(trailingOnly=TRUE)
sample <- args[1]  # Sample identifier
cell.type <- args[2]
print(sprintf("Creating CRE-gene links for %s cells in sample %s", cell.type, sample))

# Load libraries, takes between 20s and 1 min
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(harmony))
suppressMessages(library(SeuratDisk))
suppressMessages(library(ggplot2))
options(future.globals.maxSize = 4000 * 1024^2)
warnLevel <- getOption('warn')
options(warn = -1)
set.seed(1234)

# For parallelization of FeatureMatrix command
suppressMessages(library(future))
plan("multicore", workers = 16)  # Change based on availability, same for below

# Load in genomic annotations from the interwebs. Takes about 5 mins
print("Loading hg38 genome annotations from Ensembl")
annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- 'hg38'

# Load in the Seurat object
path = file.path("../indv_sample_analysis", sprintf("%s.indv.analysis.h5seurat", sample))
print(sprintf("Loading from %s", path))
adata <- LoadH5Seurat(path, assays = c("SCT"), reductions = c("pca"), graphs = FALSE, neighbors = FALSE, verbose = FALSE)
table(adata$predicted.id)

# Grab cell-type ids
ids <- colnames(adata)[adata$predicted.id == cell.type]

# Create a granges object from bed file
file = file.path(sprintf("../macs_peaks_merged/%s_peaks.narrowPeak", cell.type))
print(sprintf("Loading peaks from %s", file))
df <- read.table(file, header = FALSE)

colnames(df) <- c("chr", "start", "end", "name", "score", 
                  "strand", "fold_change", "neg_log10pvalue_summit",
                  "neg_log10qvalue_summit", "relative_summit_position")

gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
gr <- keepStandardChromosomes(gr , pruning.mode = "coarse")
gr <- subsetByOverlaps(x = gr, ranges = blacklist_hg38_unified, invert = TRUE)
gr

# # -----------------------------------------------------------------------------
# Find counts in peaks and create ChromatinAssay object
# # -----------------------------------------------------------------------------

# First create a fragments object
fragments <- CreateFragmentObject(
    path = file.path(sprintf("../%s", sample), "atac_fragments.tsv.gz"),
    cells = ids,
    validate.fragments = TRUE)

# Get the counts in fragments and store as a sparse matrix
counts <- FeatureMatrix(
    fragments = fragments,
    features = gr,
    cells = ids,
    process_n = 10000000
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
cells <- subset(x = adata, idents = cell.type)
cells[["peaks"]] <- CreateChromatinAssay(
    counts = counts,
    fragments = fragments,
    annotation = annotations,
    ranges = gr,
    genome = "hg38"
)

# # -----------------------------------------------------------------------------
# Link peaks to genes for this cell-type
# # -----------------------------------------------------------------------------

# Cells
DefaultAssay(cells) <- "peaks"
cells <- RegionStats(cells, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes. Adjust cut-offs so that all gene peak links are kept and non are discarded
cells <- LinkPeaks(
    object = cells,
    peak.assay = "peaks",
    expression.assay = "SCT",
    method="pearson",
    pvalue_cutoff = 2,
    score_cutoff = 0
)

# Save links to file
gr.links <- Links(cells)
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
gr.peaks <- granges(cells)
gr.peaks$peak <- paste0(seqnames(gr.peaks), "-", start(gr.peaks), "-", end(gr.peaks))
annnotated.links <- merge(gr.peaks, links.df, by="peak", suffixes=c("_peaks", "_links"))
write.table(annnotated.links, 
            file=file.path("../indv_sample_networks", "correlation_links_shared-peaks", sprintf("%s.%s.links.tsv", sample, cell.type)),
            quote=F, sep="\t", row.names=F, col.names=T)

# Save object to file
saveRDS(cells, file = file.path("../indv_sample_networks", "correlation_links_shared-peaks", 
                                      sprintf("%s.%s.linked.rds", sample, cell.type)))
