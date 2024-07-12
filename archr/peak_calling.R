# Usage:
# Rscript --vanilla peak_calling.R
#   -p "data/atac/archr_proj"
#   -r "data/atac/annotations.tsv"
#   -b "annotation"  # not used if -r is provided
#   -q 2
#   -l 5
#   -k 1000
#   -m "/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/bin/macs2"
#   -g "hg38"
#   -t 1
#   -s 1234


# Define options and arguments for the script.
suppressMessages(library(optparse))
option_list <- list(
  make_option(c("-p", "--proj_dir"), type="character", default=NULL, help="Path to ArchR project directory."),
  make_option(c("-r", "--annotations_path"), type="character", default=NULL, help="Two column tsv file with cell barcodes and annotations. If provided, will subset to cells with annotations and use annotations for pseudobulking. Default is NULL."),
  make_option(c("-b", "--group_by"), type="character", default="annotation", help="Group by for pseudobulk, ignored if annotations_path is provided. Default is annotation."),
  make_option(c("-q", "--min_replicates"), type="integer", default=2, help="Minimum number of pseudobulk replicates for each group to create. Default is 2."),
  make_option(c("-l", "--max_replicates"), type="integer", default=5, help="Maximum number of pseudobulk replicates for each group to create. Default is 5."),
  make_option(c("-k", "--max_cells"), type="integer", default=1000, help="Maximum number of cells per pseudobulk. Default is 1000."),
  make_option(c("-m", "--macs2_path"), type="character", default="/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/bin/macs2", help="Path to macs2 executable."),
  make_option(c("-g", "--genome"), type="character", default="hg38", help="Genome to use for ArchR project. Default is hg38."),
  make_option(c("-t", "--threads"), type="integer", default=1, help="Number of threads to use for ArchR project. Default is 1."),
  make_option(c("-s", "--seed"), type="integer", default=1234, help="Seed to use for ArchR project. Default is 1234.")
)
parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options

# Print arguments
proj_dir <- opt$proj_dir
group_by <- opt$group_by
annotations_path <- opt$annotations_path
min_replicates <- opt$min_replicates
max_replicates <- opt$max_replicates
max_cells <- opt$max_cells
macs2_path <- opt$macs2_path
genome <- opt$genome
threads <- opt$threads
seed <- opt$seed
print(paste0("ArchR project path = ", proj_dir))
print(paste0("Group by = ", group_by))
print(paste0("Annotations path = ", annotations_path))
print(paste0("Minimum replicates = ", min_replicates))
print(paste0("Maximum replicates = ", max_replicates))
print(paste0("Maximum cells = ", max_cells))
print(paste0("MACS2 path = ", macs2_path))
print(paste0("Genome = ", genome))
print(paste0("Threads = ", threads))
print(paste0("Seed = ", seed))

# Import libraries
print("Importing libraries\n")
suppressMessages(library(Seurat))
suppressMessages(library(ArchR))
suppressMessages(library(parallel))
suppressMessages(library(tidyverse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Matrix))

# Set-up
addArchRThreads(threads = threads)
set.seed(seed)
setwd(proj_dir)

# Add annotation
addArchRGenome(genome)

# Load the ArchR project
print("Loading ArchR project\n")
proj = loadArchRProject(path = "./")

# Read in annotations if any exist
if (is.null(annotations_path)) {
    if (!(group_by %in% colnames(proj@cellColData))) {
        stop(sprintf("Group by column %s not found in cell metadata.", group_by))
    }
} else {
    print("Subsetting to cells with annotations. These are expected to have matching barcodes in the first column and annotation in the second column.")
    annotations = read.csv(annotations_path, row.names = 1, sep = "\t", header = FALSE)
    annotations = as.data.frame(annotations)
    cellids = rownames(annotations)
    matched_ids = intersect(cellids, rownames(proj@cellColData))
    idxSample <- BiocGenerics::which(proj$cellNames %in% matched_ids)
    proj <- proj[idxSample,]
    annotations = annotations[proj$cellNames, ]
    proj$annotation <- annotations
    group_by = "annotation"
    print(paste0("Number of cells in subset: ", length(proj$cellNames)))
    print(paste0("Using annotations for pseudobulking. Group by = ", group_by))
}

# ArchR psuedobulks
print("Running ArchR pseudobulking\n")
proj <- addGroupCoverages(
    ArchRProj = proj, 
    groupBy = group_by,
    minReplicates = min_replicates,
    maxReplicates = max_replicates,
    maxCells = max_cells,
    force = TRUE,
)

# Run peak calling with MACS2
print("Running MACS2 peak calling with ArchR\n")
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = group_by,
    pathToMacs2 = macs2_path,
    force = TRUE,
)

# Export bed files for the peak calls
print(sprintf("Exporting bed files for the peak calls to %s/Peaks/SplitPeaks\n", proj_dir))
peaks_dir <- file.path(proj_dir, "PeakCalls", "SplitPeaks")
if (!dir.exists(peaks_dir)) {
    dir.create(peaks_dir)
}
peakset <- getPeakSet(proj)
export.bed(peakset, con=file.path(peaks_dir, "consensus_peaks.bed"))
peakset_list <- split(peakset, names(peakset))
for (i in 1:length(peakset_list)) {
    gr <- peakset_list[[i]]
    cell_type <- names(peakset_list)[i]
    export.bed(gr, con=file.path(peaks_dir, paste0(cell_type, ".bed")))
}

# Add in the peak matrix
print("Adding consensus peak matrix to ArchR project\n")
proj <- addPeakMatrix(proj, force=TRUE)

# Save peak matrix as mtx.mtx, features.tsv, and barcodes.tsv
print(sprintf("Saving peak matrix as mtx.mtx, features.tsv, and barcodes.tsv to %s/Matrices/PeakMatrix\n", proj_dir))
se <- getMatrixFromProject(
    proj, 
    useMatrix="PeakMatrix",
)
mtx_dir <- file.path(proj_dir, "Matrices", "PeakMatrix")
if (!dir.exists(mtx_dir)) {
    dir.create(mtx_dir, recursive=T)
}
peaks <- rowRanges(se)
regions <- paste0(seqnames(peaks), ":", start(peaks), "-", end(peaks))
write.table(regions, file.path(mtx_dir, "features.tsv"), row.names=F, col.names=F, quote=F)
coldata <- colData(se)
cells <- colnames(se)
write.table(cells, file.path(mtx_dir, "barcodes.tsv"), row.names=F, col.names=F, quote=F)
mtx <- assays(se)$PeakMatrix
writeMM(mtx, file.path(mtx_dir, "mtx.mtx"))

# Save per group bigwigs
print("Saving per group bigwigs\n")
getGroupBW(
    ArchRProj = proj,
    groupBy = group_by,
    normMethod = "ReadsInTSS",
    tileSize = 100,
    maxCells = 1000,
    ceiling = 4,
    verbose = TRUE,
)

# Save object with new stuff added
print("Saving ArchR project\n")
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "./",
)
