# Define options and arguments for the script.
suppressMessages(library(optparse))
option_list <- list(
  make_option(c("-i", "--input_rds"), type="character", default=NULL, help="Input RDS file."),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Output directory for ArchR project."),
  make_option(c("-g", "--group_by"), type="character", default=NULL, help="Group by variable."),
  make_option(c("-a", "--assay"), type="character", default="peaks", help="Assay to use for peak calling."),
  make_option(c("-p", "--macs2_path"), type="character", default="/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/bin/macs2", help="Path to MACS2."),
  make_option(c("-s", "--save_frags"), type="logical", default=TRUE, help="Whether to keep the pseudobulk fragments, if FALSE, the fragments will be deleted after peak calling."),
  make_option(c("-c", "--process_n"), type="integer", default=1, help="Number of processes to use for peak calling."),
  make_option(c("-t", "--store_assay"), type="character", default=NULL, help="Store the feature matrix as a new assay in the Seurat object.")
)
parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options

# Print arguments
input_rds <- opt$input_rds
output_dir <- opt$output_dir
group_by <- opt$group_by
assay <- opt$assay
macs2_path <- opt$macs2_path
save_frags <- opt$save_frags
process_n <- opt$process_n
store_assay <- opt$store_assay
cat(sprintf("Input RDS: %s\n", input_rds))
cat(sprintf("Output directory: %s\n", output_dir))
cat(sprintf("Group by: %s\n", group_by))
cat(sprintf("Assay: %s\n", assay))
cat(sprintf("MACS2 path: %s\n", macs2_path))
cat(sprintf("Save fragments: %s\n", save_frags))
cat(sprintf("Number of processes: %s\n", process_n))
cat(sprintf("Store assay: %s\n", store_assay))

# This script will not keep the intermediate fragment files!!!
cat("Importing libraries\n")
suppressMessages(library(SeuratDisk))
suppressMessages(library(SeuratData))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(Matrix))

# Read the RDS file
adata <- readRDS(input_rds)

# Call peaks, save fragments if save_frags is not NULL
cat("Calling peaks\n")
if (save_frags) {
  frag_dir <- file.path(output_dir, "fragments")
  if (!dir.exists(frag_dir)) {
    dir.create(frag_dir)
  }
} else {
  frag_dir <- NULL
}
peaks_dir <- file.path(output_dir, "peaks")
if (!dir.exists(peaks_dir)) {
  dir.create(peaks_dir)
}
peaks <- CallPeaks(
  adata,
  assay = assay,
  group.by = group_by,
  idents = NULL,
  macs2.path = macs2_path,
  broad = FALSE,
  format = "BED",
  outdir = peaks_dir,
  fragment.tempdir = frag_dir,
  combine.peaks = TRUE,
  effective.genome.size = 2.7e+09,
  extsize = 200,
  shift = -100,
  cleanup = FALSE,
  verbose = TRUE
)
cat("\n")

# Quantify counts in each peak
cat("Quantifying counts in each peak\n")
macs2_counts <- FeatureMatrix(
  fragments = Fragments(adata),
  features = peaks,
  cells = colnames(adata),
  process_n = process_n
)
cat("\n")

# Create a new assay if store_assay is not NULL
if (!is.null(store_assay)) {
  cat(sprintf("Creating new assay in %s", store_assay))
  frags <- Fragments(adata)
  adata[[store_assay]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = frags,
    annotation = Annotation(adata)
  )
  cat(sprintf("Saving Seurat object with new assay to %s", output_dir))
  saveRDS(adata, file.path(output_dir, paste0(basename(input_rds), ".rds")))
  cat("\n")
}
cat("Peak calling completed successfully!")

# Save the feature matrix as matrix exchange format
mtx_dir <- file.path(output_dir, "matrix")
if (!dir.exists(mtx_dir)) {
  dir.create(mtx_dir)
}
cat(sprintf("Saving feature matrix to %s\n", mtx_dir))
writeMM(macs2_counts, file.path(mtx_dir, "mtx.mtx"))
regions <- paste0(seqnames(peaks), ":", start(peaks), "-", end(peaks))
cells <- colnames(adata)
write.table(rownames(macs2_counts), file.path(mtx_dir, "features.tsv"), row.names=F, col.names=F, quote=F)
write.table(colnames(macs2_counts), file.path(mtx_dir, "barcodes.tsv"), row.names=F, col.names=F, quote=F)
cat("\n")
