#!/usr/bin/env Rscript

# This script takes in a processed Seurat object and a bed file and creates a feature matrix. Assumes that their exists a Fragments
# obect within the Seurat object under the input assay

# This script uses optparse to parse arguments for specifying the path to the written Seurat object and the output directory
# with optional parameters for data processing.

# Required arguments:
# rds_file: Path to RDS file containing processed Seurat object.
# bed_file: Path to BED file containing regions to serve as rows of the feature matrix
# output_dir: Path to directory to output final matrix file and all intermediates

# Optional arguments
# assay_name: the name of the assay in the Seurat object woith the Fragment object
# process_n: Number of regions to load into memory at a time, per thread. Processing more regions at once can be faster but uses more memory.
# validate: If TRUE, will run sanity checks on the fragment files using Signac functions
# store_assay: By default this is NULL, if set, it will store the new FragmentMatrix as a

# The script will:
# 1. Load in the rds file with the processed Seurat object

# Usage:
# Rscript create_pseudobulks_from_signac.R --

# Argument parsing
suppressMessages(library(optparse))
option_list <- list(
  make_option(c("-r", "--rds_file"), type="character", default=NULL,
              help="Path to RDS file containing processed Seurat object", metavar="character"),
  make_option(c("-b", "--bed_file"), type="character", default=NULL,
              help="Path to BED file containing regions for the feature matrix", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL,
              help="Path to directory to output final matrix file and intermediates", metavar="character"),
  make_option(c("-f", "--file_name"), type="character", default="matrix.mtx",
              help="Name of the matrix file to output. Default is 'matrix.mtx'", metavar="character"),
  make_option(c("-n", "--process_n"), type="integer", default=2000, 
              help="Number of regions to load into memory at a time, per thread. Processing more regions at once can be faster but uses more memory. Default is 2000.", metavar="integer"),
  make_option("--assay_name", type="character", default="ATAC",
              help="Name of the assay in the Seurat object with the Fragment object. Default is 'ATAC'", metavar="character"),
  make_option("--validate", type="logical", default=FALSE,
              help="Option to run sanity checks on the fragment files. Default is FALSE."),
  make_option("--store_assay", type="character", default=NULL,
              help="If set, will store the new FragmentMatrix as a specified assay. Default is NULL", metavar="character")
)
parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
print(opt)

# Grab the arguments
cat("Parsing arguments")
rds_file <- opt$rds_file
bed_file <- opt$bed_file
output_dir <- opt$output_dir
file_name <- opt$file_name
process_n <- opt$process_n
assay_name <- opt$assay_name
validate <- opt$validate
store_assay <- opt$store_assay
cat("\n")

# Print the arguments
cat(sprintf("rds_file: %s\n", rds_file))
cat(sprintf("bed_file: %s\n", bed_file))
cat(sprintf("output_dir: %s\n", output_dir))
cat(sprintf("file_name: %s\n", file_name))
cat(sprintf("process_n: %s\n", process_n))
cat(sprintf("assay_name: %s\n", assay_name))
cat(sprintf("validate: %s\n", validate))
cat(sprintf("store_assay: %s\n", store_assay))
cat("\n")

# Ensure all required arguments are provided
if (is.null(rds_file) | is.null(bed_file) | is.null(output_dir)) {
  stop("Missing required arguments. Exiting.")
}

# Check if output directory exists, if it does already exist, warn the user
if (dir.exists(output_dir)) {
    cat(sprintf("Warning. Output directory %s already exists.\n", output_dir))
    cat("\n")
} else {
    cat(sprintf("Creating output directory %s\n", output_dir))
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    cat("\n")
}

# Imports
cat("Loading libraries\n")
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Matrix))
source(here::here("/cellar/users/aklie/data/igvf/bin/data_wrangling", 'io.R'))
cat("\n")

# Load in the peaks from the bed file into a granges object
cat("Loading peaks from ", bed_file)
peaks <- bed_to_granges(bed_file)
cat(sprintf("Loaded %s peaks\n", length(peaks)))
cat("\n")

# Load the Seurat object from the specified RDS file.
cat("Loading Seurat object from ", rds_file)
adata <- readRDS(rds_file)
cat(sprintf("Loaded Seurat object with %s cells\n", ncol(adata)))
cat("\n")

# Verify assay in the Seurat object
if (!(assay_name %in% names(adata@assays))) {
  stop(paste("Specified assay:", assay_name, "not found in the Seurat object. Exiting."))
}

# Validate fragment files if the option is set to TRUE
if (validate) {
  cat("Validating fragment files")
  frags <- Fragments(adata)
  ValidateFragments(frags, verbose = TRUE)
  ValidateHash(frags, verbose = TRUE)
  ValidateCells(frags, verbose = TRUE)
  cat("\n")
}

# Quantify counts in each peak
cat("Quantifying counts in each peak")
macs2_counts <- FeatureMatrix(
  fragments = Fragments(adata),
  features = peaks,
  cells = colnames(adata),
  process_n = process_n
)
cat("\n")

# Save the feature matrix as matrix exchange format
mtx_file <- file.path(output_dir, file_name)
cat(sprintf("Saving feature matrix to %s", mtx_file))
writeMM(macs2_counts, mtx_file)
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
  saveRDS(adata, file.path(output_dir, paste0(basename(rds_file), ".rds")))
  cat("\n")
}
cat("Processing completed successfully!")
