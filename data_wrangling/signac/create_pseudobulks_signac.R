#!/usr/bin/env Rscript

# This script takes in a processed Seurat object that properly formatted fragment files and groupings of interest
# in the metadata and writes out fragments to BED files for each group in the grouping of interest. This is a very simple 
# script with only one real command (that takes a long time).

# This script uses optparse to parse arguments for specifying the path to the written Seurat object and the output directory
# with optional parameters for data processing.

# Required arguments:
# rds_file: Path to RDS file containing processed Seurat object.
# output_dir: Path to directory to output final object and intermediate files.

# Optional arguments
# group_by: Column in the objects metadata to split into pseudobulks. If you do not include this, it will be set to NULL and you will get either one big fragment file or just rewrite all the fragment files depending on how you set append
# assay_name: the name of the assay in the Seurat object woith the Fragment object
# append: Whether you want to keep each fragment file separate or to group all fragment files together. By default this is TRUE and you will get one pseudobulk BED file for each group. Set this to FALSE if you want a each group to have a file per fragment file (will end up with # fragment files x # groups BED files)
# validate: If TRUE, will run sanity checks on the fragment files using Signac functions

# The script will:
# 1. Load in the rds file with the processed Seurat object
# 2. Verify that the assay you passed in exists and that the group_by column is contained in the metadata
# 3. Validate the fragment file (optional)
# 4. Run SplitFragments to write out the pseudobulks

# Usage:
# Rscript create_pseudobulks_from_signac.R --

# Imports
library(optparse)
library(Seurat)
library(Signac)

# Define options and arguments for the script.
option_list <- list(
  make_option(c("-r", "--rds_file"), type="character", default=NULL,
              help="Path to RDS file containing processed Seurat object", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL,
              help="Path to directory to output final object and intermediate files", metavar="character"),
  make_option(c("-g", "--group_by"), type="character", default=NULL,
              help="Column in the object's metadata to split into pseudobulks. Default is NULL.", metavar="character"),
  make_option(c("-a", "--assay_name"), type="character", default="RNA",
              help="Name of the assay in the Seurat object with the Fragment object. Default is 'RNA'.", metavar="character"),
  make_option(c("--append"), type="logical", default=TRUE,
              help="Whether to keep each fragment file separate or to group all fragment files together. Default is TRUE."),
  make_option(c("--validate"), type="logical", default=FALSE,
              help="If TRUE, will run sanity checks on the fragment files using Signac functions. Default is FALSE.")
)

# Parse command-line arguments
args <- parse_args(OptionParser(option_list=option_list))
print(args)

# Load the Seurat object from the specified RDS file.
message("Loading Seurat object from ", args$rds_file)
adata <- readRDS(args$rds_file)

# Verify assay and group_by column
if (!(args$assay_name %in% names(adata@assays))) {
  stop(paste("Specified assay:", args$assay_name, "not found in the Seurat object. Exiting."))
}
if (!is.null(args$group_by) && !(args$group_by %in% colnames(adata@meta.data))) {
  stop(paste("Specified group_by column:", args$group_by, "not found in the metadata. Exiting."))
}

# Validate fragment files if needed
if (args$validate) {
  message("Validating fragment files...")
  frags <- Fragments(adata)
  ValidateFragments(frags, verbose = TRUE)
  ValidateHash(frags, verbose = TRUE)
  ValidateCells(frags, verbose = TRUE)
}

# Create pseudobulks using SplitFragments
message("Creating pseudobulks using SplitFragments...")
SplitFragments(
  object=adata,
  assay=args$assay_name,
  group.by=args$group_by,
  outdir=args$output_dir,
  append=args$append
)

message("Pseudobulks creation completed successfully!")
