# Script to subset a fully processed Seurat object into subset inputs
# Usage:
# Rscript script_name.R "/path/to/rds_object.rds" "/path/to/output/directory" "version_name" "column1" "column2" "column3"
# h5ad_file example: "/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/preprocess/snrna/subset/filtered.h5ad"
# out_dir example: "/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/preprocess/snrna"
# gene_arg example: "0.05"
# subset_columns example: "celltypes, timepoint"

# Load libraries
suppressMessages(library(Seurat))
suppressMessages(library(loomR))
suppressMessages(library(PISCES))
suppressMessages(library(SeuratDisk))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
theme_set(theme_cowplot())

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
rds_path <- args[1]
output_dir <- args[2]
assay_name <- args[3]
subset_cols <- args[-c(1:3)]

# Load RDS object
cat("Loading RDS object...\n")
start_time <- Sys.time()
ad <- readRDS(rds_path)
cat("RDS object loaded in", Sys.time() - start_time, "seconds\n")

# For each assay, get the counts matrix and save it as tsv
cat("Processing all cells...\n")
for (assay in assay_lst) {
  start_time <- Sys.time()
  cat("Processing", assay, "assay...\n")
  
  m <- GetAssayData(object = ad, assay = assay, slot = "data")
  print("Saving ARACNe tsv for full matrix...")
  mat.name <- file.path(out_dir, 'normalized.tsv')
  ARACNeTable(norm.mat, mat.name, subset = FALSE)

  m <- GetAssayData(object = ad, assay = assay, slot = "scale.data") 
  cat(assay, "assay processed in", Sys.time() - start_time, "seconds\n")
}

# Grab the normalized data
norm.mat <- dat@assays$SCT@scale.data

# Save as tsv

