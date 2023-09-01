#!/usr/bin/env Rscript
# Usage:
# Rscript script_name.R "/path/to/rds_object.rds" "/path/to/output/directory" "version_name" "column1" "column2" "column3"

# Load libraries
suppressMessages(library(Seurat))
suppressMessages(library(Signac))

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
rds_path <- args[1]
output_dir <- args[2]
version <- args[3]
column_names <- args[-(1:3)]

cat("Loading RDS object...\n")
start_time <- Sys.time()
ad <- readRDS(rds_path)
cat("RDS object loaded in", Sys.time() - start_time, "seconds\n")

# Get list of assays available
assay_lst <- Assays(ad)

# Function to generate output file path
generate_file_path <- function(dir, version, assay, subset = NULL, extension = "mtx") {
  file_name <- paste0(version, "_", assay, ".", extension)
  if (!is.null(subset)) {
    file_name <- paste0(version, "_", subset, "_", assay, ".", extension)
  }
  return(file.path(dir, file_name))
}

cat("Processing all cells...\n")
# For each assay, get the counts matrix and save it as mtx
for (assay in assay_lst) {
  start_time <- Sys.time()
  cat("Processing", assay, "assay...\n")
  
  m <- GetAssayData(object = ad, assay = assay, slot = "counts")
  rownames(m) <- sub("-", ":", rownames(m))
  
  write.table(colnames(m), file = generate_file_path(output_dir, version, assay, extension = "obs.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(rownames(m), file = generate_file_path(output_dir, version, assay, extension = "var.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  Matrix::writeMM(m, file = generate_file_path(output_dir, version, assay))
  
  cat(assay, "assay processed in", Sys.time() - start_time, "seconds\n")
}

cat("Processing subsets...\n")
# Get unique combinations of values in the specified columns
unique_combinations <- as.data.frame(unique(ad@meta.data[, column_names]))

# For each unique combination
for (i in seq_len(nrow(unique_combinations))) {
  subset_values <- unique_combinations[i, ]
  subset_name <- paste0(names(subset_values), subset_values, collapse = "_")
  
  cat("Processing", subset_name, "subset...\n")
  
  # For each assay, get the counts matrix and save it as mtx
  for (assay in assay_lst) {
    start_time <- Sys.time()
    cat("Processing", assay, "assay for", subset_name, "subset...\n")
    
    ad_subset <- ad
    for (column in names(subset_values)) {
      ad_subset <- subset(ad_subset, subset = ad_subset@meta.data[[column]] == subset_values[[column]])
    }
    
    m <- GetAssayData(object = ad_subset, assay = assay, slot = "counts")
    rownames(m) <- sub("-", ":", rownames(m))
    
    write.table(colnames(m), file = generate_file_path(output_dir, version, assay, subset_name, "obs.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(rownames(m), file = generate_file_path(output_dir, version, assay, subset_name, "var.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    Matrix::writeMM(m, file = generate_file_path(output_dir, version, assay, subset_name))

    cat("Processed", assay, "assay for", subset_name, "subset in", Sys.time() - start_time, "seconds\n")
  }

  cat("Processed", subset_name, "subset in", Sys.time() - start_time, "seconds\n")
}
