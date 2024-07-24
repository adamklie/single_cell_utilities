# Script to convert a Seurat object to h5ad
# Subsets the object to only the "RNA" assay within the Seurat object (assumes there is one)
# The assumption here is that one is going to want the raw counts from the Seurat object
# in order to do different analysis in Python

# Package imports
suppressMessages(library(Seurat))
suppressMessages(library(loomR))
suppressMessages(library(PISCES))
suppressMessages(library(SeuratDisk))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
theme_set(theme_cowplot())

# Get arguments from command line
args = commandArgs(trailingOnly=TRUE)
rds_file <- args[1] # /cellar/users/aklie/projects/igvf/topic_grn_links/data/mouse_heart/auxiliary_data/snrna/heart_Parse_10x_integrated.rds
out_dir <- args[2] # /cellar/users/aklie/projects/igvf/topic_grn_links/data/mouse_heart/prepare_inputs/snrna/
out_name <- args[3] # heart_Parse_10x_integrated.h5ad

# Read in the R object
print(sprintf("Reading in %s...", rds_file))
dat <- readRDS(rds_file)

# Grab only the RNA counts
print("Putting the object on a diet...")
DefaultAssay(dat) <- "RNA"
rna_dat <- DietSeurat(dat, assays=c("RNA"))

# Save as h5seuat
print("Saving h5seurat...")
SaveH5Seurat(rna_dat, filename = file.path(out_dir, "tmp", "Converted.h5Seurat"))

# Convert to h5ad
print("Saving h5ad...")
Convert(file.path(out_dir, "tmp", "Converted.h5Seurat"), dest = "h5ad")

# Move this out of tmp
print("Moving h5ad out of tmp...")
file.rename(file.path(out_dir, "tmp", "Converted.h5ad"), file.path(out_dir, out_name))
