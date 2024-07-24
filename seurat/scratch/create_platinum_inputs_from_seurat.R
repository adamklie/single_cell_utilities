# This script is meant for reading in a Seurat object from an RDS file and
# writing out a platinum file set for that object at a passed in location
# The platinum output follows the 10x CellRanger v3 output standard. 

# This script uses optparse to parse arguments and takes in the following
# rds_path - required, path to rds object containing Seurat object to process
# out_dir - required, path to dump platinum outputs, must not exist
# assay - optional, the assay of the Seurat object to get the counts matrix from, default is "RNA"
# slot - optional, the slot of the assay to grab the counts matrix from, default is "counts"
# version - optional, CellRanger version to right outputs in ("2" or "3"), default is "3"

# Argument parsing
library(optparse)
option_list <- list(
  make_option(c("-r", "--rds_path"), type="character", default=NULL, help="Path to RDS object containing Seurat object to process", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default=NULL, help="Path to dump platinum outputs", metavar="character"),
  make_option(c("-a", "--assay"), type="character", default="RNA", help="The assay of the Seurat object to get the counts matrix from", metavar="character"),
  make_option(c("-s", "--slot"), type="character", default="counts", help="The slot of the assay to grab the counts matrix from", metavar="character"),
  make_option(c("-v", "--version"), type="character", default="3", help="CellRanger version to write outputs in ('2' or '3')", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Ensure required arguments are provided
if(is.null(opt$rds_path) || is.null(opt$out_dir)){
  stop("Both rds_path and out_dir are required arguments.")
}

# Conversion libraries and Seurat
library(SeuratDisk)
library(SeuratData)
library(Seurat)
library(Signac)
library(Matrix)
library(R.utils)
source(here::here("/cellar/users/aklie/data/igvf/bin/data_wrangling", 'io.R'))

# Actually read it
adata <- readRDS(opt$rds_path)

# Using the specified CellRanger version as output, create the inputs for platinum
write10xCounts(x=adata@assays[[opt$assay]]@opt$slot, path=file.path(opt$out_dir), version=opt$version)
write.csv(x=adata@meta.data, file=file.path(opt$out_dir, "metadata.csv"), quote=FALSE)

# Test round trip
mtx <- Read10X(file.path(opt$out_dir)
metadata <- read.csv(file.path(opt$out_dir, "metadata.csv"), row.names = 1)

# Load the data and report a detailed error if this fails
adata_load <- CreateSeuratObject(counts = mtx, meta.data = metadata)
