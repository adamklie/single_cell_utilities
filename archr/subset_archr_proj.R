# Usage:
# Rscript --vanilla 2B_subset_archr_proj.R \
#   -i "data/atac/archr_proj" \
#   -n "atac_1, atac_2" \
#   -o "data/atac/archr_proj_subset" \
#   -g "hg38" \
#   -t 1 \
#   -s 1234 \
#   -f "TRUE"

# Define options and arguments for the script.
suppressMessages(library(optparse))
option_list <- list(
  make_option(c("-i", "--proj_dir"), type="character", default=NULL, help="ArchR project directory"),
  make_option(c("-n", "--sample_names"), type="character", default=NULL, help="Comma separated list of sample names to subset to."),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Output directory for ArchR project."),
  make_option(c("-g", "--genome"), type="character", default="hg38", help="Genome to use for ArchR project. Default is hg38."),
  make_option(c("-t", "--threads"), type="integer", default=1, help="Number of threads to use for ArchR project. Default is 1."),
  make_option(c("-s", "--seed"), type="integer", default=1234, help="Seed to use for ArchR project. Default is 1234."),
  make_option(c("-f", "--force"), type="character", default="FALSE", help="Force overwrite of output directory. Default is FALSE.")
)
parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options

# Print arguments
print("Arguments:")
proj_dir <- opt$proj_dir
sample_names <- strsplit(opt$sample_names, ",")[[1]]
output_dir <- opt$output_dir
genome <- opt$genome
threads <- opt$threads
seed <- opt$seed
force <- as.logical(opt$force)
print(paste0("proj_dir = ", proj_dir))
print(paste0("sample_names = ", sample_names))
print(paste0("output_dir = ", output_dir))
print(paste0("genome = ", genome))
print(paste0("threads = ", threads))
print(paste0("seed = ", seed))
print(paste0("force = ", force))

# Import libraries
print("Importing libraries\n")
suppressMessages(library(Seurat))
suppressMessages(library(ArchR))
suppressMessages(library(parallel))
suppressMessages(library(tidyverse))

# Set-up
addArchRThreads(threads = threads)
set.seed(seed)
setwd(proj_dir)

# Add annotation
addArchRGenome(genome)

# Load the ArchR project
print("Loading ArchR project\n")
proj = loadArchRProject(path = "./")

# Subset the ArchR project
print(paste0("Subsetting to ", sample_names))
idxSample <- BiocGenerics::which(proj$Sample %in% sample_names)
cellsSample <- proj$cellNames[idxSample]
print(paste0("Number of cells in subset: ", length(cellsSample)))
sub_proj = subsetArchRProject(
    ArchRProj = proj,
    outputDirectory = output_dir,
    cells = cellsSample,
    dropCells = TRUE,
    force = force
)
print("Subsetting complete")
