#!/usr/bin/env Rscript

# This script takes in BED file and uses the CallPeak functionality from Signac to call peaks on it. Often used in parallel
# after running the create_pseudobulks_signac.R

# This script uses optparse to parse arguments for specifying the path to the written Seurat object and the output directory
# with optional parameters for data processing.

# Required arguments:
# bed_file: Path to directory that includes BED files
# output_dir: Path to directory to output peak calls

# Optional arguments:
# --macs2_path: Specify the path to the macs2 executable. 
#               Default is '/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/bin/macs2'.
# --broad: Option for broad peak calling. Default is FALSE.
# --format: Define the format of the output. Default is 'BED'.
# --effective_genome_size: Set the effective genome size for peak calling. Default is 2.7e+09.
# --extsize: Set the extension size for peak calling. Default is 200.
# --shift: Set the shift size for peak calling. Default is -100.
# --cleanup: Decide whether to clean up intermediate files. Default is FALSE.
# --verbose: Enable verbose output during peak calling. Default is TRUE.

# The script will:
# 1. Parse the command-line arguments and validate the input BED file path.
# 2. Load the specified BED file for peak calling.
# 3. Utilize Signac's CallPeaks function to call peaks on the provided BED file using the specified (or default) parameters.
# 4. Save the peak calling results to the specified output directory.

# Usage:
# Rscript call_peaks_signac.R --bed_file=/path/to/your/input.bed --output_dir=/path/to/your/output/directory/

# Imports
library(optparse)
library(Signac)

# Define options and arguments for the script.
option_list <- list(
  make_option(c("-b", "--bed_file"), type="character", default=NULL,
              help="Path to directory that includes BED files", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL,
              help="Path to directory to output peak calls", metavar="character"),
  make_option("--macs2_path", type="character", default="/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/bin/macs2",
              help="Path to the macs2 executable. Default is '/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/bin/macs2'.", metavar="character"),
  make_option("--broad", type="logical", default=FALSE,
              help="Option for broad peak calling. Default is FALSE."),
  make_option("--format", type="character", default="BED",
              help="Format of the output. Default is 'BED'.", metavar="character"),
  make_option("--effective_genome_size", type="numeric", default=2.7e+09,
              help="Effective genome size. Default is 2.7e+09.", metavar="numeric"),
  make_option("--extsize", type="integer", default=200,
              help="Extension size. Default is 200.", metavar="integer"),
  make_option("--shift", type="integer", default=-100,
              help="Shift size. Default is -100.", metavar="integer"),
  make_option("--cleanup", type="logical", default=FALSE,
              help="Option to clean up intermediate files. Default is FALSE."),
  make_option("--verbose", type="logical", default=TRUE,
              help="Option for verbose output. Default is TRUE.")
)

# Parse command-line arguments
args <- parse_args(OptionParser(option_list=option_list))
print(args)

# Validate that the BED file exists.
if (!file.exists(args$bed_file)) {
  stop(paste("Specified BED file:", args$bed_file, "does not exist. Exiting."))
}

# Load in the BED file.
fragpath <- args$bed_file

# Call peaks
message("Calling peaks using Signac's CallPeaks function...")
grlist <- list()
CallPeaks(
  object = fragpath,
  macs2.path = args$macs2_path,
  outdir = args$output_dir,
  broad = args$broad,
  format = args$format,
  effective.genome.size = args$effective_genome_size,
  extsize = args$extsize,
  shift = args$shift,
  name = basename(args$bed_file), # Use the name of the BED file as the identifier.
  cleanup = args$cleanup,
  verbose = args$verbose
)

message("Peak calling completed successfully!")
