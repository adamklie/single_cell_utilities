# Usage:
# Rscript --vanilla make_archr_proj.R \
#   -i "data/atac/frag_file.tsv.gz, data/atac/frag_file2.tsv.gz, data/atac/frag_file3.tsv.gz" \
#   -n "atac_1, atac_2, atac_3" \
#   -o "data/atac/archr_proj" \
#   -m "data/atac/atac_metadata.tsv" \
#   -e 4 \
#   -f 1000 \
#   -g "hg38" \
#   -t 1 \
#   -s 1234


# Define options and arguments for the script.
suppressMessages(library(optparse))
option_list <- list(
  make_option(c("-i", "--input_file_list"), type="character", default=NULL, help="Comma separated list of input fragment or BAM files."),
  make_option(c("-n", "--sample_names"), type="character", default=NULL, help="Comma separated list of sample names."),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Output directory for ArchR project."),
  make_option(c("-m", "--project_metadata_path"), type="character", default=NULL, help="Project metadata file"),
  make_option(c("-e", "--min_tsse_enrichment"), type="numeric", default=4, help="Minimum TSS enrichment for ArchR project. Default is 4."),
  make_option(c("-f", "--min_frags"), type="numeric", default=1000, help="Minimum number of fragments for ArchR project. Default is 1000."),
  make_option(c("-g", "--genome"), type="character", default="hg38", help="Genome to use for ArchR project. Default is hg38."),
  make_option(c("-t", "--threads"), type="integer", default=1, help="Number of threads to use for ArchR project. Default is 1."),
  make_option(c("-s", "--seed"), type="integer", default=1234, help="Seed to use for ArchR project. Default is 1234.")
)
parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options

# Print arguments
input_file_list <- strsplit(opt$input_file_list, ",")[[1]]
sample_names <- strsplit(opt$sample_names, ",")[[1]]
output_dir <- opt$output_dir
project_metadata_path <- opt$project_metadata_path
min_tsse_enrichment <- opt$min_tsse_enrichment
min_frags <- opt$min_frags
genome <- opt$genome
threads <- opt$threads
seed <- opt$seed
print(paste0("input_file_list = ", input_file_list))
print(paste0("sample_names = ", sample_names))
print(paste0("output_dir = ", output_dir))
print(paste0("project_metadata_path = ", project_metadata_path))
print(paste0("min_tsse_enrichment = ", min_tsse_enrichment))
print(paste0("min_frags = ", min_frags))
print(paste0("genome = ", genome))
print(paste0("threads = ", threads))
print(paste0("seed = ", seed))

# Import libraries
print("Importing libraries\n")
suppressMessages(library(Seurat))
suppressMessages(library(ArchR))
suppressMessages(library(parallel))
suppressMessages(library(tidyverse))

# Set-up
addArchRThreads(threads = threads)
set.seed(seed)
setwd(output_dir)

# Add annotation
addArchRGenome(genome)

# Make arrow files 
print("Making arrow files and computing QCs\n")
ArrowFiles <- createArrowFiles(
  inputFiles = input_file_list,
  sampleNames = sample_names,
  minTSS = min_tsse_enrichment,
  minFrags = min_frags,
  excludeChr = c("chrM"),
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# Make archr project
print("Making ArchR project\n")
proj <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = "./"
)
print(paste0("Memory Size = ", round(object.size(proj) / 10^6, 3), " MB\n"))

# Read in project metadata if it exists
if (is.null(project_metadata_path)) {
    print("No project metadata provided, skipping...")
} else {
    print("Reading in project metadata\n")
    project_metadata <- read.csv(project_metadata_path, sep = "\t")
    
    # Clean up the archr proj metadata
    archr_metadata = as.data.frame(proj@cellColData)
    archr_metadata$sample_id = archr_metadata$Sample
    
    # Merge and add batch, timepoint, and timecourse (TODO: make this just a join that gets added to dataframe)
    archr_metadata = dplyr::left_join(archr_metadata, project_metadata, by = "sample_id")
    proj$batch = archr_metadata$batch
    proj$timepoint = archr_metadata$timepoint
    proj$condition = archr_metadata$condition
    proj$timecourse = paste(proj$batch, proj$condition, sep = "_")
}

# Save initial project metadata
proj_meta = as.data.frame(proj@cellColData)
write.table(proj_meta, file = "initial_archr_proj_meta.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

# Add doublet scores and filter
print("Adding doublet scores and filtering\n")
proj = addDoubletScores(proj, k = 10, knnMethod = "LSI")
proj = filterDoublets(proj)

# Save post-filtering project metadata
proj_meta = as.data.frame(proj@cellColData)
write.table(proj_meta, file = "post_filtering_archr_proj_meta.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Save the project
print("Saving ArchR project\n")
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "./"
)
