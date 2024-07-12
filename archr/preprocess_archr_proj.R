# Usage:
# Rscript --vanilla 3B_preprocess_archr_proj.R


# Define options and arguments for the script.
suppressMessages(library(optparse))
option_list <- list(
  make_option(c("-p", "--proj_dir"), type="character", default=NULL, help="Path to ArchR project."),
  make_option(c("-c", "--clustering_resolution"), type="numeric", default=0.5, help="Clustering resolution for ArchR project. Default is 0.5."),
  make_option(c("-b", "--run_harmony"), type="character", default="FALSE", help="Run Harmony for ArchR project. Default is FALSE."),
  make_option(c("-k", "--umap_neighbors"), type="integer", default=30, help="Number of UMAP neighbors for ArchR project. Default is 30."),
  make_option(c("-m", "--umap_min_dist"), type="numeric", default=0.5, help="Minimum UMAP distance for ArchR project. Default is 0.3."),
  make_option(c("-u", "--umap_metric"), type="character", default="cosine", help="UMAP metric for ArchR project. Default is cosine."),
  make_option(c("-g", "--genome"), type="character", default="hg38", help="Genome to use for ArchR project. Default is hg38."),
  make_option(c("-t", "--threads"), type="integer", default=1, help="Number of threads to use for ArchR project. Default is 1."),
  make_option(c("-s", "--seed"), type="integer", default=1234, help="Seed to use for ArchR project. Default is 1234.")
)
parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options

# Print arguments
proj_dir <- opt$proj_dir
clustering_resolution <- opt$clustering_resolution
run_harmony <- as.logical(opt$run_harmony)
umap_neighbors <- opt$umap_neighbors
umap_min_dist <- opt$umap_min_dist
umap_metric <- opt$umap_metric
genome <- opt$genome
threads <- opt$threads
seed <- opt$seed
print(paste0("ArchR project path = ", proj_dir))
print(paste0("Clustering resolution = ", clustering_resolution))
print(paste0("UMAP neighbors = ", umap_neighbors))
print(paste0("UMAP min dist = ", umap_min_dist))
print(paste0("UMAP metric = ", umap_metric))
print(paste0("Genome = ", genome))
print(paste0("Threads = ", threads))
print(paste0("Seed = ", seed))

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
addArchRGenome("hg38")

# Load the ArchR project
print("Loading ArchR project\n")
proj = loadArchRProject(path = "./")

# ArchR dim reduection
print("Running ArchR dimensionality reduction with IterativeLSI\n")
proj = addIterativeLSI(
    proj,  
    useMatrix = "TileMatrix",
    name = "IterativeLSI", 
    force = T
)

# Run clustering
print("Running clustering using Seurat method\n")
proj = addClusters(
    proj, 
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = clustering_resolution,
    force = TRUE
)

# Optionally run Harmony
if (run_harmony) {
    print("Running Harmony")
    proj = addHarmony(
        ArchRProj = proj, 
        reducedDims = "IterativeLSI", 
        name = "Harmony", 
        groupBy = "SampleID"
    )
    reduction = "Harmony"
} else {
    print("Using IterativeLSI")
    reduction = "IterativeLSI"
}

# Run UMAP
print("Running UMAP\n")
proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = reduction,
    name = "UMAP", 
    nNeighbors = umap_neighbors,
    minDist = umap_min_dist,
    metric = umap_metric,
)

# UMAP plots
umap_plot_vars = c("Sample", "Clusters")
p <- c()
for (i in umap_plot_vars) {
    p[[i]] <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = i, embedding = "UMAP")
}
plotPDF(
    plotList = p, 
    name = "Plot-UMAP-Sample-Clusters.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, 
    width = 5, 
    height = 5
)

# Save cell metadata to cell_metadata.tsv
print("Saving cell metadata to cell_metadata.tsv\n")
proj_meta = as.data.frame(proj@cellColData)
write.table(proj_meta, file = "cell_metadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Save object with new stuff added
print("Saving ArchR project\n")
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "./",
)
