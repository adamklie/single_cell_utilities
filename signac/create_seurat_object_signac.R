# Define options and arguments for the script.
suppressMessages(library(optparse))
option_list <- list(
  make_option(c("-i", "--input_cellranger_dirs"), type="character", default=NULL, help="Input cellranger directories."),
  make_option(c("-s", "--sample_ids"), type="character", default=NULL, help="Sample IDs."),
  make_option(c("-b", "--annotation_path"), type="character", default=NULL, help="Barcode annotation to add."),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Output directory for Seurat object."),
)
parser <- OptionParser(option_list=option_list)

# Parse arguments
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options

# Print arguments
input_cellranger_dirs <- opt$input_cellranger_dirs
sample_ids <- opt$sample_ids
annotation_path <- opt$annotation_path
output_dir <- opt$output_dir
cat(sprintf("Input cellranger directories: %s\n", input_cellranger_dirs))
cat(sprintf("Sample IDs: %s\n", sample_ids))
cat(sprintf("Barcode metadata: %s\n", annotation_path))
cat(sprintf("Output directory: %s\n", output_dir))

# Imports
cat("Importing libraries\n")
suppressMessages(library(SeuratDisk))
suppressMessages(library(SeuratData))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

# Get gene annotations for hg3a8
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# Create chromatin assays for each input cellranger dir
cat("Creating Seurat objects\n")
seurat_objs <- list()
for (i in 1:length(input_cellranger_dirs)) {
  cellranger_dir <- input_cellranger_dirs[i]
  sample_id <- sample_ids[i]
  cat(sprintf("Creating Seurat object for %s\n", sample_id))
  counts <- Read10X_h5(file.path(cellranger_dir, "filtered_feature_bc_matrix.h5"))
  assay <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = file.path(cellranger_dir, "atac_fragments.tsv.gz"),
    annotation = annotation
  )
  seurat_obj <- CreateSeuratObject(assay, assay = "peaks")
  seurat_obj$dataset <- cellranger_dir
  saveRDS(seurat_obj, file.path(output_dir, paste0(sample_id, "_seurat.rds")))
  seurat_objs[[sample_id]] <- seurat_obj
}

# Merge first with rest
combined <- merge(
  x = seurat_objs[[1]],
  y = seurat_objs[2:length(seurat_objs)],
  add.cell.ids = sample_ids,
)

# If barcode metadata is provided, filter based on it and add it to the Seurat object
annotations = read.csv(annotations_path, row.names = 1, sep = "\t", header = FALSE)
annotations = as.data.frame(annotations)
rownames(annotations) = gsub("#", "_", rownames(annotations))
cellids = rownames(annotations)
matched_ids = intersect(cellids, rownames(adata@meta.data))
idxSample <- BiocGenerics::which(rownames(adata@meta.data) %in% matched_ids)
adata <- subset(adata, cells = matched_ids)
adata$annotation <- annotations[matched_ids,]
cat(table(adata$annotation))

# Save as RDS
saveRDS(adata, file.path(output_dir, "signac.rds"))
