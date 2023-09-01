library(harmony)
library(Matrix)
library(Seurat)  
library(ggplot2)  
library(glmGamPoi)
library(tidyverse)
library(optparse)
 
option_list = list(
  make_option(c("-t", "--tissue"), type="character", default=NULL, 
              help="tissue name"),
  make_option(c("-g", "--genes"), type="double", default=300, 
              help="minimum number of genes expressed per nucleus"),
  make_option(c("-d", "--doublets"), type="double", default=0.25, 
              help="maximum doublet score"),
  make_option(c("-m", "--mito"), type="double", default=5, 
              help="maximum percent mitochondrial gene expression")
)

opt = parse_args(OptionParser(option_list=option_list))
tissue = opt$tissue
n_genes = opt$genes
doublets = opt$doublets
mito = opt$mito

source("plotting_functions.R")

seurat_obj = function(tissue){
    counts = t(readMM(paste0("../scanpy/",tissue,"_sparse_matrix.mtx")))
    meta = read.csv(paste0("../scanpy/",tissue,"_obs.csv"))
    genes = read.csv(paste0("../scanpy/",tissue,"_genes.csv"))
    
    rownames(counts) = genes$gene_name
    colnames(counts) = meta$kallisto_cellID
    meta = meta[meta$batch != 0,]
    counts = counts[,colnames(counts) %in% meta$kallisto_cellID]
    counts = counts[,match(meta$kallisto_cellID,colnames(counts))]
    
    obj = CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
    obj@meta.data = cbind(obj@meta.data,meta)
    obj[["percent.mt"]] = PercentageFeatureSet(obj, pattern = "^mt-")
    obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa")
    out = obj
    
    }

obj = seurat_obj(tissue)
obj = subset(obj,
             subset = nCount_RNA > 500 &
             nFeature_RNA > n_genes &
             doublet_scores < doublets &
             percent.mt < mito)

# Normalize
obj <- SCTransform(obj, method = "glmGamPoi", 
                   vars.to.regress = c("percent.mt","nFeature_RNA"), verbose = F)
# PCA
obj <- RunPCA(obj, verbose = T, npcs = 50)

# Harmony integration across sublibraries
obj = RunHarmony(obj, "sublibrary", plot_convergence = F, assay.use = "SCT")

# UMAP and clustering
obj = RunUMAP(obj, reduction = "harmony", dims = 1:30,verbose = F)
obj = FindNeighbors(obj,reduction = "harmony", dims = 1:30,verbose = F) 
obj = FindClusters(obj,verbose = F)

# predict celltypes from reference
obj = predict_celltypes(obj, tissue)

# add cell cycle scores 
load("../../ref/mouse_cellcycle_genes.rda")
obj<- CellCycleScoring(obj, s.features = m.s.genes, g2m.features = m.g2m.genes)

# add marker gene scores
markers = read.csv("../../ref/IGVF_curated_markers.csv")
markers = markers[markers$Tissue == tissue,]
markers = markers[markers$Gene %in% rownames(obj),]

for (celltype in unique(markers$Subtype)){
    label = paste0(celltype,"_score")
    genes = unique(markers$Gene[markers$Subtype == celltype])
    obj[[label]] = PercentageFeatureSet(obj,features = genes)
    
}

saveRDS(obj,file=paste0("../seurat/",tissue,"_seurat_novaseq.rds"))