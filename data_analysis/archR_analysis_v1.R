#ArchR parameters
atac_frag = "" # Input file  
genome = "" #either hg38 or mm10

#ArchR QC
min_tss = 4 #The minimum numeric transcription start site (TSS) enrichment score required for a cell to pass filtering
min_frags = 1000 #The minimum number of mapped ATAC-seq fragments required per cell to pass filtering for use
add_tile_mat = TRUE #A boolean value indicating whether to add a "Tile Matrix" to each ArrowFile. 
add_gene_score_mat = TRUE #A boolean value indicating whether to add a Gene-Score Matrix to each ArrowFile.

#ArchR Doublet paramaters
find_doublets = FALSE
doublet_k = 10 #The number of cells neighboring a simulated doublet to be considered as putative doublets.
doublet_knn_method = "UMAP" #Refers to the embedding to use for nearest neighbor search.
lsi_method = 1 #A number or string indicating the order of operations in the TF-IDF normalization. Possible values are: 1 or "tf-logidf", 2 or "log(tf-idf)", and 3 or "logtf-logidf". 

copy_arrow_files = TRUE #save a copy of arrow files in the ArchR project (recommended)
iter_LSI_matrix = "TileMatrix" #The name of the data matrix to retrieve from the ArrowFiles associated with the ArchRProject. Valid options are "TileMatrix" or "PeakMatrix".
threads = 1
prefix = "prefix" #project name

#ArchR Plots parameters
marker_features_test = "wilcoxon" #The name of the pairwise test method to use in comparing cell groupings to the null cell grouping during marker feature identification.
heatmap_transpose = TRUE #Boolean to transpose heatmap
heatmap_label_n = 5 #Top n genes to label per cluster in heatmap
heatmap_cutoff = "FDR <= 0.01 & Log2FC >= 0.5" #Cut-off applied to genes in heatmap

#jupyter notebook plot sizes
options(repr.plot.width=20, repr.plot.height=15)

suppressMessages(library(parallel))
suppressMessages(library(ArchR))
suppressMessages(library(magick))
suppressMessages(library(logr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpointdensity))

set.seed(1)
addArchRThreads(threads = threads)

addArchRGenome(genome)

ArrowFiles <- createArrowFiles(
      inputFiles = atac_frag,
      sampleNames = prefix,
      minTSS = min_tss, 
      minFrags = min_frags,
      addTileMat = add_tile_mat,
      addGeneScoreMat = add_gene_score_mat
    )

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = prefix,
  copyArrows = copy_arrow_files, #This is recommened so that you maintain an unaltered copy for later usage.
  showLogo = FALSE
)

df <- as.data.frame(getCellColData(proj, select = c("nFrags", "TSSEnrichment")))
                                                                                                      
ggplot(data=df, aes(x = nFrags, y = TSSEnrichment)) + 
    geom_pointdensity(method = "default") +
    scale_colour_gradientn(colors = paletteContinuous(set = "sambaNight")) +
    scale_x_continuous(trans = "log10", breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    expand_limits(x = min(df$nFrags), y = 0) + 
    geom_hline(yintercept = min_tss, lty = "dashed") + 
    geom_vline(xintercept = min_frags, lty = "dashed") +
    xlab(label ="Unique Fragments") + 
    ylab(label = "TSS Enrichment") + 
    annotation_logticks(sides = "b") +
    labs(fill="density") +
    theme(axis.title=element_text(size=14), axis.text=element_text(size=10), legend.title=element_text(size=14), legend.text=element_text(size=8.5)

plotFragmentSizes(ArchRProj = proj) + theme_gray()

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = iter_LSI_matrix, name = "IterativeLSI")

proj <- addClusters(input = proj, reducedDims = "IterativeLSI")

proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI") 

plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP") + 
                geom_point(size = 0.2)+
                theme_gray()

proj
nFrags)
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "lognFrags", embedding = "UMAP", 
                     pal = ArchRPalettes$purpleOrange) + 
        geom_point(size = 0.2)+
        theme_gray()

plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP",
                    pal = ArchRPalettes$purpleOrange) + 
        geom_point(size = 0.2)+
        theme_gray()

plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "FRIP", embedding = "UMAP",
                    pal = ArchRPalettes$purpleOrange) + 
        geom_point(size = 0.2)+
        theme_gray()

markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = marker_features_test
)

plotMarkerHeatmap(markersGS, transpose = heatmap_transpose, nLabel = heatmap_label_n, 
                         cutOff = heatmap_cutoff, plotLog2FC = TRUE)