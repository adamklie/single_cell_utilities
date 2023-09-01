# This script will not keep the intermediate fragment files!!!
suppressMessages(library(SeuratDisk))
suppressMessages(library(SeuratData))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))

# Actually read it, takes about 2 min
rds_path <- # TODO
adata <- readRDS(rds_path)

peaks <- CallPeaks(
  adata,
  assay = "mpeak",
  group.by = "predicted.cell.type",
  idents = NULL,
  macs2.path = "/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/bin/macs2",
  broad = FALSE,
  format = "BED",
  outdir = "/cellar/users/aklie/data/igvf/beta_cell_networks/aligned/igvf_sc-islet_10X-Multiome/25Aug23/CallPeaks",
  fragment.tempdir = "/cellar/users/aklie/data/igvf/beta_cell_networks/aligned/igvf_sc-islet_10X-Multiome/25Aug23/CallPeaks",
  combine.peaks = TRUE,
  effective.genome.size = 2.7e+09,
  extsize = 200,
  shift = -100,
  cleanup = FALSE,
  verbose = TRUE
)