# Wrangling scRNA-seq with ScanPy
1. Building an AnnData object
- From 10x outputs: Easiest thing to do is just to get your data in this format and then read in with the sc.read_10x function

# Wrangling scRNA-seq with Seurat
1. Building a Seurat object
- From 10x outputs: Easiest thing to do is just to get your data in this format and then read in with the Read10X function

scATAC-seq requires a lot more wrangling

# Wrangling scATAC-seq with Signac
1. Building a Seurat object from scATAC-seq data
- Fragment files are a must
- Frag object
- ChromAssay object
- Seurat object
- Pseudobulking
- Peak calling
- Consensus peaks
- Fragment matrix

# Wrangling scATAC-seq with SnapATAC2
- Fragment files are a must
- Per sample AnnDatas

# Wrangling scATAC-seq with pycisTopic
- Fragment files are a must

# Wrangling scATAC-seq with ArchR