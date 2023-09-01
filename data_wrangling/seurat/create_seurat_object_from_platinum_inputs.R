# This script builds and saves an Seurat object from a set of platinum standard inputs

# Required inputs
#matrix.mtx.gz - raw counts in sparse MEX format
#barcodes.tsv.gz - list of barcodes corresponding to column indices of the matrix
#features.tsv.gz - list of features corresponding to row indices of the matrix
#metadata.csv.gz - csv file with rows corresponding to barcodes.tsv.gz and any covariates of interest for this data

