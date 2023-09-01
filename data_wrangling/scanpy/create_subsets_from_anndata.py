# Script to subset an h5ad expression dataset and save different file formats
# Usage: python create_h5ad_subsets.py <h5ad_file> <out_dir> <gene_arg> <subset_columns>
# h5ad_file example: /cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/preprocess/snrna/subset/filtered.h5ad
# out_dir example: /cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/preprocess/snrna
# gene_arg example: "0.05"
# subset_columns example: "celltypes, timepoint"

import argparse
import os
import sys

import loompy as lp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from utils import make_dirs, qc, subset_genes, save_h5ad, save_loom, save_tsv

sc.settings.verbosity = 3
np.random.seed(13)

def main(args):
    # Main script function

    # Get args
    h5ad_file = args.h5ad_file
    out_dir = args.out_dir
    dataset_name = args.dataset_name
    gene_arg = args.gene_arg
    save_tsv_arg = args.save_tsv
    normalize_tsv = args.normalize_tsv
    subset_columns = args.subset_columns
    column_values = args.column_values
    skip_all_cells = args.skip_all_cells

    # Parse args
    print("Parsing command line args..")
    h5ad_file = h5ad_file
    if "." in gene_arg:
        cell_frac = float(gene_arg)
        n_genes = None
        use_variable_genes = False
    else:
        cell_frac = None
        n_genes = int(gene_arg)
        use_variable_genes = True
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    qc_dir = make_dirs(os.path.join(out_dir, "qc", dataset_name, gene_arg))
    h5ad_dir = make_dirs(os.path.join(out_dir, "h5ad", dataset_name, gene_arg))
    loom_dir = make_dirs(os.path.join(out_dir, "loom", dataset_name, gene_arg))
    tables_dir = make_dirs(os.path.join(out_dir, "tables", dataset_name, gene_arg))
    save_tsv_flag = save_tsv_arg
    subset_cols = [column.strip(",") for column in subset_columns]
    column_values = [value.strip(",") for value in column_values]
    if len(column_values) == 0:
        column_values = None
    all_cells_flag = skip_all_cells
    normalize_flag = normalize_tsv
    if column_values is not None:
        print(f"Will only keep cells with {column_values} in {subset_cols}...")
    else:
        print(f"Will keep all cells in {subset_cols}...")
    print(f"Will save tsv: {save_tsv_flag} with normalization: {normalize_flag}...")
    print(f"Will skip all cells: {all_cells_flag}...")
    
    # Load in anndata
    print(f"Loading {h5ad_file}...")
    adata = sc.read_h5ad(h5ad_file)

    # Prelim filtters, takes a bout a minute
    print("Prelim filters using ScanPy min_genes=200 and min_cells=3...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.raw = adata

    # Plot QC
    qc(adata, qc_dir)
    
    # Save object for All cells
    if not all_cells_flag:
        curr_ad = subset_genes(adata, use_variable_genes=use_variable_genes, n_genes=n_genes, cell_frac=cell_frac, name="All")
        save_loom(curr_ad, loom_dir, filename="All")
        save_h5ad(curr_ad, h5ad_dir, filename="All")
        if save_tsv_flag:
            save_tsv(curr_ad, tables_dir, filename="All", normalize=normalize_flag) 

    # Subset -- 500 cells, 2000 genes
    print("Subsetting genes for Subset to 500 cells, 2000 genes...")
    curr_ad = adata[adata.obs.index.isin(np.random.choice(adata.obs_names, 500)), adata.var.index.isin(np.random.choice(adata.var_names, 2000))]
    save_loom(curr_ad, loom_dir, filename="Subset")
    save_h5ad(curr_ad, h5ad_dir, filename="Subset")
    if save_tsv_flag:
        save_tsv(curr_ad, tables_dir, filename="Subset", normalize=normalize_flag)

    # Column subsets
    print(f"Found {subset_cols} columns to subset on")
    for column in subset_cols:
        if column not in adata.obs.columns:
            print(f"Column {column} not found in {adata.obs.columns}")
            continue
        options = adata.obs[column].unique()
        print(f"Found {options} values for {column}")
        for value in options:
            print(f"Processing {column}={value}...")
            if column_values is not None and value not in column_values:
                print(f"{column}=value not in {column_values}, skipping...")
                continue
            curr_ad = adata[adata.obs[column] == value]
            curr_ad = subset_genes(curr_ad, use_variable_genes=use_variable_genes, n_genes=n_genes, cell_frac=cell_frac, name=value)
            save_loom(curr_ad, out_dir=loom_dir, filename=value)
            save_h5ad(curr_ad, out_dir=h5ad_dir, filename=value)
            if save_tsv_flag:
                save_tsv(curr_ad, tables_dir, filename=value, normalize=normalize_flag)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad_file", help="h5ad file to subset", required=True)
    parser.add_argument("--out_dir", help="output directory to store all the prepared files, files will be saved in a subdirectory named after the gene_arg", required=True)
    parser.add_argument("--dataset_name", help="name of the dataset, will be used to name the output files", required=True)
    parser.add_argument("--gene_arg", help="if an int, use that many variable genes, if a float, genes that are greater than that fraction of cells", required=True)
    parser.add_argument("--save_tsv", help="save tsv format of expression matrix, the reason for not doing this by default is that it takes a long time", action="store_true")
    parser.add_argument("--normalize_tsv", help="normalize the tsv expression matrix", action="store_true")
    parser.add_argument("--subset_columns", help="columns to slice the object by and save as separate files, if multiple they should be comma separated", nargs="+", default=[])
    parser.add_argument("--column_values", help="values to subset the columns by, if multiple they should be comma separated", nargs="+", default=[])
    parser.add_argument("--skip_all_cells", help="skip the All cells subset", action="store_true")
    args = parser.parse_args()
    main(args)