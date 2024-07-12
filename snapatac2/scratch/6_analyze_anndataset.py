#!/usr/bin/env python3
"""
6_analyze_anndataset.py

This script analyzees a AnnDataset using SnapATAC2.
Expects that you have generated the AnnDataset using
5_create_processed_anndataset.ipynb or some other way.

The script performs the following steps:
    1. Reading the h5ads file.
    2. Feature selection.
    3. Spectral embedding.
    4. UMAP.
    5. Clustering.
    6. Plotting.

Usage:
    python 6_analyze_anndata.py -i <sinput_h5ads> -o <output_directory> [options]

Arguments:
    -i/--input       Path to the input h5ads file.
    -o/--output      Path to the output directory.
    [options]        Additional optional parameters.

Author: Adam Klie
"""

# Imports
import os
import sys
import time
import logging
import argparse
import random


def main(args):

    # Parse args
    input_h5ads = args.input_h5ads
    output_dir = args.output_dir
    n_features = args.n_features
    make_gene_matrix = args.make_gene_matrix
    filter_genes = args.filter_genes

    # Set up logging
    sys.path.append("/cellar/users/aklie/opt/igvf-ucsd/single_cell_utilities")
    from utils import make_dirs
    make_dirs(output_dir)
    make_dirs(os.path.join(output_dir, "logs"))
    time_id = time.strftime("%Y%m%d-%H%M%S")
    run =  random.getrandbits(128)
    run_id = time_id + "_" + str(run) +  "_analyze_anndataset"
    log_file = os.path.join(output_dir, "logs", run_id + '.log')
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info(f"Arguments: {args}")
    import snapatac2 as snap
    logging.info("SnapATAC version: " + snap.__version__)

    # Check and create output directory
    logging.info("Making output directory if doesn't already exist: " + output_dir)
    make_dirs(output_dir)

    # Load in the AnnDataset
    logging.info("Processing file: " + input_h5ads)
    adata_atac_processed = snap.read_dataset(input_h5ads)

    # Select features on the merged dataset
    logging.info(f"Selecting {n_features} features.")
    time_in = time.time()
    snap.pp.select_features(adata_atac_processed, n_features=n_features)
    time_out = time.time()
    logging.info("Feature selection took " + str(time_out - time_in) + " seconds")

    # Spectral embedding
    logging.info("Performing spectral embedding")
    time_in = time.time()
    snap.tl.spectral(adata_atac_processed)
    time_out = time.time()
    logging.info("Spectral embedding took " + str(time_out - time_in) + " seconds")

    # UMAP
    logging.info("Performing UMAP")
    time_in = time.time()
    snap.tl.umap(adata_atac_processed)
    time_out = time.time()
    logging.info("UMAP took " + str(time_out - time_in) + " seconds")

    # Clustering
    logging.info("Performing clustering")
    time_in = time.time()
    snap.pp.knn(adata_atac_processed, use_rep="X_spectral")
    snap.tl.leiden(adata_atac_processed)
    time_out = time.time()
    logging.info("Clustering performed in " + str(time_out - time_in) + " seconds")

    # Plotting
    snap.pl.umap(adata_atac_processed, color="sample", interactive=False, out_file=input_h5ads.replace(".h5ads", "_sample_umap.png"))
    snap.pl.umap(adata_atac_processed, color="leiden", interactive=False, out_file=input_h5ads.replace(".h5ads", "_leiden_umap.png"))

    if make_gene_matrix:
        import scanpy as sc

        # Create gene matrix
        logging.info("Creating gene matrix")
        time_in = time.time()
        gene_matrix = snap.pp.make_gene_matrix(
            adata=adata_atac_processed,
            gene_anno=snap.genome.hg38
        )
        time_out = time.time()
        logging.info("Gene matrix created in " + str(time_out - time_in) + " seconds")

        # Clean up the gene matrix
        sc.pp.filter_genes(gene_matrix, min_cells=filter_genes)
        sc.pp.normalize_total(gene_matrix)
        sc.pp.log1p(gene_matrix)
        
        # Perform imputation with MAGIC
        logging.info("Performing MAGIC imputation")
        time_in = time.time()
        sc.external.pp.magic(gene_matrix, solver="approximate")
        time_out = time.time()
        logging.info("MAGIC imputation took " + str(time_out - time_in) + " seconds")

        # Copy over UMAP embedding
        gene_matrix.obsm["X_umap"] = adata_atac_processed.obsm["X_umap"]

        # Write the gene matrix
        gene_matrix.write(input_h5ads.replace(".h5ads", "_gene_matrix.h5ad"))

    # Close the file
    adata_atac_processed.close()


if __name__ == "__main__":
    # Setting up argparse
    parser = argparse.ArgumentParser(description="analyze a single h5ads sample using snapatac2.")
    parser.add_argument('-i', '--input_h5ads', required=True, help="Path to the input h5ads file.")
    parser.add_argument('-o', '--output_dir', required=True, help="Path to the output directory.")

    # Additional optional arguments
    parser.add_argument('-n', '--n_features', type=int, default=10000, help="Number of features to select.")
    parser.add_argument('--make_gene_matrix', action="store_true", help="Whether to make a gene matrix.")
    parser.add_argument('--filter_genes', type=int, default=3, help="Minimum number of cells expressing a gene to keep it.")

    # Parse args
    args = parser.parse_args()

    # Call main
    main(args)
