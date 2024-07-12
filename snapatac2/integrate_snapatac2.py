#!/usr/bin/env python3
"""
integrate_snapatac2.py

This script creates AnnDatasets from multiple h5ad files and performs
feature selection, spectral embedding, clustering, and plotting using
SnapATAC2.

Usage:
    python integrate_snapatac2.py -i <input_h5ad_paths> -o <outdir_path> [options]

Author: Adam Klie
"""

# Imports
import os
import time
import logging
import argparse
import random
import hashlib


def main(args):

    # Parse args
    input_h5ad_paths = args.input_h5ad_paths
    sample_ids = args.sample_ids
    outdir_path = args.outdir_path
    output_prefix = args.output_prefix
    barcodes_path = args.barcodes_path
    n_features = args.n_features
    clustering_resolution = args.clustering_resolution
    make_gene_matrix = args.make_gene_matrix
    filter_genes = args.filter_genes
    plot_vars = args.plot_vars

    # Make output directory if doesn't already exist
    if not os.path.exists(outdir_path):
        os.makedirs(outdir_path)

    # Set up logging
    log_file = os.path.join(outdir_path, "integrate.log")
    if os.path.exists(log_file):
        os.remove(log_file)
    random_id = hashlib.md5(str(random.getrandbits(128)).encode()).hexdigest()
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    import snapatac2 as snap
    import scanpy as sc
    import numpy as np
    import pandas as pd
    logging.info("Run hash: " + random_id)
    logging.info("SnapATAC version: " + snap.__version__)
    logging.info("ScanPy version: " + sc.__version__)
    logging.info(f"Arguments: {args}")

    # If sample ids are not provided, use the file names
    if sample_ids is None:
        logging.info("Sample ids not provided. Using file names.")
        sample_ids = [os.path.basename(file).split(".")[0] for file in input_h5ad_paths]
    
    # Read in each h5ad with ScanPy, delete the X_spectral and X_umap from the obsm, and resave with new name
    logging.info("Deleting X_spectral and X_umap from obsm of AnnData and resaving.")
    cell_bcs = []
    plot_var_dict = {}
    for path in input_h5ad_paths:
        adata = sc.read_h5ad(path)
        del adata.obsm["X_spectral"]
        del adata.obsm["X_umap"]
        adata.write_h5ad(path.replace(".h5ad", "_obsm_delete.h5ad"))
        cell_bcs.extend(adata.obs_names.tolist())
        for var in plot_vars:
            if var in adata.obs.columns:
                if var in plot_var_dict:
                    plot_var_dict[var].extend(adata.obs[var].tolist())
                else:
                    plot_var_dict[var] = adata.obs[var].tolist()
            else:
                logging.warning(f"Variable {var} not found in obs. Setting values for these cells to NaN.")
                if var in plot_var_dict:
                    plot_var_dict[var].extend([np.nan] * adata.shape[0])
                else:
                    plot_var_dict[var] = [np.nan] * adata.shape[0]
    for key in plot_var_dict:
        plot_var_dict[key] = pd.Series(plot_var_dict[key], index=cell_bcs)

    # Update the input h5ad paths to the new ones
    input_h5ad_paths = [file.replace(".h5ad", "_obsm_delete.h5ad") for file in input_h5ad_paths]
    logging.info(f"Updated input h5ad paths: {input_h5ad_paths}")
    
    # Read in barcodes file if provided
    if barcodes_path is not None:
        logging.info(f"Reading in passed in barcodes file from {barcodes_path}")
        if barcodes_path.endswith(".csv"):
            barcodes = pd.read_csv(barcodes_path, header=None, index_col=0)
        elif barcodes_path.endswith(".tsv") | barcodes_path.endswith(".txt"):
            barcodes = pd.read_csv(barcodes_path, header=None, index_col=0, sep="\t")
        else:
            raise ValueError("Barcodes file must be a .csv or .tsv file.")
        barcodes = barcodes.index.tolist()
        logging.info(f"Barcodes file read in with {len(barcodes)} barcodes.")
        logging.info(f"First few barcodes: {barcodes[:5]}")

    # Create the AnnDataset
    adata_atac_merged_list = []
    for _, h5ad_file in enumerate(input_h5ad_paths):
        adata_atac = snap.read(h5ad_file)
        if barcodes_path is not None:
            logging.info(f"Subsetting AnnData to barcodes in {barcodes_path}")
            adata_atac.subset(obs_indices=np.where(pd.Index(adata_atac.obs_names).isin(barcodes))[0])
            logging.info(f"Number of cells after subset: {adata_atac.shape[0]}")
        adata_atac_merged_list.append(adata_atac)
    adatas = [(name, adata) for name, adata in zip(sample_ids, adata_atac_merged_list)]

    # Merge into one object
    logging.info(f"Creating AnnDataset from {adatas} samples.")
    adata_atac_merged = snap.AnnDataSet(
        adatas=adatas,
        filename=os.path.join(outdir_path, f"{output_prefix}.h5ads")
    )
    logging.info(f"AnnDataset created at {os.path.join(outdir_path, f'{output_prefix}.h5ads')}")

    # Close all the backed anndatas
    logging.info("Closing all the backed anndatas and the AnnDataset.")
    for adata_atac in adata_atac_merged_list:
        adata_atac.close()
    adata_atac_merged.close()

    # Read in the merged AnnDataset
    logging.info(f"Reading in the merged AnnDataset from {os.path.join(outdir_path, f'{output_prefix}.h5ads')}")
    adata_atac_merged = snap.read_dataset(os.path.join(outdir_path, f"{output_prefix}.h5ads"))
    logging.info(f"First few merged adata cell ids: {adata_atac_merged.obs_names[:5]}")
    
    # Add the plot vars to the merged AnnDataset
    logging.info("Adding plot vars to the merged AnnDataset.")
    for key in plot_var_dict:
        adata_atac_merged.obs[key] = plot_var_dict[key][adata_atac_merged.obs_names].values.tolist()
        logging.info(f"Added {key} to obs.")
    
    # Select features on the merged dataset
    logging.info(f"Selecting {n_features} features.")
    time_in = time.time()
    snap.pp.select_features(adata_atac_merged, n_features=n_features)
    time_out = time.time()
    logging.info("Feature selection took " + str(time_out - time_in) + " seconds")

    # Spectral embedding
    logging.info("Performing spectral embedding")
    time_in = time.time()
    snap.tl.spectral(adata_atac_merged)
    time_out = time.time()
    logging.info("Spectral embedding took " + str(time_out - time_in) + " seconds")

    # Run UMAP
    logging.info("Running UMAP")
    time_in = time.time()
    snap.tl.umap(adata_atac_merged, use_dims=list(range(1, adata_atac_merged.obsm["X_spectral"].shape[1])))
    time_out = time.time()
    logging.info("UMAP took " + str(time_out - time_in) + " seconds")

    # Clustering
    logging.info("Performing clustering")
    time_in = time.time()
    snap.pp.knn(adata_atac_merged, use_rep="X_spectral", use_dims=list(range(1, adata_atac_merged.obsm["X_spectral"].shape[1])))
    snap.tl.leiden(adata_atac_merged, resolution=clustering_resolution, key_added=f"leiden_{clustering_resolution}")
    time_out = time.time()
    logging.info("Clustering performed in " + str(time_out - time_in) + " seconds")

    # Make gene matrix if requested
    if make_gene_matrix:

        # Create gene matrix
        logging.info("Creating gene matrix")
        time_in = time.time()
        gene_matrix = snap.pp.make_gene_matrix(
            adata=adata_atac_merged,
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
        gene_matrix.obsm["X_umap"] = adata_atac_merged.obsm["X_umap"]

        # Write the gene matrix
        logging.info("Writing gene matrix")
        time_in = time.time()
        gene_matrix.write(os.path.join(outdir_path, f"{output_prefix}_gene_matrix.h5ad"))
        time_out = time.time()
        logging.info("Gene matrix written in " + str(time_out - time_in) + " seconds")

    # Turn into in memory AnnData
    logging.info("Turning into in memory AnnData")
    adata_mem = adata_atac_merged.to_adata()
    
    # Plot first spectral embedding against log_n_fragment
    logging.info("Plotting spectral embedding")
    import matplotlib.pyplot as plt
    if "log_n_fragment" in adata_mem.obs.columns:
        with plt.rc_context({"figure.figsize": (5, 5)}):
            sc.pl.embedding(basis="X_spectral", adata=adata_mem, color="log_n_fragment", show=False)
            plt.savefig(os.path.join(outdir_path, "spectral_embedding.png"))
            plt.close()

    # Plot the UMAP with clusters
    logging.info("Plotting UMAP")
    umap_plot_vars = ["sample", f"leiden_{clustering_resolution}"]
    umap_plot_vars.extend(list(plot_var_dict.keys()))
    with plt.rc_context({"figure.figsize": (5, 5)}):
        sc.pl.umap(adata_mem, color=umap_plot_vars, show=False, ncols=2, wspace=0.2)
        plt.savefig(os.path.join(outdir_path, "umap.png"))
        plt.close()

    # Close the file
    adata_atac_merged.close()

    # Print completion message
    logging.info("Successfully completed integration of SnapATAC2 samples.")


if __name__ == "__main__":
    # Setting up argparse
    parser = argparse.ArgumentParser(description="analyze a single h5ads sample using snapatac2.")
    parser.add_argument('--input_h5ad_paths', nargs="+", required=True, help="Paths to input h5ad files. Space separated. Must be in the same order as sample_ids.")
    parser.add_argument('--sample_ids', nargs="+", required=False, default=None, help="Sample ids for each h5ad file. Space separated. If not provided, will use the file prefixes.")
    parser.add_argument('--outdir_path', required=True, help="Path to the output directory to save all results.")
    parser.add_argument('--output_prefix', required=False, default="merged", help="Prefix for the output h5ads.")
    parser.add_argument('--barcodes_path', required=False, default=None, help="Path to a file with barcodes in the first column. Used to subset the AnnDataset to cells of interest.")
    parser.add_argument('--n_features', type=int, default=50000, help="Number of features to select prior to dimensionality reduction with spectral embedding.")
    parser.add_argument('--clustering_resolution', type=float, default=1, help="Clustering resolution for Leiden clustering.")
    parser.add_argument('--make_gene_matrix', action="store_true", help="Whether to make a gene matrix.")
    parser.add_argument('--filter_genes', type=int, default=3, help="Minimum number of cells expressing a gene to keep it in gene matrix.")
    parser.add_argument('--plot_vars', nargs="+", required=False, default=[], help="Variables to plot on UMAP. Space separated. Must be in the obs of each AnnData.")

    # Parse args
    args = parser.parse_args()

    # Call main
    main(args)
