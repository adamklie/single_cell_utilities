#!/usr/bin/env python3
"""
scRNA-seq Analysis Script

This script performs

Usage:
    python script_name.py -i <input_h5ad_path>

Arguments:
    -i/--input_h5ad_path       Path to the input h5ad file.
    [options]                  Additional optional parameters.

Author: [Your Name] (last modified: [Date])
"""

import argparse
import logging
import random
import hashlib
import os
import sys
import yaml


def main(args):
    # Parse args
    input_h5ad_path = args.input_h5ad_path
    metadata_path = args.metadata_path
    h5ad_sample_key = args.h5ad_sample_key
    metadata_sample_key = args.metadata_sample_key
    outdir_path = args.outdir_path
    groupby_keys = args.groupby_keys
    cellid_key = args.cellid_key
    compare_key = args.compare_key
    target_max_cells_per_pb = args.target_max_cells_per_pb
    obs_key = args.obs_key
    mode = args.mode
    min_cells = args.min_cells
    min_counts = args.min_counts    

    # Make output directory if it doesn't exist
    if not os.path.exists(outdir_path):
        os.makedirs(outdir_path)

    # Make a random hexidecimanl workflow hash
    workflow_hash = hashlib.sha1()
    workflow_hash.update(str(random.getrandbits(128)).encode("utf-8"))
    
    # Log file
    if os.path.exists(os.path.join(outdir_path, "pseudobulk.log")):
        os.remove(os.path.join(outdir_path, "pseudobulk.log"))
    logging.basicConfig(
        filename=os.path.join(outdir_path, "pseudobulk.log"),
        level=logging.INFO, 
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logging.info(f"Worflow hash: {workflow_hash.hexdigest()}")

    # Make and log params
    import scanpy as sc
    import decoupler as dc
    data_params = {
        "input_h5ad_path": input_h5ad_path,
        "outdir_path": outdir_path,
        "metadata_path": metadata_path,
    }
    method_params = {
        "groupby_keys": groupby_keys,
        "cellid_key": cellid_key,
        "compare_key": compare_key,
        "target_max_cells_per_pb": target_max_cells_per_pb,
        "obs_key": obs_key,
        "mode": mode,
        "min_cells": min_cells,
        "min_counts": min_counts,
    }
    version_params = {
        "workflow_hash": workflow_hash.hexdigest(),
        "Python": sys.version[:5],
        "Scanpy": sc.__version__,
        "decoupler": dc.__version__,
    }

    params = {"data": data_params, "run": method_params, "versions": version_params}
    if not os.path.exists(os.path.join(outdir_path, "pseudobulk_params.yaml")):
        logging.info("Writing params to {}".format(os.path.join(outdir_path, "pseudobulk_params.yaml")))
        with open(
            os.path.join(outdir_path, "pseudobulk_params.yaml"), "w"
        ) as outfile:
            yaml.dump(params, outfile, default_flow_style=False)
    else:
        logging.info("params.yaml already exists, will not overwrite")

    # The data to load in is a scanpy object
    logging.info("Loading in data from {}".format(input_h5ad_path))
    adata = sc.read_h5ad(input_h5ad_path)
    logging.info("Data shape: {}".format(adata.shape))

    # Load and add sample metadata
    import pandas as pd
    if metadata_path:
        # Add sample metadata to adata
        sample_metadata = pd.read_csv(metadata_path, sep="\t")
        new_obs = adata.obs.merge(sample_metadata[[metadata_sample_key, "batch", "timepoint", "condition"]], left_on=h5ad_sample_key, right_on=metadata_sample_key, how="left")
        new_obs.index = adata.obs.index
        new_obs = new_obs.drop(columns=[metadata_sample_key])
        adata.obs = new_obs

    # Get a copy of the adata_pp
    adata_pp = adata.copy()

    # Basic filtering
    logging.info("Running basic gene filtering")
    logging.info(f"Gene count before: {adata_pp.n_vars}")
    sc.pp.filter_cells(adata_pp, min_genes=200)
    sc.pp.filter_genes(adata_pp, min_cells=3)
    logging.info(f"Gene count after: {adata_pp.n_vars}")

    # Store raw counts in layers
    adata_pp.layers['counts'] = adata_pp.X.copy()

    # Verify that this is counts data
    import numpy as np
    test_data = adata_pp.layers["counts"][:10, :10].todense()
    if np.all(test_data >= 0) and np.all(test_data.astype(int) == test_data):
        logging.info("The matrix contains count data.")
    else:
        logging.info("The matrix does not contain count data.")

    from utils import get_pseudobulk_groups
    get_pseudobulk_groups(
        adata_pp, 
        groupby_cols=groupby_keys,
        target_max_cells_per_pb=target_max_cells_per_pb,
        obs_key=obs_key,
        copy=False,
    )
    print(adata_pp.obs[obs_key].value_counts())
    # No column should have multiple non-zero values across all rows
    crosstab = pd.crosstab(adata_pp.obs["condition"], adata_pp.obs[obs_key])
    if (crosstab > 0).sum(axis=0).max() > 1:
        logging.info("There are pseudobulks with cells from multiple conditions.")

    # Create a dir for the pseudobulk h5ads
    if not os.path.exists(os.path.join(outdir_path, "pdatas")):
        os.makedirs(os.path.join(outdir_path, "pdatas"))

    # Get pseudo-bulk profile WITHOUT any filtering
    logging.info("Getting pseudo-bulk profile WITHOUT any filtering")
    pdata = dc.get_pseudobulk(
        adata_pp,
        sample_col=obs_key,
        groups_col=None,
        layer='counts',
        mode=mode,
        min_cells=0,
        min_counts=0
    )

    # Plot some QCs on the pseudo-bulk data
    dc.plot_psbulk_samples(
        pdata, 
        groupby=groupby_keys,
        figsize=(14, 6), 
        save=os.path.join(outdir_path, "pseudobulk_samples_plot.png")
    )

    # Save the non-filtered pseudo-bulk adata
    logging.info("Saving the non-filtered pseudo-bulk adata")
    from utils import convert_object_in_obs
    convert_object_in_obs(pdata)
    pdata.write(os.path.join(outdir_path, "pdatas", "pseudobulk_no_filter.h5ad"))
    pdata_df = pdata.to_df()
    pdata_df.to_csv(os.path.join(outdir_path, "pdatas", "pseudobulk_no_filter.tsv"), sep="\t", index=True)

    # Get pseudo-bulk profile WITH filtering
    logging.info("Getting pseudo-bulk profile WITH filtering")
    pdata = dc.get_pseudobulk(
        adata_pp,
        sample_col=obs_key,
        groups_col=None,
        layer='counts',
        mode=mode,
        min_cells=min_cells,
        min_counts=min_counts
    )

    # Plot some QCs on the pseudo-bulk data
    dc.plot_psbulk_samples(
        pdata, 
        groupby=groupby_keys,
        figsize=(14, 6),
        save=os.path.join(outdir_path, "pseudobulk_filtered_samples_plot.png")
    )

    # Save the filtered pseudo-bulk adata
    logging.info("Saving the filtered pseudo-bulk adata")
    convert_object_in_obs(pdata)
    pdata.write(os.path.join(outdir_path, "pdatas", "pseudobulk_filter.h5ad"))
    pdata_df = pdata.to_df()
    pdata_df.to_csv(os.path.join(outdir_path, "pdatas", "pseudobulk_filter.tsv"), sep="\t", index=True)

    # Get a copy to process
    pp_pdata = pdata.copy()

    # CPM normalize, log(1+p) transformation, scale and PCA
    sc.pp.normalize_total(pp_pdata, target_sum=1e6)
    sc.pp.log1p(pp_pdata)
    sc.pp.scale(pp_pdata, max_value=10)
    sc.tl.pca(pp_pdata, n_comps=10)

    # Plot the PCA with covariates
    import matplotlib.pyplot as plt
    with plt.rc_context():
        sc.pl.pca(pp_pdata, color=groupby_keys + ['psbulk_counts'], ncols=2, size=50, show=False)
        plt.savefig(os.path.join(outdir_path, "pca_w_covariates.png"))
        plt.tight_layout()
        plt.show()

    # Get correlations with PCs
    for groupby_key in groupby_keys:
        pp_pdata.obs[groupby_key] = pp_pdata.obs[groupby_key].str.replace(".", "_")
    try:
        dc.get_metadata_associations(
            pp_pdata,
            obs_keys =  groupby_keys + ['psbulk_n_cells', 'psbulk_counts'],  # metadata columns to associate to PCs
            obsm_key='X_pca',  # where the PCs are stored
            uns_key='pca_anova',  # where the results are stored
            inplace=True
        )

        plt.figure(figsize=(7,10))
        _, _ = dc.plot_associations(
            pp_pdata,
            uns_key='pca_anova',  # summary statistics from the anova tests
            obsm_key='X_pca',  # where the PCs are stored
            stat_col='p_adj',  # which summary statistic to plot
            obs_annotation_cols = groupby_keys,
            titles=['Adjusted p-values from ANOVA', 'Principle component scores'])  # which sample annotations to plot
        plt.savefig(os.path.join(outdir_path, "ps_pc_assocation_anova.png"), bbox_inches='tight')
        plt.show()
    except:
        logging.info("Could not get metadata associations")

    if cellid_key and compare_key:
        assert cellid_key in groupby_keys, f"cellid_key {cellid_key} must be in groupby_keys {groupby_keys}"
        assert compare_key in groupby_keys, f"compare_key {compare_key} must be in groupby_keys {groupby_keys}"
        logging.info(f"Creating pseudobulk h5ads for each {cellid_key} to make DE comparisons for {compare_key}")
        
        for cellid, cell_barcodes in pdata.obs.groupby(cellid_key).groups.items():
            logging.info(f"Selecting pseudobulks {cellid}")
            # Select profiles of interest
            pdata_deseq = pdata[cell_barcodes, :].copy()

            # Obtain genes that pass the thresholds
            genes = dc.filter_by_expr(pdata_deseq, group=compare_key, min_count=10, min_total_count=15)

            # Filter by these genes
            pdata_deseq = pdata_deseq[:, genes].copy()

            # Print the number of genes and samples
            logging.info(f"Group: {cellid}")
            logging.info(f"Genes: {pdata_deseq.n_vars}")
            logging.info(f"Samples: {pdata_deseq.n_obs}")
            logging.info("Samples in each group:")
            logging.info(pdata_deseq.obs[compare_key].value_counts().to_string())

            # Save the object
            pdata_deseq.write(os.path.join(outdir_path, "pdatas", f"pseudobulk_{cellid}_filtered_genes.h5ad"))
            pdata_deseq_df = pdata_deseq.to_df()
            pdata_deseq_df.to_csv(os.path.join(outdir_path, "pdatas", f"pseudobulk_{cellid}_filtered_genes.tsv"), sep="\t", index=True)

    logging.info("Completed")


if __name__ == "__main__":
    # Parse args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input_h5ad_path", type=str, required=True, help="Path to the input h5ad file."
    )
    parser.add_argument(
        "-o", "--outdir_path", type=str, required=True, help="Path to the output directory."
    )
    parser.add_argument(
        "--groupby_keys", type=str, required=True, nargs="+", help="Keys to groupby."
    )
    parser.add_argument(
        "--metadata_path", type=str, default=None, help="Path to the metadata file."
    )
    parser.add_argument(
        "--h5ad_sample_key", type=str, required=False, default="sample", help="Key to sample column in h5ad."
    )
    parser.add_argument(
        "--metadata_sample_key", type=str, required=False, default="sample_id", help="Key to sample column in metadata."
    )

    parser.add_argument(
        "--cellid_key", type=str, required=False, default=None, help="Key to groupby for cell id."
    )
    parser.add_argument(
        "--compare_key", type=str, required=False, default=None, help="Key to groupby for comparison."
    )
    parser.add_argument(
        "--target_max_cells_per_pb", type=int, required=False, default=None, help="Target max cells per pseudobulk."
    )
    parser.add_argument(
        "--obs_key", type=str, required=False, default="pseudobulk", help="Key to groupby for comparison."
    )
    parser.add_argument(
        "--mode", type=str, required=False, default="sum", help="Mode for pseudobulk creation."
    )
    parser.add_argument(
        "--min_cells", type=int, required=False, default=10, help="Minimum number of cells for filtering."
    )
    parser.add_argument(
        "--min_counts", type=int, required=False, default=1000, help="Minimum number of counts for filtering."
    )
    args = parser.parse_args()

    # Run main
    main(args)
