#!/usr/bin/env python3
import argparse
import logging
import random
import hashlib
import os
import sys
import yaml
sys.path.append("/cellar/users/aklie/projects/igvf/single_cell_utilities")
from process_data.normalize import normalize_tpm


def main(args):
    
    # Parse args
    path_h5ad = args.path_h5ad
    layer = args.layer
    if layer == "None":
        layer = None
    path_gene_lengths = args.path_gene_lengths
    if path_gene_lengths == "None":
        path_gene_lengths = None
    path_out = args.path_out
    groupby_keys = args.groupby_keys
    cellid_key = args.cellid_key
    if cellid_key == "None":
        cellid_key = None
    compare_key = args.compare_key
    if compare_key == "None":
        compare_key = None
    target_max_cells_per_pb = args.target_max_cells_per_pb
    if target_max_cells_per_pb == "None":
        target_max_cells_per_pb = None
    else:
        target_max_cells_per_pb = int(target_max_cells_per_pb)
    obs_key = args.obs_key
    mode = args.mode
    min_cells = args.min_cells
    min_counts = args.min_counts
    random_state = args.random_state

    # Make output directory if it doesn't exist
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    # Make a random hexidecimanl workflow hash
    workflow_hash = hashlib.sha1()
    workflow_hash.update(str(random.getrandbits(128)).encode("utf-8"))
    
    # Log file
    if os.path.exists(os.path.join(path_out, "pseudobulk.log")):
        os.remove(os.path.join(path_out, "pseudobulk.log"))
    logging.basicConfig(
        filename=os.path.join(path_out, "pseudobulk.log"),
        level=logging.INFO, 
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logging.info(f"Worflow hash: {workflow_hash.hexdigest()}")

    # Make and log params
    import scanpy as sc
    import decoupler as dc
    data_params = {
        "path_h5ad": path_h5ad,
        "layer": layer,
        "path_gene_lengths": path_gene_lengths,
        "path_out": path_out,
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
    logging.info("Writing params to {}".format(os.path.join(path_out, "pseudobulk_params.yaml")))
    with open(os.path.join(path_out, "pseudobulk_params.yaml"), "w") as outfile:
        yaml.dump(params, outfile, default_flow_style=False)

    # The data to load in is a scanpy object
    logging.info("Loading in data from {}".format(path_h5ad))
    adata = sc.read_h5ad(path_h5ad)
    logging.info("Data shape: {}".format(adata.shape))

    # Get a copy of the adata_pp
    adata_pp = adata.copy()

    # Basic filtering
    logging.info("Running basic feature filtering")
    logging.info(f"Feature count before: {adata_pp.n_vars}")
    sc.pp.filter_cells(adata_pp, min_genes=200)
    sc.pp.filter_genes(adata_pp, min_cells=3)
    logging.info(f"Feature count after: {adata_pp.n_vars}")

    # Verify that this is counts data
    import numpy as np
    if layer is not None:
        test_data = adata_pp.layers[layer][:10, :10].todense()
    else:
        test_data = adata_pp.X[:10, :10].todense()
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
        random_state=random_state
    )
    print(adata_pp.obs[obs_key].value_counts())

    # No column should have multiple non-zero values across all rows
    import pandas as pd
    crosstab = pd.crosstab(adata_pp.obs["condition"], adata_pp.obs[obs_key])
    if (crosstab > 0).sum(axis=0).max() > 1:
        logging.info("There are pseudobulks with cells from multiple conditions.")

    # Create a dir for the pseudobulk h5ads
    if not os.path.exists(os.path.join(path_out, "pdatas")):
        os.makedirs(os.path.join(path_out, "pdatas"))

    # Get pseudo-bulk profile WITHOUT any filtering
    logging.info("Getting pseudo-bulk profile WITHOUT any filtering")
    pdata = dc.get_pseudobulk(
        adata_pp,
        sample_col=obs_key,
        groups_col=None,
        layer=layer,
        mode=mode,
        min_cells=0,
        min_counts=0
    )

    # Plot some QCs on the pseudo-bulk data
    dc.plot_psbulk_samples(
        pdata, 
        groupby=groupby_keys,
        figsize=(14, 6), 
        save=os.path.join(path_out, "pseudobulk_samples_plot.png")
    )

    # Save the non-filtered pseudo-bulk adata
    logging.info("Saving the non-filtered pseudo-bulk adata")
    from utils import convert_object_in_obs
    convert_object_in_obs(pdata)
    pdata.write(os.path.join(path_out, "pdatas", "pseudobulk_no_filter.h5ad"))
    pdata_df = pdata.to_df().T
    pdata_df.to_csv(os.path.join(path_out, "pdatas", "pseudobulk_no_filter.tsv"), sep="\t", index=True)

    # Normalize the pseudo-bulk data with TPM and save
    if path_gene_lengths:
        gene_lengths = pd.read_csv(path_gene_lengths, sep="\t", index_col=0)
        gene_lengths = gene_lengths.loc[[gene_name for gene_name in pdata.var_names if gene_name in gene_lengths.index]]
        pdata_norm_df = pdata_df.loc[gene_lengths.index]
        pdata_norm_df = normalize_tpm(pdata_norm_df, gene_lengths["length"])
        pdata_norm_df.to_csv(os.path.join(path_out, "pdatas", "pseudobulk_no_filter_tpm.tsv"), sep="\t", index=True)

    # Get pseudo-bulk profile WITH filtering
    logging.info("Getting pseudo-bulk profile WITH filtering")
    pdata = dc.get_pseudobulk(
        adata_pp,
        sample_col=obs_key,
        groups_col=None,
        layer=layer,
        mode=mode,
        min_cells=min_cells,
        min_counts=min_counts
    )

    # Plot some QCs on the pseudo-bulk data
    dc.plot_psbulk_samples(
        pdata, 
        groupby=groupby_keys,
        figsize=(14, 6),
        save=os.path.join(path_out, "pseudobulk_filtered_samples_plot.png")
    )

    # Save the filtered pseudo-bulk adata
    logging.info("Saving the filtered pseudo-bulk adata")
    convert_object_in_obs(pdata)
    pdata.write(os.path.join(path_out, "pdatas", "pseudobulk_filter.h5ad"))
    pdata_df = pdata.to_df().T
    pdata_df.to_csv(os.path.join(path_out, "pdatas", "pseudobulk_filter.tsv"), sep="\t", index=True)

    # Normalize the pseudo-bulk data with TPM and save
    if path_gene_lengths:
        gene_lengths = pd.read_csv(path_gene_lengths, sep="\t", index_col=0)
        gene_lengths = gene_lengths.loc[[gene_name for gene_name in pdata.var_names if gene_name in gene_lengths.index]]
        pdata_norm_df = pdata_df.loc[gene_lengths.index]
        pdata_norm_df = normalize_tpm(pdata_norm_df, gene_lengths["length"])
        pdata_norm_df.to_csv(os.path.join(path_out, "pdatas", "pseudobulk_filter_tpm.tsv"), sep="\t", index=True)

    # Get a copy to process
    pp_pdata = pdata.copy()

    # CPM normalize, log(1+p) transformation, scale and PCA
    sc.pp.normalize_total(pp_pdata, target_sum=1e6)
    sc.pp.log1p(pp_pdata)
    sc.pp.scale(pp_pdata, max_value=10)
    sc.tl.pca(pp_pdata, n_comps=10)

    # Plot the PCA with covariates
    import matplotlib.pyplot as plt
    with plt.rc_context({"figure.figsize": (4, 4)}):
        sc.pl.pca(pp_pdata, color=groupby_keys + ['psbulk_counts'], ncols=2, size=50, show=False)
        plt.savefig(os.path.join(path_out, "pca_w_covariates.png"))
        plt.tight_layout()
        plt.show()

    # Get correlations with PCs
    for groupby_key in groupby_keys:
        pp_pdata.obs[groupby_key] = pp_pdata.obs[groupby_key].str.replace(".", "_")
    try:
        dc.get_metadata_associations(
            pp_pdata,
            obs_keys = groupby_keys + ['psbulk_n_cells', 'psbulk_counts'],  # metadata columns to associate to PCs
            obsm_key='X_pca',  # where the PCs are stored
            uns_key='pca_anova',  # where the results are stored
            inplace=True,
        )
        plt.figure(figsize=(7,10))
        dc.plot_associations(
            pp_pdata,
            uns_key='pca_anova',  # summary statistics from the anova tests
            obsm_key='X_pca',  # where the PCs are stored
            stat_col='p_adj',  # which summary statistic to plot
            obs_annotation_cols = groupby_keys,
            titles=['Adjusted p-values from ANOVA', 'Principle component scores'],  # which sample annotations to plot
            cmap_cats='tab20',
        )
        plt.savefig(os.path.join(path_out, "ps_pc_assocation_anova.png"), bbox_inches='tight')
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
            logging.info(f"Features: {pdata_deseq.n_vars}")
            logging.info(f"Samples: {pdata_deseq.n_obs}")
            logging.info("Samples in each group:")
            logging.info(pdata_deseq.obs[compare_key].value_counts().to_string())

            # Save the object
            pdata_deseq.write(os.path.join(path_out, "pdatas", f"pseudobulk_{cellid}_filtered_genes.h5ad"))
            pdata_deseq_df = pdata_deseq.to_df().T
            pdata_deseq_df.to_csv(os.path.join(path_out, "pdatas", f"pseudobulk_{cellid}_filtered_genes.tsv"), sep="\t", index=True)

            # Normalize the pseudo-bulk data with TPM and save
            if path_gene_lengths:
                gene_lengths = gene_lengths.loc[[gene_name for gene_name in pdata_deseq.var_names if gene_name in gene_lengths.index]]
                pdata_deseq_norm_df = pdata_deseq_df.loc[gene_lengths.index]
                pdata_deseq_norm_df = normalize_tpm(pdata_deseq_norm_df, gene_lengths["length"])
                pdata_deseq_norm_df.to_csv(os.path.join(path_out, "pdatas", f"pseudobulk_{cellid}_filtered_genes_tpm.tsv"), sep="\t", index=True)

    logging.info("Completed pseudobulking")


if __name__ == "__main__":
    # Parse args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--path_h5ad", type=str, required=True, help="Path to the input h5ad file."
    )
    parser.add_argument(
        "--layer", type=str, required=False, default="counts", help="Layer to use."
    )
    parser.add_argument(
        "--path_gene_lengths", type=str, required=False, default=None, help="Path to the gene lengths file."
    )
    parser.add_argument(
        "-o", "--path_out", type=str, required=True, help="Path to the output directory."
    )
    parser.add_argument(
        "--groupby_keys", type=str, required=True, nargs="+", help="Keys to groupby."
    )
    parser.add_argument(
        "--cellid_key", type=str, required=False, default=None, help="Key to groupby for cell id."
    )
    parser.add_argument(
        "--compare_key", type=str, required=False, default=None, help="Key to groupby for comparison."
    )
    parser.add_argument(
        "--target_max_cells_per_pb", required=False, default=None, help="Target max cells per pseudobulk."
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
    parser.add_argument(
        "--random_state", type=int, required=False, default=0, help="Random state for reproducibility."
    )
    args = parser.parse_args()

    # Run main
    main(args)
