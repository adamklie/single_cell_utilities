#!/usr/bin/env python3
"""
scRNA-seq Analysis Script

This script performs cell type annotation on scRNA-seq data using decoupler.

Usage:
    python script_name.py -i <input_h5ad_path>

Arguments:
    -i/--input_h5ad_path       Path to the input h5ad file.
    [options]                  Additional optional parameters.

Author: [Your Name] (last modified: [Date])
"""

import argparse
import logging
import os
import sys
import yaml


def main(args):
    # Parse args
    input_h5ad_path = args.input_h5ad_path
    outdir_path = args.outdir_path
    annotation_resource = args.annotation_resource
    groupby_key = args.groupby_key
    normalize_data = args.normalize_data
    method = args.method
    obsm_key = args.obsm_key
    plot_cell_types = args.plot_cell_types
    n_dotplot = args.n_dotplot

    # Make 
    if not os.path.exists(outdir_path):
        os.makedirs(outdir_path)

    # Set up logging
    logging.basicConfig(
        filename=os.path.join(outdir_path, "cell_type_annotation.log"),
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    # Make and log params
    import scanpy as sc
    import decoupler as dc
    data_params = {
        "input_h5ad_path": input_h5ad_path,
        "outdir_path": outdir_path,
    }
    decoupler_params = {
        "annotation_resource": annotation_resource,
        "normalize_data": normalize_data,
        "method": method,
    }
    version_params = {
        "Python": sys.version[:5],
        "scanpy": sc.__version__,
        "decoupler": dc.__version__,
    }
    params = {"data": data_params, "run": decoupler_params, "versions": version_params}
    if not os.path.exists(os.path.join(outdir_path, "cell_type_annotation_params.yaml")):
        logging.info("Writing params to {}".format(outdir_path))
        with open(
            os.path.join(outdir_path, "cell_type_annotation_params.yaml"), "w"
        ) as f:
            yaml.dump(params, f)
    else:
        logging.info("params.yaml already exists, will not overwrite")

    # The data to load in is a scanpy object
    logging.info("Loading in data from {}".format(input_h5ad_path))
    adata = sc.read_h5ad(input_h5ad_path)
    logging.info("Data shape: {}".format(adata.shape))

    # Query Omnipath and get PanglaoDB
    logging.info("Querying Omnipath for {}".format(annotation_resource))
    markers = dc.get_resource("PanglaoDB")
    if markers:
        logging.info("Resource found, filtering to canonical human markers and saving")
        markers = markers[
            (markers["human"] == "True") & (markers["canonical_marker"] == "True")
        ]
        markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]
        markers.to_csv(os.path.join(outdir_path, "markers_used.tsv"), sep="\t")

    # Get a copy of the adata_pp
    adata_pp = adata_pp.copy()

    # Basic filtering
    logging.info("Running basic gene filtering")
    logging.info(f"Gene count before: {adata_pp.n_vars}")
    sc.pp.filter_cells(adata_pp, min_genes=200)
    sc.pp.filter_genes(adata_pp, min_cells=3)
    logging.info(f"Gene count after: {adata_pp.n_vars}")

    if normalize_data:
        logging.info("Normalizing data")
        sc.pp.normalize_total(adata_pp, target_sum=1e4)
        sc.pp.log1p(adata_pp)
        adata_pp.raw = adata_pp

    # Run ORA on top 5% of genes in each cell
    logging.info("Running ORA on top 5% of genes in each cell")
    dc.run_ora(
        mat=adata_pp,
        net=markers,
        source="cell_type",
        target="genesymbol",
        min_n=3,
        use_raw=True,
        verbose=True,
    )

    # Basically just pulls out the above obsm into a new AnnDatas.X and copies the rest of it
    import numpy as np
    logging.info("Pulling out ORA activity results")
    acts = dc.get_acts(adata_pp, obsm_key="ora_estimate")
    acts_v = acts.X.ravel()
    max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
    acts.X[~np.isfinite(acts.X)] = max_e

    # We can scale the obtained activities for better visualizations
    logging.info("Scaling the obtained activities for better visualizations")
    sc.pp.scale(acts)

    # Plot a few familiar cell types
    import matplotlib.pyplot as plt
    with plt.rc_context():
        sc.pl.embedding(
            acts,
            basis=obsm_key,
            color=plot_cell_types,
            cmap="RdBu_r",
            vmin=-1,
            vmax=3,
            show=False,
        )
        plt.savefig(os.path.join(outdir_path, "selected_cell_type_annotation_umap.png"))
        plt.show()
        plt.close()


    # Plot the violin plots
    if args.groupby_key is not None:
        with plt.rc_context():
            sc.pl.violin(
                acts,
                keys=plot_cell_types,
                groupby=groupby_key,
                multi_panel=True,
                rotation=90,
                show=False,
            )
            plt.savefig(os.path.join(outdir_path, "selected_cell_type_annotation_violin.png"))
            plt.show()
            plt.close()

        # Does just that
        df = dc.rank_sources_groups(
            acts,
            groupby=groupby_key,
            reference="rest",
            method="t-test_overestim_var",
        )

    # Save the annotation results
    df.to_csv(os.path.join(outdir_path, "cell_type_annotation_by_groups.tsv"), sep="\t")

    # We can grab a dictionary of the top X for each group
    n_ctypes = n_dotplot
    ctypes_dict = (
        df.groupby("group")
        .head(n_ctypes)
        .groupby("group")["names"]
        .apply(lambda x: list(x))
        .to_dict()
    )
    ctypes_dict

    # And plot that as a matrixplot
    with plt.rc_context():
        sc.pl.matrixplot(
            acts,
            ctypes_dict,
            groupby_key,
            dendrogram=True,
            colorbar_title="Z-scaled scores",
            vmin=-1,
            vmax=1,
            cmap="RdBu_r",
            show=False,
        )
        plt.savefig(
            os.path.join(outdir_path, f"{n_ctypes}_cell_type_annotation_matrixplot.png")
        )
        plt.show()



if __name__ == "__main__":
    # Setting up argparse
    parser = argparse.ArgumentParser(
        description="Cell type annotation on scRNA-seq data using decoupler."
    )
    parser.add_argument(
        "-i", "--input_h5ad_path", required=True, help="Path to the input h5ad file."
    )
    parser.add_argument(
        "-o", "--outdir_path", required=True, help="Path to the output directory."
    )
    parser.add_argument(
        "--annotation_resource",
        type=str,
        default="PanglaoDB",
        help="The resource used for annotation.",
    )
    parser.add_argument(
        "--normalize_data",
        type=bool,
        default=True,
        help="Whether to normalize the data.",
    )
    parser.add_argument(
        "--method", type=str, default="ora", help="The method used for analysis."
    )
    parser.add_argument(
        "--obsm_key",
        type=str,
        default="X_umap",
        help="The key for the obsm to use for UMAP plotting.",
    )
    parser.add_argument(
        "--groupby_key",
        type=str,
        default="integrated_manual_cellid_annotation",
        help="The key for the annotation to use for grouping for plotting "
        "and differential analysis by groups."
    )
    parser.add_argument(
        "--plot_cell_types",
        nargs="+",
        type=str,
        default="Alpha cells Beta cells Delta cells Enterochromaffin cells",
        help="The cell types to plot.",
    )
    parser.add_argument(
        "--n_dotplot",
        type=int,
        default=5,
        help="The number of cell types to plot in the dotplot.",
    )

    # Parse the arguments
    args = parser.parse_args()

    # Call the main function
    main(args)
