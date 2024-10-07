# This script is meant for running a DESeq2 analysis on an AnnData object

# This script uses argparse to parse arguments
# for specifying h5ad file, DESeq2 parameters, and output directory

# The script will:
# 1. Parse arguments and write them to a log file
# 1. Read in the h5ad file as an AnnData object
# 2. Build the DESeq2 objects
# 3. Run DESeq2
# 4. For each contrast
#  - Extract and save the results table
#  - Optionally shrink LFCs and save the results table
#  - Plot a volcano plot
# Usage:
# python 

import os
import logging
import hashlib
import random
import argparse

def main(args):
    
    # Parse args
    input_h5ad_path = args.input_h5ad_path
    outdir_path = args.outdir_path
    design_factors = args.design_factors  # Comma-separated list of factors
    continuous_factors = args.continuous_factors  # Comma-separated list of factors
    reference_factor = args.reference_factor  # String
    reference_value = args.reference_value # String
    refit_cooks = args.refit_cooks  # True or False
    n_cpus = int(args.n_cpus)  # Integer

    # Make output directory if it doesn't exist
    if not os.path.exists(outdir_path):
        os.makedirs(outdir_path)

    # Log file
    if os.path.exists(os.path.join(outdir_path, "pydeseq2.log")):
        os.remove(os.path.join(outdir_path, "pydeseq2.log"))
    
    # Make a random workflow hash
    workflow_hash = hashlib.sha1()
    workflow_hash.update(str(random.getrandbits(128)).encode("utf-8"))

    # Log file
    if os.path.exists(os.path.join(outdir_path, "pydeseq2.log")):
        os.remove(os.path.join(outdir_path, "pydeseq2.log"))
    logging.basicConfig(
        filename=os.path.join(outdir_path, "pydeseq2.log"),
        level=logging.INFO, 
        format='%(asctime)s - %(levelname)s - %(message)s'
        
    )
    logging.info(f"Workflow hash: {workflow_hash.hexdigest()}")

    # Log versions
    import pydeseq2
    logging.info("pydeseq2 version: %s", pydeseq2.__version__)

    # Write arguments to log file
    logging.info("input_h5ad_path: %s", input_h5ad_path)
    logging.info("outdir_path: %s", outdir_path)
    logging.info("design_factors: %s", design_factors)
    logging.info("continuous_factors: %s", continuous_factors)
    logging.info("reference_factor: %s", reference_factor)
    logging.info("reference_value: %s", reference_value)
    logging.info("refit_cooks: %s", refit_cooks)
    logging.info("n_cpus: %s", n_cpus)

    # Read in the h5ad file as an AnnData object
    import scanpy as sc
    logging.info("Reading in h5ad file as AnnData object")
    pdata = sc.read_h5ad(input_h5ad_path)

    # Build the DESeq2 objects
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference  # this is for the inference engine we will use
    logging.info("Building DESeq2 objects")
    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(
        adata=pdata,
        design_factors=design_factors,
        continuous_factors=continuous_factors,
        ref_level=[reference_factor, reference_value],
        refit_cooks=refit_cooks,
        inference=inference
    )

    # Run DESeq2
    logging.info("Running DESeq2")
    dds.deseq2()

    # Check available contrasts
    contrast_factors = pdata.obs[reference_factor].unique()
    logging.info(f"Available contrasts: {contrast_factors}")

    ## For each contrast, get DESeq2 statistics and generate a few plots
    from pydeseq2.ds import DeseqStats
    import decoupler as dc
    for curr_contrast in contrast_factors:
        if curr_contrast == reference_value:
            continue
        else:
            curr_contrast = curr_contrast.replace("_", "-")
            logging.info(f"Running DESeq2 for {curr_contrast} vs {reference_value}")
            stat_res = DeseqStats(dds, contrast=[reference_factor, curr_contrast, reference_value], inference=inference)
            stat_res.summary()
            results_df = stat_res.results_df
            results_df.sort_values("padj").to_csv(os.path.join(outdir_path, f"{curr_contrast}_vs_{reference_value}.tsv"), sep="\t")
            dc.plot_volcano_df(
                results_df, 
                x='log2FoldChange', 
                y='padj', 
                top=20, 
                save=os.path.join(outdir_path, f"volcano_plot_{curr_contrast}_vs_{reference_value}.png")
            )
            stat_res.lfc_shrink(coeff=f"{reference_factor}_{curr_contrast}_vs_{reference_value}")
            shrunk_results_df = stat_res.results_df
            shrunk_results_df.sort_values("padj").to_csv(os.path.join(outdir_path, f"{curr_contrast}_vs_{reference_value}_shrunkLFC.tsv"), sep="\t")
            dc.plot_volcano_df(
                shrunk_results_df, 
                x='log2FoldChange', 
                y='padj', 
                top=20, 
                save=os.path.join(outdir_path, f"volcano_plot_{curr_contrast}_vs_{reference_value}_shrunkLFCs.png")
            )


if __name__ == "__main__":
    # Setting up argparse
    parser = argparse.ArgumentParser(
        description="Run DESeq2 on an AnnData object."
    )
    parser.add_argument(
        "-i", "--input_h5ad_path", required=True, help="Path to the input h5ad file."
    )
    parser.add_argument(
        "-o", "--outdir_path", required=True, help="Path to the output directory."
    )
    parser.add_argument(
        "-d", "--design_factors", required=True, nargs="+", help="Comma-separated list of factors to include in the design formula."
    )
    parser.add_argument(
        "-c", "--continuous_factors", required=False, default=None, help="Comma-separated list of factors to include in the design formula."
    )
    parser.add_argument(
        "-r", "--reference_factor", required=True, help="String specifying the reference factor."
    )
    parser.add_argument(
        "-v", "--reference_value", required=True, help="String specifying the reference value."
    )
    parser.add_argument(
        "-f", "--refit_cooks", required=False, action="store_true", help="Boolean specifying whether to refit cooks distance."
    )
    parser.add_argument(
        "-n", "--n_cpus", required=False, default=1, type=int, help="Integer specifying the number of CPUs to use."
    )
    args = parser.parse_args()
    main(args)
    