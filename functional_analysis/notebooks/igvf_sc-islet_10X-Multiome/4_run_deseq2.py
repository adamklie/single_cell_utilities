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
import sys
import time
import random
import logging
import argparse
from tqdm.auto import tqdm

def main(args):
    
    # Parse args
    h5ad_file = args.h5ad_file  # Simple path like
    design_factors = args.design_factors  # Comma-separated list of factors
    continuous_factors = args.continuous_factors  # Comma-separated list of factors
    reference_factor = args.reference_factor  # String
    reference_value = args.reference_value # String
    refit_cooks = args.refit_cooks  # True or False
    shrink = args.shrink  # True or False
    n_cpus = args.n_cpus  # Integer
    output_dir = args.output_dir  # Simple path like

    # Set up logging file in output directory
    sys.path.append("/cellar/users/aklie/data/igvf/bin")
    from utils import make_dirs
    make_dirs(output_dir)
    time_id = time.strftime("%Y%m%d-%H%M%S")
    run =  random.getrandbits(128)
    run_id = time_id + "_" + str(run) +  "_deseq2"
    log_file = os.path.join(output_dir, run_id + '.log')
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    import pydeseq2
    logging.info("pydeseq2 version: %s", pydeseq2.__version__)

    # Write arguments to log file
    logging.info("h5ad_file: %s", h5ad_file)
    logging.info("design_factors: %s", design_factors)
    logging.info("continuous_factors: %s", continuous_factors)
    logging.info("ref_level: %s", ref_level)
    logging.info("refit_cooks: %s", refit_cooks)
    logging.info("shrink: %s", shrink)
    logging.info("n_cpus: %s", n_cpus)
    logging.info("output_dir: %s", output_dir)

    # Read in the h5ad file as an AnnData object
    import scanpy as sc
    logging.info("Reading in h5ad file as AnnData object")
    adata = sc.read_h5ad(h5ad_file)

    # Build the DESeq2 objects
    from pydeseq2.dds import DeseqDataSet
    logging.info("Building DESeq2 objects")
    dds = DeseqDataSet(
        adata=adata,
        design_factors=design_factors,
        continuous_factors=continuous_factors,
        ref_level=ref_level,
        refit_cooks=refit_cooks,
        n_cpus=n_cpus,
    )

    # Run DESeq2
    logging.info("Running DESeq2")
    dds.deseq2()

    # Extract contrasts
    logging.info("Extracting contrasts")
    from pydeseq2.ds import DeseqStats
    stat_res = DeseqStats(dds, contrast=["condition", 'IFNg', 'control'], n_cpus=3)

    # Compute Wald test
    logging.info("Computing Wald test")
    stat_res.summary()

    # Shrink LFCs
    if shrink:
        logging.info("Shrinking LFCs")
        stat_res.lfc_shrink(coeff='condition_IFNg_vs_control')

    # Extract results
    results_df = stat_res.results_df
    
    # 
    import decoupler as dc
    dc.plot_volcano_df(results_df, x='log2FoldChange', y='padj', top=20)