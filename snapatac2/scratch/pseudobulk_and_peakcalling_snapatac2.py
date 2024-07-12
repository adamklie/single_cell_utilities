#!/usr/bin/env python3
"""
integrate_snapatac2.py
Run with SnapATAC2 version 2.3

This script creates AnnDatasets
and performs feature selection, spectral embedding, clustering, and plotting.

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
    input_h5ads_path = args.input_h5ads_path
    input_h5ad_path = args.input_h5ad_path
    outdir_path = args.outdir_path
    save_coverage = args.save_coverage
    save_peak_matrix = args.save_peak_matrix
    groupby_key = args.groupby_key

    # Make output directory if doesn't already exist
    if not os.path.exists(outdir_path):
        os.makedirs(outdir_path)

    # Set up logging
    log_file = os.path.join(outdir_path, "pseudobulk_and_peakcalling.log")
    if os.path.exists(log_file):
        os.remove(log_file)
    random_id = hashlib.md5(str(random.getrandbits(128)).encode()).hexdigest()
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    import snapatac2 as snap
    logging.info("Run hash: " + random_id)
    logging.info("SnapATAC version: " + snap.__version__)
    logging.info(f"Arguments: {args}")
    
    # Load anndataset
    if input_h5ads_path:
        logging.info(f"Loading AnnDataset from {input_h5ads_path}")
        adata = snap.read_dataset(input_h5ads_path)
    elif input_h5ad_path:
        logging.info(f"Loading AnnData from {input_h5ad_path}")
        adata = snap.read(input_h5ad_path)

    # Exporting bigwig files
    if save_coverage:
        # Make directory if doesn't already exist
        if not os.path.exists(save_coverage):
            os.makedirs(save_coverage)
        
        # Export bigWig files for visualization
        logging.info("Exporting coverage in bigwig format file.")
        time_in = time.time()
        snap.ex.export_coverage(
            adata=adata,
            groupby=groupby_key,
            out_dir=save_coverage,
        )
        time_out = time.time()
        logging.info(f"Coverage exported to {save_coverage} in {time_out - time_in} seconds.")

        logging.info("Exporting fragments to BED format.")
        snap.ex.export_bigwig(
            adata=adata,
            groupby=groupby_key,
            out_dir=save_coverage,
        )
        logging.info(f"Pseudobulk bigwig files exported to {save_coverage} in {time_out - time_in} seconds.")

    # Run peak calling
    logging.info("Running peak calling with MACS2.")
    time_in = time.time()
    snap.tl.call_peaks(
        adata=adata,
        groupby=groupby_key,
        out_dir=outdir_path,
    )
    time_out = time.time()
    logging.info(f"Peak calling completed in {time_out - time_in} seconds.")

    # Export peak matrix
    if save_peak_matrix:
        # Make directory if doesn't already exist
        if not os.path.exists(save_peak_matrix):
            os.makedirs(save_peak_matrix)
        peak_mat = snap.pp.make_peak_matrix(adata)
        peak_mat.write(os.path.join(save_peak_matrix, "peak_mat.h5ad"))


if __name__ == "__main__":
    # Setting up argparse
    parser = argparse.ArgumentParser(description="analyze a single h5ads sample using snapatac2.")
    parser.add_argument('--input_h5ads_path', type=str, required=False, help='Path to input h5ads.')
    parser.add_argument('--input_h5ad_path', type=str, required=False, help='Path to input h5ad.')
    parser.add_argument('--groupby_key', type=str, required=True, help='Key to groupby.')
    parser.add_argument('--outdir_path', type=str, required=True, help='Path to output directory.')
    parser.add_argument('--output_prefix', type=str, required=False, help='Prefix for output files.')
    parser.add_argument('--save_coverage', type=str, required=False, default=None, help='Path to save pseudobulk.')
    parser.add_argument('--save_peak_matrix', type=str, required=False, default=None, help='Path to save peak matrix.')

    # Parse args
    args = parser.parse_args()

    # Call main
    main(args)
