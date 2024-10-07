#!/usr/bin/env python3
"""
3_preprocess_anndata.py

This script preprocesses a single h5ad sample using SnapATAC2.
Needs to be run after 1_create_anndatas_from_frag_files.sh or
creating AnnDatas from fragment files some other way.

The script performs the following steps:
    1. Reading and processing of the h5ad file.
    2. Plotting tsse.
    3. Filtering of cells.
    4. Generating tile matrix.
    5. Feature selection.
    6. Doublet detection and filtering.

Usage:
    python 3_preprocess_anndata.py -i <input_h5ad> -o <output_dir> [options]
s
Arguments:
    -i/--input       Path to the input h5ad file.
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
    input_h5ad = args.input_h5ad
    output_dir = args.output_dir
    min_counts = args.min_counts
    min_tsse = args.min_tsse
    max_counts = args.max_counts
    bin_size = args.bin_size
    n_features = args.n_features

    # Set up logging
    sys.path.append("/cellar/users/aklie/data/igvf/bin")
    from utils import make_dirs
    make_dirs(output_dir)
    make_dirs(os.path.join(output_dir, "logs"))
    time_id = time.strftime("%Y%m%d-%H%M%S")
    run =  random.getrandbits(128)
    run_id = time_id + "_" + str(run) +  "_preprocess_anndata"
    log_file = os.path.join(output_dir, "logs", run_id + '.log')
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    import snapatac2 as snap
    logging.info("SnapATAC version: " + snap.__version__)
    
    # Check and create output directory
    logging.info("Making output directory if doesn't already exist: " + output_dir)
    make_dirs(output_dir)

    # Print the file name
    logging.info("Processing file: " + input_h5ad)

    # Read in the h5ad file
    adata_atac = snap.read(input_h5ad)

    # Create processed h5ad file path
    logging.info("Creating processed h5ad")
    time_in = time.time()
    base_name = os.path.basename(input_h5ad).split('.')[0] + "_" + "processed.h5ad"
    out_path = os.path.join(output_dir, base_name)
    if os.path.exists(out_path):
        logging.info("File already exists, exiting")
        sys.exit()
    adata_atac_processed = adata_atac.copy(filename=out_path)
    adata_atac.close()
    time_out = time.time()
    logging.info("Copy of h5ad created in " + str(time_out - time_in) + " seconds")

    # Plotting tsse
    logging.info("Plotting tsse")
    time_in = time.time()
    snap.pl.tsse(adata_atac_processed, interactive=False, out_file=out_path.replace(".h5ad", "_tsse.png"))
    time_out = time.time()
    logging.info("tsse plotted in " + str(time_out - time_in) + " seconds")

    # Filtering
    logging.info("Filtering cells")
    time_in = time.time()
    snap.pp.filter_cells(adata_atac_processed, min_counts=min_counts, min_tsse=min_tsse, max_counts=max_counts)
    time_out = time.time()
    logging.info("Cells filtered in " + str(time_out - time_in) + " seconds")

    # Tile matrix
    logging.info("Generating tile matrix")
    time_in = time.time()
    snap.pp.add_tile_matrix(adata_atac_processed, bin_size=bin_size)
    time_out = time.time()
    logging.info("Tile matrix generated in " + str(time_out - time_in) + " seconds")

    # Feature selection
    logging.info("Selecting features")
    time_in = time.time()
    snap.pp.select_features(adata_atac_processed, n_features=n_features)
    time_out = time.time()
    logging.info("Features selected in " + str(time_out - time_in) + " seconds")

    # Doublet detection
    logging.info("Detecting doublets")
    time_in = time.time()
    snap.pp.scrublet(adata_atac_processed)
    time_out = time.time()
    logging.info("Doublets detected in " + str(time_out - time_in) + " seconds")

    # Doublet filtering
    logging.info("Filtering doublets")
    time_in = time.time()
    snap.pp.filter_doublets(adata_atac_processed)
    time_out = time.time()
    logging.info("Doublets filtered in " + str(time_out - time_in) + " seconds")

    # Close file
    adata_atac_processed.close()


if __name__ == "__main__":
    # Setting up argparse
    parser = argparse.ArgumentParser(description="Preprocess a single h5ad sample using snapatac2.")
    parser.add_argument('-i', '--input_h5ad', required=True, help="Path to the input h5ad file.")
    parser.add_argument('-o', '--output_dir', required=True, help="Path to the output directory.")
    
    # Additional optional arguments
    parser.add_argument('--min_counts', type=int, default=5000, help="Minimum counts for cell filtering.")
    parser.add_argument('--min_tsse', type=int, default=10, help="Minimum tsse for cell filtering.")
    parser.add_argument('--max_counts', type=int, default=100000, help="Maximum counts for cell filtering.")
    parser.add_argument('--bin_size', type=int, default=5000, help="Bin size for tile matrix generation.")
    parser.add_argument('--n_features', type=int, default=50000, help="Number of features for feature selection.")

    # Parse the arguments
    args = parser.parse_args()

    # Call the main function
    main(args)
