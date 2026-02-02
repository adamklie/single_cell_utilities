#!/usr/bin/env python3
"""
preprocess.py

Preprocesses a single h5ad sample using SnapATAC2.

This script performs the following steps:
    1. Reading and processing of the h5ad file.
    2. Plotting TSSE.
    3. Filtering of cells.
    4. Generating tile matrix.
    5. Feature selection.
    6. Doublet detection and filtering.

Usage:
    python preprocess.py -i <input_h5ad> -o <output_dir> [options]

Arguments:
    -i/--input_h5ad   Path to the input h5ad file.
    -o/--output_dir   Path to the output directory.
    [options]         Additional optional parameters.

Author: Adam Klie
"""

import os
import sys
import time
import logging
import argparse
import random


def make_dirs(path: str) -> None:
    """Create directory if it doesn't exist."""
    if not os.path.exists(path):
        os.makedirs(path)


def main(args):
    # Parse args
    input_h5ad = args.input_h5ad
    output_dir = args.output_dir
    min_counts = args.min_counts
    min_tsse = args.min_tsse
    max_counts = args.max_counts
    bin_size = args.bin_size
    n_features = args.n_features

    # Set up output directory and logging
    make_dirs(output_dir)
    make_dirs(os.path.join(output_dir, "logs"))
    time_id = time.strftime("%Y%m%d-%H%M%S")
    run = random.getrandbits(128)
    run_id = f"{time_id}_{run}_preprocess_anndata"
    log_file = os.path.join(output_dir, "logs", f"{run_id}.log")
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    # Import snapatac2 after logging is set up
    import snapatac2 as snap

    logging.info(f"SnapATAC version: {snap.__version__}")
    logging.info(f"Output directory: {output_dir}")
    logging.info(f"Processing file: {input_h5ad}")

    # Read in the h5ad file
    adata_atac = snap.read(input_h5ad)

    # Create processed h5ad file path
    logging.info("Creating processed h5ad")
    time_in = time.time()
    base_name = os.path.basename(input_h5ad).split(".")[0] + "_processed.h5ad"
    out_path = os.path.join(output_dir, base_name)
    if os.path.exists(out_path):
        logging.info("File already exists, exiting")
        sys.exit()
    adata_atac_processed = adata_atac.copy(filename=out_path)
    adata_atac.close()
    time_out = time.time()
    logging.info(f"Copy of h5ad created in {time_out - time_in:.2f} seconds")

    # Plotting TSSE
    logging.info("Plotting TSSE")
    time_in = time.time()
    snap.pl.tsse(
        adata_atac_processed,
        interactive=False,
        out_file=out_path.replace(".h5ad", "_tsse.png"),
    )
    time_out = time.time()
    logging.info(f"TSSE plotted in {time_out - time_in:.2f} seconds")

    # Filtering
    logging.info(
        f"Filtering cells (min_counts={min_counts}, min_tsse={min_tsse}, max_counts={max_counts})"
    )
    time_in = time.time()
    snap.pp.filter_cells(
        adata_atac_processed,
        min_counts=min_counts,
        min_tsse=min_tsse,
        max_counts=max_counts,
    )
    time_out = time.time()
    logging.info(f"Cells filtered in {time_out - time_in:.2f} seconds")

    # Tile matrix
    logging.info(f"Generating tile matrix (bin_size={bin_size})")
    time_in = time.time()
    snap.pp.add_tile_matrix(adata_atac_processed, bin_size=bin_size)
    time_out = time.time()
    logging.info(f"Tile matrix generated in {time_out - time_in:.2f} seconds")

    # Feature selection
    logging.info(f"Selecting features (n_features={n_features})")
    time_in = time.time()
    snap.pp.select_features(adata_atac_processed, n_features=n_features)
    time_out = time.time()
    logging.info(f"Features selected in {time_out - time_in:.2f} seconds")

    # Doublet detection
    logging.info("Detecting doublets")
    time_in = time.time()
    snap.pp.scrublet(adata_atac_processed)
    time_out = time.time()
    logging.info(f"Doublets detected in {time_out - time_in:.2f} seconds")

    # Doublet filtering
    logging.info("Filtering doublets")
    time_in = time.time()
    snap.pp.filter_doublets(adata_atac_processed)
    time_out = time.time()
    logging.info(f"Doublets filtered in {time_out - time_in:.2f} seconds")

    # Close file
    adata_atac_processed.close()
    logging.info("Preprocessing complete")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Preprocess a single h5ad sample using SnapATAC2."
    )
    parser.add_argument(
        "-i", "--input_h5ad", required=True, help="Path to the input h5ad file."
    )
    parser.add_argument(
        "-o", "--output_dir", required=True, help="Path to the output directory."
    )

    # Additional optional arguments
    parser.add_argument(
        "--min_counts",
        type=int,
        default=5000,
        help="Minimum counts for cell filtering. Default: 5000",
    )
    parser.add_argument(
        "--min_tsse",
        type=int,
        default=10,
        help="Minimum TSSE for cell filtering. Default: 10",
    )
    parser.add_argument(
        "--max_counts",
        type=int,
        default=100000,
        help="Maximum counts for cell filtering. Default: 100000",
    )
    parser.add_argument(
        "--bin_size",
        type=int,
        default=5000,
        help="Bin size for tile matrix generation. Default: 5000",
    )
    parser.add_argument(
        "--n_features",
        type=int,
        default=50000,
        help="Number of features for feature selection. Default: 50000",
    )

    args = parser.parse_args()
    main(args)
