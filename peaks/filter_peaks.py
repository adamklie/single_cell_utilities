#!/usr/bin/env python3
"""
filter_peaks.py

This script filters an input peak set based on specified criteria.

Usage:
    python filter_peaks.py --path_input_peaks <path_to_peak_file> [options]

Author: Adam Klie
"""

# Imports
import os
import sys
import logging
import argparse
import hashlib
import random
import pandas as pd
import pyranges as pr

sys.path.append("/cellar/users/aklie/projects/igvf/single_cell_utilities")
from peaks.granges import filter_peaks


def main(args):

    # Parse args
    path_input_peaks = args.path_input_peaks
    n_peaks = args.n_peaks
    min_q_value = args.min_q_value
    score_col = args.score_col
    path_blacklist = args.path_blacklist
    peak_schema = args.peak_schema
    path_filtered_peaks = args.path_filtered_peaks

    # Make output directory if it doesn't already exist
    path_outdir = os.path.dirname(path_filtered_peaks)
    if not os.path.exists(path_outdir):
        os.makedirs(path_outdir)

    # Set up logging
    log_file = os.path.join(path_outdir, "filter_peaks.log")
    if os.path.exists(log_file):
        os.remove(log_file)
    random_id = hashlib.md5(str(random.getrandbits(128)).encode()).hexdigest()
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info("Run hash: " + random_id)
    logging.info(f"Arguments: {args}")

    # Read the blacklist file
    logging.info(f"Reading blacklist from {path_blacklist}")
    blacklist = pr.read_bed(path_blacklist)

    # Read peaks
    logging.info(f"Reading peaks from {path_input_peaks}")
    peaks = pd.read_csv(path_input_peaks, sep="\t", header=None, names=peak_schema)

    # Filter peaks
    logging.info(f"Filtering {path_input_peaks}...")
    filtered_peaks = filter_peaks(
        pr.PyRanges(peaks),
        blacklist,
        score_col=score_col,
        n_peaks=n_peaks,
        min_q_value=min_q_value
    )

    # Save filtered peaks
    logging.info(f"Saving filtered peaks to {path_filtered_peaks}")
    filtered_peaks.to_bed(path_filtered_peaks)

    # Done
    logging.info("Successfully completed processing peak files.")

if __name__ == "__main__":
    # Setting up argparse
    parser = argparse.ArgumentParser(description="Filter peak files based on specified criteria.")
    parser.add_argument('--path_input_peaks', type=str, required=True, help="Path to the input peak file.")
    parser.add_argument('--path_filtered_peaks', type=str, required=True, help="Path to save the filtered peak file.")
    parser.add_argument('--path_blacklist', type=str, default="/cellar/users/aklie/data/ref/genomes/hg38/blacklist/blacklist.bed.gz", help="Path to the blacklist BED file.")
    parser.add_argument('--n_peaks', type=int, default=None, help="Number of top peaks to retain after filtering.")
    parser.add_argument('--min_q_value', type=float, default=None, help="Minimum q-value to retain a peak.")
    parser.add_argument('--score_col', type=str, default="q_value", help="Column name of the score column in the peak file.")
    parser.add_argument('--peak_schema', nargs="+", default=['Chromosome', 'Start', 'End', 'name', 'score', 'strand', 'signal_value', 'p_value', 'q_value', 'peak'], help="List of column names for the peak file schema.")

    # Parse args
    args = parser.parse_args()

    # Call main
    main(args)
