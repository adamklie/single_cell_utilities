#!/usr/bin/env python3
"""
annotate_peaks.py

This script processes narrowPeak files by filtering peaks based on specified
criteria and annotates them using HOMER.

Usage:
    python annotate_peaks.py --path_peaks <list_of_peak_paths> [options]

Author: Adam Klie
"""

# Imports
import os
import sys
import logging
import argparse
import glob
import hashlib
import random
import pandas as pd
import pyranges as pr
import tqdm.auto as tqdm

sys.path.append("/cellar/users/aklie/projects/igvf/single_cell_utilities")
from utils import run_shell_cmd
from peaks.granges import filter_peaks

def main(args):

    # Parse args
    path_peaks = args.path_peaks
    n_peaks = args.n_peaks
    min_q_value = args.min_q_value
    score_col = args.score_col
    blacklist_path = args.blacklist_path
    path_annotate_pl = args.path_annotate_pl
    narrowpeak_schema = args.narrowpeak_schema
    outdir_path = args.outdir_path

    # Make output directory if it doesn't already exist
    if not os.path.exists(outdir_path):
        os.makedirs(outdir_path)

    # Set up logging
    log_file = os.path.join(outdir_path, "annotate_peaks.log")
    if os.path.exists(log_file):
        os.remove(log_file)
    random_id = hashlib.md5(str(random.getrandbits(128)).encode()).hexdigest()
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info("Run hash: " + random_id)
    logging.info(f"Arguments: {args}")

    # Read the blacklist file
    logging.info(f"Reading blacklist from {blacklist_path}")
    blacklist = pr.read_bed(blacklist_path)

    # Process each narrowPeak file
    for peak_file in tqdm.tqdm(path_peaks):
        logging.info(f"Processing {peak_file}")

        # Read peaks
        peaks = pd.read_csv(peak_file, sep="\t", header=None, names=narrowpeak_schema)

        # Filter peaks
        logging.info(f"Filtering peaks")
        filtered_peaks = filter_peaks(
            pr.PyRanges(peaks),
            blacklist,
            score_col=score_col,
            n_peaks=n_peaks,
            min_q_value=min_q_value
        )

        # Save filtered peaks
        filtered_peaks_file = peak_file.replace('.narrowPeak', '.filt.narrowPeak')
        logging.info(f"Saving filtered peaks to {filtered_peaks_file}")
        filtered_peaks.to_bed(filtered_peaks_file)

        # Annotate peaks using HOMER
        annotated_peaks_file = peak_file.replace('.narrowPeak', '.filt.annot.bed')
        cmd = f"{path_annotate_pl} {filtered_peaks_file} hg38 > {annotated_peaks_file}"
        logging.info(f"Annotating peaks with command: {cmd}")
        run_shell_cmd(cmd)
        logging.info(f"Annotated peaks saved to {annotated_peaks_file}")

        # Move files to output directory
        logging.info(f"Moving files to output directory")
        os.rename(filtered_peaks_file, os.path.join(outdir_path, os.path.basename(filtered_peaks_file)))
        os.rename(annotated_peaks_file, os.path.join(outdir_path, os.path.basename(annotated_peaks_file)))

    logging.info("Successfully completed processing narrowPeak files.")

if __name__ == "__main__":
    # Setting up argparse
    parser = argparse.ArgumentParser(description="Process narrowPeak files and annotate peaks using HOMER.")
    parser.add_argument('--path_peaks', nargs="+", required=True, help="List of paths to narrowPeak files. Space separated.")
    parser.add_argument('--n_peaks', type=int, default=None, help="Number of top peaks to retain after filtering.")
    parser.add_argument('--min_q_value', type=float, default=None, help="Minimum q-value to retain a peak.")
    parser.add_argument('--score_col', type=str, default="q_value", help="Column name of the score column in the narrowPeak file.")
    parser.add_argument('--blacklist_path', type=str, default="/cellar/users/aklie/data/ref/genomes/hg38/blacklist/blacklist.bed.gz", help="Path to the blacklist BED file.")
    parser.add_argument('--path_annotate_pl', type=str, default="/cellar/users/aklie/opt/homer/bin/annotatePeaks.pl", help="Full path to HOMER's annotatePeaks.pl script.")
    parser.add_argument('--outdir_path', type=str, default=".", help="Path to the output directory to save all results.")
    parser.add_argument('--narrowpeak_schema', nargs="+", default=['Chromosome', 'Start', 'End', 'name', 'score', 'strand', 'signal_value', 'p_value', 'q_value', 'peak'], help="List of column names for the narrowPeak file schema.")

    # Parse args
    args = parser.parse_args()

    # Call main
    main(args)
