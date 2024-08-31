#!/usr/bin/env python3
"""
annotate_peaks.py

This script annotates a peak set using HOMER's annotatePeaks.pl script.

Usage:
    python annotate_peaks.py --path_input_peaks <path_to_peak_file>

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
from utils import run_shell_cmd
from peaks.granges import filter_peaks


def main(args):

    # Parse args
    path_input_peaks = args.path_input_peaks
    path_annotate_pl = args.path_annotate_pl
    path_annotated_peaks = args.path_annotated_peaks

    # Make output directory if it doesn't already exist
    path_outdir = os.path.dirname(path_annotated_peaks)
    if not os.path.exists(path_outdir):
        os.makedirs(os.path.dirname(path_annotated_peaks))

    # Set up logging
    log_file = os.path.join(path_outdir, "annotate_peaks.log")
    if os.path.exists(log_file):
        os.remove(log_file)
    random_id = hashlib.md5(str(random.getrandbits(128)).encode()).hexdigest()
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info("Run hash: " + random_id)
    logging.info(f"Arguments: {args}")

    # Run preflight checks
    if not os.path.exists(path_input_peaks):
        logging.error(f"Could not find input peak file at {path_input_peaks}. Exiting.")
        sys.exit(1)
    logging.info(f"Annotating peaks in {path_input_peaks} using HOMER...")
    if not os.path.exists(path_annotate_pl):
        logging.error(f"Could not find annotatePeaks.pl script at {path_annotate_pl}. Exiting.")
        sys.exit(1)

    cmd = f"{path_annotate_pl} {path_input_peaks} hg38 > {path_annotated_peaks}"
    logging.info(f"Annotating peaks with command: {cmd}")
    run_shell_cmd(cmd)
    logging.info(f"Annotated peaks saved to {path_annotated_peaks}")

    # DONE
    logging.info("Successfully completed peak annotation.")

if __name__ == "__main__":
    
    # Setting up argparse
    parser = argparse.ArgumentParser(description="Annotate peaks using HOMER's annotatePeaks.pl script.")
    parser.add_argument('--path_input_peaks', type=str, required=True, help="Path to the input peak file.")
    parser.add_argument('--path_annotate_pl', type=str, default="/cellar/users/aklie/opt/homer/bin/annotatePeaks.pl", help="Full path to HOMER's annotatePeaks.pl script.")
    parser.add_argument('--path_outdir', type=str, default="annotated_peaks.bed", help="Path to the output directory to save all results.")

    # Parse args
    args = parser.parse_args()

    # Call main
    main(args)
