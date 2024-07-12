#!/usr/bin/env python3
"""
peak_analysis_snapatac2.py
Run with SnapATAC2 version 2.5.3

Usage:
    python peak_analysis_snapatac2.py -i <input_h5ad_paths> -o <outdir_path> [options]

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
    input_path = args.input_path
    groupby_key = args.groupby_key
    outdir_path = args.outdir_path
    blacklist_path = args.blacklist_path

    # Make output directory if doesn't already exist
    if not os.path.exists(outdir_path):
        os.makedirs(outdir_path)

    # Set up logging
    log_file = os.path.join(outdir_path, "peak_analysis.log")
    if os.path.exists(log_file):
        os.remove(log_file)
    random_id = hashlib.md5(str(random.getrandbits(128)).encode()).hexdigest()
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    import snapatac2 as snap
    logging.info("Run hash: " + random_id)
    logging.info("SnapATAC version: " + snap.__version__)
    logging.info(f"Arguments: {args}")
    
    # Load anndata or andataset
    if input_path.endswith(".h5ad"):
        adata = snap.read(input_path)
        logging.info(f"Loading AnnData from {input_path}")
    elif input_path.endswith(".h5ads"):
        adata = snap.read_dataset(input_path)
        logging.info(f"Loading AnnDataset from {input_path}")

    # Filter peaks
    import pyranges as pr
    from utils import filter_peaks
    blacklist = pr.read_bed(blacklist_path)
    peaksets = adata.uns["macs3"].copy()
    for k, v in peaksets.items():
        v.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"}, inplace=True)
    peaksets_pr = {k: pr.PyRanges(v) for k, v in peaksets.items()}
    for k, v in peaksets_pr.items():
        peaksets_pr[k] = filter_peaks(v, blacklist, score_col="score")
    for k, v in peaksets_pr.items():
        v.to_bed(os.path.join(outdir_path, f"{k}.filt.narrowPeak"))

    # Annotate peak
    from utils import run_shell_cmd
    for k, v in peaksets_pr.items():
        infile = f"/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/annotation/2024_01_15/sample/mo1/snapatac2/peak_calls/{k}.filt.narrowPeak"
        outfile = f"/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/annotation/2024_01_15/sample/mo1/snapatac2/peak_calls/{k}.filt.annot.bed"
        cmd = f"annotatePeaks.pl {infile} hg38 > {outfile}"
        run_shell_cmd(cmd)

    # Fraction of reads in peaks
    regions = {}
    for peakset in peaksets:
        regions[peakset] = peaksets[peakset].apply(lambda x: f"{x[0]}:{x[1]}-{x[2]}", axis=1)
        regions[peakset] = regions[peakset].values
    snap.metrics.frip(adata, regions=regions, n_jobs=-1)

    # Find marker regions based on passed in obs column name
    peaks = snap.tl.merge_peaks(adata.uns['macs3'], snap.genome.hg38)
    peak_mat = snap.pp.make_peak_matrix(adata, use_rep=peaks['Peaks'])
    marker_peaks = snap.tl.marker_regions(peak_mat, groupby=groupby_key, pvalue=0.05)
    snap.pl.regions(peak_mat, groupby=groupby_key, peaks=marker_peaks, interactive=False)

    # Look for motif enrichment in marker regions
    motifs = snap.tl.motif_enrichment(
        motifs=snap.datasets.cis_bp(unique=True),
        regions=marker_peaks,
        genome_fasta=snap.genome.hg38,
    )
    snap.pl.motif_enrichment(motifs, max_fdr=0.0001, height=1600, interactive=False)

    # Perform differential test between conditions if column specified
    group1 = "Naive B"
    group2 = "Memory B"
    naive_B = data.obs['cell_type'] == group1
    memory_B = data.obs['cell_type'] == group2
    peaks_selected = np.logical_or(
        peaks[group1].to_numpy(),
        peaks[group2].to_numpy(),
    )
    diff_peaks = snap.tl.diff_test(
        peak_mat,
        cell_group1=naive_B,
        cell_group2=memory_B,
        features=peaks_selected,
    )
    diff_peaks = diff_peaks.filter(pl.col('adjusted p-value') < 0.01)
    diff_peaks.head()

if __name__ == "__main__":
    # Setting up argparse
    parser = argparse.ArgumentParser(description="analyze a single h5ads sample using snapatac2.")
    parser.add_argument('--input_path', type=str, required=True, help='Path to input AnnData or AnnDataset, will infer from file extension.')
    parser.add_argument('--groupby_key', type=str, required=True, help='Key to groupby.')
    parser.add_argument('--outdir_path', type=str, required=True, help='Path to output directory.')
    parser.add_argument('--blacklist_path', type=str, required=False, default="/cellar/users/aklie/data/ref/blacklists/hg38/temp.bed", help='Path to blacklist bed file.')
    parser
    
    
    # Parse args
    args = parser.parse_args()

    # Call main
    main(args)