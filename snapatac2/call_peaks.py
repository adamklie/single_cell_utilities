#!/usr/bin/env python3
"""
call_peaks.py

This script calls peaks on an AnnDataSet or AnnData object (assuming it has fragments or insertions stored)
using SnapATAC2's MACS3 wrapper.

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
    path_input = args.path_input
    path_outdir = args.path_outdir
    path_annotations = args.path_annotations
    annotations_key = args.annotations_key
    groupby_key = args.groupby_key
    selections = args.selections
    replicate_key = args.replicate_key
    replicate_qvalue = args.replicate_qvalue
    min_len = args.min_len
    save_peaks = args.save_peaks
    save_peak_matrix = args.save_peak_matrix
    chromsizes_path = args.chromsizes_path
    n_jobs = args.n_jobs

    # Make output directory if doesn't already exist
    if not os.path.exists(path_outdir):
        os.makedirs(path_outdir)

    # Set up logging
    log_file = os.path.join(path_outdir, "call_peaks.log")
    if os.path.exists(log_file):
        os.remove(log_file)
    random_id = hashlib.md5(str(random.getrandbits(128)).encode()).hexdigest()
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    import snapatac2 as snap
    import anndata as ad
    logging.info("Run hash: " + random_id)
    logging.info("SnapATAC version: " + snap.__version__)
    logging.info(f"Arguments: {args}")
    
    # Load anndata or andataset
    if path_input.endswith(".h5ad"):
        adata = snap.read(path_input)
        logging.info(f"Loading AnnData from {path_input}")
    elif path_input.endswith(".h5ads"):
        adata = snap.read_dataset(path_input)
        logging.info(f"Loading AnnDataset from {path_input}")

    # If an annotation file is passed in, we need to subset
    import pandas as pd
    import numpy as np
    if path_annotations:
        logging.info(f"Subsetting AnnData based on cell ids in first column in annotation file {path_annotations} and adding annotation to obs.")
        if path_annotations.endswith(".tsv") | path_annotations.endswith(".txt"):
            cell_annotations = pd.read_csv(path_annotations, sep="\t", index_col=0, header=None)[1]
        elif path_annotations.endswith(".csv"):
            cell_annotations = pd.read_csv(path_annotations, index_col=0, header=None)[1]
        else:
            raise ValueError("Annotation file must be in .tsv, .txt, or .csv format.")
        cell_annotations = cell_annotations[~cell_annotations.isna()]
        logging.info(f"Annotation file contains {len(cell_annotations)} cell ids.")
        if isinstance(adata, snap.AnnDataSet):
            subset_out = os.path.join(path_outdir, "annotated.h5ads")
            adata = adata.subset(obs_indices=np.where(pd.Index(adata.obs_names).isin(cell_annotations.index))[0], out=subset_out)[0]
            adata.obs[annotations_key] = cell_annotations[adata.obs_names].values.tolist()
        elif isinstance(adata, ad.AnnData):
            adata = adata[pd.Index(adata.obs_names).isin(cell_annotations.index)]
            adata.obs[annotations_key] = cell_annotations[adata.obs_names].values.tolist()
        groupby_key = annotations_key
    else:
        assert groupby_key, "Groupby key must be provided if annotation file is not passed in."

    # Run peak calling
    if 'macs3' in adata.uns and not args.overwrite:
        logging.info("Peak calling already run on this AnnData object. Use --overwrite to rerun.")
    else:
        logging.info(f"Running peak calling with SnapATAC2 version {snap.__version__}.")
        time_in = time.time()
        snap.tl.macs3(
            adata=adata, 
            groupby=groupby_key,
            replicate=replicate_key,
            replicate_qvalue=replicate_qvalue if replicate_key is not None else None,
            selections=selections,
            n_jobs=n_jobs,
            min_len=min_len,
        )
        time_out = time.time()
        logging.info(f"Peak calling completed in {time_out - time_in} seconds.")

    # Save peaks
    if save_peaks:
        
        # Make directory if doesn't already exist
        if not os.path.exists(save_peaks):
            os.makedirs(save_peaks)
            
        # Get consensus peaks
        if chromsizes_path:
            logging.info(f"Using chromsizes file at {chromsizes_path}.")
            chromsizes = pd.read_csv(chromsizes_path, sep="\t", header=None, names=["chrom", "size"])
            chromsizes = chromsizes.set_index("chrom").to_dict()["size"]
            chromsizes = {k: int(v) for k, v in chromsizes.items()}
            peaks = snap.tl.merge_peaks(adata.uns['macs3'], chromsizes)
        else:
            peaks = snap.tl.merge_peaks(adata.uns['macs3'], snap.genome.hg38)
        
        # Save each peakset to narrowPeak format in output directory
        logging.info("Saving peaksets to narrowPeak format.")
        time_in = time.time()
        import polars as pl
        peaks.to_pandas().to_csv(os.path.join(save_peaks, "peak_sources.tsv"), sep="\t", index=False, header=True)
        peak_ids = peaks.with_columns([
            pl.col('Peaks').str.split(":").list.get(0).alias('chrom'),
            pl.col('Peaks').str.split(":").list.get(1).str.split("-").list.get(0).cast(pl.Int64).alias('start'),
            pl.col('Peaks').str.split(":").list.get(1).str.split("-").list.get(1).cast(pl.Int64).alias('end')
        ]).select(['chrom', 'start', 'end']).to_pandas()
        peak_ids.to_csv(
            os.path.join(save_peaks, "consensus_peaks.bed"),
            sep="\t",
            header=False,
            index=False,
        )
        for cell_id in adata.uns['macs3']:
            adata.uns['macs3'][cell_id].write_csv(os.path.join(save_peaks, f"{cell_id}.narrowPeak"), separator="\t", include_header=False)
            columns = adata.uns['macs3'][cell_id].columns
            with open(os.path.join(save_peaks, f"{cell_id}.narrowPeak.schema"), "w") as f:
                f.write("\n".join(columns))
        time_out = time.time()
        logging.info(f"Peaksets saved in {time_out - time_in} seconds.")

    # Export peak matrix
    if save_peak_matrix:
        
        # Make directory if doesn't already exist
        if not os.path.exists(save_peak_matrix):
            os.makedirs(save_peak_matrix)

        # Build peak matrix
        peak_mat = snap.pp.make_peak_matrix(adata, use_rep=peaks['Peaks'])

        # Save as h5ad
        peak_mat.write(os.path.join(save_peak_matrix, "peak_mat.h5ad"))

        # Save as mtx, features.tsv.gz and barcodes.tsv.gz
        from scipy import io
        io.mmwrite(os.path.join(save_peak_matrix, "mtx.mtx"), peak_mat.X)
        peak_mat.var.index.to_series().to_csv(os.path.join(save_peak_matrix, "features.tsv.gz"), sep="\t", header=False, index=False, compression="gzip")
        peak_mat.obs.index.to_series().to_csv(os.path.join(save_peak_matrix, "barcodes.tsv.gz"), sep="\t", header=False, index=False, compression="gzip")


if __name__ == "__main__":
    # Setting up argparse
    parser = argparse.ArgumentParser(description="Call peaks on an AnnDataSet or AnnData object using SnapATAC2's MACS3 wrapper.")
    parser.add_argument('--path_input', type=str, required=True, help='Path to input AnnData or AnnDataset, will infer from file extension.')
    parser.add_argument('--path_outdir', type=str, required=True, help='Path to output directory.')
    parser.add_argument('--path_annotations', type=str, required=False, help='Path to annotation file.')
    parser.add_argument('--annotations_key', type=str, required=False, default="annotations", help='Key to add to .obs for passed in annotations.')
    parser.add_argument('--groupby_key', type=str, required=False, help='Key to pseudobulk cells prior to peak calling.')
    parser.add_argument('--selections', type=str, required=False, default=None, help='Which particular groups to run peak calling on.')
    parser.add_argument('--replicate_key', type=str, required=False, default=None, help='Key to split cells by replicates prior to peak calling in each group.')
    parser.add_argument('--replicate_qvalue', type=float, required=False, default=0.1, help='Q-value cutoff for replicates.')
    parser.add_argument('--min_len', type=int, required=False, default=None, help='Minimum length of peaks to call.')
    parser.add_argument('--save_peaks', type=str, required=False, default=None, help='Path to save peak calls to.')
    parser.add_argument('--save_peak_matrix', type=str, required=False, default=None, help='Path to save peak matrix to.')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing peak calls.')
    parser.add_argument("--chromsizes_path", type=str, required=False, default=None, help="Path to chromsizes file for genome of interest")
    parser.add_argument("--n_jobs", type=int, required=False, default=1, help="Number of jobs to run in parallel.")

    # Parse args
    args = parser.parse_args()

    # Call main
    main(args)
