#!/usr/bin/env python3
"""
integrate_snapatac2.py
Run with SnapATAC2 version 2.5.3

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
    input_path = args.input_path
    outdir_path = args.outdir_path
    groupby_key = args.groupby_key
    replicate_key = args.replicate_key
    annotations_path = args.annotations_path
    chromsizes_path = args.chromsizes_path
    bg2bw_path = args.bg2bw_path
    save_peaks = args.save_peaks
    save_fragments = args.save_fragments
    save_coverage = args.save_coverage
    save_peak_matrix = args.save_peak_matrix
    n_jobs = args.n_jobs

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
    
    # Load anndata or andataset
    if input_path.endswith(".h5ad"):
        adata = snap.read(input_path)
        logging.info(f"Loading AnnData from {input_path}")
    elif input_path.endswith(".h5ads"):
        adata = snap.read_dataset(input_path)
        logging.info(f"Loading AnnDataset from {input_path}")

    # If an annotation file is passed in, we need to subset
    import pandas as pd
    import numpy as np
    if annotations_path:
        logging.info(f"Subsetting AnnData based on cell ids in first column in annotation file {annotations_path} and adding annotation to obs.")
        if annotations_path.endswith(".tsv") | annotations_path.endswith(".txt"):
            cell_annotations = pd.read_csv(annotations_path, sep="\t", index_col=0, header=None)[1]
        elif annotations_path.endswith(".csv"):
            cell_annotations = pd.read_csv(annotations_path, index_col=0, header=None)[1]
        else:
            raise ValueError("Annotation file must be in .tsv, .txt, or .csv format.")
        if isinstance(adata, snap.AnnDataSet):
            subset_out = os.path.join(outdir_path, "annotated.h5ads")
            adata = adata.subset(obs_indices=np.where(pd.Index(adata.obs_names).isin(cell_annotations.index))[0], out=subset_out)[0]
            adata.obs["annotation"] = cell_annotations[adata.obs_names].values.tolist()
        elif isinstance(adata, snap.AnnData):
            adata = adata[pd.Index(adata.obs_names).isin(cell_annotations.index)]
            adata.obs["annotation"] = cell_annotations[adata.obs_names].values.tolist()
        groupby_key = "annotation"
    else:
        assert groupby_key, "Groupby key must be provided if annotation file is passed in."

    # Exporting bigwig files
    if save_coverage:

        # Make directory if doesn't already exist
        if not os.path.exists(save_coverage):
            os.makedirs(save_coverage)
        
        # Export bedgraph files for visualization
        logging.info("Exporting coverage in bedgraph format file.")
        time_in = time.time()
        snap.ex.export_coverage(
            adata=adata,
            groupby=groupby_key,
            suffix=".bedgraph",
            output_format="bedgraph",
            out_dir=save_coverage,
        )
        time_out = time.time()
        logging.info(f"Coverage exported to {save_coverage} in {time_out - time_in} seconds.")
        
        # Converting to bigwig
        logging.info("Creating browser extensible data files.")
        time_in = time.time()
        for file in os.listdir(save_coverage):
            if file.endswith(".bedgraph"):
                os.system(f"sort -k1,1 -k2,2n {os.path.join(save_coverage, file)} > {os.path.join(save_coverage, file).replace('.bedgraph', '.sorted.bedgraph')}")
                os.system(f"{bg2bw_path} {os.path.join(save_coverage, file).replace('.bedgraph', '.sorted.bedgraph')} {chromsizes_path} {os.path.join(save_coverage, file).replace('.bedgraph', '.bw')}")
                os.system(f"bgzip -f {os.path.join(save_coverage, file).replace('.bedgraph', '.sorted.bedgraph')}")
                os.system(f"tabix -p bed {os.path.join(save_coverage, file).replace('.bedgraph', '.sorted.bedgraph')}.gz")                
                os.remove(os.path.join(save_coverage, file))
        time_out = time.time()
        logging.info(f"Coverage exported to {save_coverage} in {time_out - time_in} seconds.")

    # Export fragments to BED format
    if save_fragments:

        # Make directory if doesn't already exist
        if not os.path.exists(save_fragments):
            os.makedirs(save_fragments)

        # Export fragments to BED format
        logging.info("Exporting fragments to BED format.")
        time_in = time.time()
        snap.ex.export_fragments(
            adata=adata,
            groupby=groupby_key,
            suffix=".bed.gz",
            out_dir=save_fragments,
        )
        time_out = time.time()
        logging.info(f"Fragments exported to {save_fragments} in {time_out - time_in} seconds.")

    # Run peak calling
    logging.info(f"Running peak calling with SnapATAC2 version {snap.__version__}.")
    time_in = time.time()
    snap.tl.macs3(
        adata=adata, 
        groupby=groupby_key,
        replicate=replicate_key,
        replicate_qvalue=0.1 if replicate_key is not None else None,
        n_jobs=n_jobs,
        #min_len=200  for snapatac2 v2.6 upgrade
    )
    time_out = time.time()
    logging.info(f"Peak calling completed in {time_out - time_in} seconds.")

    # Save peaks
    if save_peaks:
        
        # Make directory if doesn't already exist
        if not os.path.exists(save_peaks):
            os.makedirs(save_peaks)
            
        # Get consensus peaks
        peaks = snap.tl.merge_peaks(adata.uns['macs3'], snap.genome.hg38)
        
        # Save each peakset to narrowPeak format in output directory
        logging.info("Saving peaksets to narrowPeak format.")
        time_in = time.time()
        peak_ids = peaks.with_columns(
            chrom=peaks['Peaks'].str.split(":").apply(lambda x: x[0]),
            start=peaks['Peaks'].str.split(":").apply(lambda x: int(x[1].split("-")[0])),
            end=peaks['Peaks'].str.split(":").apply(lambda x: int(x[1].split("-")[1])),
        )[['chrom', 'start', 'end']].to_pandas()
        peak_ids.to_csv(
            os.path.join(save_peaks, "consensus_peaks.bed"),
            sep="\t",
            header=False,
            index=False,
        )
        for cell_id in adata.uns['macs3']:
            adata.uns['macs3'][cell_id].write_csv(os.path.join(save_peaks, f"{cell_id}.narrowPeak"), separator="\t", has_header=False)
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
    parser = argparse.ArgumentParser(description="analyze a single h5ads sample using snapatac2.")
    parser.add_argument('--input_path', type=str, required=False, help='Path to input AnnData or AnnDataset, will infer from file extension.')
    parser.add_argument('--outdir_path', type=str, required=True, help='Path to output directory.')
    parser.add_argument('--groupby_key', type=str, required=False, help='Key to groupby.')
    parser.add_argument('--replicate_key', type=str, required=False, default=None, help='Key to replicate.')
    parser.add_argument('--annotations_path', type=str, required=False, help='Path to annotation file.')
    parser.add_argument('--bg2bw_path', type=str, required=False, help='Path to bedGraphToBigWig.', default="/cellar/users/aklie/opt/bedGraphToBigWig")
    parser.add_argument("--chromsizes_path", type=str, required=False, default="/cellar/users/aklie/data/ref/genomes/hg38/GRCh38_EBV.chrom.sizes", help="Path to chromsizes file for genome of interest")
    parser.add_argument('--output_prefix', type=str, required=False, help='Prefix for output files.')
    parser.add_argument('--save_peaks', type=str, required=False, default=None, help='Path to save peaks.')
    parser.add_argument('--save_fragments', type=str, required=False, default=None, help='Path to save fragments.')
    parser.add_argument('--save_coverage', type=str, required=False, default=None, help='Path to save bigWig coverage.')
    parser.add_argument('--save_peak_matrix', type=str, required=False, default=None, help='Path to save peak matrix.')
    parser.add_argument("--n_jobs", type=int, required=False, default=1, help="Number of jobs to run in parallel.")

    # Parse args
    args = parser.parse_args()

    # Call main
    main(args)
