#!/usr/bin/env python3
"""
pseudobulk.py

This script creates pseudobulk fragment and bigwig files from an AnnDataSet or AnnData object.

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
    chromsizes_path = args.chromsizes_path
    bg2bw_path = args.bg2bw_path
    save_fragments = args.save_fragments
    save_coverage = args.save_coverage
    n_jobs = args.n_jobs
    log_file_name = args.log_file_name

    # Make output directory if doesn't already exist
    if not os.path.exists(path_outdir):
        os.makedirs(path_outdir)

    # Set up logging
    log_file = os.path.join(path_outdir, log_file_name)
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
        adata = ad.read_h5ad(path_input)
        logging.info(f"Loading AnnData from {path_input}")
        logging.info(f"AnnData contains {adata.shape[0]} cells and {adata.shape[1]} peaks.")
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
            adata = adata[pd.Index(adata.obs_names).isin(cell_annotations.index), :].copy()
            adata.obs[annotations_key] = cell_annotations[adata.obs_names].values.tolist()
            adata.obs[annotations_key] = adata.obs[annotations_key].astype(str)
            logging.info(f"Writing annotated AnnData to {os.path.join(path_outdir, 'annotated.h5ad')}")
            adata.write(os.path.join(path_outdir, "annotated.h5ad"))
        groupby_key = annotations_key
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
            n_jobs=n_jobs,
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


if __name__ == "__main__":
    # Setting up argparse
    parser = argparse.ArgumentParser(description="Create pseudobulk fragment and bigwig files from an AnnDataSet or AnnData object.")
    parser.add_argument('--path_input', type=str, required=True, help='Path to input AnnData or AnnDataset, will infer from file extension.')
    parser.add_argument('--path_outdir', type=str, required=True, help='Path to output directory.')
    parser.add_argument('--path_annotations', type=str, required=False, help='Path to annotation file.')
    parser.add_argument('--annotations_key', type=str, required=False, help='Key to store annotations.')
    parser.add_argument('--groupby_key', type=str, required=False, help='Key to groupby.')
    parser.add_argument('--bg2bw_path', type=str, required=False, help='Path to bedGraphToBigWig.', default="/cellar/users/aklie/opt/bedGraphToBigWig")
    parser.add_argument("--chromsizes_path", type=str, required=False, default="/cellar/users/aklie/data/ref/genomes/hg38/GRCh38_EBV.chrom.sizes", help="Path to chromsizes file for genome of interest")
    parser.add_argument('--save_fragments', type=str, required=False, default=None, help='Path to save fragments.')
    parser.add_argument('--save_coverage', type=str, required=False, default=None, help='Path to save bigWig coverage.')
    parser.add_argument('--n_jobs', type=int, required=False, default=1, help='Number of jobs to run in parallel.')
    parser.add_argument('--log_file_name', type=str, required=False, default="pseudobulk.log", help='Name of log file.')

    # Parse args
    args = parser.parse_args()

    # Call main
    main(args)
