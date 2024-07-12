# This script is meant for reading in a set of passed-in fragment files along with
# cell barcodes to generate pseudobulked bed and bigwig files using pycisTopic.

# This script uses argparse to parse arguments for specifying fragment files, 
# output directory, and optional parameters for data processing.

# Required arguments:
# fragment_files: List of paths to ATAC fragment files.
# sample_ids: List of sample IDs.
# metadata_file: File containing metadata information. The file must have a "barcode" column 
#                corresponding to the barcodes found in the fragment files or the barcode must be 
#                the index. Your metadata file must contain a sample ID column that matches the 
#                sample_ids passed in
# output_dir: Path to directory to output bed and bigwig files

# The script will:
# 1. Retrieve the list of fragment files.
# 3. Read in the metadata.
# 4. Check for required columns in metadata.
# 5. Create a dictionary linking sample IDs to fragment files.
# 6. Fetch chromosome sizes.
# 7. Export pseudobulk data.

# Usage:
# python script_name.py --fragment_files [list_of_fragment_files] --sample_ids [list_of_sample_ids] 
#                       --metadata_file [metadata_file_path] --out_dir [output_directory] [optional_arguments]

import os
import glob
import time
import logging
import argparse

def main(args):
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    
    # Log arguments
    logging.info("Arguments: %s", args)
    
    # Get list of fragment files and sample IDs from argparse
    fragment_files = args.fragment_files
    sample_ids = args.sample_ids
    
    # Read in metadata 
    import pandas as pd
    if args.metadata_file.endswith('.csv'):
        metadata = pd.read_csv(args.metadata_file, low_memory=False, index_col=0)
    elif args.metadata_file.endswith('.tsv'):
        metadata = pd.read_csv(args.metadata_file, sep='\t', low_memory=False, index_col=0)
    else:
        logging.error("Metadata file must be a .csv or .tsv file.")
        return
    
    # Check to make sure pseudobulk_column exists in metadata columns
    if args.pseudobulk_column not in metadata.columns:
        logging.error("pseudobulk_column not found in metadata columns.")
        return
    
    # Check to make sure sample_id_column exists in metadata columns
    if args.sample_id_column not in metadata.columns:
        logging.error("sample_id_column not found in metadata columns.")
        return
    
    # Keep only fragment files for samples that are in metadata
    import numpy as np
    metadata_samples = metadata[args.sample_id_column].unique()
    sample_inds = [np.where(np.array(sample_ids) == x)[0][0] for x in metadata_samples]
    fragment_files = [fragment_files[x] for x in sample_inds]
    sample_ids = [sample_ids[x] for x in sample_inds]
    logging.info(f"Fragment file remaining: {fragment_files}")
    if len(fragment_files) == 0:
        logging.error("No fragment files remaining after filtering.")
        return
    
    # Create fragments_dict and log it
    fragments_dict = dict(zip(sample_ids, fragment_files))
    logging.info(f"Fragments dictionary: {fragments_dict}")
    
    # Get the number of categories and log it
    if args.num_cpus is None:
        num_pseudobulks = metadata[args.pseudobulk_column].nunique()
    else:
        num_pseudobulks = args.num_cpus
    logging.info(f"Number of pseudobulks: {num_pseudobulks}")
    
    # Get chromosome sizes
    import requests
    import pyranges as pr
    target_url = args.chrom_sizes  # default: 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
    chromsizes = pd.read_csv(target_url, sep='\t', header=None)
    chromsizes.columns = ['Chromosome', 'End']
    chromsizes['Start'] = [0]*chromsizes.shape[0]
    chromsizes = chromsizes.loc[:,['Chromosome', 'Start', 'End']]
    if args.cellranger_arc_annot:  # default is true
        # Exceptionally in this case, to agree with CellRangerARC annotations
        chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('v', '.') for x in range(len(chromsizes['Chromosome']))]
        chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].split('_')[1] if len(chromsizes['Chromosome'][x].split('_')) > 1 else chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))]
    chromsizes = pr.PyRanges(chromsizes)
    
    # Run the pseudobulk export
    from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
    bw_paths, bed_paths = export_pseudobulk(
        input_data = metadata,
        variable = args.pseudobulk_column,
        sample_id_col = args.sample_id_column,
        chromsizes = chromsizes,
        bed_path = args.out_dir,
        bigwig_path = args.out_dir,
        path_to_fragments = fragments_dict,
        n_cpu = num_pseudobulks,
        normalize_bigwig = True,
        remove_duplicates = True,
        _temp_dir = args.temp_dir,
        split_pattern = args.split_pattern,
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process and merge ATAC fragment files.")
    
    # Required args
    parser.add_argument("--fragment_files", type=str, required=True, nargs='+', help="List of paths to ATAC fragment files.")
    parser.add_argument("--sample_ids", type=str, required=True, nargs='+', help="List of sample IDs.")
    parser.add_argument("--metadata_file", type=str, required=True, help="Path to the metadata file.")
    parser.add_argument("--out_dir", type=str, required=True, help="Directory to save output files.")
    
    # Optional args
    parser.add_argument("--pseudobulk_column", type=str, default="Pseudobulk", help="Column name in the metadata file representing pseudobulk categories.")
    parser.add_argument("--num_cpus", type=int, default=None, help="Number of CPUs to use for Ray parallelization. Default is number of groups in pseudobulk_column")
    parser.add_argument("--sample_id_column", type=str, default="SampleID", help="Column name in the metadata file representing sample IDs.")
    parser.add_argument("--chrom_sizes", type=str, default='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes', help="URL to fetch chromosome sizes.")
    parser.add_argument("--split_pattern", type=str, default="___", help="Pattern to split the fragment file names on to get the sample ID. Default is '___'.")
    parser.add_argument("--cellranger_arc_annot", action="store_true", help="Whether to use CellRangerARC annotations for chromosomes.")
    parser.add_argument("--temp_dir", type=str, default="/tmp", help="Temporary directory for intermediate files.")
    
    # Parse args
    args = parser.parse_args()    
    main(args)
