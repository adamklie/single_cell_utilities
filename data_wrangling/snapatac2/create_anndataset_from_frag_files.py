# This script is meant for reading in a a set of passde in fragment files
# and generating a SnapATAC2 AnnDataset object written to disk

# This script uses argparse to parse arguments
# for specifying fragment files, output directory, and optional parameters for data processing.

# The script will:
# 1. Retrieve the list of fragment files.
# 2. For each fragment file, it checks if the processed AnnData object exists in the output directory.
#    - If it exists, it reads the AnnData object directly.
#    - If not, it processes the fragment file to create the AnnData object and saves it to the output directory.
# 3. After processing all fragment files, it merges all the AnnData objects into a single merged AnnData object.
# 4. The merged AnnData object is saved to the output directory.
# 
# Usage:
# python script_name.py --fragment_files [list_of_fragment_files] --out_dir [output_directory] [optional_arguments]


import os
import glob
import time
import logging
import argparse
from tqdm.auto import tqdm

def main(args):
    
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    
    # Log arguments
    logging.info("Arguments: %s", args)

    # Get list of fragment files from argparse
    fragment_files = args.fragment_files
    
    # Read in all the fragment files into separate anndatas
    import snapatac2 as snap
    adata_atac_list = []
    for i, frag_file in tqdm(enumerate(fragment_files), total=len(fragment_files)):
        logging.info("Current fragment file: %s", frag_file)
        if os.path.exists(os.path.join(args.out_dir, f"adata_atac_{i+1}.h5ad")):
            logging.info("Reading: %s", os.path.join(args.out_dir, f"adata_atac_{i+1}.h5ad"))
            adata_atac = snap.read(os.path.join(args.out_dir, f"adata_atac_{i+1}.h5ad"))
        elif os.path.exists(os.path.join(frag_file)):
            logging.info("Processing: %s", frag_file)
            time_in = time.time()
            adata_atac = snap.pp.import_data(
                fragment_file=frag_file,
                genome=snap.genome.hg38,
                min_tsse=args.min_tsse,
                min_num_fragments=args.min_num_fragments,
                file=os.path.join(args.out_dir, f"adata_atac_{i+1}.h5ad"),
                sorted_by_barcode=args.sorted_by_barcode,
                chunk_size=args.chunk_size,
                low_memory=args.low_memory
            )
            time_out = time.time()
            logging.info("Processing time: %s", time_out - time_in)
        adata_atac_list.append(adata_atac)

    # Merge
    time_in = time.time()
    adata_atac_merged = snap.AnnDataSet(
        adatas=[(name, adata) for name, adata in zip(fragment_files, adata_atac_list)],
        filename=os.path.join(args.out_dir, "adata_atac_merged.h5ads")
    )
    time_out = time.time()
    logging.info("Time elapsed for merging: %s", time_out - time_in)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process and merge ATAC fragment files.")
    
    # Required args
    parser.add_argument("--fragment_files", type=str, required=True, nargs='+', help="List of paths to ATAC fragment files.")
    parser.add_argument("--out_dir", type=str, required=True, help="Directory to save h5ad files.")

    # Optional args for import_data function
    parser.add_argument("--min_tsse", type=int, default=7, help="Minimum tsse value.")
    parser.add_argument("--min_num_fragments", type=int, default=1000, help="Minimum number of fragments.")
    parser.add_argument("--sorted_by_barcode", action="store_true", help="If the data is sorted by barcode.")
    parser.add_argument("--chunk_size", type=int, default=100000, help="Size of data chunks.")
    parser.add_argument("--low_memory", action="store_true", help="Low memory mode.")
    args = parser.parse_args()
    
    main(args)