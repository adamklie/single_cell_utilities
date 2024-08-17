# This script is meant for reading in a a set of passde in fragment files
# and generating a SnapATAC2 adatas

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
# python script_sample_id.py --fragment_files [list_of_fragment_files] --sample_ids [sample_ids] --out_dir [output_directory] [optional_arguments]

import os
import time
import logging
import argparse
from tqdm.auto import tqdm

def main(args):
    
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    
    # Log arguments
    logging.info("Arguments: %s", args)

    # Parse args
    frag_file = args.fragment_file
    out_file = args.out_file
    
    # Read in all the fragment file into anndata
    import snapatac2 as snap
    loaded = False
    if os.path.exists(out_file):
        logging.info("%s already exists. Trying to load it", out_file)
        try:
            adata_atac = snap.read(out_file)
            logging.info("Successfully loaded %s file will not reprocess it. Delete it if you want to reprocess", out_file)
            loaded = True
        except:
            logging.info("Failed to load %s. Deleting file and reprocessing it", out_file)
            os.remove(out_file)
    if not loaded:
        logging.info("Processing: %s to %s", frag_file, out_file)
        time_in = time.time()
        adata_atac = snap.pp.import_data(
            fragment_file=frag_file,
            genome=snap.genome.hg38,
            min_tsse=args.min_tsse,
            min_num_fragments=args.min_num_fragments,
            file=out_file,
            sorted_by_barcode=args.sorted_by_barcode,
            chunk_size=args.chunk_size,
            low_memory=args.low_memory
        )
        adata_atac.close()
        time_out = time.time()
        logging.info("Processing time: %s", time_out - time_in)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process and merge a single ATAC fragment file to ann AnnData object.")
    
    # Required args
    parser.add_argument("--fragment_file", type=str, required=True, help="Path to fragment file.")
    parser.add_argument("--out_file", type=str, required=True, help="Directory to save h5ad files.")

    # Optional args for import_data function
    parser.add_argument("--min_tsse", type=int, default=7, help="Minimum tsse value.")
    parser.add_argument("--min_num_fragments", type=int, default=1000, help="Minimum number of fragments.")
    parser.add_argument("--sorted_by_barcode", action="store_true", help="If the data is sorted by barcode.")
    parser.add_argument("--chunk_size", type=int, default=100000, help="Size of data chunks.")
    parser.add_argument("--low_memory", action="store_true", help="Low memory mode.")
    args = parser.parse_args()
    
    main(args)
