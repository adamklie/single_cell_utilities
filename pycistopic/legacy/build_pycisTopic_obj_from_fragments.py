import argparse
import logging
import os
import pickle

import numpy as np
import pandas as pd
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments


def configure_logger():
    logging.basicConfig(
        format="%(asctime)s %(levelname)-8s %(message)s",
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger(__name__)
    return logger


def main(args):
    # Define arguments
    fragments_dir = args.fragments_dir
    fragment_files = args.fragment_files
    out_dir = args.out_dir
    barcodes_file = args.barcode_list
    regions_file = args.regions_list
    cell_metadata_file = args.cell_metadata_file
    blacklist_file = args.blacklist_file
    split_pattern = args.split_pattern
    project_name = args.project_name
    n_cpu = args.n_cpu

    # Log arguments
    logger = configure_logger()
    logger.info(
        "fragments_dir: %s", fragments_dir
    ) if fragments_dir is not None else None
    logger.info(
        "fragment_files: %s", fragment_files
    ) if fragment_files is not None else None
    logger.info("out_dir: %s", out_dir)
    logger.info("barcodes_file: %s", barcodes_file)
    logger.info("regions_file: %s", regions_file)
    logger.info("cell_metadata_file: %s", cell_metadata_file)
    logger.info("blacklist_file: %s", blacklist_file)
    logger.info("split_pattern: %s", split_pattern)
    logger.info("project_name: %s", project_name)
    logger.info("n_cpu: %s", n_cpu)

    # Create a fragment file dictionary for cisTopic
    logger.info("Creating fragment file dictionary")
    if fragments_dir is not None:
        fragment_files = [
            os.path.join(fragments_dir, f) for f in os.listdir(fragments_dir)
        ]
        fragment_file_dict = {
            os.path.basename(f).split(split_pattern)[0]: f for f in fragment_files
        }
    elif fragment_files is not None:
        fragment_file_dict = {
            os.path.basename(f).split(split_pattern)[0]: f for f in fragment_files
        }
    else:
        raise ValueError("Must specify either fragments_dir or fragment_files")
    logger.info("Fragment file dictionary created: %s", fragment_file_dict)

    # Load the path to regions, if more than one fragment file, use the same regions for all
    path_to_regions = dict.fromkeys(fragment_file_dict.keys(), regions_file)

    logger.info("Loading barcodes...")
    bcs = pd.read_csv(barcodes_file, sep="\t", header=None)[0].values
    logger.info("{} barcodes loaded: {}...".format(len(bcs), bcs[:5]))

    logger.info("Creating cisTopic object...")
    cistopic_obj_list = [
        create_cistopic_object_from_fragments(
            path_to_fragments=fragment_file_dict[key],
            path_to_regions=path_to_regions,
            path_to_blacklist=blacklist_file,
            valid_bc=bcs,
            n_cpu=n_cpu,
            project=key,
            split_pattern=split_pattern,
        )
        for key in fragment_file_dict.keys()
    ]
    if len(cistopic_obj_list) > 1:
        logger.info("Merging cisTopic objects...")
        cistopic_obj = merge_cistopic_objects(cistopic_obj_list)
    else:
        logger.info("Only one fragment file, no need to merge")
        cistopic_obj = cistopic_obj_list[0]

    if cell_data is not None:
        logger.info("Loading cell metadata...")
        cell_data = pd.read_csv(cell_metadata_path, index_col=0)
        logger.info(f"Cells metadata loaded with shape {cell_data.shape}")
        cistopic_obj.add_cell_data(cell_data)

    output_path = os.path.join(out_dir, cistopic_obj.project + ".pkl")
    logger.info("Saving cisTopic object to %s", output_path)
    pickle.dump(cistopic_obj, open(output_path, "wb"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create cisTopic object from fragment files"
    )
    parser.add_argument(
        "--fragments_dir",
        type=str,
        help="Directory containing a set of fragment files. Must be specified if not using --fragment_files",
        required=False,
    )
    parser.add_argument(
        "--fragment_files",
        type=str,
        nargs="+",
        help="List of fragment files. Must be specified if not using --fragments_dir",
        required=False,
    )
    parser.add_argument(
        "--barcodes_file",
        help="Path to barcodes list as single column .csv file with no header",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--regions_file",
        type=str,
        help="Path to regions list in csv or bed format. If csv, this should be a single column with one region per line and no header. If bed format, this should be a bed file with no header",
        required=True,
    )
    parser.add_argument(
        "--cell_metadata_file",
        type=str,
        help="Path to cell metadata file in csv format. This should be a table with cells as rows and metadata as columns. The first column should be the cell barcodes",
        required=True,
    )
    parser.add_argument(
        "--out_dir",
        type=str,
        help="Path to output directory, will save as <project_name>.pkl",
        required=True,
    )
    parser.add_argument(
        "--blacklist_file",
        type=str,
        help="Path to blacklist regions in bed format",
        required=False,
        default="/cellar/users/aklie/opt/pycisTopic/blacklist/hg38-blacklist.v2.bed",
    )
    parser.add_argument(
        "--split_pattern",
        type=str,
        help="Pattern to split fragment file names on. Default is '.'",
        required=False,
        default=".",
    )
    parser.add_argument(
        "--project_name",
        type=str,
        help="Name of project",
        required=False,
        default="cistopic",
    )
    parser.add_argument(
        "--n_cpu",
        type=int,
        help="Number of CPUs to use. Use this if you are loading from multiple fragment files at a time",
        required=False,
        default=1,
    )
    args = parser.parse_args()
    main(args)
