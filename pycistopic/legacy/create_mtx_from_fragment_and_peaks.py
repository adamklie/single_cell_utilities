import argparse
import logging
import os
import sys
import time
from typing import List, Tuple

import numpy as np
import pandas as pd
import pyranges as pr
from pycisTopic.utils import (collapse_duplicates, get_position_index,
                              read_fragments_from_file)
from scipy import sparse
from scipy.io import mmwrite

sys.path.append("/cellar/users/aklie/data/igvf/bin")
from utils import make_dirs


def create_fragment_dataframe(
    fragment_file: str,
    regions: pr.PyRanges,
    check_for_duplicates: bool,
    barcodes: np.ndarray = None,
    ncpu: int = 1,
) -> sparse.csr_matrix:
    """
    Create a fragment matrix from the provided fragment file and regions.

    Parameters
    ----------
    fragment_file : str
        Path to the fragment file.
    regions : pr.PyRanges
        PyRanges object of regions.
    check_for_duplicates : bool
        Whether to check for duplicate fragments.
    cell_data : pd.DataFrame
        Cell data.

    Returns
    -------
    csr_matrix
        Sparse fragment matrix.
    """

    # Read in fragments as pyranges
    print(f"Reading fragments from {fragment_file}")
    fragments = read_fragments_from_file(fragment_file, use_polars=True)

    # Keep only fragments that map to those barcodes
    print(f"Keeping only fragments that map to barcodes")
    if barcodes:
        fragments = fragments.subset(lambda df: df.Name.isin(set(barcodes)))

    # Check for duplicates
    print(f"Checking for duplicates")
    if "Score" not in fragments.df:
        fragments_df = fragments.df
        if check_for_duplicates:
            fragments_df = pd.concat(
                [
                    collapse_duplicates(fragments_df[fragments_df.Chromosome == x])
                    for x in fragments_df.Chromosome.cat.categories.values
                ]
            )
        else:
            fragments_df["Score"] = 1
        fragments = pr.PyRanges(fragments_df)

    # Joing fragments and peaks
    print(f"Joining fragments and peaks")
    fragments_in_regions = regions.join(fragments, nb_cpu=ncpu)

    # Create a dataframe of counts (long dataframe)
    print(f"Creating a long dataframe of counts")
    counts_df = pd.concat(
        [
            fragments_in_regions.regionID.astype("category"),
            fragments_in_regions.Name.astype("category"),
            fragments_in_regions.Score.astype(np.int32),
        ],
        axis=1,
        sort=False,
    )

    # Create a matrix
    print(f"Creating dataframe of counts (wannabe matrix)")
    fragment_matrix = (
        counts_df.groupby(["Name", "regionID"], sort=False, observed=True)
        .size()
        .unstack(level="Name", fill_value=0)
        .astype(np.int32)
    )
    fragment_matrix.columns.names = [None]

    # Return
    return fragment_matrix


def merge_sparse_matrices(
    matrix_list: List[sparse.csr_matrix],
    column_names_list: List[np.ndarray],
    row_names_list: List[np.ndarray],
):
    # Grab the first one as a starting point
    matrix = matrix_list[0]
    row_names = row_names_list[0]
    column_names = column_names_list[0]

    # Loop through the rest
    for i in range(1, len(row_names_list)):
        # Deal with rows
        row_names_to_add = row_names_list[i]
        common_rows = list(set(row_names) & set(row_names_to_add))
        diff_rows = list(set(row_names) ^ set(row_names_to_add))
        common_index_fm = get_position_index(common_rows, row_names)
        common_index_fm_to_add = get_position_index(common_rows, row_names_to_add)

        # Deal with columns
        column_names_to_add = column_names_list[i]
        column_names = column_names + column_names_to_add

        # Deal with fragment matrix
        matrix_to_add = matrix_list[i]

        # Add common rows
        matrix_common = sparse.hstack(
            [matrix[common_index_fm,], matrix_to_add[common_index_fm_to_add,]]
        )

        # Add different rows
        if len(diff_rows) > 0:
            # Ones that are in the original fragment matrix but not in the new one
            diff_rows_1 = list(np.setdiff1d(row_names, row_names_to_add))
            diff_index_fm_1 = get_position_index(diff_rows_1, row_names)
            matrix_diff_1 = sparse.hstack(
                [
                    matrix[diff_index_fm_1,],
                    np.zeros((len(diff_rows_1), matrix_to_add.shape[1])),
                ]
            )

            # Ones that are in the new fragment matrix but not in the original one
            diff_rows_2 = list(np.setdiff1d(row_names_to_add, row_names))
            diff_index_fm_2 = get_position_index(diff_rows_2, row_names_to_add)
            matrix_diff_2 = sparse.hstack(
                [
                    np.zeros((len(diff_rows_2), matrix.shape[1])),
                    matrix_to_add[diff_index_fm_2,],
                ]
            )

            # Put everything together
            matrix = sparse.vstack(
                [
                    matrix_common,
                    matrix_diff_1,
                    matrix_diff_2,
                ]
            )
            row_names = common_rows + diff_rows_1 + diff_rows_2

        # Otherwise, just add the common rows
        else:
            matrix = matrix_common
            row_names = common_rows

        # Convert to CSR matrix
        matrix = sparse.csr_matrix(matrix, dtype=np.int32)

    return matrix, row_names, column_names


def main(args) -> None:
    """
    Main function to execute the script.
    """
    logging.basicConfig(level=logging.INFO)
    logging.info(f"Arguments: {args}")

    # Parse arguments
    fragment_files = args.fragment_files
    peak_file = args.peak_file
    out_dir = args.out_dir
    blacklist_file = args.blacklist_file
    barcodes_file = args.barcodes_file
    check_for_duplicates = args.check_for_duplicates
    write_intermediates = args.write_intermediates
    ncpu = args.ncpu

    # Grab some ids
    sample_ids = [
        os.path.basename(file).replace(".tsv.gz", "") for file in fragment_files
    ]

    # Make output directory if it doesn't exist, otherwise exit
    if os.path.exists(out_dir):
        logging.error(f"Output directory already exists: {out_dir}")
        sys.exit(1)
    else:
        make_dirs(out_dir)

    logging.info(f"Reading regions from {peak_file}")
    regions = pr.read_bed(peak_file)
    regions = regions[["Chromosome", "Start", "End"]]

    # Read in the blacklist file as a PyRanges object
    if args.blacklist_file:
        logging.info(f"Reading blacklist file from {blacklist_file}")
        blacklist = pr.read_bed(blacklist_file)
        pre_regions = len(regions)
        regions = regions.overlap(blacklist, invert=True)
        logging.info(f"Removed {pre_regions - len(regions)} regions from blacklist.")
    selected_regions = (
        regions.Chromosome.astype(str)
        + ":"
        + regions.Start.astype(str)
        + "-"
        + regions.End.astype(str)
    ).to_list()
    regions.regionID = selected_regions
    region_names = selected_regions
    logging.info(f"Number of regions: {len(regions)}")

    # Read in barcodes
    if args.barcodes_file:
        logging.info(f"Reading barcodes from {barcodes_file}")
        barcodes = pd.read_csv(barcodes_file, header=None)[0].to_list()
    else:
        barcodes = None

    # Sequentially read in fragments and create matrices
    fragment_matrix_list = []
    cell_names_list = []
    region_names_list = []
    for i, frag_file in enumerate(fragment_files):
        # Create fragment dataframe
        logging.info(f"Creating fragment dataframe from {frag_file}")
        curr_frag_matrix = create_fragment_dataframe(
            fragment_file=frag_file,
            regions=regions,
            check_for_duplicates=check_for_duplicates,
            barcodes=barcodes,
            ncpu=ncpu,
        )

        # Grab cell barcodes and region names
        curr_cell_names = curr_frag_matrix.columns.values.to_list()
        curr_region_names = curr_frag_matrix.index.values.to_list()
        print(len(curr_cell_names), len(curr_region_names))

        # Make it sparse
        curr_fragment_matrix = sparse.csr_matrix(
            curr_frag_matrix.to_numpy(), dtype=np.int32
        )
        print(curr_fragment_matrix.shape)

        # Add to growing lists
        cell_names_list.append(curr_cell_names)
        region_names_list.append(curr_region_names)
        fragment_matrix_list.append(curr_fragment_matrix)

        # Write intermediate files
        if write_intermediates:
            logging.info(f"Writing intermediate for {sample_ids[i]}")
            mmwrite(
                target=os.path.join(out_dir, f"{sample_ids[i]}.mtx"),
                a=curr_fragment_matrix,
            )
            with open(os.path.join(out_dir, f"{sample_ids[i]}_barcodes.tsv"), "w") as f:
                f.write("\n".join(curr_cell_names))
            with open(os.path.join(out_dir, f"{sample_ids[i]}_features.tsv"), "w") as f:
                f.write("\n".join(curr_region_names))

    # Concatenate matrices
    logging.info(f"Concatenating matrices")
    merged_matrix, region_names, cell_names = merge_sparse_matrices(
        matrix_list=fragment_matrix_list,
        column_names_list=cell_names_list,
        row_names_list=region_names_list,
    )

    # Write final matrix
    logging.info(f"Writing final matrix to {out_dir}")
    mmwrite(
        target=os.path.join(out_dir, "matrix.mtx"),
        a=merged_matrix,
    )
    with open(os.path.join(out_dir, "barcodes.tsv"), "w") as f:
        f.write("\n".join(cell_names))
    with open(os.path.join(out_dir, "features.tsv"), "w") as f:
        f.write("\n".join(region_names))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a matrix from peaks and fragments."
    )
    parser.add_argument(
        "--fragment_files",
        type=str,
        required=True,
        nargs="+",
        help="List of paths to ATAC fragment files.",
    )
    parser.add_argument("--peak_file", required=True, help="Path to the peak file.")
    parser.add_argument(
        "--out_dir", type=str, required=True, help="Directory to final mtx file."
    )
    parser.add_argument(
        "--blacklist_file", required=False, help="Path to the blacklist file."
    )
    parser.add_argument(
        "--barcodes_file", required=False, help="Path to the barcodes to keep."
    )
    parser.add_argument(
        "--check_for_duplicates",
        action="store_true",
        help="Check for duplicate fragments.",
    )
    parser.add_argument(
        "--write_intermediates", action="store_true", help="Write intermediate files."
    )
    parser.add_argument("--ncpu", type=int, default=1, help="Number of CPUs to use.")
    args = parser.parse_args()

    main(args)
