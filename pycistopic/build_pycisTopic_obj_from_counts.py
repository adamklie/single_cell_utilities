import argparse
import logging
import os
import pickle


def configure_logger():
    logging.basicConfig(
        format="%(asctime)s %(levelname)-8s %(message)s",
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger(__name__)
    return logger


def main(args):
    # Get arguments
    counts_matrix_path = args.counts_matrix_path
    barcodes_path = args.barcodes_path
    regions_path = args.regions_path
    cell_metadata_path = args.cell_metadata_path
    blacklist_path = args.blacklist_path
    out_dir = args.out_dir
    project_name = args.project_name

    # Log arguments
    logger = configure_logger()
    logger.info("counts_matrix_path: %s", args.counts_matrix_path)
    logger.info("barcodes_path: %s", args.barcodes_path)
    logger.info("regions_path: %s", args.regions_path)
    logger.info("cell_metadata_path: %s", args.cell_metadata_path)
    logger.info("blacklist_path: %s", args.blacklist_path)
    logger.info("out_dir: %s", args.out_dir)
    logger.info("project_name: %s", args.project_name)

    logger.info("Loading barcodes...")
    import pandas as pd

    bcs = pd.read_csv(barcodes_path, sep="\t", header=None)[0].values
    logger.info("{} barcodes loaded: {}...".format(len(bcs), bcs[:5]))

    logger.info("Loading regions...")
    regions = pd.read_csv(regions_path, sep="\t", header=None)[0].values
    logger.info("{} regions loaded: {}...".format(len(regions), regions[:5]))

    logger.info("Loading counts matrix...")
    from scipy.io import mmread

    cnt_mtx = mmread(counts_matrix_path).tocsr()
    logger.info("Counts matrix loaded: %s...", cnt_mtx.shape)

    if cell_metadata_path is not None:
        logger.info("Loading cell metadata...")
        if cell_metadata_path.endswith(".tsv") or cell_metadata_path.endswith(".txt"):
            cell_data = pd.read_csv(cell_metadata_path, sep="\t", index_col=0)
        elif cell_metadata_path.endswith(".csv"):
            cell_data = pd.read_csv(cell_metadata_path, index_col=0)
        else:
            raise ValueError("cell_metadata_path must be .tsv, .txt, or .csv")
        logger.info(f"Cells metadata loaded with shape {cell_data.shape}")

    logger.info("Creating cisTopic object...")
    from pycisTopic.cistopic_class import create_cistopic_object
    logger.info("Checking if counts matrix needs to be transposed...")
    if cnt_mtx.shape[0] == len(regions) and cnt_mtx.shape[1] == len(bcs):
        logger.info("Counts matrix is already in the correct shape")
    elif cnt_mtx.shape[1] == len(regions) and cnt_mtx.shape[0] == len(bcs):
        logger.info("Transposing counts matrix to match barcodes and regions numbers...")
        cnt_mtx = cnt_mtx.T
    else:
        raise ValueError(
            "Counts matrix shape does not match either barcodes or regions list"
        )
    df = pd.DataFrame.sparse.from_spmatrix(cnt_mtx)
    df.columns = bcs
    df.index = regions
    if cell_metadata_path is not None:
        bcs = df.columns.intersection(cell_data.index)
        logger.info(f"{len(bcs)} cells have metadata...")
        df = df[bcs]
        logger.info(f"Counts matrix shape after filtering: {df.shape}")
    cistopic_obj = create_cistopic_object(
        fragment_matrix=df,
        cell_names=bcs,
        region_names=regions,
        path_to_blacklist=blacklist_path,
        project=project_name,
    )

    if cell_metadata_path is not None:
        logger.info("Adding cell metadata to cisTopic object...")
        cell_data = cell_data.loc[bcs]
        cistopic_obj.add_cell_data(cell_data)

    output_path = os.path.join(out_dir, cistopic_obj.project + ".pkl")
    logger.info("Saving cisTopic object to %s", output_path)
    pickle.dump(cistopic_obj, open(output_path, "wb"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build cisTopic object from a prespecified counts matrix"
    )
    parser.add_argument(
        "--counts_matrix_path",
        help="Path to counts matrix file in .mtx format",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--barcodes_path",
        help="Path to barcodes list as single column .csv file with no header",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--regions_path",
        help="Path to regions list as single column .csv file with no header",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--out_dir",
        help="Path to output directory, output will be saved as <project_name>.pkl",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--cell_metadata_path",
        help="Optional path to cell metadata to add to cisTopic object",
        type=str,
        required=False,
        default=None,
    )
    parser.add_argument(
        "--blacklist_path",
        help="Path to blacklist regions bed file",
        type=str,
        required=False,
        default=None,
    )
    parser.add_argument(
        "--project_name",
        help="Name of project, will be used as name of output file",
        type=str,
        default="cistopic",
    )
    args = parser.parse_args()
    main(args)
