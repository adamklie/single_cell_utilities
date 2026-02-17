#!/usr/bin/env python
"""Build a cisTopic object from a peak matrix (MTX or H5AD).

Outputs:
    {output_dir}/{project_name}.pkl  — pickled CistopicObject
"""

import argparse
import logging
import os
import pickle
import re
import sys

import numpy as np
import pandas as pd
from scipy.io import mmread

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

REGION_PATTERN = re.compile(r"^chr[\dXYMUn_]+[:\-][\d]+[:\-][\d]+$")


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------

def validate_mtx_args(args):
    """Ensure all three MTX files are provided and exist."""
    for name, path in [
        ("--matrix_path", args.matrix_path),
        ("--barcodes_path", args.barcodes_path),
        ("--regions_path", args.regions_path),
    ]:
        if path is None:
            logger.error("MTX mode requires %s", name)
            sys.exit(1)
        if not os.path.isfile(path):
            logger.error("File not found: %s (%s)", path, name)
            sys.exit(1)


def validate_h5ad_args(args):
    """Ensure H5AD file exists."""
    if not os.path.isfile(args.h5ad_path):
        logger.error("File not found: %s (--h5ad_path)", args.h5ad_path)
        sys.exit(1)


def check_region_format(regions, separator):
    """Warn if regions don't match expected chr:start-end pattern."""
    sample = regions[:min(20, len(regions))]
    n_match = sum(1 for r in sample if REGION_PATTERN.match(str(r)))
    if n_match < len(sample) * 0.5:
        logger.warning(
            "Most regions don't match the expected chr:start-end format "
            "(e.g. '%s'). If your separator is different, use --region_separator. "
            "Current separator: '%s'",
            sample[0],
            separator,
        )


def check_duplicates(values, label):
    """Warn about duplicate entries."""
    n_dup = len(values) - len(set(values))
    if n_dup > 0:
        logger.warning("%d duplicate %s found", n_dup, label)


def validate_blacklist(path):
    """Basic sanity check that a blacklist looks like a BED file."""
    if path is None:
        return
    if not os.path.isfile(path):
        logger.error("Blacklist file not found: %s", path)
        sys.exit(1)
    try:
        sample = pd.read_csv(
            path, sep="\t", header=None, nrows=5, comment="#"
        )
    except Exception as e:
        logger.error("Could not parse blacklist file: %s", e)
        sys.exit(1)
    if sample.shape[1] < 3:
        logger.error(
            "Blacklist file has %d columns; expected >= 3 (BED format)",
            sample.shape[1],
        )
        sys.exit(1)
    chrom_sample = sample.iloc[:, 0].astype(str)
    if not chrom_sample.str.startswith("chr").any():
        logger.warning(
            "Blacklist first column doesn't look like chromosome names: %s",
            chrom_sample.tolist(),
        )


def log_matrix_summary(label, n_regions, n_cells, nnz):
    """Log a summary of matrix dimensions."""
    total = n_regions * n_cells
    sparsity = (1 - nnz / total) * 100 if total > 0 else 0
    logger.info(
        "%s: n_regions=%d, n_cells=%d, n_nonzero=%d, sparsity=%.2f%%",
        label,
        n_regions,
        n_cells,
        nnz,
        sparsity,
    )


# ---------------------------------------------------------------------------
# Loading functions
# ---------------------------------------------------------------------------

def load_from_mtx(args):
    """Load matrix, barcodes, and regions from MTX files."""
    logger.info("Loading barcodes from %s", args.barcodes_path)
    bcs = pd.read_csv(args.barcodes_path, sep="\t", header=None)[0].values
    logger.info("  %d barcodes loaded", len(bcs))

    logger.info("Loading regions from %s", args.regions_path)
    regions = pd.read_csv(args.regions_path, sep="\t", header=None)[0].values
    logger.info("  %d regions loaded", len(regions))

    logger.info("Loading counts matrix from %s", args.matrix_path)
    cnt_mtx = mmread(args.matrix_path).tocsr()
    logger.info("  Raw matrix shape: %s", cnt_mtx.shape)

    # Auto-detect orientation
    if cnt_mtx.shape[0] == len(regions) and cnt_mtx.shape[1] == len(bcs):
        logger.info("  Matrix is regions x cells — no transpose needed")
    elif cnt_mtx.shape[0] == len(bcs) and cnt_mtx.shape[1] == len(regions):
        logger.info("  Matrix is cells x regions — transposing to regions x cells")
        cnt_mtx = cnt_mtx.T
    else:
        logger.error(
            "Matrix shape %s does not match barcodes (%d) x regions (%d) "
            "in either orientation",
            cnt_mtx.shape,
            len(bcs),
            len(regions),
        )
        sys.exit(1)

    return cnt_mtx, bcs, regions


def load_from_h5ad(args):
    """Load matrix, barcodes, and regions from an H5AD file."""
    import anndata

    logger.info("Loading AnnData from %s", args.h5ad_path)
    adata = anndata.read_h5ad(args.h5ad_path)
    logger.info("  AnnData shape: %s (obs x var)", adata.shape)

    bcs = adata.obs_names.values
    regions = adata.var_names.values

    # .X is cells x regions; we need regions x cells
    from scipy import sparse

    X = adata.X
    if sparse.issparse(X):
        cnt_mtx = X.T.tocsr()
    else:
        cnt_mtx = sparse.csr_matrix(X.T)
    logger.info("  Transposed to regions x cells: %s", cnt_mtx.shape)

    return cnt_mtx, bcs, regions


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(args):
    # --- Validate mutual exclusivity ---
    has_mtx = args.matrix_path is not None
    has_h5ad = args.h5ad_path is not None

    if has_mtx and has_h5ad:
        logger.error(
            "Provide EITHER MTX files (--matrix_path, --barcodes_path, "
            "--regions_path) OR --h5ad_path, not both"
        )
        sys.exit(1)
    if not has_mtx and not has_h5ad:
        logger.error(
            "Must provide either MTX files (--matrix_path + --barcodes_path + "
            "--regions_path) or --h5ad_path"
        )
        sys.exit(1)

    # --- Load data ---
    if has_mtx:
        validate_mtx_args(args)
        cnt_mtx, bcs, regions = load_from_mtx(args)
    else:
        validate_h5ad_args(args)
        cnt_mtx, bcs, regions = load_from_h5ad(args)

    # --- Guard rails ---
    check_duplicates(bcs, "barcodes")
    check_duplicates(regions, "regions")
    check_region_format(regions, args.region_separator)
    validate_blacklist(args.blacklist_path)

    log_matrix_summary("Before filtering", len(regions), len(bcs), cnt_mtx.nnz)

    # --- Build DataFrame for pycisTopic ---
    logger.info("Building sparse DataFrame...")
    df = pd.DataFrame.sparse.from_spmatrix(cnt_mtx)
    df.columns = bcs
    df.index = regions

    # --- Cell metadata ---
    if args.cell_metadata_path is not None:
        logger.info("Loading cell metadata from %s", args.cell_metadata_path)
        if not os.path.isfile(args.cell_metadata_path):
            logger.error(
                "Cell metadata file not found: %s", args.cell_metadata_path
            )
            sys.exit(1)
        ext = os.path.splitext(args.cell_metadata_path)[1].lower()
        if ext in (".tsv", ".txt"):
            cell_data = pd.read_csv(
                args.cell_metadata_path, sep="\t", index_col=0
            )
        elif ext == ".csv":
            cell_data = pd.read_csv(args.cell_metadata_path, index_col=0)
        else:
            logger.error("Cell metadata must be .tsv, .txt, or .csv (got %s)", ext)
            sys.exit(1)
        logger.info("  Loaded cell metadata: %s", cell_data.shape)

        overlap = df.columns.intersection(cell_data.index)
        n_overlap = len(overlap)
        n_total = len(df.columns)
        pct = n_overlap / n_total * 100 if n_total > 0 else 0
        logger.info(
            "  %d / %d barcodes (%.1f%%) have metadata",
            n_overlap,
            n_total,
            pct,
        )
        if pct < 50:
            logger.warning(
                "Less than 50%% of barcodes match cell metadata. "
                "Check that barcode formats match."
            )
        bcs = overlap
        df = df[bcs]
        logger.info("  Matrix shape after metadata filtering: %s", df.shape)

    # --- Create cisTopic object ---
    logger.info("Creating cisTopic object...")
    from pycisTopic.cistopic_class import create_cistopic_object

    cistopic_obj = create_cistopic_object(
        fragment_matrix=df,
        cell_names=bcs,
        region_names=regions,
        path_to_blacklist=args.blacklist_path,
        project=args.project_name,
    )

    if args.cell_metadata_path is not None:
        logger.info("Adding cell metadata to cisTopic object...")
        cell_data = cell_data.loc[bcs]
        cistopic_obj.add_cell_data(cell_data)

    # --- Final summary ---
    from scipy import sparse

    final_mtx = cistopic_obj.fragment_matrix
    if sparse.issparse(final_mtx):
        nnz = final_mtx.nnz
    else:
        nnz = np.count_nonzero(final_mtx)
    n_regions_final = final_mtx.shape[0]
    n_cells_final = final_mtx.shape[1]
    log_matrix_summary("After creation", n_regions_final, n_cells_final, nnz)

    # --- Save ---
    os.makedirs(args.output_dir, exist_ok=True)
    output_path = os.path.join(args.output_dir, cistopic_obj.project + ".pkl")
    logger.info("Saving cisTopic object to %s", output_path)
    with open(output_path, "wb") as f:
        pickle.dump(cistopic_obj, f)
    logger.info("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build a cisTopic object from a peak matrix (MTX or H5AD).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # From MTX files
  python build_cistopic_obj.py \\
      --matrix_path counts.mtx --barcodes_path barcodes.tsv --regions_path regions.tsv \\
      --output_dir ./results --project_name my_sample

  # From H5AD
  python build_cistopic_obj.py \\
      --h5ad_path peaks.h5ad \\
      --output_dir ./results --project_name my_sample
""",
    )

    # --- MTX mode ---
    mtx_group = parser.add_argument_group("MTX mode (mutually exclusive with H5AD)")
    mtx_group.add_argument(
        "--matrix_path",
        type=str,
        default=None,
        help="Path to counts matrix in .mtx format",
    )
    mtx_group.add_argument(
        "--barcodes_path",
        type=str,
        default=None,
        help="Path to barcodes file (one per line, no header)",
    )
    mtx_group.add_argument(
        "--regions_path",
        type=str,
        default=None,
        help="Path to regions file (one per line, no header)",
    )

    # --- H5AD mode ---
    h5ad_group = parser.add_argument_group("H5AD mode (mutually exclusive with MTX)")
    h5ad_group.add_argument(
        "--h5ad_path",
        type=str,
        default=None,
        help="Path to .h5ad AnnData file (cells x regions)",
    )

    # --- Required ---
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory (created if needed). Object saved as {project_name}.pkl",
    )

    # --- Optional ---
    parser.add_argument(
        "--project_name",
        type=str,
        default="cistopic",
        help="Project name, used as output filename (default: cistopic)",
    )
    parser.add_argument(
        "--blacklist_path",
        type=str,
        default=None,
        help="Path to blacklist regions BED file",
    )
    parser.add_argument(
        "--cell_metadata_path",
        type=str,
        default=None,
        help="Path to cell metadata (.tsv, .txt, or .csv with index column)",
    )
    parser.add_argument(
        "--region_separator",
        type=str,
        default=":-",
        help="Expected separator pattern in region names (default: ':' and '-')",
    )

    args = parser.parse_args()
    main(args)
