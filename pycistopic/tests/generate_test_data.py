#!/usr/bin/env python
"""Generate minimal test data for the pycisTopic pipeline.

Creates a small synthetic dataset with ~100 regions x 50 cells that can be used
to test build_cistopic_obj.py, run_models.py, and evaluate_models.py locally.

Usage:
    python tests/generate_test_data.py [--output_dir tests/test_data]

Only requires numpy, scipy, and pandas (no pycisTopic needed).
"""

import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse


def generate_test_data(output_dir: str) -> None:
    os.makedirs(output_dir, exist_ok=True)

    n_cells = 50
    n_regions = 100
    rng = np.random.default_rng(42)

    # --- Regions ---
    # Realistic region names: chr:start-end
    regions = []
    chrom = "chr1"
    pos = 1000
    for i in range(n_regions):
        start = pos + i * 2000
        end = start + 1000
        regions.append(f"{chrom}:{start}-{end}")
    regions_path = os.path.join(output_dir, "regions.tsv")
    with open(regions_path, "w") as f:
        for r in regions:
            f.write(r + "\n")

    # --- Barcodes ---
    barcodes = [f"CELL{i+1:03d}" for i in range(n_cells)]
    barcodes_path = os.path.join(output_dir, "barcodes.tsv")
    with open(barcodes_path, "w") as f:
        for b in barcodes:
            f.write(b + "\n")

    # --- Sparse count matrix (regions x cells) ---
    # ~5% density — realistic for scATAC-seq
    density = 0.05
    data = scipy.sparse.random(
        n_regions, n_cells, density=density, format="coo", random_state=42
    )
    # Scale to small integer counts (0-5)
    data.data = np.ceil(data.data * 5).astype(np.int32)
    data = data.tocsc()

    mtx_path = os.path.join(output_dir, "counts.mtx")
    scipy.io.mmwrite(mtx_path, data)

    # --- Cell metadata ---
    cell_types = rng.choice(["TypeA", "TypeB", "TypeC"], size=n_cells)
    meta = pd.DataFrame(
        {
            "barcode": barcodes,
            "cell_type": cell_types,
            "sample": ["test_sample"] * n_cells,
        }
    )
    meta = meta.set_index("barcode")
    meta_path = os.path.join(output_dir, "cell_metadata.tsv")
    meta.to_csv(meta_path, sep="\t")

    # --- Blacklist BED ---
    blacklist_regions = [
        "chr1:5000-6000",
        "chr1:15000-16000",
        "chr1:25000-26000",
        "chr1:35000-36000",
        "chr1:45000-46000",
    ]
    blacklist_path = os.path.join(output_dir, "blacklist.bed")
    with open(blacklist_path, "w") as f:
        for r in blacklist_regions:
            chrom_bl, coords = r.split(":")
            start_bl, end_bl = coords.split("-")
            f.write(f"{chrom_bl}\t{start_bl}\t{end_bl}\n")

    # --- H5AD (cells x regions) ---
    # AnnData expects cells-as-rows, regions-as-columns
    adata = ad.AnnData(
        X=data.T.tocsr(),  # cells x regions
        obs=meta.copy(),
        var=pd.DataFrame(index=regions),
    )
    h5ad_path = os.path.join(output_dir, "peaks.h5ad")
    adata.write_h5ad(h5ad_path)

    print(f"Test data generated in {output_dir}/")
    print(f"  counts.mtx       : {n_regions} regions x {n_cells} cells")
    print(f"  barcodes.tsv     : {n_cells} barcodes")
    print(f"  regions.tsv      : {n_regions} regions")
    print(f"  cell_metadata.tsv: {n_cells} rows, columns: cell_type, sample")
    print(f"  blacklist.bed    : {len(blacklist_regions)} regions")
    print(f"  peaks.h5ad       : {adata.shape[0]} cells x {adata.shape[1]} regions")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate test data for pycisTopic pipeline")
    parser.add_argument(
        "--output_dir",
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data"),
        help="Output directory for test data (default: tests/test_data/)",
    )
    args = parser.parse_args()
    generate_test_data(args.output_dir)
