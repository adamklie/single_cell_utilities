import argparse
import logging
import os
import pickle


def configure_logger():
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)
    return logger

def main(args):
    input_path = args.input_path
    input_files = args.input_files
    out_dir = args.output_dir
    tmp_dir = args.tmp_dir
    ctx_db = args.ctx_db
    motif_annot = args.motif_annotation
    annotation_version = args.annotation_version
    species = args.species
    num_cpu = int(args.num_cpu)

    # Configure logger
    logger = configure_logger()
    logger.info("Using the following arguments:")
    logger.info(f"input_path: {input_path}")
    logger.info(f"input_files: {input_files}")
    logger.info(f"out_dir: {out_dir}")
    logger.info(f"tmp_dir: {tmp_dir}")
    logger.info(f"ctx_db: {ctx_db}")
    logger.info(f"motif_annot: {motif_annot}")
    logger.info(f"annotation_version: {annotation_version}")
    logger.info(f"species: {species}")
    logger.info(f"num_cpu: {num_cpu}")

    # Import packages
    logger.info("Importing packages")
    import pandas as pd
    import pyranges as pr
    from pycistarget.motif_enrichment_cistarget import run_cistarget

    # Read input files
    if input_path:
        logger.info(f"Reading all .narrowPeak files from input directory: {input_path}")
        region_set_files = [os.path.join(input_path, f) for f in os.listdir(input_path) if f.endswith(".narrowPeak")]
        input_file_exts = [".narrowPeak"] * len(region_set_files)
    elif input_files:
        logger.info(f"Reading input files: {input_files}")
        region_set_files = input_files
        input_file_exts = [os.path.splitext(x)[1] for x in input_files]
    else:
        raise ValueError("Either input_path or input_files must be set.")
    region_sets = {x.split("/")[-1].replace(input_file_exts[i], ""): pr.read_bed(x) for i, x in enumerate(region_set_files)}
    print(region_sets.keys(), region_sets.values())
    
    # Create output directory
    if not os.path.exists(out_dir):
        logger.info(f"Creating output directory: {out_dir}")
        os.makedirs(out_dir)

    # Run cisTarget
    logger.info("Running cisTarget")
    cistarget_dict = run_cistarget(
        ctx_db=ctx_db,
        region_sets=region_sets,
        specie=species,
        annotation=['Direct_annot', 'Orthology_annot'],
        annotation_version=annotation_version,
        path_to_motif_annotations=motif_annot,
        n_cpu=num_cpu,
        _temp_dir=tmp_dir
    )
    logger.info("Finished running cisTarget")

    # Save to output directory
    pkl_path = os.path.join(out_dir, "cisTarget.pkl")
    logger.info(f"Saving to {pkl_path}")
    with open(pkl_path, 'wb') as f:
        pickle.dump(cistarget_dict, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pycisTarget with provided input paths and output directory.")
    parser.add_argument("--input_path", help="Path to the directory containing region set files. If not specified, input files must be set", nargs='?', default=None)
    parser.add_argument("--input_files", help="List of input file paths.", nargs='*', default=None)
    parser.add_argument("--output_dir", help="Output directory for saving the results.", default="./")
    parser.add_argument("--ctx_db", help="Path to the cisTarget database file.", default="/cellar/users/aklie/opt/shared/SCENIC+/databases/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather")
    parser.add_argument("--motif_annotation", help="Path to the motif annotation file.", default="/cellar/users/aklie/opt/shared/SCENIC+/motif_annotation/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
    parser.add_argument("--annotation_version", help="Motif annotation version to use", default="v10nr_clust")
    parser.add_argument("--species", help="Species to use", default="homo_sapiens")
    parser.add_argument("--tmp_dir", help="Temporary directory.", default="/cellar/users/aklie/tmp")
    parser.add_argument("--num_cpu", help="Number of CPUs to use.", default=4)
    parser.add_argument("--verbosity", default=1)
    
    args = parser.parse_args()
    main(args)
