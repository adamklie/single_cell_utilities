import os
import pickle
import warnings
import argparse
import logging
import time
from pycisTopic.cistopic_class import *
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *


def main(args):
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    
    # Suppress warnings and stderr messages
    warnings.simplefilter(action='ignore', category=FutureWarning)
    
    # Set directories
    output_dir = args.output_dir
    tmp_dir = args.tmp_dir

    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    logging.info(f"Output will be saved to {output_dir}.")

    # Load the object
    logging.info("Loading cisTopic object...")
    start_time = time.time()
    cistopic_obj = pickle.load(open(args.cistopic_obj_path, 'rb'))
    logging.info(f"Loaded cisTopic object in {time.time() - start_time:.2f} seconds.")

    # Binarize topics to get regions
    logging.info("Binarizing regions to create sets of candidate enhancers...")
    start_time = time.time()
    logging.info("Binarizing regions using Otsu's method...")
    region_bin_topics_otsu = binarize_topics(
        cistopic_obj, 
        method='otsu',
        save=os.path.join(output_dir, 'region_bin_topics_otsu.png')
    )
    logging.info(f"Binarizing regions using top {args.ntop} topics...")
    region_bin_topics_top3k = binarize_topics(
        cistopic_obj, 
        method='ntop', 
        ntop=args.ntop,
        save=os.path.join(output_dir, 'region_bin_topics_top3k.png')
    )
    logging.info(f"Binarized regions in {time.time() - start_time:.2f} seconds.")

    # Save these
    pickle.dump(region_bin_topics_otsu, open(os.path.join(output_dir, 'region_bin_topics_otsu.pkl'), 'wb'))
    pickle.dump(region_bin_topics_top3k, open(os.path.join(output_dir, 'region_bin_topics_top3k.pkl'), 'wb'))

    # Impute accessibility
    logging.info("Imputing accessibility...")
    start_time = time.time()
    imputed_acc_obj = impute_accessibility(
        cistopic_obj,  # cisTopic object
        selected_cells=None,  # certain cells to use for imputation
        selected_regions=None,  # certain regions to use for imputation
        scale_factor=args.scale_factor_impute,  # A number to multiply the imputed values for.
        project=cistopic_obj.project + "impute"  # A string to add to the name of the imputed object.
    )
    logging.info(f"Imputed accessibility in {time.time() - start_time:.2f} seconds.")

    # Find variable regions
    logging.info("Get DARs as candidate enhancers...")
    start_time = time.time()
    normalized_imputed_acc_obj = normalize_scores(
        imputed_acc_obj, 
        scale_factor=args.scale_factor_normalize  # Scale factor for cell-level normalization
    )
    variable_regions = find_highly_variable_features(
        normalized_imputed_acc_obj, 
        plot=False
        save=os.path.join(output_dir, 'variable_regions.png')
    )
    del normalized_imputed_acc_obj
    logging.info(f"Found DARs in {time.time() - start_time:.2f} seconds.")

    # Get marker regions and save
    logging.info("Getting marker regions...")
    start_time = time.time()
    markers_dict = find_diff_features(
        cistopic_obj, 
        imputed_acc_obj, 
        variable=args.dar_column, 
        var_features=variable_regions, 
        split_pattern=args.split_pattern
    )
    pickle.dump(markers_dict, open(os.path.join(output_dir, 'markers_dict.pkl'), 'wb'))
    logging.info(f"Got marker regions in {time.time() - start_time:.2f} seconds.")

if __name__ == "__main__":
    # Set up argparse
    parser = argparse.ArgumentParser(description="Process cisTopic objects and generate candidate enhancers.")
    parser.add_argument("--cistopic_obj_path", required=True, help="Path to the cistopic object.")
    parser.add_argument("--dar_column", required=True, help="Column name to use for finding DARs.")
    parser.add_argument("--output_dir", required=True, help="Directory to save the results.")
    parser.add_argument("--tmp_dir", required=False, default="/cellar/users/aklie/tmp/", help="Temporary directory for intermediate files.")
    parser.add_argument("--ntop", required=False, type=int, default=3000, help="Number of top topics for binarization.")
    parser.add_argument("--scale_factor_impute", required=False, type=float, default=10**6, help="Scale factor for imputing accessibility.")
    parser.add_argument("--scale_factor_normalize", required=False, type=float, default=10**4, help="Scale factor for normalizing scores.")
    parser.add_argument("--split_pattern", required=False, default="-", help="Pattern to split cell types for finding differential features.")
    
    args = parser.parse_args()
    main(args)
