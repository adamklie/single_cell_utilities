import os
import pickle
import sys
import warnings

from pycisTopic.cistopic_class import *
from pycisTopic.diff_features import *
from pycisTopic.topic_binarization import *

warnings.simplefilter(action="ignore", category=FutureWarning)
_stderr = sys.stderr
null = open(os.devnull, "wb")

# Dirs and paths
work_dir = "../mouse_adrenal"
tmp_dir = "/cellar/users/aklie/tmp/"

# Load the object
print("Loading cisTopic object...")
cistopic_obj = pickle.load(
    open(os.path.join(work_dir, "scATAC/cistopic_obj_50_iter_model.pkl"), "rb")
)

# Binarize topics to get regions
print("Binarizing regions to create sets of candidate enhancers...")
region_bin_topics_otsu = binarize_topics(cistopic_obj, method="otsu")
region_bin_topics_top3k = binarize_topics(cistopic_obj, method="ntop", ntop=3000)

# Save these
if not os.path.exists(os.path.join(work_dir, "scATAC/candidate_enhancers")):
    os.makedirs(os.path.join(work_dir, "scATAC/candidate_enhancers"))
pickle.dump(
    region_bin_topics_otsu,
    open(
        os.path.join(work_dir, "scATAC/candidate_enhancers/region_bin_topics_otsu.pkl"),
        "wb",
    ),
)
pickle.dump(
    region_bin_topics_top3k,
    open(
        os.path.join(
            work_dir, "scATAC/candidate_enhancers/region_bin_topics_top3k.pkl"
        ),
        "wb",
    ),
)

# Impute accessibility
print("Imputing accessibility...")
imputed_acc_obj = impute_accessibility(
    cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6
)

# Find variable regions
print("Get DARs as candidate enhancers...")
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot=False)
del normalized_imputed_acc_obj

# Get marker regions and save
markers_dict = find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable="celltype",
    var_features=variable_regions,
    split_pattern="-",
)
pickle.dump(
    markers_dict,
    open(os.path.join(work_dir, "scATAC/candidate_enhancers/markers_dict.pkl"), "wb"),
)
