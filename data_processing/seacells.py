# following this tutorial: https://github.com/dpeerlab/SEACells/blob/main/notebooks/SEACell_computation.ipynb
import numpy as np
import pandas as pd
import scanpy as sc
import SEACells

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

import scipy
from scipy import io

# load the dataset downloaded above into Python
adata = sc.read_h5ad('cd34_multiome_rna.h5ad')

# retain the unprocessed UMI counts matrix in the .raw slot
raw_ad = sc.AnnData(adata.X)
raw_ad.obs_names, raw_ad.var_names = adata.obs_names, adata.var_names
adata.raw = raw_ad

# process data with SCANPY
# note that we don't scale the data matrix before PCA. this is how
# they do it in the SEACells tutorial so we do it that way here.
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1500)
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)

##################################################################################
# Running SEACells
##################################################################################

# they recommend one metacell for every 75 real cells
n_SEACells = int(np.floor(adata.obs.shape[0] / 75))

build_kernel_on = 'X_pca' # key in ad.obsm to use for computing metacells
                          # This would be replaced by 'X_svd' for ATAC data

## Additional parameters
n_waypoint_eigs = 10 # Number of eigenvalues to consider when initializing metacells
waypoint_proportion = 0.9 # Proportion of metacells to initialize using waypoint analysis,
                        # the remainder of cells are selected by greedy selection

# set up the model
model = SEACells.core.SEACells(adata,
                  build_kernel_on=build_kernel_on,
                  n_SEACells=n_SEACells,
                  n_waypoint_eigs=n_waypoint_eigs,
                  waypt_proportion=waypoint_proportion,
                  convergence_epsilon = 1e-5)

# Initialize archetypes
model.initialize_archetypes()

# fit the model
model.fit(n_iter=20)

# create aggregated metacell expression dataset
SEACell_ad = SEACells.core.summarize_by_SEACell(adata, SEACells_label='SEACell', summarize_layer='raw')