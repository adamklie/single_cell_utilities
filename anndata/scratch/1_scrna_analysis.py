import os
import sys
import warnings
import scanpy as sc

# Set-up
warnings.simplefilter(action='ignore', category=FutureWarning)
_stderr = sys.stderr
null = open(os.devnull,'wb')

# Dirs and paths
work_dir = 'mouse_adrenal'

# Make a directory for to store the processed scRNA-seq data.
if not os.path.exists(os.path.join(work_dir, 'scRNA')):
    os.makedirs(os.path.join(work_dir, 'scRNA'))
    
# Read in the anndata
adata = sc.read_h5ad(os.path.join(work_dir, 'data/adata.h5ad'))
adata.var_names_make_unique()

# If you don't processed RNA data, process it here. Remember to get cell type labels
# TODO

# Save
adata.write(os.path.join(work_dir, 'scRNA/adata.h5ad'), compression='gzip')