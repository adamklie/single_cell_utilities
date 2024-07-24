import numpy as np
import scanpy as sc
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
import pandas as pd
from typing import List, Union, Optional, Any


def construct_metacells(
    adata: sc.AnnData,
    name: str = 'agg',
    k: int = 25,
    reduction: str = 'umap',
    layer: str = 'counts',
    return_metacell: bool = False,
    mode: str = 'average',
    max_shared: int = 15,
    target_metacells: int = 1000,
    max_iter: int = 5000,
    verbose: bool = False,
    cells_use: Optional[np.array] = None
) -> sc.AnnData:
    
    if reduction not in adata.obsm.keys():
        raise ValueError(f"Invalid reduction ({reduction}). Reductions in adata: {list(adata.obsm.keys())}")

    if layer not in adata.layers.keys():
        raise ValueError(f"Invalid layer ({layer}). Layers in adata: {list(adata.layers.keys())}")

    if cells_use is not None:
        adata = adata[cells_use]

    reduced_coordinates = adata.obsm[reduction]
    tree = cKDTree(reduced_coordinates)
    _, nn_map = tree.query(reduced_coordinates, k=k+1)
    nn_map = nn_map[:, 1:]  # remove self index

    choices = list(range(len(adata)))
    np.random.shuffle(choices)
    chosen = [choices.pop()]

    it = 0
    k2 = k * 2

    while len(choices) > 0 and len(chosen) < target_metacells and it < max_iter:
        it += 1
        new_chosen = chosen + [choices.pop()]
        cell_sample = nn_map[new_chosen]

        shared = [k2 - len(np.unique(cell_sample[-1, :].tolist() + cell_sample[i, :].tolist())) for i in range(len(cell_sample) - 1)]
        if max(shared) <= max_shared:
            chosen = new_chosen

    cell_sample = nn_map[chosen]
    cells_merged = [",".join(adata.obs_names[cell_sample[i]] for i in row) for row in cell_sample]

    exprs = adata.layers[layer].todense()
    mask = np.zeros((adata.shape[1], len(chosen)), dtype=bool)
    for idx, ch in enumerate(chosen):
        mask[ch, idx] = True
    mask = csr_matrix(mask)

    new_exprs = exprs @ mask
    if mode == 'average':
        new_exprs /= k

    new_adata = sc.AnnData(new_exprs.T)
    new_adata.obs_names = [f"{name}_{i}" for i in range(1, len(chosen)+1)]
    new_adata.var_names = adata.var_names
    new_adata.obs['cells_merged'] = cells_merged

    return new_adata if return_metacell else adata


def metacells_by_groups(
    adata: sc.AnnData, 
    group_by: List[str] = ['seurat_clusters'],
    k: int = 25, 
    reduction: str = 'pca', 
    assay: Optional[str] = None,
    cells_use: Optional[List[str]] = None,
    slot: str = 'counts', 
    mode: str = 'average',
    min_cells: int = 100,
    max_shared: int = 15,
    target_metacells: int = 1000,
    max_iter: int = 5000, 
    verbose: bool = False,
    wgcna_name: Optional[str] = None
) -> sc.AnnData:

    # Placeholder: Validate inputs
    def validate_inputs():
        pass

    # Placeholder: Get default assay if none provided
    def get_default_assay(data):
        # Implement default assay logic
        return "default_assay_name"
    
    # Placeholder: Set WGCNA params
    def set_wgcna_params(data, params):
        # Implement logic for setting WGCNA parameters
        pass
    
    # Placeholder: Set metacell object
    def set_metacell_object(data, metacell, wgcna_name):
        # Implement logic to save the metacell data in a structured manner
        pass

    validate_inputs()

    # If assay is None, get the default assay
    if assay is None:
        assay = get_default_assay(adata)

    # Grouping
    if len(group_by) > 1:
        adata.obs['metacell_grouping'] = adata.obs[group_by].astype(str).agg('#'.join, axis=1)
    else:
        adata.obs['metacell_grouping'] = adata.obs[group_by[0]]

    groupings = adata.obs['metacell_grouping'].unique().tolist()
    group_counts = adata.obs['metacell_grouping'].value_counts()

    # Remove groups that don't meet min_cells
    groupings = [group for group in groupings if group_counts[group] >= min_cells]
    if not groupings:
        raise ValueError("No groups met the min_cells requirement.")

    # Subset Data
    if cells_use:
        adata_full = adata
        adata = adata[cells_use, :]

    # Metacell Computation
    metacell_list = {}
    for group in groupings:
        subset_adata = adata[adata.obs['metacell_grouping'] == group].copy()
        metacell = construct_metacells(subset_adata, k=k)
        metacell_list[group] = metacell
    
    # Merge metacells
    merged_adata = metacell_list[groupings[0]].concatenate(*[metacell_list[group] for group in groupings[1:]], join='outer')

    # Set idents for metacell object
    # Placeholder: This is Seurat-specific, need equivalent action in Scanpy
    def set_idents(data, idents):
        pass
    set_idents(merged_adata, "idents_name_here")

    # Revert to full AnnData object if we subsetted earlier
    if cells_use:
        adata = adata_full

    # Placeholder: Add metacell AnnData object to the main AnnData object
    set_metacell_object(adata, merged_adata, wgcna_name)

    # Placeholder: Add other info such as WGCNA parameters
    wgcna_params = {
        'metacell_k': k,
        'metacell_reduction': reduction,
        'metacell_slot': slot,
        'metacell_assay': assay,
        'metacell_stats': "stats_data_here"  # Placeholder
    }
    set_wgcna_params(adata, wgcna_params)

    return adata


def get_metacell_object(adata: sc.AnnData, wgcna_name: Optional[str] = None) -> sc.AnnData:
    """Placeholder to retrieve metacell object."""
    # Implement retrieval logic
    return adata

def set_metacell_object(adata: sc.AnnData, metacell_obj: sc.AnnData, wgcna_name: Optional[str] = None) -> None:
    """Placeholder to set metacell object."""
    # Implement set logic
    pass

def normalize_metacells(adata: sc.AnnData, wgcna_name: Optional[str] = None, **kwargs) -> None:
    metacell_obj = get_metacell_object(adata, wgcna_name)
    sc.pp.normalize_total(metacell_obj, **kwargs)
    set_metacell_object(adata, metacell_obj, wgcna_name)

def scale_metacells(adata: sc.AnnData, wgcna_name: Optional[str] = None, **kwargs) -> None:
    if 'features' not in locals():
        features = adata.var_names
    metacell_obj = get_metacell_object(adata, wgcna_name)
    sc.pp.scale(metacell_obj, **kwargs)
    set_metacell_object(adata, metacell_obj, wgcna_name)

def run_pca_metacells(adata: sc.AnnData, wgcna_name: Optional[str] = None, **kwargs) -> None:
    metacell_obj = get_metacell_object(adata, wgcna_name)
    sc.tl.pca(metacell_obj, **kwargs)
    set_metacell_object(adata, metacell_obj, wgcna_name)

# Assuming you've installed Scanpy's Harmony (scanpy.external.pp.harmony_integrate)
def run_harmony_metacells(adata: sc.AnnData, **kwargs) -> None:
    from scanpy.external.pp import harmony_integrate
    wgcna_name = adata.uns['active_wgcna']  # Placeholder for active WGCNA name
    harmony_integrate(adata, basis=wgcna_name, **kwargs)

def run_umap_metacells(adata: sc.AnnData, **kwargs) -> None:
    wgcna_name = adata.uns['active_wgcna']  # Placeholder for active WGCNA name
    sc.tl.umap(adata, **kwargs)

def dimplot_metacells(adata: sc.AnnData, **kwargs) -> None:
    wgcna_name = adata.uns['active_wgcna']  # Placeholder for active WGCNA name
    sc.pl.embedding(adata, basis=wgcna_name, **kwargs)
