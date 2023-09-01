import scanpy as sc

def dowsample_cells():
    """Downsample adata to target_cells overall
    
    We want to be able to randomly select a number of cells from our AnnData
    regardless of covariates
    
    Pretty sure ScanPy already has this in their API
    """
    pass


def dowsample_counts(
    adata,
    counts_per_cell,
    total_counts,
    random_state,
    replace,
    layer_in=None,
    layer_out=None,
):
    """Downsample the counts table uniformly to some fraction of the original counts"""
    if layer_in is None:
        if layer_out is None:
            sc.pp.downsample_counts(
                adata,
                counts_per_cell=counts_per_cell,
                total_counts=total_counts,
                random_state=random_state,
                replace=replace,
            )
        else:
            adata_copy = sc.pp.downsample_counts(
                adata,
                counts_per_cell=counts_per_cell,
                total_counts=total_counts,
                random_state=random_state,
                replace=replace,
                copy=True,
            )
            adata.layers[layer_out] = adata_copy.X
    else:
        if layer_out is None:
            adata_copy = adata.copy()
            adata_copy.X = adata_copy.layers[layer_in]
            sc.pp.downsample_counts(
                adata_copy,
                counts_per_cell=counts_per_cell,
                total_counts=total_counts,
                random_state=random_state,
                replace=replace,
            )
            adata.layers[layer_in] = adata_copy.X
        else:
            adata_copy = adata.copy()
            adata_copy.X = adata_copy.layers[layer_in]
            sc.pp.downsample_counts(
                adata_copy,
                counts_per_cell=counts_per_cell,
                total_counts=total_counts,
                random_state=random_state,
                replace=replace,
            )
            adata.layers[layer_out] = adata_copy.X


def stratified_downsample(
    adata, 
    group_by, 
    target_cells_per_group=1000
):
    """Downsample adata to target_cells per group in group_by.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    group_by : str
        Key for observation grouping
    target_cells_per_group : int, optional
        Target number of cells per group, by default 1000

    Returns
    -------
    AnnData
        Downsampled AnnData
    """
    adatas = [adata[adata.obs[group_by] == clust] for clust in adata.obs[group_by].cat.categories]
    for dat in adatas:
        if dat.n_obs > target_cells_per_group:
            sc.pp.subsample(dat, n_obs=target_cells_per_group)
    adata_downsampled = adatas[0].concatenate(*adatas[1:])
    return adata_downsampled
    