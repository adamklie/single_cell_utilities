import scanpy as sc


def select_variable_genes(
    adata,
    layer=None,
    n_genes=2000,
):
    sc.pp.highly_variable_genes(adata, layer=layer, n_top_genes=n_genes, flavor='seurat_v3')
    adata = adata[:, adata.var.highly_variable].copy()
    return adata


def select_cell_fraction_genes(
    adata,
    cell_frac=0.1,
):
    """Filter genes by cell fraction
    
    Parameters
    ----------
    adata : AnnData
        AnnData object
    cell_frac : float
        Fraction of cells to use
    
    Returns
    -------
    AnnData
        AnnData object with subsetted genes
    """
    print(f"Subsetting genes with cell_frac={cell_frac}...")
    keep_mask = sc.pp.filter_genes(adata, min_cells=adata.shape[0]*cell_frac, inplace=False)[0]
    adata = adata[:, keep_mask].copy()
    return adata


# Deprecated
def subset_genes_deprecated(adata, use_variable_genes=True, n_genes=None, cell_frac=None, name="All"):
    """Subset genes based on use_variable_genes or cell_frac
    
    Parameters
    ----------
    adata : AnnData
        AnnData object
    use_variable_genes : bool
        Whether to use variable genes, if False, use cell_frac
    n_genes : int
        Number of genes to use if use_variable_genes is True
    cell_frac : float
        Fraction of cells to use if use_variable_genes is False
    name : str
        Name of subset for logging purposes
    
    Returns
    -------
    AnnData
        AnnData object with subsetted genes
    """
    print(f"Subsetting genes for {name} with use_variable_genes={use_variable_genes}, n_genes={n_genes}, cell_frac={cell_frac}...")
    if use_variable_genes:
        assert n_genes is not None
        adata_copy = sc.pp.normalize_total(adata, target_sum=1e4, copy=True)
        sc.pp.log1p(adata_copy)
        sc.pp.highly_variable_genes(adata_copy, n_top_genes=n_genes, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata = adata[:, adata_copy.var.highly_variable]
    else:
        assert cell_frac is not None
        sc.pp.filter_genes(adata, min_cells=adata.shape[0]*cell_frac)
    return adata