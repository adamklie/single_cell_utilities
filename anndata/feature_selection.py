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
