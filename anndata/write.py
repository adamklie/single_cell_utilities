import os
import numpy as np
import pandas as pd
import loompy as lp
from scipy import io
import scipy.sparse as sp


def write_h5ad(adata, filename, use_raw=False):
    """Write AnnData object as h5ad

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape `n_obs` x `n_vars`.
        Rows correspond to cells and columns to genes.
    filename : str
        Filename to save to

    Returns
    -------
    None
    """
    print(f"Saving as h5ad to {filename}")
    if use_raw:
        adata = adata.raw.to_adata()  #only if adata has RAW saved and thats what you want!!
    adata.write(filename)


def write_loom(adata, filename, layer=None, use_raw=False):
    """Write AnnData object as loom that is compatible with SCENIC

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape `n_obs` x `n_vars`.
        Rows correspond to cells and columns to genes.
    filename : str
        Filename to save to
    layer : str
        Layer to save instead of `X`. If `None`, `X` is saved.

    Returns
    -------
    None
    """
    print(f"Saving as loom to {filename}")
    if use_raw:
        adata = adata.raw.to_adata()  #only if adata has RAW saved and thats what you want!!
    row_attrs = dict(zip(adata.var.reset_index().columns, adata.var.reset_index().values.T))
    col_attrs = dict(zip(adata.obs.reset_index().columns, adata.obs.reset_index().values.T))
    row_attrs["Gene"] = np.array(adata.var_names)
    col_attrs["CellID"] = np.array(adata.obs_names)
    if layer is not None:
        X = adata.layers[layer]
    else:
        X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    col_attrs["nGene"] = np.array(np.sum(X.transpose() > 0, axis=0)).flatten()
    col_attrs["nUMI"] = np.array(np.sum(X.transpose(), axis=0)).flatten()
    lp.create(filename, X.transpose(), row_attrs, col_attrs)


def write_tsv(adata, filename="expr.tsv", layer=None, use_raw=False):
    """Write AnnData matrix as tsv

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape `n_obs` x `n_vars`.
        Rows correspond to cells and columns to genes.
    filename : str
        Filename to save to
    layer : str
        Layer to save instead of `X`. If `None`, `X` is saved.

    Returns
    -------
    None
    """
    print(f"Saving tsv to {filename}")
    if use_raw:
        adata = adata.raw.to_adata()  #only if adata has RAW saved and thats what you want!!
    if layer is not None:
        X = adata.layers[layer]
    else:
        X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    expr = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names).T
    expr.index.name = "gene"
    expr.to_csv(filename, sep="\t")
    
    
def write_10x(adata, dirname, layer=None, use_raw=False):
    """Write AnnData object to 10x CellRanger v3 output

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape `n_obs` x `n_vars`.
        Rows correspond to cells and columns to genes.
    filename : str
        Filename to save to
    layer : str
        Layer to save instead of `X`. If `None`, `X` is saved.

    Returns
    -------
    None
    """
    if use_raw:
        adata = adata.raw.to_adata()  #only if adata has RAW saved and thats what you want!!
    with open(os.path.join(dirname, 'barcodes.tsv'), 'w') as f:
        for item in adata.obs_names:
            f.write(item + '\n')
    with open(os.path.join(dirname, 'features.tsv'), 'w') as f:
        for item in ['\t'.join([x, x,'Gene Expression']) for x in adata.var_names]:
            f.write(item + '\n')
    io.mmwrite(os.path.join(dirname, 'matrix.mtx'), adata.X.T)
    #command = ['gzip', os.path.join(dirname, "*")]
    #subprocess.run(command, check=True)
    adata.obs.to_csv(os.path.join(dirname, 'metadata.csv'))
    