# Good resources
# https://github.com/pachterlab/BHGP_2022/blob/main/analysis/angelidis_normalize.ipynb
# https://github.com/pachterlab/BHGP_2022/blob/main/scripts/utils.py

import scanpy as sc
import numpy as np
from sklearn.preprocessing import normalize, scale
from scipy.stats import rankdata
try:
    from scipy.stats import median_absolute_deviation
except ImportError:
    from scipy.stats import median_abs_deviation as median_absolute_deviation


def do_pf(
    mtx, 
    target_sum=None
):
    pf = mtx.sum(1).flatten()
    mtx_pf = mtx / (pf / pf.mean())[:, None].reshape(-1, 1)
    if target_sum:
        mtx_pf = normalize(mtx, norm="l1") * target_sum
    return mtx_pf


def total_normalize(
    adata,
    target_sum=1e4,
    layer_in=None,
    layer_out=None,
):
    """Counts per million (CPM) transform data matrix
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    layer : str, optional
        Layer to add the Counts per million transform to. If None, will add to X
    Returns
    -------
    adata : AnnData
        Annotated data matrix with CPM transformed data matrix in .X if layer is None
        Otherwise, CPM transformed data matrix in .layers[layer]
    """
    if layer_in is None:
        norm_mtx = sc.pp.normalize_total(
            adata=adata,
            target_sum=target_sum,
            inplace=False
        )["X"]
    else:
        norm_mtx = sc.pp.normalize_total(
            adata=adata,
            target_sum=target_sum,
            layer=layer_in,
            inplace=False
        )["X"]
    if layer_out is None:
        adata.X = norm_mtx
    else:
        adata.layers[layer_out] = norm_mtx


def cpm_normalize(
    adata,
    layer_in=None,
    layer_out=None,
):
    """Counts per million (CPM) transform data matrix
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    layer : str, optional
        Layer to add the Counts per million transform to. If None, will add to X
    Returns
    -------
    adata : AnnData
        Annotated data matrix with CPM transformed data matrix in .X if layer is None
        Otherwise, CPM transformed data matrix in .layers[layer]
    """
    total_normalize(
        adata=adata,
        target_sum=1e6,
        layer_in=layer_in,
        layer_out=layer_out,
    )


def rank_normalize(
    adata,
    layer_in=None,
    layer_out=None
):
    """Rank transform data matrix a la PISCES tutorial: 
    https://github.com/califano-lab/single-cell-pipeline/blob/master/functions/process-utils.R
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    layer : str, optional
        Layer to add the rank transform to. If None, will add to X
    Returns
    -------
    adata : AnnData
        Annotated data matrix with rank transformed data matrix in .X if layer is None
        Otherwise, rank transformed data matrix in .layers[layer]
    """
    if layer_in is None:
        mtx = adata.X.toarray()
    else:
        mtx = adata.layers[layer_in].toarray()
    rank_mat = np.apply_along_axis(rankdata, 0, mtx)
    median = np.apply_along_axis(np.median, 1, rank_mat)
    mad = np.apply_along_axis(median_absolute_deviation, 1, rank_mat)
    rank_mat = (rank_mat - median[:, np.newaxis]) / mad[:, np.newaxis]
    if layer_out is None:
        adata.X = rank_mat
    else:
        adata.layers[layer_out] = rank_mat


def log1p_normalize(
    adata,
    base=None,
    layer_in=None,
    layer_out=None,
):
    """Log1p normalize data matrix.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    base : float, optional
        Base of logarithm. If None, defaults to e if `layer_in` is None, 
        otherwise defaults to 2.
    layer_in : str, optional
        Name of layer to use as input. If None, uses `adata.X`.
    layer_out : str, optional
        Name of layer to use as output. If None, uses `adata.layers['log1p']`.
    
    Returns
    -------
    AnnData
        Annotated data matrix with log1p normalized data matrix.
    """
    if layer_in is None:
        mtx = adata.X
    else:
        mtx = adata.layers[layer_in]
    if base is None:
        base = np.e
    if layer_out is None:
        adata.X = np.log1p(mtx) / np.log(base)
    else:
        adata.layers[layer_out] = np.log1p(mtx) / np.log(base)


def scale_normalize(
    adata,
    zero_center=True,
    max_value=None,
    layer_in=None,
    layer_out=None
):
    """Scale normalize data matrix.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    zero_center : bool, optional
        Whether to zero center the data, by default True
    max_value : int, optional
        Max value to clip data, by default None
    layer_in : str, optional
        Key for input layer, by default None
    layer_out : str, optional
        Key for output layer, by default None
    
    Returns
    -------
    AnnData
        Annotated data matrix with scaled data matrix in layer_out.
    """
    if layer_in is None:
        if layer_out is None:
            sc.pp.scale(
                adata, 
                zero_center=zero_center, 
                max_value=max_value
            )
        else:
            adata_copy = sc.pp.scale(
                adata, 
                zero_center=zero_center, 
                max_value=max_value,
                copy=True
            )
            adata.layers[layer_out] = adata_copy.X
    else:
        if layer_out is None:
            sc.pp.scale(
                adata, 
                zero_center=zero_center, 
                max_value=max_value,
                layer=layer_in
            )
        else:
            adata_copy = sc.pp.scale(
                adata, 
                zero_center=zero_center, 
                max_value=max_value,
                layer=layer_in,
                copy=True
            )
            adata.layers[layer_out] = adata_copy.layers[layer_in]
    

def proportional_filtering_normalize(
    adata,
    target_sum=None,
    layer_in=None,
    layer_out=None
):
    """Proportional filtering normalization.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    target_sum : float, optional
        Target sum of data matrix, by default None
    layer_in : str, optional
        Key for input layer, by default None
    layer_out : str, optional
        Key for output layer, by default None
    
    Returns
    -------
    AnnData
        Annotated data matrix with proportional filtered data matrix in layer_out.
    """
    if layer_in is None:
        if layer_out is None:
            adata.X = do_pf(adata.X, target_sum=target_sum)
            if isinstance(adata.X, np.matrix):
                adata.X = adata.X.A
        else:
            adata.layers[layer_out] = do_pf(adata.X, target_sum=target_sum)
            if isinstance(adata.layers[layer_out], np.matrix):
                adata.layers[layer_out] = adata.layers[layer_out].A
    else:
        if layer_out is None:
            adata.layers[layer_in] = do_pf(adata.layers[layer_in], target_sum=target_sum)
            if isinstance(adata.layers[layer_in], np.matrix):
                adata.layers[layer_in] = adata.layers[layer_in].A
        else:
            adata.layers[layer_out] = do_pf(adata.layers[layer_in], target_sum=target_sum)
            if isinstance(adata.layers[layer_out], np.matrix):
                adata.layers[layer_out] = adata.layers[layer_out].A


def log_proportional_filtering_normalize(
    adata,
    pc=0.5,
    iter=1,
    target_sum=None,
    layer_in=None,
    layer_out=None
):
    """Log proportional filtering normalization.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    pc : float, optional
        Pseudocount, by default 0.5
    iter : int, optional
        Number of iterations, by default 1
    target_sum : float, optional
        Target sum of data matrix, by default None
    layer_in : str, optional
        Key for input layer, by default None
    layer_out : str, optional
        Key for output layer, by default None
    
    Returns
    -------
    AnnData
        Annotated data matrix with log proportional filtered data matrix in layer_out.
    """
    if layer_in is None:
        if layer_out is None:
            adata.X = do_log_pf(adata.X, pc=pc, iter=iter)
        else:
            adata.layers[layer_out] = do_log_pf(adata.X, pc=pc, iter=iter)
    else:
        if layer_out is None:
            adata.layers[layer_in] = do_log_pf(adata.layers[layer_in], pc=pc, iter=iter)
        else:
            adata.layers[layer_out] = do_log_pf(adata.layers[layer_in], pc=pc, iter=iter
)


def sqrt_normalize(
    adata,
    layer_in=None,
    layer_out=None
):
    """Square root normalize data matrix.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    layer_in : str, optional
        Key for input layer, by default None
    layer_out : str, optional
        Key for output layer, by default None
    
    Returns
    -------
    AnnData
        Annotated data matrix with square root normalized data matrix in layer_out.
    """
    if layer_in is None:
        if layer_out is None:
            adata.X = np.sqrt(adata.X)
        else:
            adata.layers[layer_out] = np.sqrt(adata.X)
    else:
        if layer_out is None:
            adata.layers[layer_in] = np.sqrt(adata.layers[layer_in])
        else:
            adata.layers[layer_out] = np.sqrt(adata.layers[layer_in])


def sct_normalize(
    adata,
    var_features_n=None,
    vst_flavor="v2",
    keep_genes=False,
    layer_in=None, # TODO: add layer_in
    layer_out="sct_normalized_counts",
):
    """Normalize AnnData object using SCTransform

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix, cell by gene
    var_features_n : int, optional
        Number of variable features to use, by default None
    vst_flavor : str, optional
        Flavor of SCTransform to use, by default "v2"
    layer_in : str, optional
        Key for input layer, by default None (TODO: add layer_in)
    layer_out : str, optional
        Key for output layer, by default "sct_normalized_counts"
    """
    try:
        from pysctransform import SCTransform
    except ImportError:
        raise ImportError("Please install pysctransform: https://github.com/saketkc/pySCTransform/tree/develop")
    if var_features_n is None:
        var_features_n = adata.n_vars
    residuals = SCTransform(adata, var_features_n=var_features_n, vst_flavor=vst_flavor)
    sctgenes = residuals.columns.values
    reorder_gidx = np.array([adata.var.index.get_loc(i) for i in sctgenes])
    if keep_genes:
        adata.layers[layer_out] = sp.sparse.csr_matrix(adata.X.shape, dtype=np.float32)
        adata.layers[layer_out][:, reorder_gidx] = residuals.values
        adata.var["sctgenes"] = adata.var.index.isin(sctgenes)
    else:
        adata = adata[:, reorder_gidx]
        adata.layers[layer_out] = residuals.values

# Experimental

import scipy as sp
def pf_log1p_pf_normalization(mtx):
    pf = mtx.sum(axis=1).A.ravel()
    log1p_pf = np.log1p(sp.sparse.diags(pf.mean()/pf) @ mtx)
    
    pf = log1p_pf.sum(axis=1).A.ravel()
    pf_log1p_pf = sp.sparse.diags(pf.mean()/pf) @ log1p_pf
    
    return pf_log1p_pf


# TODO:
def sctransform(mtx):
    #https://github.com/saketkc/pySCTransform
    raise NotImplementedError
    

def proportional_filtering(
    mtx, 
    target_sum=None
):
    """proportional_filtering(mxt)"""
    pf = mtx.sum(1).flatten()
    mtx_pf = mtx / (pf / pf.mean())[:, None]
    if target_sum:
        mtx_pf = normalize(mtx, norm="l1") * target_sum
    return mtx_pf


def proportional_filtering_log(mtx, pc=0.5, iter=1):
    """PC is pseudocount"""
    """proportional_filtering(mtx)"""
    print(f"iter: {iter}")
    log = np.log(mtx + pc)
    pf = do_pf(log)

    iter -= 1
    if iter == 0:
        return pf
    pf_up = do_pf(np.exp(pf) - 1)
    return do_log_pf(pf_up, iter)


# Deprecated
def rank_normalize_deprecated(
        mtx):
    """Rank transform data matrix a la PISCES tutorial: 
    https://github.com/califano-lab/single-cell-pipeline/blob/master/functions/process-utils.R
    
    Parameters
    ----------
    mtx : np.ndarray
        Data matrix
    
    Returns
    -------
    np.ndarray
        Rank transformed data matrix
    """
    rank_mat = np.apply_along_axis(rankdata, 0, dat_mat)
    median = np.apply_along_axis(np.median, 1, rank_mat)
    mad = np.apply_along_axis(median_absolute_deviation, 1, rank_mat)
    rank_mat = (rank_mat - median[:, np.newaxis]) / mad[:, np.newaxis]
    return rank_mat
