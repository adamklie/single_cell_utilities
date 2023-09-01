# Good resources
# https://github.com/scverse/scanpy/issues/189
# https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pp.magic.html#scanpy.external.pp.magic

import scanpy as sc
import numpy as np
        
        
def magic_impute(
    adata,
    normalized=False,
    layer_in=None,
    layer_out=None,
):
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
        
        
    sc.pp.normalize_total(gene_matrix)
    sc.pp.log1p(gene_matrix)

    # Perform imputation with MAGIC
    logging.info("Performing MAGIC imputation")
    time_in = time.time()
    sc.external.pp.magic(gene_matrix, solver="approximate")
    
    
def knn_imputation(
    self, k: int=None, 
    metric: str="euclidean", 
    diag: float=1,
    n_pca_dims: 
    int=None, 
    maximum: bool=False,
    balanced: bool=False, 
    b_sight: int=None, 
    b_maxl: int=None,
    group_constraint: Union[str, np.ndarray]=None, 
    n_jobs: int=8
) -> None:
    """Performs k-nn smoothing of the data matrix
    https://github.com/morris-lab/CellOracle/blob/8fdb96fbdf4cb6559db7a323dd7baf6d2f164d13/celloracle/trajectory/modified_VelocytoLoom_class.py#L183

    Arguments
    ---------
    k: int
        number of neighbors. If None the default it is chosen to be `0.025 * Ncells`
    metric: str
        "euclidean" or "correlation"
    diag: int, default=1
        before smoothing this value is substituted in the diagonal of the knn contiguity matrix
        Resulting in a reduction of the smoothing effect.
        E.g. if diag=8 and k=10 value of Si = (8 * S_i + sum(S_n, with n in 5nn of i)) / (8+5)
    maximum: bool, default=False
        If True the maximum value of the smoothing and the original matrix entry is taken.
    n_pca_dims: int, default=None
        number of pca to use for the knn distance metric. If None all pcs will be used. (used only if pca_space == True)
    balanced: bool
        whether to use BalancedKNN version
    b_sight: int
        the sight parameter of BalancedKNN (used only if balanced == True)
    b_maxl: int
        the maxl parameter of BalancedKNN (used only if balanced == True)

    n_jobs: int, default 8
        number of parallel jobs in knn calculation

    Returns
    -------
    Nothing but it creates the attributes:
    knn: scipy.sparse.csr_matrix
        knn contiguity matrix
    knn_smoothing_w: scipy.sparse.lil_matrix
        the weights used for the smoothing
    Sx: np.ndarray
        smoothed spliced
    Ux: np.ndarray
        smoothed unspliced

    """
    X = _adata_to_matrix(self.adata, "normalized_count")


    N = self.adata.shape[0] # cell number

    if k is None:
        k = int(N * 0.025)
    if b_sight is None and balanced:
        b_sight = int(k * 8)
    if b_maxl is None and balanced:
        b_maxl = int(k * 4)

    space = self.pcs[:, :n_pca_dims]

    if balanced:
        bknn = BalancedKNN(k=k, sight_k=b_sight, maxl=b_maxl,
                           metric=metric, mode="distance", n_jobs=n_jobs)
        bknn.fit(space)
        self.knn = bknn.kneighbors_graph(mode="distance")
    else:

        self.knn = knn_distance_matrix(space, metric=metric, k=k,
                                       mode="distance", n_jobs=n_jobs)
    connectivity = (self.knn > 0).astype(float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # SparseEfficiencyWarning: Changing the sparsity structure of a csr_matrix is expensive. lil_matrix is more efficient.
        connectivity.setdiag(diag)
    self.knn_smoothing_w = connectivity_to_weights(connectivity)

    ###
    Xx = convolve_by_sparse_weights(X, self.knn_smoothing_w)
    self.adata.layers["imputed_count"] = Xx.transpose().copy()

    self.k_knn_imputation = k