import random
import numpy as np
import pandas as pd
from anndata import AnnData
import scanpy as sc


def create_pseudobulks(
    adata: AnnData,
    groupby: str,
    layer: str = None,
    use_raw: bool = False,
    n_pseudoreplicates: int = 1,
    groupby_metadata: list = [],
):
    if use_raw:
        adata = adata.raw.to_adata()

    # Check if groupby in obs
    if groupby not in adata.obs.columns:
        raise ValueError(f"{groupby} not found in adata.obs.columns")
    
    pseudobulk_adatas = []
    for group in adata.obs[groupby].unique():
        group_adata = adata[adata.obs[groupby] == group].copy()
        group_adata.X = group_adata.layers[layer] if layer else group_adata.X
        
        indices = list(group_adata.obs_names)
        random.shuffle(indices)
        indices = np.array_split(np.array(indices), n_pseudoreplicates)
        for i, pseudo_rep in enumerate(indices):
            group_rep_adata = group_adata[group_adata.obs.index.isin(pseudo_rep)].copy()
            rep_adata = sc.AnnData(
                X = group_rep_adata.X.sum(axis = 0),
                var = group_adata.var[[]]
            )

            rep_adata.obs_names = [group + '_rep' + str(i)]
            for col in groupby_metadata + [groupby]:
                # check to see if the column is in the metadata
                if col in group_adata.obs.columns:
                    rep_adata.obs[col] = group_adata.obs[col].iloc[0]
            rep_adata.obs['replicate'] = i

            pseudobulk_adatas.append(rep_adata)
    return sc.concat(pseudobulk_adatas)
