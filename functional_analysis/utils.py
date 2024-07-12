import random
import decoupler as dc
import numpy as np
from anndata import AnnData
from typing import Union, List, Tuple


def get_pseudobulk_groups(
    adata: AnnData, 
    groupby_cols: Union[str, List[str]],
    target_max_cells_per_pb: int = None,
    obs_key: str = "pseudobulk",
    copy: bool = False,
) -> Union[AnnData, None]:
    """Add a pseudobulk column to adata.obs, which assigns cells to pseudobulks.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    groupby_cols : str or list of str
        Column(s) in adata.obs to group cells by.
        No pseudobulks will contain cells from different combinations of groupby_cols.
    target_max_target_max_cells_per_pb : int, optional
        The maximum number of cells per pseudobulk. 
        If None, then psuedobulks will be every cell with the same combination of groupby_cols.
        If not None, then pseudobulks will be randomly assigned to cells with the same combination of groupby_cols.
        The default is None.
    obs_key : str, optional
        The name of the pseudobulk column in adata.obs. The default is "pseudobulk".
    copy : bool, optional
        If True, return a copy of adata. If False, perform operation inplace. The default is False.
    
    Returns
    -------
    adata : AnnData
        Annotated data matrix including a pseudobulk column in adata.obs if copy=True, otherwise None.
    """
    # Set the seed for reproducibility
    random.seed(1234)

    # If copy=True, make a copy of adata
    if copy:
        adata = adata.copy()
    
    # If groupby_cols is a string, convert to a list
    if isinstance(groupby_cols, str):
        groupby_cols = [groupby_cols]

    # If target_max_target_max_cells_per_pb is None, then the pseudobulk column should be the same as the groupby_col
    if target_max_cells_per_pb is not None:

        # Create a new column in adata.obs for pseudobulk assignment
        adata.obs[obs_key] = np.nan
        
        # For each combination of groupby_cols, assign all cells to the pseudobulks
        for group, cell_indices in adata.obs.groupby(groupby_cols).groups.items():
            group_long = "_".join([str(x) for x in group])
            random.shuffle(cell_indices.values)
            for i, start in enumerate(range(0, len(cell_indices), target_max_cells_per_pb)):
                adata.obs.loc[cell_indices[start:start+target_max_cells_per_pb], obs_key] = f"{group_long}_{i}"
    
    # If target_max_target_max_cells_per_pb is not None, then there should be 1 pseudobulk per groupby_cols combination
    else:
        for group, cells in adata.obs.groupby(groupby_cols).groups.items():
            group_long = "_".join([str(x) for x in group])
            adata.obs.loc[cells, obs_key] = f"{group_long}"
            
    # Check that every cell has a pseudobulk assignment
    nas = adata.obs.pseudobulk.isna().sum()
    if nas == 0:
        print("All cells have a pseudobulk assignment.")
    else:
        print(f"{nas} cells do not have a pseudobulk assignment.")


# Convert all object columns to category
def convert_object_in_obs(adata):
    for col in adata.obs.columns:
        if adata.obs[col].dtype == "object":
            adata.obs[col] = adata.obs[col].astype("category")
        if adata.obs[col].dtype == "category":
            adata.obs[col] = adata.obs[col].astype("str")