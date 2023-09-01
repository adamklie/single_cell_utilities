import os
import sys
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sys.path.append("/cellar/users/aklie/data/igvf/bin")
from utils import show_values


def filtering_plot(adata, lo_counts=500, hi_counts=20000, genes=300, mito=25, doublet_score=0.2):
    """Plot number of cells passing each filter as a barplot

    Parameters
    ----------
    adata : AnnData
        AnnData object
    lo_counts : int
        Low threshold for total counts
    hi_counts : int
        High threshold for total counts
    genes : int
        Threshold for number of genes
    mito : int
        Threshold for mitochondrial counts percentage
    doublet_score : float
        Threshold for doublet score
    
    Returns
    -------
    None
    """
    conditions = [
        (adata.obs['doublet_score'] > doublet_score),
        (adata.obs['total_counts'] < lo_counts),
        (adata.obs['total_counts'] > hi_counts),
        (adata.obs['n_genes'] < genes),
        (adata.obs['pct_counts_mt'] > mito),
        ((adata.obs['doublet_score'] <= doublet_score) & (adata.obs['total_counts'] >= lo_counts) & (adata.obs['total_counts'] <= hi_counts) & (adata.obs['n_genes'] >= genes) & adata.obs['pct_counts_mt'] <= mito)
    ]
    values = ['Doublet', 'Low_nFeature', 'High_nFeature', 'Low_ngenes', 'High_MT', 'Passing']
    adata.obs['QC'] = np.select(conditions, values)
    adata.obs['QC'] = adata.obs['QC'].astype('category')
    new_df1 = pd.DataFrame(adata.obs['QC'].value_counts()).reset_index()
    p = sns.barplot(x='index', y='QC', data=new_df1, color='sandybrown')
    show_values(p)
    

def basic_scatter_plots(
    adata, 
    mito_filter=15, 
    n_counts_filter=6000,
    save=None
):
    """Perform QC on adata and save scatter plot to out_dir/qc/qc.png
    
    Parameters
    ----------
    adata : AnnData
        AnnData object
    out_dir : str
        Output directory
    mito_filter : int
        Draw horizontal red lines indicating threshold for mitochondrial counts percentage
    n_counts_filter : int
        Draw horizontal red lines indicating threshold for total counts
    
    Returns
    -------
    None
    """
    # Another way of doing the exact same thing
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Total counts and mito filters, #draw horizontal red lines indicating thresholds.
    fig, axs = plt.subplots(ncols = 2, figsize = (8,4))
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax = axs[0], show=False)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax = axs[1], show = False)
    axs[0].hlines(y = mito_filter, xmin = 0, xmax = max(adata.obs['total_counts']), color = 'red', ls = 'dashed')
    axs[1].hlines(y = n_counts_filter, xmin = 0, xmax = max(adata.obs['total_counts']), color = 'red', ls = 'dashed')
    fig.tight_layout()
    if save:
        print(f"Saving plot to {save}...")
        plt.savefig(save)


# Deprecated
def qc_deprecated(adata, out_dir, mito_filter=15, n_counts_filter=6000):
    """Perform QC on adata and save scatter plot to out_dir/qc/qc.png
    
    Parameters
    ----------
    adata : AnnData
        AnnData object
    out_dir : str
        Output directory
    mito_filter : int
        Draw horizontal red lines indicating threshold for mitochondrial counts percentage
    n_counts_filter : int
        Draw horizontal red lines indicating threshold for total counts
    
    Returns
    -------
    None
    """
    print(f"Calculating metrics and saving plot to {out_dir}/qc.png...")

    # Another way of doing the exact same thing
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Total counts and mito filters, #draw horizontal red lines indicating thresholds.
    fig, axs = plt.subplots(ncols = 2, figsize = (8,4))
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax = axs[0], show=False)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax = axs[1], show = False)
    axs[0].hlines(y = mito_filter, xmin = 0, xmax = max(adata.obs['total_counts']), color = 'red', ls = 'dashed')
    axs[1].hlines(y = n_counts_filter, xmin = 0, xmax = max(adata.obs['total_counts']), color = 'red', ls = 'dashed')
    fig.tight_layout()
    plt.savefig(os.path.join(out_dir, "qc.png"))





