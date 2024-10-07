import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc


def grouped_proportion_barplot(
    adata,
    groupby_key,
    proportion_key,
    groupby_order=None,
    colors=None,
    save=None
):
    # Stacked barplot of groupby_key proportions in grouped groups
    grouped = adata.obs[groupby_key].unique().tolist()
    proportions_lst = []
    for group in grouped:
        proportions = adata.obs[adata.obs[groupby_key] == group][proportion_key].value_counts(normalize=True)
        proportions.name = group
        proportions_lst.append(proportions)
    proportions_df = pd.concat(proportions_lst, axis=1).T
    if groupby_order is not None:
        proportions_df = proportions_df.loc[groupby_order]
    if colors is not None:
        proportions_df = proportions_df[colors.keys()]
        colors = colors.values()

    # Plot horizontal barplot, legend outside to right
    _, ax = plt.subplots(figsize=(7.5, 10))
    proportions_df.plot.barh(stacked=True, ax=ax, color=colors)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xlabel("Proportion")
    plt.ylabel(f"{groupby_key}")
    plt.title(f"{proportion_key} proportions in {groupby_key} groups")
    if save:
        plt.savefig(save, bbox_inches="tight")


def prettier_grouped_proportion_barplot(
    adata,
    groupby_key,
    proportion_key,
    groupby_order=None,
    colors=None,
    figsize=(2, 7),
    fontscale=1.5,
    sns_style="white",
    save=None,
):
    # Set font scale and style
    sns.set(font_scale=1.5)
    sns.set_style("white")

    # Stacked barplot of groupby_key proportions in grouped groups
    grouped = adata.obs[groupby_key].unique().tolist()
    proportions_lst = []
    for group in grouped:
        proportions = adata.obs[adata.obs[groupby_key] == group][proportion_key].value_counts(normalize=True)
        proportions.name = group
        proportions_lst.append(proportions)
    proportions_df = pd.concat(proportions_lst, axis=1).T
    if groupby_order is not None:
        proportions_df = proportions_df.loc[groupby_order]
    if colors is not None:
        proportions_df = proportions_df[colors.keys()]
        colors = colors.values()

    # Plot horizontal barplot, legend outside to right
    _, ax = plt.subplots(figsize=figsize)
    proportions_df.plot.barh(stacked=True, ax=ax, color=colors)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xlabel("Proportion")
    plt.ylabel(f"{groupby_key}")
    plt.title(f"{proportion_key} proportions in {groupby_key} groups")
    if save:
        plt.savefig(save, bbox_inches="tight")
    plt.show()


def plot_cluster_counts(
    adata,
    clustering_key,
    cluster_order=None,
    save=None
):
    # Make fontsizes larger
    sns.set(font_scale=1.5)
    sns.set_style("white")

    # Plot horizontal barplot, legend outside to right
    _, ax = plt.subplots(figsize=(2, 7))
    if cluster_order is None:
        s = adata.obs[clustering_key].value_counts(ascending=True)
    else:
        s = adata.obs[clustering_key].value_counts(ascending=True).reindex(cluster_order)
    s.plot.barh(stacked=True, ax=ax, color="#86cefb")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xlabel("# nuclei")
    plt.ylabel("Cluster")
    plt.title("")
    if save:
        plt.savefig(save, bbox_inches="tight")
    plt.show()


def combined_cluster_proportions(
    adata,
    clustering_key,
    proportion_keys,
    cluster_order=None,
    colors_dicts=None,
    figsize=(20, 7),
    fontscale=1.5,
    sns_style="white",
    title=None,
    save=None
):
    sns.set(font_scale=fontscale)
    sns.set_style(sns_style)

    n_plots = len(proportion_keys) + 1
    fig, axes = plt.subplots(1, n_plots, figsize=figsize, sharey=True)

    # Plot the cluster counts
    sns.set(font_scale=fontscale)
    sns.set_style(sns_style)
    if cluster_order is None:
        s = adata.obs[clustering_key].value_counts(ascending=True)
    else:
        s = adata.obs[clustering_key].value_counts(ascending=True).reindex(cluster_order)
    s.plot.barh(stacked=True, ax=axes[0], color="#86cefb")
    axes[0].set_xlabel("# nuclei")
    axes[0].set_ylabel("Cluster")
    axes[0].set_title("")

    # Plot each proportion barplot
    legend_handles = []
    legend_labels = []
    legend_titles = []

    for i, (proportion_key, colors) in enumerate(zip(proportion_keys, colors_dicts)):
        grouped = adata.obs[clustering_key].unique().tolist()
        proportions_lst = []
        for group in grouped:
            proportions = adata.obs[adata.obs[clustering_key] == group][proportion_key].value_counts(normalize=True)
            proportions.name = group
            proportions_lst.append(proportions)
        proportions_df = pd.concat(proportions_lst, axis=1).T
        if cluster_order is not None:
            proportions_df = proportions_df.loc[cluster_order]
        if colors is not None:
            proportions_df = proportions_df[colors.keys()]
            colors = colors.values()

        ax = axes[i + 1]
        proportions_df.plot.barh(stacked=True, ax=ax, color=colors)
        ax.set_xlabel("Proportion")
        ax.set_title(f"")

        # Collect legend handles and labels
        handles, labels = ax.get_legend_handles_labels()
        legend_handles.append(handles)
        legend_labels.append(labels)
        legend_titles.append(proportion_key)
        ax.get_legend().remove()

    # Set the title for the entire figure
    if title:
        plt.suptitle(title, y=1.02, fontsize=16)

    # Adjust layout to make space for the legends
    plt.tight_layout(rect=[0, 0, 0.8, 1])

    # Create a new axis for each legend
    for i, (handles, labels, title) in enumerate(zip(legend_handles, legend_labels, legend_titles)):
        legend_ax = fig.add_axes([0.82, 0.8 - 0.2 * i, 0.1, 0.2])  # Adjust the position of each legend
        legend_ax.axis('off')
        legend_ax.legend(handles, labels, title=title, loc='upper left')

    if save:
        plt.savefig(save, bbox_inches="tight")
    plt.show()


def pretty_umap(
    adata, 
    color_key, 
    figsize=(6,5), 
    s=25, 
    legend_loc="on data", 
    legend_fontsize=12, 
    legend_fontoutline=2,
    save=None
):
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    sc.pl.umap(
        adata, 
        color=[color_key],
        s=s, 
        show=False, 
        frameon=False, 
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        legend_fontoutline=legend_fontoutline,
        ax=ax
    )
    ax.axis("on")
    ax.tick_params(
        top="off",
        bottom="on",
        left="on",
        right="off",
        labelleft="on",
        labelbottom="off",
    )
    _ = ax.set_xlabel("umap1", fontsize=16)
    _ = ax.set_ylabel("umap2", fontsize=16)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches="tight")
        plt.close()
        

def countplot(
    data: pd.DataFrame,
    y_column: str,
    ax: plt.Axes = None,
    figsize: tuple = (4, 6),
    order: list = None,
    palette: dict = None,
    save: str = None,
    xlim: int = None,
    **kwargs,
):
    """
    Horizontal barplot showing number of cells (x) in each sample (y) colored by updated annotation
    """
    # if 
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    # if order is not provided, order by count
    if order is None:
        order = data[y_column].value_counts().index
        
    # Plot
    sns.countplot(
        y=y_column,
        data=data,
        ax=ax,
        palette=palette,
        order=order,
        **kwargs,
    )

    # Add count numbers to right of each bar
    for p in ax.patches:
        width = p.get_width()
        ax.text(width, p.get_y() + p.get_height() / 2, f"{int(width)}", ha="left", va="center")

    # Set xlim
    if xlim:
        ax.set_xlim(0, xlim)

    # Save
    if save:
        plt.savefig(save, bbox_inches="tight")
        