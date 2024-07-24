import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


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