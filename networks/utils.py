import os
import pandas as pd
import numpy as np
import pyranges
import pickle
from pyarrow import feather


# Helper function for creating weight vector for adjacency list
def create_adjacency_weights(x, sample):
    smpls = np.array(x["sample"])
    ind = np.array(np.where(smpls==sample, True, False))
    scores = np.array(x["score"])
    if len(smpls) == 1:
        #print("HERE")
        if ind:
            return scores[0]
        else:
            return 0
    else:
        if ind.sum() == 0:
            return 0
        else:
            return scores[ind].mean()
        
        
def celltype_clusters(links, celltype, slack_bp=10000):
    print("Collapsing {} peaks across samples: {}".format(celltype, links["sample"].value_counts().index.values))
    
    # Create pyranges object to aid with manipulation
    celltype_links = links[links["celltype"] == celltype]
    pr = pyranges.PyRanges(celltype_links, strands=celltype_links.strand.values)

    # Collapse peaks in same sample that are within some slack bp of each other
    clusters = \
        pr.cluster(strand=None, by=["gene"], slack=slack_bp) \
        .df.groupby(['Cluster']).agg(
            {'Chromosome':'first', 
             'Start':'min', 
             'End':'max',
             'strand': lambda x: ", ".join(sorted(list(set(x)))),
             'sample': lambda x: x.tolist(),
             'celltype': lambda x: ", ".join(sorted(list(set(x)))),
             'gene': lambda x: ", ".join(sorted(list(set(x)))),
             'score': lambda x: x.tolist(),
             'qvalue': lambda x: x.tolist()})

    print("Merged {} across {} for all samples for final total merged peak count of {}".format(len(celltype_links)-len(clusters), len(celltype_links), len(clusters)))
    return clusters


# Function to create a long format list of links between peaks and genes and return it
def celltype_links(clusters, celltype, samples, return_weight_mtx=False, save_pickle=None):
    sample_links = dict()
    genes = clusters["gene"].values
    peaks = ["-".join(peak) for peak in clusters[["Chromosome", "Start", "End"]].astype(str).values]
    if return_weight_mtx:
        weight_mtx = np.empty((len(clusters), len(samples)))
        print("Returning weight matrix with dims: [{} X {}]".format(weight_mtx.shape[0], weight_mtx.shape[1]))
    
    for i, sample in enumerate(samples):
        print("Creating list of links for {}".format(sample))
        weights = clusters.apply(create_adjacency_weights, axis=1, sample=sample).values
        if return_weight_mtx:
            weight_mtx[:, i] = weights
        sample_links[sample] = pd.DataFrame(data={"gene": genes, "peak": peaks, "weight": weights})

    # Save as pickle
    if save_pickle != None:
        with open(os.path.join(save_pickle, "multi.{}.links.long.pickle".format(celltype)), "wb") as handle:
            pickle.dump(sample_links, handle)
    if return_weight_mtx:
        return sample_links, weight_mtx
    else:
        return sample_links


# Generate gene x peak adjecency matrices from a links
def celltype_adj_mtx(links, celltype, path=None):
    sample_mtxs = dict()
    for key in links.keys():  
        print("Creating adjacency matrix for {}".format(key))
        sample_mtxs[key] = links[key].pivot(index="gene", columns="peak", values="weight").fillna(0)
        if path != None:
            full_path = os.path.join(path, "{}.{}.adj.mtx.feather".format(key, celltype))
            print("Saving {}".format(full_path))
            #sample_mtxs.to_csv(os.path.join(path, "{}.{}.adj.mtx.tsv".format(key, celltype)), sep="\t")
            feather.write_feather(sample_mtxs[key].reset_index(), full_path)
    return sample_mtxs
