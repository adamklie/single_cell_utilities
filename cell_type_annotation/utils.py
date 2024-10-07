import sys
import os
import re


def check_marker_genes(adata, marker_genes):
    if isinstance(marker_genes, dict):
        # Check markers against adata.var
        marker_genes_in_data = dict()
        for ct, markers in marker_genes.items():
            markers_found = list()
            for marker in markers:
                if marker in adata.var.index:
                    markers_found.append(marker)
            marker_genes_in_data[ct] = markers_found

        # Remove any keys with lists of length 0
        marker_genes_in_data = {k: v for k, v in marker_genes_in_data.items() if len(v) > 0}
        
        return marker_genes_in_data
    