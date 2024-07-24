# For each reduction in the list save a tsv file for the cell embeddings and the feature loadings
write_embeddings <- function(adata, dirname) {
    reduction_lst <- Reductions(adata)
    for(reduction in reduction_lst) {
        cell_embeddings <- Embeddings(adata, reduction = reduction)
        write.table(
            cell_embeddings, 
            file = file.path(dirname, sprintf("%s.tsv", reduction)),
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE
        )
    }
}

# For each reduction in the list save a tsv file for the cell embeddings and the feature loadings
write_loadings <- function(adata, dirname) {
    reduction_lst <- Reductions(adata)
    for(reduction in reduction_lst) {
        feature_loadings <- Loadings(adata, reduction = reduction)
        write.table(
            feature_loadings, 
            file = file.path(dirname, sprintf("%s.tsv", reduction)),
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE
        )
    }
}