subset_seurat_object <- function() {
}
    # For each unique combination
for (i in seq_len(nrow(unique_combinations))) {
  subset_values <- unique_combinations[i, ]
  subset_name <- paste0(subset_values, collapse = "_")
  
  cat("Processing", subset_name, "subset...\n")
  
  # For each assay, get the counts matrix and save it as mtx
  for (assay in assay_lst) {
    start_time <- Sys.time()
    cat("Processing", assay, "assay for", subset_name, "subset...\n")
    
    #ad_subset <- ad
    for (column in names(subset_values)) {
      #ad_subset <- subset(ad_subset, subset = ad_subset@meta.data[[column]] == subset_values[[column]])
      print(column)
    }
    
    #m <- GetAssayData(object = ad_subset, assay = assay, slot = "counts")
    #rownames(m) <- sub("-", ":", rownames(m))
    
    #write.table(colnames(m), file = generate_file_path(output_dir, version, assay, subset_name, "obs.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    #write.table(rownames(m), file = generate_file_path(output_dir, version, assay, subset_name, "var.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    #Matrix::writeMM(m, file = generate_file_path(output_dir, version, assay, subset_name))
    print(generate_file_path(output_dir, version, assay, subset_name, "obs.tsv"))
    print(generate_file_path(output_dir, version, assay, subset_name, "var.tsv"))
    print(generate_file_path(output_dir, version, assay, subset_name))

    cat("Processed", assay, "assay for", subset_name, "subset in", Sys.time() - start_time, "seconds\n")
  }

  cat("Processed", subset_name, "subset in", Sys.time() - start_time, "seconds\n")
}


# For each cell type and assay, save the variable names as a single column tsv file
for (cell_type in cell_types) { 
    for (assay in assay_lst) {
        # Get and save subset
        ad_subset <- subset(ad, subset = cell.type.1 == cell_type)
        m <- GetAssayData(object = ad_subset, assay=assay, slot = "counts")
        rownames(m) <- sub("-", ":", rownames(m))
        write.table(rownames(m), file=paste0("/cellar/users/aklie/data/igvf/beta_cell_networks/multiome_stimulated_sc/matrix/dm023_palmitate/dm023_palmitate_endocrine_", cell_type, "_", assay, ".var.tsv"), row.names=F, col.names=F, quote=F)
    }
}

# For each cell type and assay, get the counts matrix and save it as mtx
for (cell_type in cell_types) { 
    for (assay in assay_lst) {
        # Get and save subset
        ad_subset <- subset(ad, subset = cell.type.1 == cell_type)
        m <- GetAssayData(object = ad_subset, assay=assay, slot = "counts")
        rownames(m) <- sub("-", ":", rownames(m))
        write.table(colnames(m), file=paste0("/cellar/users/aklie/data/igvf/beta_cell_networks/multiome_stimulated_sc/matrix/dm023_palmitate/dm023_palmitate_endocrine_", cell_type, "_", assay, ".obs.tsv"), row.names=F, col.names=F, quote=F)
        write.table(rownames(m), file=paste0("/cellar/users/aklie/data/igvf/beta_cell_networks/multiome_stimulated_sc/matrix/dm023_palmitate/dm023_palmitate_endocrine_", cell_type, "_", assay, ".var.tsv"), row.names=F, col.names=F, quote=F)
        Matrix::writeMM(m, file=paste0("/cellar/users/aklie/data/igvf/beta_cell_networks/multiome_stimulated_sc/matrix/dm023_palmitate/dm023_palmitate_endocrine_", cell_type, "_", assay, ".count.mtx"))
    }
}

# For each cell type and assay, get the counts matrix and save it as mtx
for (cell_type in cell_types) { 
    for (assay in assay_lst) {
        # Get and save subset
        ad_subset <- subset(ad, subset = cell.type.1 == cell_type)
        m <- GetAssayData(object = ad_subset, assay=assay, slot = "counts")
        rownames(m) <- sub("-", ":", rownames(m))
        write.table(colnames(m), file=paste0("/cellar/users/aklie/data/igvf/beta_cell_networks/multiome_stimulated_sc/matrix/dm023_palmitate/dm023_palmitate_endocrine_", cell_type, "_", assay, ".obs.tsv"), row.names=F, col.names=F, quote=F)
        write.table(rownames(m), file=paste0("/cellar/users/aklie/data/igvf/beta_cell_networks/multiome_stimulated_sc/matrix/dm023_palmitate/dm023_palmitate_endocrine_", cell_type, "_", assay, ".var.tsv"), row.names=F, col.names=F, quote=F)
        Matrix::writeMM(m, file=paste0("/cellar/users/aklie/data/igvf/beta_cell_networks/multiome_stimulated_sc/matrix/dm023_palmitate/dm023_palmitate_endocrine_", cell_type, "_", assay, ".count.mtx"))
    }
}
for (subset in subset_combinations) {
  subset_name <- paste0(subset, collapse = "_")
  cat("Processing", subset_name, "subset...\n")
  
  # For each assay, get the counts matrix and save it as mtx
  for (assay in assay_lst) {
    start_time <- Sys.time()
    cat("Processing", assay, "assay for", subset_name, "subset...\n")
    
    ad_subset <- ad
    for (column in subset) {
      ad_subset <- subset(ad_subset, subset = ad_subset@meta.data[[column]] == 1)
    }
    m <- GetAssayData(object = ad_subset, assay = assay, slot = "counts")
    rownames(m) <- sub("-", ":", rownames(m))
    
    write.table(colnames(m), file = generate_file_path(output_dir, version, assay, subset_name, "obs.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(rownames(m), file = generate_file_path(output_dir, version, assay, subset_name, "var.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    Matrix::writeMM(m, file