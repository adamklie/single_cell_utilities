for (cell_type in cell_types) {
  
  # Get and save subset
  ad_subset <- subset(ad, subset = cell.type.1 == cell_type)
  
  # Save cell metadata as csv
  write.csv(ad_subset@meta.data, file = paste0(out_dir, "/metadata/dm023_palmitate/dm023_palmitate_endocrine_", cell_type, "_metadata.csv"), row.names = TRUE, quote = FALSE)

  # Save bcs as single column
  write.table(rownames(ad_subset@meta.data), file = paste0(out_dir, "/barcodes/dm023_palmitate/dm023_palmitate_endocrine_", cell_type, "_barcodes.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Save subset as rds
  saveRDS(ad_subset, file = paste0(out_dir, "/rds/dm023_palmitate/dm023_palmitate_endocrine_", cell_type, ".rds"))  
}