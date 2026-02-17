# Run cisTopic in R using a preprocessed Seurat object saved in RDS format

library(cisTopic)

input_rds <- "/cellar/users/aklie/projects/igvf/beta_cell_networks/data/multiome_stimulated_sc/igvf/data/multiome/DM023_palmitate/rds/01sep22_igvf_mo1.qc.rds"
output_directory <- "/cellar/users/aklie/projects/igvf/beta_cell_networks/bin/qc_preprocess/cistopic/results/multiome_stimulated_sc"

# create output directory if it doesn't exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory, showWarnings = FALSE)
}

###################
# get ATAC matrix #
###################
ad <- readRDS(input_rds)
m <- GetAssayData(object = adata, assay="mpeak", slot = "counts")
m@x <- rep(1, length(m@i))
colnames(m) <- ad$cells$cell_id
rownames(m) <- sub("-", ":", rownames(m))

#############
# run model #
#############
cisTopicObject <- createcisTopicObject(m, project.name='multiome_stimulated_sc')
cisTopicObject <- runCGSModels(cisTopicObject, topic=c(2,3,4), seed=123, nCores=3, burnin = 5, iterations = 10, addModels=FALSE)
#cisTopicObject <- runCGSModels(cisTopicObject, topic=c(2, 5, 10, 20, 30, 40, 50, 60, 80, 100), seed=987, nCores=10, burnin = 120, iterations = 200, addModels=FALSE)
saveRDS(cisTopicObject, file=file.path(output_directory, "cistopicmodel.rds"))

###################
# some evaluations
###################
pdf(file.path(output_directory, "select_models.pdf"), 10, 5)
par(mfrow=c(1,2))
cisTopicObject <- selectModel(cisTopicObject)
dev.off()

#n_topics <- c(2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 80, 100)
n_topics <- c(2, 3, 4)
pdf(file.path(output_directory, "convergence.pdf"), 5, 5)
logLikelihoodByIter(cisTopicObject, select=n_topics)
dev.off()

# save all models
dir.create(file.path(output_directory, 'models'))
for (i in 1:length(n_topics)) {
  a <- selectModel(cisTopicObject, select=n_topics[i])
  cellassign <- modelMatSelection(a, 'cell', 'Probability')
  write.csv(t(cellassign), file=file.path(output_directory, sprintf("models/topics_%d.csv", n_topics[i])))
  print(i)
}
dev.off()
