#####
# USAGE:
# bash cellcommander_annotate.sh <input_path> <outdir_path> <marker_gene_path>
#####

# Date
date

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/cellcommander
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# Inputs
input_path=$1
marker_gene_path=$2
outdir_path=$3

# Run the script
cmd="cellcommander annotate \
--input_path $input_path \
--outdir_path $outdir_path \
--layer log1p_norm \
--method manual \
--annotation-key manual_annotation \
--marker-gene-list $marker_gene_path \
--dim-reduction X_seurat_default_umap \
--cluster-key leiden_1 \
--skip-dea"
echo -e "Running command:\n$cmd"
eval $cmd

date
