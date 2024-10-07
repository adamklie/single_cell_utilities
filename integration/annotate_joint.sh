#####
# USAGE:
# bash cellcommadner_annotate.sh <input_path> <outdir_path>
#####

# Date
date

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/cellcommander
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# Inputs
input_path=$1
outdir_path=$2
marker_gene_path=$3

# Run the script
cmd="cellcommander annotate \
--input_path $input_path \
--outdir_path $outdir_path \
--modality_key rna \
--layer log1p_norm \
--method manual \
--annotation-key sample_annotation \
--marker-gene-list $marker_gene_path \
--dim-reduction X_umap_wnn \
--cluster-key wnn_leiden_1 \
--skip-dea"
echo -e "Running command:\n$cmd"
eval $cmd

date