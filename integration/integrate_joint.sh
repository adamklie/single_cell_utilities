#!/bin/bash

#####
# USAGE:
# sbatch integrate_joint.sh --SLURM_SETINGS ... <input_tsv> <output_dir>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/cellcommander
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# Inputs
input_tsv=$1
rna_h5ad_paths=($(cut -f1 $input_tsv))
atac_h5ad_paths=($(cut -f2 $input_tsv))
sample_ids=($(cut -f3 $input_tsv))
rna_h5ad_path=${rna_h5ad_paths[$SLURM_ARRAY_TASK_ID-1]}
atac_h5ad_path=${atac_h5ad_paths[$SLURM_ARRAY_TASK_ID-1]}
sample_id=${sample_ids[$SLURM_ARRAY_TASK_ID-1]}
outdir_path=$2/${sample_id}/joint

# Echo inputs and number of inputs
echo -e "rna_h5ad_path: $rna_h5ad_path"
echo -e "atac_h5ad_path: $atac_h5ad_path"
echo -e "sample_id: $sample_id"
echo -e "outdir_path: $outdir_path\n"

# Run the command
cmd="cellcommander joint-integrate \
--rna_h5ad_path $rna_h5ad_path \
--atac_h5ad_path $atac_h5ad_path \
--outdir_path $outdir_path \
--method wnn \
--rna-dim-reduction X_seurat_default \
--atac-dim-reduction X_spectral \
--clust-resolution 1 \
--cluster-key wnn_leiden_1 \
--random-state 1234"
echo -e "Running:\n $cmd\n"
eval $cmd

date
