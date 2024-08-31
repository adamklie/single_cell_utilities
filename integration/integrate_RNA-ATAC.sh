#!/bin/bash

#####
# USAGE:
# sbatch integrate_RNA-ATAC.sh --SLURM_SETINGS ... <rna_h5ad_paths> <atac_h5ad_paths> <output_dir>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/cellcommander
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# Inputs
rna_h5ad_path=$1
rna_dim_reduction=$2
atac_h5ad_path=$3
atac_dim_reduction=$4
outdir_path=$5

# Echo inputs and number of inputs
echo -e "rna_h5ad_path: $rna_h5ad_path"
echo -e "rna_dim_reduction: $rna_dim_reduction"
echo -e "atac_h5ad_path: $atac_h5ad_path"
echo -e "atac_dim_reduction: $atac_dim_reduction"
echo -e "outdir_path: $outdir_path\n"

# Make the output directory if it doesn't exist
mkdir -p $outdir_path

# Run the command
cmd="cellcommander joint-integrate \
--rna_h5ad_path $rna_h5ad_path \
--atac_h5ad_path $atac_h5ad_path \
--outdir_path $outdir_path \
--method wnn \
--rna-dim-reduction $rna_dim_reduction \
--atac-dim-reduction $atac_dim_reduction \
--clust-resolution 1.0 \
--cluster-key wnn_leiden_1.0 \
--random-state 1234"
echo -e "Running:\n $cmd\n"
eval $cmd

date
