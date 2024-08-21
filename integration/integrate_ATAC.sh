#!/bin/bash

#####
# USAGE:
# sbatch integrate_ATAC.sh --SLURM_SETINGS ... <input_tsv> <output_dir>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py39
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/
script_path=/cellar/users/aklie/projects/igvf/single_cell_utilities/snapatac2/integrate_snapatac2.py

# Inputs
input_tsv=$1
outdir_path=$2
input_h5ad_paths=($(cut -f1 $input_tsv))
sample_ids=($(cut -f2 $input_tsv))

# Echo inputs and number of inputs
echo -e "total inputs: ${#input_h5ad_paths[@]}"
echo -e "input_h5ad_paths: ${input_h5ad_paths[@]}"
echo -e "sample_ids: ${sample_ids[@]}"
echo -e "outdir_path: $outdir_path"

# Make outdir path if it doesn't exist
if [ ! -d $outdir_path ]; then
    mkdir -p $outdir_path
fi 

cmd="python $script_path \
--input_h5ad_paths ${input_h5ad_paths[@]} \
--sample_ids ${sample_ids[@]} \
--outdir_path $outdir_path \
--plot_vars log_n_fragment tsse"
echo -e "Running:\n $cmd\n"
eval $cmd

date