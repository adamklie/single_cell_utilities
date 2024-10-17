#!/bin/bash

#####
# USAGE:
# sbatch call_peaks.sh --SLURM_SETINGS ... <input_h5ad_path> <groupby_key> <outdir_path>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py39
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/
script_path=/cellar/users/aklie/projects/igvf/single_cell_utilities/snapatac2/call_peaks.py

# Inputs
input_h5ad_path=$1
groupy_key=$2
outdir_path=$3

# Echo inputs and number of inputs
echo -e "input_h5ad_path: $input_h5ad_path"
echo -e "groupy_key: $groupy_key"
echo -e "outdir_path: $outdir_path"

# Make outdir path if it doesn't exist
if [ ! -d $outdir_path ]; then
    mkdir -p $outdir_path
fi 

cmd="python $script_path \
--input_path $input_h5ad_path \
--outdir_path $outdir_path \
--groupby_key $groupby_key \
--save_peaks $outdir_path/peak_calls \
--save_peak_matrix $outdir_path/peak_matrices \
--n_jobs $SLURM_CPUS_PER_TASK"
echo -e "Running:\n $cmd\n"
eval $cmd

date
