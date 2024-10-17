#!/bin/bash

#####
# USAGE:
# sbatch pseudobulk.sh --SLURM_SETINGS ... <input_h5ad_path> <annotations_path> <outdir_path>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py39
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/
script_path=/cellar/users/aklie/projects/igvf/single_cell_utilities/snapatac2/pseudobulk.py

# Inputs
input_h5ad_path=$1
annotations_path=$2
outdir_path=$3

# Echo inputs and number of inputs
echo -e "input_h5ad_path: $input_h5ad_path"
echo -e "annotations_path: $annotations_path"
echo -e "outdir_path: $outdir_path"

# Make outdir path if it doesn't exist
if [ ! -d $outdir_path ]; then
    mkdir -p $outdir_path
fi 

cmd="python $script_path \
--path_input $input_h5ad_path \
--path_outdir $outdir_path \
--path_annotations $annotations_path \
--annotations_key pseudobulk \
--save_fragments $outdir_path/fragments \
--n_jobs $SLURM_CPUS_PER_TASK"
echo -e "Running:\n $cmd\n"
eval $cmd

date
