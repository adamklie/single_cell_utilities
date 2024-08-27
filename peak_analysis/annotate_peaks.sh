#!/bin/bash

#####
# USAGE:
# sbatch annotate_peaks.sh --SLURM_SETINGS <peaks_dir> <outdir_path>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py39
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/
path_script=/cellar/users/aklie/projects/igvf/single_cell_utilities/peak_analysis/annotate_peaks.py

# Inputs
peaks_dir=$1
narrowPeak_files=$(ls $peaks_dir/*.narrowPeak)
outdir_path=$2

# Echo inputs 
echo -e "peaks_dir: $peaks_dir"
echo -e "narrowPeak_files: $narrowPeak_files"
echo -e "outdir_path: $outdir_path"

# Make outdir path if it doesn't exist
if [ ! -d $outdir_path ]; then
    mkdir -p $outdir_path
fi 

cmd="python $path_script \
--path_peaks $narrowPeak_files \
--outdir_path $outdir_path"
echo -e "Running:\n $cmd\n"
eval $cmd

date