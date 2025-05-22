#!/bin/bash

#####
# USAGE:
# sbatch annotate_peaks.sh --SLURM_SETINGS <input_tsv>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py39
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/
path_script=/cellar/users/aklie/projects/igvf/single_cell_utilities/peaks/annotate_peaks.py

# Inputs
input_tsv=$1
outdir_path=$2
input_peak_paths=($(cut -f1 $input_tsv))
names=($(cut -f2 $input_tsv))
input_peak_path=${input_peak_paths[$SLURM_ARRAY_TASK_ID-1]}
name=${names[$SLURM_ARRAY_TASK_ID-1]}
output_annot_path=$outdir_path/${name}.annotated.bed

# Echo inputs 
echo -e "input_peak_path: $input_peak_path"
echo -e "name: $name"
echo -e "outdir_path: $outdir_path"

# Make outdir path if it doesn't exist
if [ ! -d $outdir_path ]; then
    mkdir -p $outdir_path
fi 

# Step 2: Annotate peaks
cmd="python $path_script \
--path_input_peaks $input_peak_path \
--path_annotated_peaks $output_annot_path"
echo -e "Running:\n $cmd\n"
eval $cmd

date
