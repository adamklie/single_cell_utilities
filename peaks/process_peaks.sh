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
path_filter_script=/cellar/users/aklie/projects/igvf/single_cell_utilities/peak_analysis/filter_peaks.py
path_annotate_script=/cellar/users/aklie/projects/igvf/single_cell_utilities/peak_analysis/annotate_peaks.py

# Inputs
input_tsv=$1
outdir_path=$2
input_peak_paths=($(cut -f1 $input_tsv))
names=($(cut -f2 $input_tsv))
input_peak_path=${input_peak_paths[$SLURM_ARRAY_TASK_ID-1]}
name=${names[$SLURM_ARRAY_TASK_ID-1]}

# TODOs
path_blacklist=/cellar/users/aklie/data/ref/genomes/hg38/blacklist/blacklist.bed
n_peaks=null
min_q_value=null

# Echo inputs 
echo -e "input_peak_path: $input_peak_path"
echo -e "name: $name"
echo -e "outdir_path: $outdir_path"

# Make outdir path if it doesn't exist
if [ ! -d $outdir_path ]; then
    mkdir -p $outdir_path
fi 

# Step 1: Filter peaks
cmd="python $path_filter_script \
--path_input_peaks $input_peak_path \
--path_filtered_peaks $outdir_path/${name}.filt.bed"
echo -e "Running:\n $cmd\n"
eval $cmd

# Step 2: Annotate peaks
cmd="python $path_annotate_script \
--path_input_peaks $outdir_path/${name}.filt.bed \
--path_annotate_peaks $outdir_path/${name}.annotated.bed"
echo -e "Running:\n $cmd\n"
eval $cmd

date
