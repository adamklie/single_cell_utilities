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
overwrite=$4

# Echo inputs and number of inputs
echo -e "input_h5ad_path: $input_h5ad_path"
echo -e "annotations_path: $annotations_path"
echo -e "outdir_path: $outdir_path"
echo -e "overwrite: $overwrite\n"

# Make outdir path if it doesn't exist
if [ ! -d $outdir_path ]; then
    mkdir -p $outdir_path
fi 

# Cmd for making initial AnnDataset, if annotated.h5ads exists and overwrite is False, overwrite it
if [ -d "$outdir_path/annotated.h5ads" ] && [ "$overwrite" = "False" ]; then
    cmd="echo -e '$outdir_path/annotated.h5ads already exists with overwrite set to False, skipping'"
else
    echo -e "Writing $outdir_path/annotated.h5ads"
    cmd="python $script_path \
    --path_input $input_h5ad_path \
    --path_outdir $outdir_path \
    --path_annotations $annotations_path \
    --annotations_key pseudobulk \
    --n_jobs $SLURM_CPUS_PER_TASK \
    --log_file_name pseudobulk.log"
    echo -e "Running $cmd\n"
fi
eval $cmd

# Check if .bed.gz files exist in $outdir_path/fragments and whether to overwrite
if [ -d "$outdir_path/fragments" ] && find "$outdir_path/fragments" -name "*.bed.gz" | grep -q . && [ "$overwrite" = "False" ]; then
    cmd="echo -e '$outdir_path/fragments already exists with .bed.gz files and overwrite set to False, skipping'"
else
    echo -e "Either no .bed.gz files found or overwrite set to True, proceeding with $outdir_path/fragments"
    cmd="python $script_path \
    --path_input $outdir_path/annotated.h5ads/_dataset.h5ads \
    --path_outdir $outdir_path \
    --groupby_key pseudobulk \
    --save_fragments $outdir_path/fragments \
    --n_jobs $SLURM_CPUS_PER_TASK \
    --log_file_name write_fragments.log"
    echo -e "Running $cmd\n"
fi
eval $cmd

date
