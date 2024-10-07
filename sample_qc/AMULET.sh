#!/bin/bash

#####
# USAGE:
# sbatch qc_ATAC.sh --SLURM_SETINGS ... <input_tsv> <output_dir>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/cellcommander
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# Inputs
input_tsv=$1
input_h5_paths=($(cut -f1 $input_tsv))
input_frag_paths=($(cut -f2 $input_tsv))
sample_ids=($(cut -f3 $input_tsv))
input_h5_path=${input_h5_paths[$SLURM_ARRAY_TASK_ID-1]}
input_frag_path=${input_frag_paths[$SLURM_ARRAY_TASK_ID-1]}
sample_id=${sample_ids[$SLURM_ARRAY_TASK_ID-1]}
outdir_path=$2/${sample_id}/atac

# Echo inputs and number of inputs
echo -e "total inputs: ${#input_h5_paths[@]}"
echo -e "input_h5_path: $input_h5_path"
echo -e "input_frag_path: $input_frag_path"
echo -e "sample_id: $sample_id"
echo -e "outdir_path: $outdir_path\n"

# If output dir does not exist, create it
if [ ! -d $outdir_path ]; then
    mkdir -p $outdir_path
fi

# Run AMULET
cmd="cellcommander detect-doublets \
--input_h5ad_path $input_h5_path \
--outdir_path $outdir_path/amulet \
--fragments_file $input_frag_path \
--method amulet \
--random-state 1234"

echo -e "Running command:\n$cmd\n"
eval $cmd

date
