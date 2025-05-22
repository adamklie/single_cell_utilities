#!/bin/bash

#####
# SLURM USAGE:
# sbatch --partition=carter-compute --job-name=factor_analysis --array=1-10%10 --output=%x.%A.%a.out --cpus-per-task=1 --mem=16G --time=01-00:00:00 run_factor_analysis.sh
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py39
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/
path_script=/cellar/users/aklie/opt/gene_program_evaluation/src/inference/program_models/factor_analysis/factor_analysis.py

# Params
n_components=(5 10 15 20 25 30 35 40 45 50)
n_components=${n_components[$SLURM_ARRAY_TASK_ID-1]}
path_data=$1
path_config_dir=$2
path_config=$path_config_dir/K${n_components}.gin
layer=$3
path_out=$4

# Echo inputs and parameters
echo -e "n_components: $n_components"
echo -e "path_script: $path_script"
echo -e "path_data: $path_data"
echo -e "path_config: $path_config"
echo -e "layer: $layer"
echo -e "path_out: $path_out\n"

# Script
cmd="python $path_script \
$path_data \
--config_path $path_config \
--prog_key factor_analysis_K${n_components} \
--data_key rna \
--layer $layer \
--output \
--path_out $path_out"
echo "Running command: $cmd"
eval $cmd

date
echo "Job completed."
