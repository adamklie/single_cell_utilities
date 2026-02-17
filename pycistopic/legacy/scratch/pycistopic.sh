#!/bin/bash

# Run cisTopic models
script_path=/cellar/users/aklie/projects/igvf/single_cell_utilities/pycistopic/runModels_lda_cgs.sh

# Slurm settings
job_name=2025_12_05_sc-islet-differentiation_10X-Multiome_pycistopic_runModels
partition=carter-compute
cpus_per_task=12
mem=600G
time="14-00:00:00"
output_path="/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/scratch/2025_12_05/pycistopic/%x.%A.out"

# Inputs
inputcisTopic_obj=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/scratch/2025_11_30/cistopic_obj.pkl
save_path=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/scratch/2025_12_05/pycistopic
output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/scratch/2025_12_05/pycistopic/models.pkl

# Cmd
cmd="sbatch \
--job-name=$job_name \
--partition=$partition \
--cpus-per-task=$cpus_per_task \
--mem=$mem \
--time=$time \
--output=$output_path \
$script_path $inputcisTopic_obj $save_path $output"
echo -e "Running command:\n$cmd\n"
eval $cmd
