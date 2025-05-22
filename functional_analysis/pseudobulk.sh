#!/bin/bash

#####
# USAGE:
# sbatch pseudobulk.sh --SLURM_SETINGS <config_yaml>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py39
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/
script_path=/cellar/users/aklie/projects/igvf/single_cell_utilities/functional_analysis/pseudobulk.py

# Load YAML parameters
config_yaml=$1
parse_script_path=/cellar/users/aklie/projects/igvf/single_cell_utilities/parse_yaml.py
source <(python3 $parse_script_path $config_yaml)

#cmd
cmd="python $script_path \
--path_h5ad $io_path_h5ad \
--layer $pseudobulk_layer \
--path_gene_lengths $pseudobulk_path_gene_lengths \
--path_out $io_path_out \
--groupby_keys $pseudobulk_groupby_keys \
--cellid_key $pseudobulk_cellid_key \
--target_max_cells_per_pb $pseudobulk_max_cells_per_pb \
--mode $pseudobulk_mode \
--random_state $pseudobulk_random_state"
echo -e "Running:\n$cmd\n"
eval $cmd

# date
date
