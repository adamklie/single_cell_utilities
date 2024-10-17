#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/bin/slurm_logs/%x.%A.out
#SBATCH --time=14-00:00:00

#####
# USAGE:
# sbatch --job-name=igvf_sc-islet_10X-Multiome_snapatac2_peak_calling snapatac2_peak_calling.sh
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py39
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/
script_path=/cellar/users/aklie/opt/igvf-ucsd/single_cell_utilities/snapatac2/pseudobulk_and_peakcalling_snapatac2.py

# Inputs
input_h5ad_path=$1
outdir_path=$2
annotations_path=$3

# Echo inputs and number of inputs
echo -e "input_h5ad_path: $input_h5ad_path"
echo -e "outdir_path: $outdir_path"

# Make outdir path if it doesn't exist
if [ ! -d $outdir_path ]; then
    mkdir -p $outdir_path
fi 

cmd="python $script_path \
--input_path $input_h5ad_path \
--outdir_path $outdir_path \
--annotations_path $annotations_path \
--save_peaks $outdir_path/peak_calls \
--save_fragments $outdir_path/fragments \
--save_coverage $outdir_path/coverage \
--save_peak_matrix $outdir_path/peak_matrices \
--n_jobs $SLURM_CPUS_PER_TASK"
echo -e "Running:\n $cmd\n"
eval $cmd

date
