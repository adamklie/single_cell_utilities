#!/bin/bash
#SBATCH --partition=PARTITION_NAME
#SBATCH --output=logs/%x.%A.out
#SBATCH --error=logs/%x.%A.err
#SBATCH --time=14-00:00:00

#####
# SLURM Wrapper Job Template
#
# This template takes positional arguments for flexible command-line usage.
# Useful for peak calling, pseudobulk generation, and other workflows.
#
# USAGE:
#   sbatch --job-name=my_job wrapper_job.sh /path/to/input.h5ad /path/to/output /path/to/annotations.tsv
#
# NOTES:
#   - 14-day time limit for long-running jobs
#   - Creates output directory if it doesn't exist
#   - Uses $SLURM_CPUS_PER_TASK for parallelization (request CPUs with --cpus-per-task)
#####

# Date and job info
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Activate environment
source activate your_conda_env
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# Script path
script=/path/to/your_script.py

# Inputs from positional arguments
input_path=$1
output_dir=$2
annotations_path=$3

# Echo inputs
echo -e "input_path: $input_path"
echo -e "output_dir: $output_dir"
echo -e "annotations_path: $annotations_path"

# Make output directory if it doesn't exist
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

# Run the script
cmd="python $script \
--input_path $input_path \
--outdir_path $output_dir \
--annotations_path $annotations_path \
--save_peaks $output_dir/peak_calls \
--save_fragments $output_dir/fragments \
--save_coverage $output_dir/coverage \
--n_jobs ${SLURM_CPUS_PER_TASK:-1}"

echo -e "Running:\n $cmd\n"
eval $cmd

date
