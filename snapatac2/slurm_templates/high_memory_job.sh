#!/bin/bash
#SBATCH --partition=PARTITION_NAME
#SBATCH --output=logs/%x.%A.out
#SBATCH --error=logs/%x.%A.err
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G
#SBATCH --time=14:00:00

#####
# SLURM High-Memory Job Template
#
# This template is for single jobs that require high memory, such as:
#   - Merging multiple AnnData objects into AnnDataSet
#   - Running analysis on merged datasets
#   - Integration/embedding of large datasets
#
# USAGE:
#   1. Update PARTITION_NAME to your cluster partition
#   2. Adjust mem-per-cpu as needed (128G is typical for merged datasets)
#   3. Submit: sbatch --job-name=my_job_name high_memory_job.sh
#
# NOTES:
#   - For merged dataset analysis, 128G is often needed
#   - Time limit of 14 hours is suitable for spectral embedding + clustering
#####

# Date and job info
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Activate environment
source activate your_conda_env

# Input/output paths
input_h5ads=/path/to/merged_dataset.h5ads
output_dir=/path/to/output

# Script path
script=/path/to/analyze_script.py

# Parameters
n_features=50000

# Run the script
CMD="python $script \
--input_h5ads $input_h5ads \
--output_dir $output_dir \
--n_features $n_features"

echo -e "Running:\n $CMD\n"
eval $CMD

# Date
date
