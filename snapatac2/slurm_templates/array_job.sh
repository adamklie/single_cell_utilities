#!/bin/bash
#SBATCH --partition=PARTITION_NAME
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64G
#SBATCH --time=02:00:00
#SBATCH --array=1-N%N

#####
# SLURM Array Job Template
#
# This template demonstrates how to run parallel jobs across multiple samples.
# Each array task processes one sample independently.
#
# USAGE:
#   1. Update PARTITION_NAME to your cluster partition
#   2. Update N in --array to match number of samples
#   3. Fill in sample arrays below
#   4. Submit: sbatch --job-name=my_job_name array_job.sh
#
# NOTES:
#   - %A = master job ID, %a = array task ID
#   - Array index is 1-based, bash arrays are 0-based (hence $SLURM_ARRAY_TASK_ID-1)
#   - mem-per-cpu=64G is suitable for fragment file import; reduce to 32G for lighter tasks
#####

# Date and job info
date
echo -e "Job ID: $SLURM_JOB_ID"
echo -e "Array Task ID: $SLURM_ARRAY_TASK_ID\n"

# Activate environment
source activate your_conda_env

# Define input files (one per array task)
input_files=(
    '/path/to/sample1/input.file'
    '/path/to/sample2/input.file'
    # ... add more samples
)

# Define sample IDs (parallel arrays)
sample_ids=(
    'sample1'
    'sample2'
    # ... add more samples
)

# Get current task's input and sample ID
input_file=${input_files[$SLURM_ARRAY_TASK_ID-1]}
sample_id=${sample_ids[$SLURM_ARRAY_TASK_ID-1]}

# Output directory
output_dir=/path/to/output

# Your script path
script=/path/to/your_script.py

# Run the script
CMD="python $script \
--input $input_file \
--output $output_dir/${sample_id}_output.h5ad"

echo -e "Running:\n $CMD\n"
eval $CMD

# Date
date
