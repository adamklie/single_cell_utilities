#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/path/to/output/directory/out/%x.%A.out
#SBATCH --error=/path/to/output/directory/err/%x.%A.err
#SBATCH --mem-per-cpu=32G
#SBATCH -n 1
#SBATCH -t 14-00:00:00

#####
# INFO:
# Script to run create_pseudobulks_from_signac.R on a user-defined input.
#####

#####
# USAGE:
# sbatch --job-name=your_job_name create_pseudobulks.slurm.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# ======================
# BEGIN USER CONFIGURATION
# ======================

# User must specify the path to the input RDS file
RDS_FILE="/path/to/your/input/file.rds" # Modify with your actual path

# User must specify the directory to output results from the R script
OUTPUT_DIR="/path/to/your/output/directory" # Modify with your actual path

# Parameters: Modify these as per your requirements
GROUP_BY="your_group_by_column"
ASSAY_NAME="your_assay_name"
APPEND=TRUE
VALIDATE=FALSE

# If you have a specific R environment, specify it here
# source activate R_env

# ======================
# END USER CONFIGURATION
# ======================

echo -e "Processing: $RDS_FILE\n"

# Run the R script
cmd="Rscript create_pseudobulks_from_signac.R \
--rds_file $RDS_FILE \
--output_dir $OUTPUT_DIR \
--group_by $GROUP_BY \
--assay_name $ASSAY_NAME \
--append $APPEND \
--validate $VALIDATE"
echo -e $cmd
echo
eval $cmd

date
