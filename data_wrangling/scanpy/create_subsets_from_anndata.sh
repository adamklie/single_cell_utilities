#!/bin/bash
#SBATCH --partition=carter-compute

##############################################
# USAGE:
# sbatch \
#   --job-name=subset_h5ad \
#   --output=%x.out \
#   --error=%x.err \
#   --cpus-per-task=1 \
#   --mem=128G \
#   $script \
#   $H5AD_FILE $OUT_DIR $DATASET_NAME $GENE_ARG $SUBSET_COLUMNS
##############################################

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scenicplus

# Set-up dirs
H5AD_FILE=$1
OUT_DIR=$2
DATASET_NAME=$3
GENE_ARG=$4 # 0.05
SUBSET_COLUMNS=$5 # "celltypes"
SCRIPT=/cellar/users/aklie/data/igvf/scripts/create_h5ad_subsets.py

# Print messages
echo -e "Subsetting h5ad at: $H5AD_FILE"
echo -e "Saving inputs to folders in $OUT_DIR"
echo -e "Using dataset name: $DATASET_NAME"
echo -e "Using gene argument: $GENE_ARG"
echo -e "No saving csv files in this version of the script"
echo -e "Subsetting columns: $SUBSET_COLUMNS"

# Create h5ad and loom files in Python
CMD="python $SCRIPT \
    --h5ad_file $H5AD_FILE \
    --out_dir $OUT_DIR \
    --dataset_name $DATASET_NAME \
    --gene_arg $GENE_ARG
    --subset_columns $SUBSET_COLUMNS"

echo -e "Running:\n $CMD\n"
$CMD
date
