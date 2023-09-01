#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/igvf_sc-islet_10X-Multiome/out/%x.%A_%a.out
#SBATCH --error=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/igvf_sc-islet_10X-Multiome/err/%x.%A_%a.err
#SBATCH --mem-per-cpu=32G
#SBATCH -n 1
#SBATCH -t 14-00:00:00
#SBATCH --array=1-6%6

#####
# INFO:
# Script to run peak calling using Signac's CallPeaks function on pseudobulk outputs.
#####

#####
# USAGE:
# sbatch --job-name=igvf_sc-ilset_10X-Multiome_peakcalling call_peaks.slurm.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring environment (if needed)
# source activate R_env

# ======================
# BEGIN USER CONFIGURATION
# ======================

# Get the file name based on the SLURM_ARRAY_TASK_ID
FILE=$(ls /cellar/users/aklie/data/igvf/beta_cell_networks/aligned/igvf_sc-islet_10X-Multiome/10Aug23/pycistopic/*bed.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")
BASENAME=$(basename $FILE .bed.gz)

# Configuring output
bed=$FILE
out_dir=/cellar/users/aklie/data/igvf/beta_cell_networks/peaks/igvf_sc-islet_10X-Multiome/10Aug23/pycistopic

# Make the output directory if it doesn't exist
if [ ! -d $out_dir ]
then
    mkdir -p $out_dir
fi

# ======================
# END USER CONFIGURATION
# ======================

echo -e "Input file: $bed\n"
echo -e "Output directory: $out_dir\n"

# Check if the output file already exists
if [ -f $out_dir/$BASENAME*.narrowPeak ]
then
    echo -e "Output narrowPeak file for $BASENAME already exists. Skipping...\n"
    exit 0
fi

# Run the R script for peak calling
cmd="Rscript call_peaks_signac.R \
--bed_file=$bed \
--output_dir=$out_dir"
echo -e $cmd
echo
eval $cmd

date
