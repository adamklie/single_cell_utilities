#! /bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/mo_EndoC-bH1_ATAC-seq/out/%x.%A.out
#SBATCH --error=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/mo_EndoC-bH1_ATAC-seq/err/%x.%A.err
#SBATCH --mem=64G
#SBATCH -n 4
#SBATCH -t 14-00:00:00

#####
# INFO:
# Script to run 
#####

#####
# USAGE: sbatch --job-name=mo_EndoC-bH1_ATAC-seq_prep_nonpeaks chrombpnet_prep_nonpeaks.sh $genome $peaks $chromsizes $fold $blacklist $out_dir
#####
date
echo -e "Job ID: $SLURM_JOB_ID\n"

source activate chrombpnet
genome=$1
peaks=$2
chromsizes=$3
fold=$4
blacklist=$5
out_dir=$6
cmd="chrombpnet prep nonpeaks \
    -g $genome \
    -p $peaks \
    -c $chromsizes \
    -fl $fold \
    -br $blacklist \
    -o $out_dir"
echo -e $cmd
eval $cmd

date