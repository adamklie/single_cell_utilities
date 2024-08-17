#! /bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/mo_EndoC-bH1_ATAC-seq/out/%x.%A.out
#SBATCH --error=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/mo_EndoC-bH1_ATAC-seq/err/%x.%A.err
#SBATCH --mem=64G
#SBATCH -n 4
#SBATCH -t 14-00:00:00

#####
# INFO:
# Script to 
#####

#####
# USAGE: sbatch --job-name=${dataset_name}_remove_blacklist chrombpnet_remove_blacklist.sh $peaks $blacklist $peaks_no_blacklist
#####

source activate chrombpnet
peaks=$1
blacklist=$2
peaks_no_blacklist=$3
cmd="bedtools intersect \
    -v \
    -a $peaks \
    -b $blacklist \
    > $peaks_no_blacklist"
echo -e $cmd
eval $cmd
