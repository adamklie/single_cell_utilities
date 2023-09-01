#! /bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/igvf/beta_cell_networks/download/out/%x.%A.out
#SBATCH --error=/cellar/users/aklie/data/igvf/beta_cell_networks/download/err/%x.%A.err
#SBATCH --mem=60G
#SBATCH -n 32
#SBATCH -t 14-00:00:00

# USAGE: sbatch --job-name=Wang2023_islet_snATAC-seq_fastq-dump 1_sra_download.sh

python 1_sra_download.py