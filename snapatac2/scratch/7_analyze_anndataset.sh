#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/bin/data_annotation/snapatac2/logs/%x.%A.out
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G
#SBATCH --time=14:00:00

#####
# USAGE:
# sbatch --job-name=igvf_sc-islet_10X-Multiome_analyze_anndataset 7_analyze_anndataset.sh
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
script=/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/bin/data_annotation/snapatac2/6_analyze_anndataset.py
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py38

h5ads_file=/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/annotation/16Aug23/snapatac2/analysis/adata_atac_merged_processed.h5ads
output_dir=/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/annotation/16Aug23/snapatac2/analysis
n_features=50000

# Run the script
CMD="python $script \
--input_h5ads $h5ads_file \
--output_dir $output_dir \
--n_features $n_features \
--make_gene_matrix >> $output_dir/${SLURM_JOB_ID}_out.txt 2>>$output_dir/${SLURM_JOB_ID}_err.txt"
echo -e "Running:\n $CMD\n"
echo
eval $CMD

# Date
date
