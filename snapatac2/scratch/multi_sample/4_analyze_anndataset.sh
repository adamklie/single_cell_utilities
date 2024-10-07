#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/igvf_sc-islet_10X-Multiome/out/%x.%A.out
#SBATCH --error=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/igvf_sc-islet_10X-Multiome/err/%x.%A.err
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
script=/cellar/users/aklie/data/igvf/beta_cell_networks/notebooks/igvf_sc-islet_10X-Multiome/snapatac2/6_analyze_anndataset.py
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py38

h5ads_file=/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/processed/adata_atac_merged_processed.h5ads
output_dir=/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/analysis
n_features=50000

# Run the script
CMD="python $script \
--input_h5ads $h5ads_file \
--output_dir $output_dir \
--n_features $n_features \
--make_gene_matrix"
echo -e "Running:\n $CMD\n"
echo
eval $CMD

# Date
date
