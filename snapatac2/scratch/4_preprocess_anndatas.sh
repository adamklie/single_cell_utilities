#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/igvf_sc-islet_10X-Multiome/out/%x.%A_%a.out
#SBATCH --error=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/igvf_sc-islet_10X-Multiome/err/%x.%A_%a.err
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=02:00:00
#SBATCH --array=1-27%27

#####
# USAGE:
# sbatch --job-name=igvf_sc-islet_10X-Multiome_preprocess_anndatas 4_preprocess_anndatas.sh
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
script=/cellar/users/aklie/data/igvf/beta_cell_networks/notebooks/igvf_sc-islet_10X-Multiome/snapatac2/3_preprocess_anndata.py
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py38

h5ad_files=(
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM0B.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM11A.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM12B.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM14B.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM21A.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM23A.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM24A.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM25A.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM31A.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM32A.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM33A.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM34A.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM35A.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM42B.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM43B.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM44A.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_DM45A.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_MO1.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_MO14.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_MO18.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_MO22.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_MO26.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_MO29.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_MO3.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_MO33.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_MO38.h5ad'
    '/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/adata_atac_MO9.h5ad'
)
h5ad_file=${h5ad_files[$SLURM_ARRAY_TASK_ID-1]}
out_dir=/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23/snapatac2/processed
min_counts=5000
min_tsse=10
max_counts=100000
bin_size=5000
n_features=50000

# Run the script
CMD="python $script \
--input_h5ad $h5ad_file \
--output_dir $out_dir \
--min_counts $min_counts \
--min_tsse $min_tsse \
--max_counts $max_counts \
--bin_size $bin_size \
--n_features $n_features"
echo -e "Running:\n $CMD\n"
echo
eval $CMD

# Date
date
