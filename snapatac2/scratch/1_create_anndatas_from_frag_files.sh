#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/igvf_sc-islet_10X-Multiome/out/%x.%A_%a.out
#SBATCH --error=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/igvf_sc-islet_10X-Multiome/err/%x.%A_%a.err
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64G
#SBATCH --time=02:00:00
#SBATCH --array=1-27%27

#####
# USAGE:
# sbatch --job-name=igvf_sc-islet_10X-Multiome_create_anndatas_from_frag_files create_anndatas_from_frag_files.sh
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
script=/cellar/users/aklie/data/igvf/bin/data_wrangling/create_anndata_from_frag_file.py
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py38

fragment_files=(
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm45a_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm43b_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm44a_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_mo1_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm31a_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_mo22_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm24a_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm33a_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_mo29_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm23a_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_mo9_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_mo14_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm34a_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm11a_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_mo3_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_mo18_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm12b_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm0b_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm21a_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm35a_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_mo38_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_mo33_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm25a_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm42b_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm32a_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_dm14b_deep/outs/atac_fragments.tsv.gz'
  '/cellar/users/aklie/data/igvf/beta_cell_networks/cellranger/igvf_sc-islet_10X-Multiome/igvf_mo26_deep/outs/atac_fragments.tsv.gz'
)
frag_file=${fragment_files[$SLURM_ARRAY_TASK_ID-1]}
sample_ids=(
  'DM45A'
  'DM43B'
  'DM44A'
  'MO1'
  'DM31A'
  'MO22'
  'DM24A'
  'DM33A'
  'MO29'
  'DM23A'
  'MO9'
  'MO14'
  'DM34A'
  'DM11A'
  'MO3'
  'MO18'
  'DM12B'
  'DM0B'
  'DM21A'
  'DM35A'
  'MO38'
  'MO33'
  'DM25A'
  'DM42B'
  'DM32A'
  'DM14B'
  'MO26'
)
sample_id=${sample_ids[$SLURM_ARRAY_TASK_ID-1]}
out_dir=/cellar/users/aklie/data/igvf/beta_cell_networks/h5ad/igvf_sc-islet_10X-Multiome/16Aug23
out_file=$out_dir/adata_atac_$sample_id.h5ad

# Run the script
CMD="python $script \
--fragment_file $frag_file \
--out_file $out_file"
echo -e "Running:\n $CMD\n"
$CMD

# Date
date
