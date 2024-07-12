#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH --job-name=create_pseudobulks_from_frag_files_igvf_sc-islet_10X-Multiome
#SBATCH --output=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/igvf_sc-islet_10X-Multiome/out/%x.%A.out
#SBATCH --error=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/igvf_sc-islet_10X-Multiome/err/%x.%A.err
#SBATCH --cpus-per-task=6
#SBATCH --mem=256G
#SBATCH --time=14-00:00:00

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scenicplus

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
metadata_file='/cellar/users/aklie/data/igvf/beta_cell_networks/platinum/igvf_sc-islet_10X-Multiome/10Aug23/scATAC/metadata.csv'
out_dir='/cellar/users/aklie/data/igvf/beta_cell_networks/aligned/igvf_sc-islet_10X-Multiome/10Aug23/pycistopic'
pseudobulk_column='predicted.cell.type'
num_cpus=1
sample_id_column='sample'
temp_dir='/cellar/users/aklie/tmp/'

# Run the script
CMD="python /cellar/users/aklie/data/igvf/bin/data_wrangling/create_pseudobulks_from_frag_files.py \
--fragment_files "${fragment_files[@]}" \
--sample_ids "${sample_ids[@]}" \
--metadata_file $metadata_file \
--out_dir $out_dir \
--pseudobulk_column $pseudobulk_column \
--num_cpus $num_cpus \
--sample_id_column $sample_id_column \
--cellranger_arc_annot \
--temp_dir $temp_dir"
echo -e "Running:\n $CMD\n"
$CMD

# Date
date
