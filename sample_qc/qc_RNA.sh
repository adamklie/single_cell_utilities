#!/bin/bash

#####
# USAGE:
# sbatch qc_RNA.sh --SLURM_SETINGS ... <input_tsv> <config_yaml> <output_dir>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/cellcommander
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# Inputs
input_tsv=$1
input_h5_paths=($(cut -f1 $input_tsv))
sample_ids=($(cut -f2 $input_tsv))
input_h5_path=${input_h5_paths[$SLURM_ARRAY_TASK_ID-1]}
sample_id=${sample_ids[$SLURM_ARRAY_TASK_ID-1]}
config_yaml=$2
outdir_path=$3/${sample_id}/rna

# Load YAML parameters
parse_script_path=/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/bin/3_sample_qc/scripts/parse_yaml.py
source <(python3 $parse_script_path $config_yaml)

# Echo inputs and number of inputs
echo -e "total inputs: ${#input_h5_paths[@]}"
echo -e "input_h5_path: $input_h5_path"
echo -e "sample_id: $sample_id"
echo -e "outdir_path: $outdir_path\n"

# If output dir does not exist, create it
if [ ! -d $outdir_path ]; then
    mkdir -p $outdir_path
fi

# Step 1 -- QC and filtering
per_barcode_metrics_path=$(dirname $input_h5_path)/per_barcode_metrics.csv
echo -e "Running step 1 -- QC and filtering\n"
cmd="cellcommander qc \
--input_h5_path $input_h5_path \
--outdir_path $outdir_path/threshold_qc \
--sample_name $sample_id \
--metadata_path $per_barcode_metrics_path \
--metadata_source cellranger \
--output_prefix threshold_qc \
--mode rna \
--multimodal_input \
--filtering_strategy $qc_filtering_strategy \
--n_features_low $qc_n_features_low \
--n_features_hi $qc_n_features_hi \
--pct_counts_mt_hi $qc_pct_counts_mt_hi \
--pct_counts_ribo_hi $qc_pct_counts_ribo_hi \
--random_state $random_state"
echo -e "Running:\n $cmd\n"
eval $cmd
echo -e "Done with step 1\n"

# Step 2 -- Background removal
echo -e "Running step 2 -- Background removal\n"
cmd="cellcommander remove-background \
--input_h5ad_path $outdir_path/threshold_qc/threshold_qc.h5ad \
--outdir_path $outdir_path/remove_background \
--method $remove_background_method \
--raw-h5-path $input_h5_path \
--markers_path $remove_background_markers_path \
--layer $remove_background_layer \
--random-state $random_state"
echo -e "Running:\n $cmd\n"
eval $cmd
echo -e "Done with step 2\n"

# Step 3 -- Detect doublets
echo -e "Running step 3 -- Detect doublets\n"
cmd="cellcommander detect-doublets \
--input_h5ad_path $outdir_path/remove_background/remove_background.h5ad \
--outdir_path $outdir_path/detect_doublets \
--output_prefix consensus \
--method $detect_doublelets_method \
--consensus-methods $detect_doublelets_consensus_methods \
--consensus-strategy $detect_doublelets_consensus_strategy \
--random-state $random_state"
echo -e "Running:\n $cmd\n"
eval $cmd
echo -e "Done with step 3\n"

# Step 4 -- Normalize data
echo -e "Running step 4 -- Normalize data\n"
cmd="cellcommander normalize \
--input_h5ad_path $outdir_path/detect_doublets/consensus.h5ad \
--outdir_path $outdir_path/normalize \
--output_prefix normalize \
--save-normalized-mtx \
--methods $normalize_methods \
--vars_to_regress $normalize_vars_to_regress \
--random-state $random_state"
echo -e "Running:\n $cmd\n"
eval $cmd
echo -e "Done with step 4\n"

# Step 5 -- Reduce dimensionality
echo -e "Running step 5 -- Reduce dimensionality\n"
cmd="cellcommander reduce-dimensions \
--input_h5ad_path $outdir_path/normalize/normalize.h5ad \
--outdir_path $outdir_path/reduce_dimensions \
--output_prefix seurat_default_pca \
--method $reduce_dimensions_method \
--obsm_key $reduce_dimensions_obsm_key \
--variable-features-key $reduce_dimensions_variable_features_key \
--random-state $random_state"
echo -e "Running:\n $cmd\n"
eval $cmd
echo -e "Done with step 5\n"

date
