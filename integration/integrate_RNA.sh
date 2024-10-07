#!/bin/bash

#####
# USAGE:
# sbatch integrate_RNA.sh --SLURM_SETINGS ... <input_tsv> <output_dir> <batch_info_script> <sample_metadata_path>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/cellcommander
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# Inputs
input_tsv=$1
outdir_path=$2
input_files=($(cut -f1 $input_tsv))
sample_ids=($(cut -f2 $input_tsv))
batch_info_script=$3
sample_metadata_path=$4

# Echo inputs and number of inputs
echo -e "total inputs: ${#input_files[@]}"
echo -e "input_files: ${input_files[@]}"
echo -e "sample_ids: ${sample_ids[@]}"
echo -e "outdir_path: $outdir_path"

# Step 1 -- merge
echo -e "Running step 1 -- Merge data\n"
cmd="cellcommander merge \
--input_paths ${input_files[@]} \
--outdir_path $outdir_path/merge \
--names ${sample_ids[@]} \
--layer soupx_counts"
echo -e "Running:\n $cmd\n"
eval $cmd
echo -e "Done with step 1\n"

# Step 2 -- normalize
echo -e "Running step 2 -- Normalize data\n"
cmd="cellcommander normalize \
--input_h5ad_path $outdir_path/merge/merge.h5ad \
--outdir_path $outdir_path/normalize \
--methods log1p sctransform \
--vars_to_regress pct_counts_mt pct_counts_ribo \
--random-state 1234"
echo -e "Running:\n $cmd\n"
eval $cmd
echo -e "Done with step 2\n"

# Step 3 -- dim reduction
echo -e "Running step 3 -- Dimensionality reduction\n"
cmd="cellcommander reduce-dimensions \
--input_h5ad_path $outdir_path/normalize/normalize.h5ad \
--outdir_path $outdir_path/reduce_dimensions \
--output_prefix seurat_default_pca \
--method seurat_default \
--obsm_key sctransform_scale_data \
--variable-features-key sctransform_genes \
--random-state 1234"
echo -e "Running:\n $cmd\n"
eval $cmd
echo -e "Done with step 3\n"

# Step 4a -- make metadata file
if [ ! -z "$batch_info_script" ]; then
    echo -e "Running step 4a -- Make metadata file\n"
    cmd="python $batch_info_script $outdir_path/reduce_dimensions/cell_metadata.tsv $sample_metadata_path"
    echo -e "Running:\n $cmd\n"
    eval $cmd
    echo -e "Done with step 4a\n"
fi

# Step 4b -- correct batch effect 
cmd="cellcommander integrate \
--input_h5ad_path $outdir_path/reduce_dimensions/seurat_default_pca.h5ad \
--batch_file $outdir_path/reduce_dimensions/bc_sample_info.tsv \
--outdir_path $outdir_path/integrate \
--vars_to_correct sample \
--obsm_key X_seurat_default \
--method harmonyR \
--corrected_obsm_key X_seurat_default_harmony"
echo -e "Running command:\n $cmd\n"
eval $cmd
echo -e "Done with step 4\n"
date
