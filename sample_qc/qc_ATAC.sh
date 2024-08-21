#!/bin/bash

#####
# USAGE:
# sbatch qc_ATAC.sh --SLURM_SETINGS ... <input_tsv> <config_yaml> <output_dir>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py39
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# Inputs
input_tsv=$1
outdir_path=$2
input_frag_paths=($(cut -f1 $input_tsv))
sample_ids=($(cut -f2 $input_tsv))
input_frag_path=${input_frag_paths[$SLURM_ARRAY_TASK_ID-1]}
sample_id=${sample_ids[$SLURM_ARRAY_TASK_ID-1]}
config_yaml=$2
outdir_path=$3/${sample_id}/atac


# Echo inputs and number of inputs
echo -e "total inputs: ${#input_frag_paths[@]}"
echo -e "input_frag_path: $input_frag_path"
echo -e "sample_id: $sample_id"
echo -e "outdir_path: $outdir_path"
echo -e "config_yaml: $config_yaml"
echo -e "annot_path: $annot_path\n"

# If output dir does not exist, create it
if [ ! -d $outdir_path ]; then
    mkdir -p $outdir_path
fi

# Command
per_barcode_metrics_path=$(dirname $input_frag_path)/per_barcode_metrics.csv
cmd="cellcommander recipes \
--input_paths $input_frag_path \
--outdir_path $outdir_path \
--sample_name $sample_id \
--method snapatac2 \
--mode single-sample \
--metadata_path $per_barcode_metrics_path \
--metadata_source cellranger \
--params_path $config_yaml"
echo -e "Running command:\n$cmd\n"
eval $cmd

date
