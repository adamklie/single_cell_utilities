#!/bin/bash

#####
# USAGE:
# sbatch frag_to_count_bw.sh --SLURM_SETINGS <input_tsv> <output_dir> <genome> <chr_sizes>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/chrombpnet
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

input_tsv=$1
output_dir=$2
input_frag_paths=($(cut -f1 $input_tsv))
names=($(cut -f2 $input_tsv))
input_frag_path=${input_frag_paths[$SLURM_ARRAY_TASK_ID-1]}
name=${names[$SLURM_ARRAY_TASK_ID-1]}
output_prefix=$output_dir/${name}
genome=$3
chr_sizes=$4

# Echo inputs
echo -e "frag_file: $input_frag_path"
echo -e "name: $name"
echo -e "output_dir: $output_dir"
echo -e "genome: $genome"
echo -e "chr_sizes: $chr_sizes"
echo -e "output_prefix: $output_prefix\n"

# env
script_path=/cellar/users/aklie/opt/chrombpnet/chrombpnet/helpers/preprocessing/reads_to_bigwig.py

# If output dir does not exist, create it
if [ ! -d $output_dir ]; then
    mkdir -p $output_dir
fi

# cmd
cmd="python $script_path \
--genome $genome \
--input-fragment-file $input_frag_path \
--chrom-sizes $chr_sizes \
--output-prefix $output_prefix \
--data-type ATAC"
echo $cmd
eval $cmd

# Date
echo -e "\n"
date
