#!/bin/bash

#####
# USAGE:
# sbatch frag_to_tagAlign.sh --SLURM_SETINGS <input_tsv> <output_dir> <threads>
#####

# based on https://github.com/EngreitzLab/sc-E2G/blob/main/workflow/rules/frag_to_tagAlign.smk

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/chrombpnet
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# Define the input arguments
input_tsv=$1
output_dir=$2
input_frag_paths=($(cut -f1 $input_tsv))
names=($(cut -f2 $input_tsv))
input_frag_path=${input_frag_paths[$SLURM_ARRAY_TASK_ID-1]}
name=${names[$SLURM_ARRAY_TASK_ID-1]}
tagAlign_sort_file=$output_dir/${name}.tagAlign.sort.gz
threads=$3

# Echo inputs 
echo -e "frag_file: $input_frag_path"
echo -e "name: $name"
echo -e "output_dir: $output_dir"
echo -e "threads: $threads"

# If output dir does not exist, create it
if [ ! -d $output_dir ]; then
    mkdir -p $output_dir
fi

# Create, sort, and compress the tagAlign file
LC_ALL=C zcat "$input_frag_path" | sed '/^#/d' | \
    awk -v OFS='\t' '{mid=int(($2+$3)/2); print $1,$2,mid,"N",1000,"+"; print $1,mid+1,$3,"N",1000,"-"}' | \
    sort -k 1,1V -k 2,2n -k3,3n --parallel "$threads" | \
    bgzip -c > "$tagAlign_sort_file"

# Index the tagAlign file
tabix -p bed "$tagAlign_sort_file"
