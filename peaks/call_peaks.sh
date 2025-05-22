#!/bin/bash

#####
# USAGE:
# sbatch call_peaks.sh --SLURM_SETINGS <input_tsv> <output_dir>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/chrombpnet
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# Define the input arguments
input_tsv=$1
output_dir=$2
input_ta_paths=($(cut -f1 $input_tsv))
names=($(cut -f2 $input_tsv))
input_ta_path=${input_ta_paths[$SLURM_ARRAY_TASK_ID-1]}
name=${names[$SLURM_ARRAY_TASK_ID-1]}

# Echo inputs 
echo -e "input_ta_path: $input_ta_path"
echo -e "name: $name"
echo -e "output_dir: $output_dir"

# Make outdir path if it doesn't exist
if [ ! -d $outdir_path ]; then
    mkdir -p $outdir_path
fi 

# Call peaks
cmd="macs2 callpeak \
-t $input_ta_path \
-f BED \
-n $name \
-g "hs" \
-p 0.01 \
--outdir $output_dir \
--shift 75 \
--extsize 150 \
--nomodel \
-B \
--SPMR \
--keep-dup "all" \
--call-summits"
echo -e $cmd
eval $cmd

date
