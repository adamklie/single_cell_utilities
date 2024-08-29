#!/bin/bash

#####
# USAGE:
# sbatch motif_finding.sh --SLURM_SETINGS <input_bed> <genome> <output> <size>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)

# Inputs
input_tsv=$1
input_peak_paths=($(cut -f1 $input_tsv))
names=($(cut -f2 $input_tsv))
input_peak_path=${input_peak_paths[$SLURM_ARRAY_TASK_ID-1]}
name=${names[$SLURM_ARRAY_TASK_ID-1]}
genome=$2
size=$3
out=$4/${name}

# Echo inputs
echo -e "input_peak_path: $input_peak_path"
echo -e "name: $name"
echo -e "genome: $genome"
echo -e "size: $size"
echo -e "out: $out"

# Cmd
cmd="findMotifsGenome.pl \
$input_peak_path \
$genome \
$out \
-size $size \
-mask \
-p $SLURM_CPUS_PER_TASK"

# Run the command
echo "Running command: $cmd"
eval $cmd

# Date
date

