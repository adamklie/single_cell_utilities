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
bed_files=$1
genome=$2
size=$3
out=$4

# Get the input file corresponding to the current array task ID

# Echo inputs
echo -e "input_bed: $input_bed"
echo -e "genome: $genome"
echo -e "size: $size"
echo -e "out: $out"

# Cmd
cmd="findMotifsGenome.pl \
$input_bed \
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

