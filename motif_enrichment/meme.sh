#!/bin/bash

#####
# USAGE:
# sbatch meme.sh --SLURM_SETINGS <input_bed> <genome> <output>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/chrombpnet
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# Inputs
input_tsv=$1
input_peak_paths=($(cut -f1 $input_tsv))
names=($(cut -f2 $input_tsv))
input_peak_path=${input_peak_paths[$SLURM_ARRAY_TASK_ID-1]}
name=${names[$SLURM_ARRAY_TASK_ID-1]}
genome=$2
out=$3/${name}

# Echo inputs
echo -e "input_peak_path: $input_peak_path"
echo -e "name: $name"
echo -e "genome: $genome"
echo -e "out: $out"

# make output directory
mkdir -p $out

# use bedtools to extract temporary fasta file from input peak file
fasta=$out/temp.fa
bedtools getfasta -fi $genome -bed $input_peak_path -fo $fasta
 
# Cmd
cmd="meme \
$fasta \
-dna \
-revcomp \
-mod anr \
-nmotifs 3 \
-minw 6 \
-maxw 50 \
-objfun classic \
-markov_order 0 \
-oc $out \
-seed 1234"

# Run the command
echo "Running command: $cmd"
eval $cmd

# clean up 
rm $fasta

# Date
date
