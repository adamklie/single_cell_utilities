#!/bin/bash

#####
# USAGE:
# sbatch frag_to_norm_bigWig.sh --SLURM_SETINGS <input_tsv> <output_dir> <threads>
#####

# based on https://github.com/EngreitzLab/sc-E2G/blob/main/workflow/rules/frag_to_norm_bigWig.smk

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
bedGraph_file=$output_dir/${name}.bedGraph
bigWig_file=$output_dir/${name}.fpm.bw
chr_sizes=$3
threads=$4

# Echo inputs 
echo -e "frag_file: $input_frag_path"
echo -e "name: $name"
echo -e "output_dir: $output_dir"
echo -e "chr_sizes: $chr_sizes"
echo -e "threads: $threads\n"

# Calculate the scale factor
echo -e "Calculating scale factor"
frag_count=$(zcat $input_frag_path | awk '$1 !~ /_/' | wc -l)
scale_factor=$(echo "$frag_count / 1000000" | bc)
echo -e "Scale factor: $scale_factor"

# Save the scale factor to output_dir/$name.scale_factor.txt
echo -e "$scale_factor" > "$output_dir/$name.scale_factor.txt"

# Generate the normalized BedGraph and BigWig files
echo -e "Generating normalized BedGraph and BigWig files"
LC_ALL=C
cmd="zcat $input_frag_path | awk '\$1 !~ /_/' | sort -k 1,1 | bedtools genomecov -bg -i stdin -g $chr_sizes -scale $scale_factor | sort -k1,1 -k2,2n --parallel=$threads > $bedGraph_file"
echo -e "Running $cmd\n"
eval $cmd

# Convert BedGraph to BigWig
echo -e "Converting BedGraph to BigWig"
cmd="bedGraphToBigWig $bedGraph_file $chr_sizes $bigWig_file"
echo -e "Running $cmd\n"
eval $cmd

# Clean up BedGraph file
rm -f "$bedGraph_file"

# Date
echo -e "\n"
date
