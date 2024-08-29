#!/bin/bash

#####
# USAGE:
# sbatch --SLURM_SETINGS build_pycisTopic_obj_from_counts.sh <mtx_path> <barcodes_path> <regions_path> <cell_metadata_path> <outdir_path> <name>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/test_scenicplus_dev
script_path=/cellar/users/aklie/projects/igvf/single_cell_utilities/pycistopic/build_pycisTopic_obj_from_counts.py

# Inputs
mtx_path=$1
barcodes_path=$2
regions_path=$3
cell_metadata_path=$4
outdir_path=$5
name=$6
blacklist_path=/cellar/users/aklie/data/ref/genomes/hg38/blacklist/blacklist.bed.gz

# Make outdir path if it doesn't exist
if [ ! -d $outdir_path ]; then
    mkdir -p $outdir_path
fi

cmd="python $script_path \
--counts_matrix_path $mtx_path \
--barcodes_path $barcodes_path \
--regions_path $regions_path \
--cell_metadata_path $cell_metadata_path \
--blacklist_path $blacklist_path \
--out_dir $outdir_path \
--project_name $name"
echo -e "Running:\n $cmd\n"
eval $cmd

# Date
date
