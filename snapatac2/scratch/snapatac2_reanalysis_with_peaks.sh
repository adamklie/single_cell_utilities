#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/bin/slurm_logs/%x.%A.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=01-00:00:00

#####
# USAGE:
# sbatch --job-name=igvf_sc-islet_10X-Multiome_snapatac2_peak_calling 3_snapatac2_peak_calling.sh
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py39
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# step 1 -- normalize
echo -e "Running step 1 -- Normalize data\n"
cmd="cellcommander normalize \
--input_h5ad_path $outdir_path/merge/merge.h5ad \
--outdir_path $outdir_path/normalize \
--methods tfidf \
--random-state 1234"
echo -e "Running:\n $cmd\n"
#eval $cmd
echo -e "Done with step 1\n"

# step 2 -- Select features
echo -e "Running step 2 -- Select features\n"
cmd="cellcommander select-features \
--input_h5ad_path $outdir_path/normalize/normalize.h5ad \
--outdir_path $outdir_path/select_features \
--methods snapatac2 \
--n-top-genes 50000 \
--skip-plotting \
--random-state 1234"
echo -e "Running:\n $cmd\n"
#eval $cmd
echo -e "Done with step 2\n"

# step 3 -- Reduce dimensionality
echo -e "Running step 4 -- Reduce dimensionality\n"
cmd="cellcommander reduce-dimensions \
--input_h5ad_path $outdir_path/select_features/select_features.h5ad \
--outdir_path $outdir_path/reduce_dimensions \
--method spectral \
--variable-features-key highly_variable_snapatac2 \
--components-to-remove 0 \
--random-state 1234"
echo -e "Running:\n $cmd\n"
eval $cmd
echo -e "Done with step 4\n"