#! /bin/bash
#SBATCH --partition=carter-compute
#SBATCH -n 1
#SBATCH -t 14-00:00:00

##############################################
# USAGE: 
# sbatch \
#     --job-name=dataset_create_h5ad_from_rds \
#     --mem=64G \
#     --output=/path/to/out/%x.%A.out \
#     --error=/path/to/out/err/%x.%A.err \
#     create_h5ad_from_rds.sh $rds_file $out_dir $out_name
# Date 07/27/2023
##############################################

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-R413

# Set-up dirs
rds_file=$1 
out_dir=$2
out_name=$3
script=/cellar/users/aklie/data/igvf/scripts/create_h5ad_from_rds.R

# Print messages
echo -e "Loading rds object from $rds_file"
echo -e "Saving h5ad file $out_name to $out_dir\n"

# Convert RDS to H5AD
CMD="Rscript --vanilla $script \
    $rds_file \
    $out_dir \
    $out_name"
echo -e "Running:\n $CMD\n"
$CMD

date
