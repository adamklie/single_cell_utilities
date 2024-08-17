#! /bin/bash
#SBATCH --account carter-gpu
#SBATCH --partition carter-gpu
#SBATCH --gpus=1
#SBATCH --output=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/mo_EndoC-bH1_ATAC-seq/out/%x.%A.out
#SBATCH --error=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/mo_EndoC-bH1_ATAC-seq/err/%x.%A.err
#SBATCH --mem=128G
#SBATCH -n 2
#SBATCH -t 14-00:00:00

#####
# INFO:
# Script to run chrombpnet pipeline on single input sample in bam format
#####

#####
# USAGE: sbatch --job-name=${dataset_name}_chrombpnet --error /path/to/err --output /path/to/out chrombpnet_pipeline.sh $coverage $genome $peaks $chromsizes $negatives $fold $bias_model $out_dir
#####

# Set-up env
source activate chrombpnet
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/lib
python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"

# Set up args
ta=$1
genome=$2
peaks=$3
chromsizes=$4
negatives=$5
fold=$6
bias_model=$7
out_dir=$8

# Run cmd
cmd="chrombpnet pipeline \
    -itag $coverage \
    -d "ATAC" \
    -g $genome \
    -c $chromsizes \
    -p $peaks \
    -n $negatives \
    -fl $fold \
    -b $bias_model \
    -o $out_dir"
echo -e $cmd
eval $cmd
