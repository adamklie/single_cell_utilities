#! /bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/mo_EndoC-bH1_ATAC-seq/out/%x.%A.out
#SBATCH --error=/cellar/users/aklie/data/igvf/beta_cell_networks/scripts/mo_EndoC-bH1_ATAC-seq/err/%x.%A.err
#SBATCH --mem=64G
#SBATCH -n 4
#SBATCH -t 14-00:00:00

#####
# INFO:
# Script to run macs2 callpeak on single input sample in bam format
# More info on parameters: https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/docs/input.md
# https://manpages.ubuntu.com/manpages/xenial/man1/macs2_callpeak.1.html
#####

#####
# USAGE: sbatch --job-name=mo_EndoC-bH1_ATAC-seq_macs2 encode_atac_pipeline_macs2_command.sh $bam $name $out_dir
#####

source activate chrombpnet
bam=$1
name=$2
out_dir=$3
cmd="macs2 callpeak \
    -t $bam \
    -f BAM \
    -n $name \
    -g "hs" \
    -p 0.01 \
    --outdir $out_dir \
    --shift 75 \
    --extsize 150 \
    --nomodel \
    -B \
    --SPMR \
    --keep-dup "all" \
    --call-summits"
echo -e $cmd
eval $cmd
