#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH -o /cellar/users/aklie/projects/igvf/beta_cell_networks/infer_grns/scenicplus/bin/igvf_sc-islet_10X-Multiome/out/%x.%A.out
#SBATCH -e /cellar/users/aklie/projects/igvf/beta_cell_networks/infer_grns/scenicplus/bin/igvf_sc-islet_10X-Multiome/err/%x.%A.err
#SBATCH --time=14-00:00:00
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G

##############################################
# USAGE: sbatch --job-name=igvf_sc-islet_10X-Multiome_runModels_lda_cgs runModels_lda_cgs.sh
# Date 02/17/2022
##############################################

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scenicplus
script=/cellar/users/aklie/projects/igvf/beta_cell_networks/infer_grns/scenicplus/bin/igvf_sc-islet_10X-Multiome/runModels_lda_cgs.py

# Configure input arguments
inputcisTopic_obj=/cellar/users/aklie/projects/igvf/beta_cell_networks/infer_grns/scenicplus/results/igvf_sc-islet_10X-Multiome/10Aug23/cistopic_obj.pkl
save_path=/cellar/users/aklie/projects/igvf/beta_cell_networks/infer_grns/scenicplus/results/igvf_sc-islet_10X-Multiome/10Aug23/cgs
output=$save_path/models.pkl
n_topics="2, 4, 8, 16, 32, 48, 64, 80"
n_cpu=8
n_iter=150
alpha=50
alpha_by_topic=True
eta=0.1
eta_by_topic=False
seed=555
temp_dir=/cellar/users/aklie/tmp/

# Print input arguments
echo -e "Input cisTopic object: $inputcisTopic_obj"
echo -e "Path to save intermediate files: $save_path"
echo -e "Output file: $output"
echo -e "Number of topics: $n_topics"
echo -e "Number of CPU cores: $n_cpu"
echo -e "Number of iterations: $n_iter"
echo -e "Alpha: $alpha"
echo -e "Divide alpha by the number of topics: $alpha_by_topic"
echo -e "Eta: $eta"
echo -e "Divide eta by the number of topics: $eta_by_topic"
echo -e "Seed: $seed"
echo -e "Path to TMP dir: $temp_dir"

# Make the output directory if it doesn't exist
if [ ! -d $save_path ]
then
    mkdir -p $save_path
fi

# Save a the input arguments to the save_path under params.yaml
echo -e "lda_model: cgs\ninputcisTopic_obj: $inputcisTopic_obj\noutput: $output\nn_topics: $n_topics\nn_cpu: $n_cpu\nn_iter: $n_iter\nalpha: $alpha\nalpha_by_topic: $alpha_by_topic\neta: $eta\neta_by_topic: $eta_by_topic\nsave_path: $save_path\nseed: $seed\ntemp_dir: $temp_dir" > $save_path/params.yaml

# Run Python script from the CLI
CMD="python $script \
--inputcisTopic_obj $inputcisTopic_obj \
--output $output \
--n_topics $n_topics \
--n_cpu $n_cpu \
--n_iter $n_iter \
--alpha $alpha \
--alpha_by_topic $alpha_by_topic \
--eta $eta \
--eta_by_topic $eta_by_topic \
--save_path $save_path \
--seed $seed \
--temp_dir $temp_dir"
echo -e "Running:\n$CMD"
echo
$CMD

date
