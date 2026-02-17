#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH -o /cellar/users/aklie/projects/igvf/beta_cell_networks/infer_cellular_programs/cistopic/bin/igvf_sc-islet_10X-Multiome/out/%x.%A.%a.out
#SBATCH -e /cellar/users/aklie/projects/igvf/beta_cell_networks/infer_cellular_programs/cistopic/bin/igvf_sc-islet_10X-Multiome/err/%x.%A.%a.err
#SBATCH --time 14-00:00:00
#SBATCH --array=1-5%21
#SBATCH --cpus-per-task=8
#SBATCH --mem=102G

##############################################
# USAGE: sbatch --job-name=igvf_sc-islet_10X-Multiome_runModels_lda_mallet runModels_lda_mallet.sh
# You will need to manually create a models.pkl file from the individual models unlike the cgs run
# Date 02/17/2022
##############################################

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env
source activate /cellar/users/aklie/opt/miniconda3/envs/scenicplus
script=/cellar/users/aklie/projects/igvf/beta_cell_networks/infer_cellular_programs/cistopic/bin/igvf_sc-islet_10X-Multiome/runModels_lda_mallet.py

# Configure input arguments that likely will not change
inputcisTopic_obj=/cellar/users/aklie/projects/igvf/beta_cell_networks/infer_cellular_programs/cistopic/results/igvf_sc-islet_10X-Multiome/10Aug23/cistopic_obj.pkl
mallet_binary=/cellar/users/aklie/opt/Mallet-202108/bin/mallet
save_path=/cellar/users/aklie/projects/igvf/beta_cell_networks/infer_cellular_programs/cistopic/results/igvf_sc-islet_10X-Multiome/10Aug23/mallet
corpus_dir=$save_path/corpus

# Changeable input arguments
n_topics="2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100"
n_iter=500
n_cpu=$SLURM_CPUS_PER_TASK
alpha=50
alpha_by_topic=True
eta=0.1
eta_by_topic=False
seed=555
reuse_corpus=True

# Configure place to put outputs
name=${n_iter}iter_${alpha}alpha_${alpha_by_topic}abt_${eta}eta_${eta_by_topic}ebt_${seed}seed
save_path=$save_path/$name

# Save a the input arguments to the save_path under $save_path/params.yaml if doesn't exist
if [ ! -d $save_path/params.yaml ]
then
    mkdir -p $save_path
    echo -e "lda_model: mallet\ninputcisTopic_obj: $inputcisTopic_obj\nn_topics: $n_topics\nn_cpu: $n_cpu\nn_iter: $n_iter\nalpha: $alpha\nalpha_by_topic: $alpha_by_topic\neta: $eta\neta_by_topic: $eta_by_topic\nsave_path: $save_path\nseed: $seed\ncorpus_dir: $corpus_dir\nreuse_corpus: $reuse_corpus\nmallet_binary: $mallet_binary" > $save_path/params.yaml
fi

# Choose the number of topics to run
selected_topic=$(echo $n_topics | cut -d ',' -f $SLURM_ARRAY_TASK_ID | tr -d ' ')
temp_dir=$corpus_dir/Topic${selected_topic}

# Print input arguments
echo -e "Input cisTopic object: $inputcisTopic_obj"
echo -e "Path to Mallet binary: $mallet_binary"
echo -e "Path to current corpus dir: $temp_dir"
echo -e "Path to save intermediate files: $save_path"
echo -e "Number of iterations: $n_iter"
echo -e "Number of CPU cores: $n_cpu"
echo -e "Alpha: $alpha"
echo -e "Divide alpha by the number of topics: $alpha_by_topic"
echo -e "Eta: $eta"
echo -e "Divide eta by the number of topics: $eta_by_topic"
echo -e "Seed: $seed"
echo -e "Reuse Mallet corpus: $reuse_corpus"
echo -e "Number of topics: $selected_topic\n"

# Run Python script from the CLI
CMD="python $script \
--inputcisTopic_obj $inputcisTopic_obj \
--n_topics $n_topics \
--n_cpu $n_cpu \
--n_iter $n_iter \
--alpha $alpha \
--alpha_by_topic $alpha_by_topic \
--eta $eta \
--eta_by_topic $eta_by_topic \
--save_path $save_path \
--seed $seed \
--temp_dir $temp_dir \
--reuse_corpus $reuse_corpus"
echo -e "Running:\n$CMD"
echo
#eval $CMD

date
