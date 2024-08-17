#! /bin/bash
#SBATCH --mem=128G
#SBATCH --partition carter-compute
#SBATCH -o ./out/%x.%A.%a.out # STDOUT
#SBATCH -e ./err/%x.%A.%a.err # STDERR
#SBATCH --array=1-5%5
#SBATCH --cpus-per-task 20
#SBATCH --job-name generatePeakToGeneLinksAlpha

# Set-up timer
start=$SECONDS

source activate Renv_dev
samples=(R207 R221 R223 R226 R275)
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}
echo -e "Generating celltype specific networks for $sample"
Rscript --vanilla generatePeakToGeneLinks.R $sample alpha

# Ouptut time of command
duration=$(( SECONDS - start ))
echo -e "\nTime elapsed for this job of array in seconds: $duration"
