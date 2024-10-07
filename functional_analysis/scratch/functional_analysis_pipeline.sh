#####
# USAGE:
# sbatch functional_analysis_pipeline.sh --SLURM_SETINGS --input_h5ad_path INPUT_H5AD_PATH --condition CONDITION --outdir_path OUTDIR_PATH
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py38
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/
pseudobulk_script_path=/cellar/users/aklie/projects/igvf/single_cell_utilities/functional_analysis/pseudobulk.py
deseq2_script_path=/cellar/users/aklie/projects/igvf/single_cell_utilities/functional_analysis/pyDESeq2.py
fa_script_path=/cellar/users/aklie/projects/igvf/single_cell_utilities/functional_analysis/functional_analysis.py

# Conditions
input_h5ad_path=$1
condition=$2
outdir_path=$3

# Step 1 -- Create the psuedobulks
echo -e "Running step 1 -- Create pseudobulks\n"
for n_cells in 100 500 1000
do
    cmd="python $pseudobulk_script_path \
--input_h5ad_path $input_h5ad_path \
--outdir_path $outdir_path/pseudobulk/${n_cells}_cells \
--groupby_keys updated_annotation condition \
--cellid_key updated_annotation \
--compare_key condition \
--target_max_cells_per_pb $n_cells \
--mode sum"
    echo -e "Running:\n $cmd\n"
    eval $cmd
done

cmd="python $pseudobulk_script_path \
--input_h5ad_path $input_h5ad_path \
--outdir_path $outdir_path/pseudobulk/sample \
--groupby_keys sample_id updated_annotation condition \
--cellid_key updated_annotation \
--compare_key condition \
--mode sum"
echo -e "Running:\n $cmd\n"
eval $cmd
echo -e "Done with step 1\n"

# Step 2 -- Run DESeq2
echo -e "Running step 2 -- Running DESeq2\n"
for n_cells in 100_cells 500_cells 1000_cells sample
do
    input_h5ad_path=$outdir_path/pseudobulk/${n_cells}/pdatas/pseudobulk_SC.beta_filtered_genes.h5ad
    cmd="python $deseq2_script_path \
--input_h5ad_path $input_h5ad_path \
--outdir_path $outdir_path/deseq2/${n_cells}/SC.beta \
--design_factors condition \
--reference_factor condition \
--reference_value control \
--refit_cooks \
--n_cpus 3"
    echo -e "Running:\n $cmd\n"
    eval $cmd
done
echo -e "Done with step 2\n"

# Step 3 -- Run functional analysis
echo -e "Running step 3 -- Running functional analysis\n"
for n_cells in 100_cells 500_cells 1000_cells sample
do
    if [ $condition == "ext4" ]
    then
        condition="Ex-4-HG"
    fi
    input_tsv_path=$outdir_path/deseq2/${n_cells}/SC.beta/${condition}_vs_control_shrunkLFC.tsv
    cmd="python $fa_script_path \
--input_contrast_tsv $input_tsv_path \
--outdir_path $outdir_path/functional_analysis/${n_cells}/SC.beta \
--name SC_beta"
    echo -e "Running:\n $cmd\n"
    eval $cmd
done
echo -e "Done with step 3\n"

date
