#!/bin/bash

#####
# USAGE:
# sbatch functional_analysis_pipeline.sh --SLURM_SETINGS <config_yaml>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py39
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/
pseudobulk_script_path=/cellar/users/aklie/projects/igvf/single_cell_utilities/functional_analysis/pseudobulk.py
deseq2_script_path=/cellar/users/aklie/projects/igvf/single_cell_utilities/functional_analysis/pyDESeq2.py
fa_script_path=/cellar/users/aklie/projects/igvf/single_cell_utilities/functional_analysis/functional_analysis.py

# Load YAML parameters
config_yaml=$1
parse_script_path=/cellar/users/aklie/projects/igvf/single_cell_utilities/sample_qc/parse_yaml.py
source <(python3 $parse_script_path $config_yaml)

# Step 1 -- Create the pseudobulks
echo -e "Running step 1 -- Create pseudobulks\n"
for n_cells in $pseudobulk_max_cells_per_pb; do
    cmd="python $pseudobulk_script_path \
--path_h5ad $io_path_h5ad \
--layer $pseudobulk_layer \
--path_gene_lengths $pseudobulk_path_gene_lengths \
--path_out $io_path_out/pseudobulk/${n_cells}_cells \
--groupby_keys $pseudobulk_groupby_keys \
--cellid_key $pseudobulk_cellid_key \
--compare_key $pseudobulk_compare_key \
--target_max_cells_per_pb $n_cells \
--mode $pseudobulk_mode \
--random_state $pseudobulk_random_state"
    echo -e "Running:\n$cmd\n"
    eval $cmd
done

# Step 2 -- Run DESeq2
echo -e "Running step 2 -- Running DESeq2\n"
for n_cells in $pseudobulk_max_cells_per_pb; do
    n_cells=${n_cells}_cells
    for celltype in $deseq2_celltypes; do
        input_h5ad_path=$io_path_out/pseudobulk/${n_cells}/pdatas/pseudobulk_${celltype}_filtered_genes.h5ad
        cmd="python $deseq2_script_path \
--input_h5ad_path $input_h5ad_path \
--outdir_path $io_path_out/deseq2/${n_cells}/${celltype} \
--design_factors $deseq2_design_factors \
--reference_factor $deseq2_reference_factor \
--reference_value $deseq2_reference_value \
--refit_cooks \
--n_cpus $deseq2_n_cpus"
        echo -e "Running:\n$cmd\n"
        eval $cmd
    done
done
echo -e "Done with step 2\n"

# Step 3 -- Run functional analysis
echo -e "Running step 3 -- Running functional analysis\n"
for n_cells in $pseudobulk_max_cells_per_pb; do
    n_cells=${n_cells}_cells
    for celltype in $deseq2_celltypes; do
        for condition in $functional_analysis_conditions; do
            input_tsv_path=$io_path_out/deseq2/${n_cells}/${celltype}/${condition}_vs_control_shrunkLFC.tsv
            cmd="python $fa_script_path \
--input_contrast_tsv $input_tsv_path \
--outdir_path $io_path_out/functional_analysis/${n_cells}/${celltype} \
--name ${celltype}_${condition}_vs_control"
            echo -e "Running:\n$cmd\n"
            eval $cmd
        done
    done
done
echo -e "Done with step 3\n"

# date
date
