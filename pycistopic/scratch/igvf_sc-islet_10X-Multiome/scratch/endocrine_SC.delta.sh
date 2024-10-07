OUT_DIR="/cellar/users/aklie/projects/igvf/beta_cell_networks/scratch/infer_cellular_programs/cistopic/results"
PROJECT_NAME="SC.delta_script_test"

################################
### Create pycisTopic object ###
################################
BARCODES_PATH="/cellar/users/aklie/projects/igvf/beta_cell_networks/data/multiome_stimulated_sc/barcodes/dm023_palmitate/dm023_palmitate_endocrine_SC.delta_barcodes.tsv"
REGIONS_PATH="/cellar/users/aklie/projects/igvf/beta_cell_networks/data/multiome_stimulated_sc/matrix/dm023_palmitate/dm023_palmitate_endocrine_SC.delta_mpeak.var.tsv"
CELL_METADATA_PATH="/cellar/users/aklie/projects/igvf/beta_cell_networks/data/multiome_stimulated_sc/metadata/dm023_palmitate/dm023_palmitate_endocrine_SC.delta_metadata.csv"
COUNTS_MATRIX_PATH="/cellar/users/aklie/projects/igvf/beta_cell_networks/data/multiome_stimulated_sc/matrix/dm023_palmitate/dm023_palmitate_endocrine_SC.delta_mpeak.count.mtx"
SCRIPT="/cellar/users/aklie/projects/igvf/beta_cell_networks/bin/infer_cellular_programs/cistopic/scripts/create_pycisTopic_obj.py"
echo -e "Running script: ${SCRIPT} to create pycisTopic object"
cmd="python $SCRIPT \
    --out_dir $OUT_DIR \
    --barcodes_path $BARCODES_PATH \
    --regions_path $REGIONS_PATH \
    --cell_metadata_path $CELL_METADATA_PATH \
    --counts_matrix_path $COUNTS_MATRIX_PATH \
    --project_name $PROJECT_NAME"
#echo "Running command: ${cmd}"
#$cmd

#########################
### Run LDA with CGS ###
########################
TMP_DIR='/cellar/users/aklie/tmp/'
SCRIPT="/cellar/users/aklie/projects/igvf/beta_cell_networks/bin/infer_cellular_programs/cistopic/scripts/runModels_lda_cgs.py"
echo -e "Running script: ${SCRIPT} to run LDA with CGS"
cmd="python $SCRIPT \
    --inputcisTopic_obj ${OUT_DIR}/${PROJECT_NAME}.pkl \
    --output ${OUT_DIR}/${PROJECT_NAME}_lda_cgs/ \
    --n_topics 2,5,10,15,20,25,30,35,40,45,50 \
    --n_cpu 12 \
    --n_iter 500 \
    --seed 555 \
    --alpha 50 \
    --alpha_by_topic True \
    --eta 0.1 \
    --eta_by_topic False \
    --save_path ${OUT_DIR}/${PROJECT_NAME}_lda_cgs/intermediate_models/ \
    --temp_dir $TMP_DIR"
echo "Running command: ${cmd}"
$cmd
