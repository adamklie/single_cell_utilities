GENOME=/cellar/users/aklie/data/genomes/hg38/hg38.fa
BLACKLIST=/cellar/users/aklie/data/genomes/hg38/blacklist.bed.gz
CHROM_SIZES=/cellar/users/aklie/data/genomes/hg38/hg38.chrom.sizes
OUT_DIR=/cellar/users/aklie/projects/igvf/beta_cell_networks/bin/infer_grns/chrombpnet/pipeline
FOLD=/cellar/users/aklie/projects/igvf/beta_cell_networks/bin/infer_grns/chrombpnet/splits/fold_0.json
PEAKS=/cellar/users/aklie/data/igvf/beta_cell_networks/multiome_stimulated_sc/peaks/dm023_palmitate_endocrine_SC.beta.narrowPeak.wrangled.tmp
BAM=/cellar/users/aklie/data/igvf/beta_cell_networks/multiome_stimulated_sc/coverage/dm023_palmitate/dm023_palmitate_SC.beta.bam
NEG_REGIONS=/cellar/users/aklie/projects/igvf/beta_cell_networks/bin/infer_grns/chrombpnet/bg_regions/SC.beta_fold_0_bg.bed_negatives.bed
BIAS_MODEL=/cellar/users/aklie/projects/igvf/beta_cell_networks/bin/infer_grns/chrombpnet/bias_models/scATAC_dermal_fibroblast.h5

cmd="chrombpnet pipeline \
        -ibam $BAM \
        -d "ATAC" \
        -g $GENOME \
        -p $PEAKS \
        -c $CHROM_SIZES \
        -n $NEG_REGIONS \
        -fl $FOLD \
        -b $BIAS_MODEL \
        -o $OUT_DIR/SC.beta_fold_0"
echo $cmd
eval $cmd