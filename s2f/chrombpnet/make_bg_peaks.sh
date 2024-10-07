GENOME=/cellar/users/aklie/data/genomes/hg38/hg38.fa
BLACKLIST=/cellar/users/aklie/data/genomes/hg38/blacklist.bed.gz
CHROM_SIZES=/cellar/users/aklie/data/genomes/hg38/hg38.chrom.sizes
OUT_DIR=/cellar/users/aklie/projects/igvf/beta_cell_networks/bin/infer_grns/chrombpnet/bg_regions
FOLD_DIR=/cellar/users/aklie/projects/igvf/beta_cell_networks/bin/infer_grns/chrombpnet/splits
FOLDS=(0 1 2 3 4)
PEAKS=/cellar/users/aklie/data/igvf/beta_cell_networks/multiome_stimulated_sc/peaks/dm023_palmitate_endocrine_SC.beta.narrowPeak

# Remove columns 4 and 12 from the peak file and remove the header
awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$5"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' $PEAKS > ${PEAKS}.wrangled
tail -n +2 ${PEAKS}.wrangled > ${PEAKS}.wrangled.tmp

# Make background regions
for fold in ${FOLDS[@]}; do
    cmd="chrombpnet prep nonpeaks \
        -g $GENOME \
        -p ${PEAKS}.wrangled.tmp \
        -c $CHROM_SIZES \
        -fl $FOLD_DIR/fold_${fold}.json \
        -br $BLACKLIST \
        -o $OUT_DIR/SC.beta_fold_${fold}_bg.bed"
    echo $cmd
    eval $cmd
done

