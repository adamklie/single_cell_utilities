# define data directory
main_dir=/cellar/users/aklie/data/igvf/tutorials/chrombpnet

# download bias model
mkdir -p ${main_dir}/bias_model
wget https://storage.googleapis.com/chrombpnet_data/input_files/bias_models/ATAC/ENCSR868FGK_bias_fold_0.h5 -O ${main_dir}/bias_model/ENCSR868FGK_bias_fold_0.h5

# train a chrombpnet model
chrombpnet pipeline \
        -ibam ${main_dir}/data/downloads/merged.bam \
        -d "ATAC" \
        -g ${main_dir}/data/downloads/hg38.fa \
        -c ${main_dir}/data/downloads/hg38.chrom.sizes \
        -p ${main_dir}/data/peaks_no_blacklist.bed \
        -n ${main_dir}/data/output_negatives.bed \
        -fl ${main_dir}/data/splits/fold_0.json \
        -b ${main_dir}/bias_model/ENCSR868FGK_bias_fold_0.h5 \
        -o ${main_dir}/chrombpnet_model/
        