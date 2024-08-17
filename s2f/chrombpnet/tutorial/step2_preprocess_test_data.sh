# define data directory
data_dir=/cellar/users/aklie/data/igvf/tutorials/chrombpnet/data

# merge and index bam files
samtools merge -f ${data_dir}/downloads/merged_unsorted.bam ${data_dir}/downloads/rep1.bam  ${data_dir}/downloads/rep2.bam  ${data_dir}/downloads/rep3.bam
samtools sort -@4 ${data_dir}/downloads/merged_unsorted.bam -o ${data_dir}/downloads/merged.bam
samtools index ${data_dir}/downloads/merged.bam

# remove blacklist regions from overlap peaks
bedtools slop -i ${data_dir}/downloads/blacklist.bed.gz -g ${data_dir}/downloads/hg38.chrom.sizes -b 1057 > ${data_dir}/downloads/temp.bed
bedtools intersect -v -a ${data_dir}/downloads/overlap.bed.gz -b ${data_dir}/downloads/temp.bed  > ${data_dir}/peaks_no_blacklist.bed

# create subset of chromosome sizes
head -n 24  ${data_dir}/downloads/hg38.chrom.sizes >  ${data_dir}/downloads/hg38.chrom.subset.sizes

# create splits
mkdir -p ${data_dir}/splits
chrombpnet prep splits -c ${data_dir}/downloads/hg38.chrom.subset.sizes -tcr chr1 chr3 chr6 -vcr chr8 chr20 -op ${data_dir}/splits/fold_0

# generate background non-peaks
chrombpnet prep nonpeaks -g ${data_dir}/downloads/hg38.fa -p ${data_dir}/peaks_no_blacklist.bed -c  ${data_dir}/downloads/hg38.chrom.sizes -fl ${data_dir}/splits/fold_0.json -br ${data_dir}/downloads/blacklist.bed.gz -o ${data_dir}/output
