#!/bin/bash

# Usage: ./frag_to_norm_bigWig.sh <frag_file> <chr_sizes_file> <output_dir> <threads>
# Example: ./frag_to_norm_bigWig.sh input_fragments.tsv chrom.sizes /path/to/output 16

frag_file="$1"
chr_sizes="$2"
bigWig_file="$3"
bedGraph_file="$4"
threads="$4"

# Make sure the output directory exists
mkdir -p "$(dirname "$bigWig_file")"

# Calculate the scale factor
frag_count=$(zcat "$frag_file" | awk '$1 !~ /_/' | wc -l)
scale_factor=$(echo "$frag_count / 1000000" | bc)

# Generate the normalized BedGraph and BigWig files
LC_ALL=C
zcat "$frag_file" | awk '$1 !~ /_/' | \
    bedtools genomecov -bg -i stdin -g "$chr_sizes" -scale "$scale_factor" | \
    sort -k1,1 -k2,2n --parallel="$threads" > "$bedGraph_file"

# Convert BedGraph to BigWig
bedGraphToBigWig "$bedGraph_file" "$chr_sizes" "$bigWig_file"

# Clean up BedGraph file
rm "$bedGraph_file"
