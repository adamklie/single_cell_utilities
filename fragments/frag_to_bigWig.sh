#!/bin/bash

# Usage: ./frag_to_bigWig.sh <frag_file> <chr_sizes_file> <output_dir> <threads>
# Example: ./frag_to_bigWig.sh input_fragments.tsv chrom.sizes /path/to/output 16

frag_file="$1"
chr_sizes="$2"
bigWig_file="$3"
bedGraph_file="$4"
threads="$5"

# Make sure the output directory exists
mkdir -p "$(dirname "$bigWig_file")"

# Generate the BedGraph and BigWig files
LC_ALL=C
zcat "$frag_file" | awk '$1 !~ /_/' | \
    bedtools genomecov -bg -i stdin -g "$chr_sizes" | \
    sort -k1,1 -k2,2n --parallel="$threads" > "$bedGraph_file"

# Convert BedGraph to BigWig
bedGraphToBigWig "$bedGraph_file" "$chr_sizes" "$bigWig_file"

# Clean up BedGraph file
rm "$bedGraph_file"
