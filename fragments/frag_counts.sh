#!/bin/bash

# Usage: ./frag_counts.sh <input_tsv> <output_tsv>
# Example: ./frag_counts.sh input_fragments.tsv output_counts.tsv

input_tsv=$1
output_tsv=$2

# Get the number of fragments for each sample
frag_counts=($(cut -f1 $input_tsv | xargs -I {} wc -l {} | cut -d' ' -f1))
names=($(cut -f2 $input_tsv))

# Write the output tsv (name/tab/frag_count)
if [ -f $output_tsv ]; then
    rm -f $output_tsv
fi
echo -e "name\tfragments" > $output_tsv
for i in ${!frag_counts[@]}; do
    echo -e "${names[$i]}\t${frag_counts[$i]}" >> $output_tsv
done

# Print the output tsv path
echo -e "Output written to $output_tsv"
