#! /bin/bash

#####
# INFO:
# Script to 
#####

date
echo ""

start=$SECONDS

# Get the array of SRA IDs
sra_ids=($(cat $1))

# Print the length of the array
echo "Number of SRA IDs: ${#sra_ids[@]}"
echo ""

# Loop through the array and prefetch and fasterq-dump each SRA ID
for sra_id in "${sra_ids[@]}"
do
    echo "prefetching $sra_id"
    prefetch $sra_id --max-size 100g

    echo "fasterq-dump $sra_id"
    fasterq-dump --threads 16 --include-technical --split-files --progress $sra_id -O fastq

    # Need to gzip all the new $sra_id*.fastq files
    echo -e "gzipping $sra_id\n"
    gzip fastq/${sra_id}*.fastq

    echo -e "completed $sra_id\n"
done

duration=$(( SECONDS - start ))
echo -e "Time elapsed in seconds: $duration"
date