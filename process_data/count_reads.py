import gzip
import os
import concurrent.futures
import csv

def count_reads_in_fastq_file(file_path):
    """
    Count reads in a single gzipped FASTQ file.
    """
    count = 0
    with gzip.open(file_path, 'rt') as f:  # 'rt' mode is used to read gzip file as text file
        for line in f:
            count += 1
    return file_path, count // 4  # Each FASTQ record consists of 4 lines


def main(directory, output_tsv):
    fastq_files = [f for f in os.listdir(directory) if f.endswith('.fastq.gz')][:5]

    # Use ThreadPoolExecutor for parallel processing
    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Map file paths to the function and retrieve results
        results = list(executor.map(count_reads_in_fastq_file, [os.path.join(directory, file) for file in fastq_files]))

    # Write results to TSV
    with open(output_tsv, 'w', newline='') as f_out:
        tsv_writer = csv.writer(f_out, delimiter='\t')
        for result in results:
            tsv_writer.writerow(result)


if __name__ == "__main__":
    directory = '/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/fastq/21Jul23'  # specify your directory with fastq.gz files
    output_tsv = '/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/metadata/11Oct23/reads_count.tsv'  # specify path to the output TSV file
    main(directory, output_tsv)
