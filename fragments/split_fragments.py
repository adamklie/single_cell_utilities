import multiprocessing
import os
import tempfile
import gzip
import shutil
import argparse


def process_fragment_chunk(tmp_file_name, barcode_group_dict, output_dir, debug):
    local_output_files = {}
    for group in set(barcode_group_dict.values()):
        local_output_files[group] = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=output_dir)

    line_count = 0
    missing_barcodes = 0
    with open(tmp_file_name, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.split()
                barcode = fields[3]
                group = barcode_group_dict.get(barcode, None)
                if group:  # Process if barcode group is found
                    local_output_files[group].write(line)
                    line_count += 1
                elif debug:
                    missing_barcodes += 1

    for file in local_output_files.values():
        file.close()

    if debug:
        print(f"Debug: Processed {line_count} lines in {tmp_file_name}.")
        if missing_barcodes > 0:
            print(f"{missing_barcodes} lines in {tmp_file_name} did not have a barcode group.")
    return {group: file.name for group, file in local_output_files.items()}

def compress_file(file_name, output_dir):
    output_file = os.path.join(output_dir, os.path.basename(file_name) + '.gz')
    with open(file_name, 'rb') as f_in, gzip.open(output_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(file_name)
    return output_file

def concatenate_files(group_files, final_output_file):
    with open(final_output_file, 'w') as outfile:
        for fname in group_files:
            with open(fname) as infile:
                shutil.copyfileobj(infile, outfile)
            os.remove(fname)

def split_fragment_file(filename, num_chunks, is_gzipped, debug):
    if is_gzipped:
        with gzip.open(filename, 'rt') as f_in:
            with open(filename[:-3], 'w') as f_out:
                shutil.copyfileobj(f_in, f_out)
        filename = filename[:-3]

    total_lines = sum(1 for line in open(filename) if not line.startswith('#'))
    lines_per_chunk = total_lines // num_chunks
    tmp_files = []
    current_line = 0
    if debug:
        print(f"Debug: Splitting {filename} into {num_chunks} chunks of approximately {lines_per_chunk} lines each.")

    with open(filename, 'r') as file:
        for _ in range(num_chunks):
            tmp_file = tempfile.NamedTemporaryFile(delete=False, mode='w')
            tmp_files.append(tmp_file.name)
            for line in file:
                if not line.startswith('#'):
                    tmp_file.write(line)
                    current_line += 1
                if current_line >= lines_per_chunk * (len(tmp_files)):
                    break
            tmp_file.close()
            if current_line >= total_lines:
                break

    return tmp_files, total_lines

def main(fragment_file, barcode_group_file, output_dir, num_processes, debug):
    is_gzipped = fragment_file.endswith('.gz')
    barcode_group_dict = {}
    with open(barcode_group_file, 'r') as f:
        for line in f:
            barcode, group = line.strip().split()
            barcode_group_dict[barcode] = group

    temp_dir = tempfile.mkdtemp()
    output_temp_files = {group: os.path.join(temp_dir, f"{group}.tsv") for group in set(barcode_group_dict.values())}
    
    tmp_files, total_lines = split_fragment_file(fragment_file, num_processes, is_gzipped, debug)
    
    pool = multiprocessing.Pool(processes=num_processes)
    results = []
    for tmp_file_name in tmp_files:
        result = pool.apply_async(process_fragment_chunk, args=(tmp_file_name, barcode_group_dict, temp_dir, debug))
        results.append(result)

    pool.close()
    pool.join()

    group_files = {group: [] for group in set(barcode_group_dict.values())}
    for result in results:
        local_files = result.get()
        for group, file_name in local_files.items():
            group_files[group].append(file_name)

    for group, files in group_files.items():
        concatenate_files(files, output_temp_files[group])

    os.makedirs(output_dir, exist_ok=True)
    compressed_files = {}
    for group, temp_file in output_temp_files.items():
        compressed_file = compress_file(temp_file, output_dir)
        compressed_files[group] = compressed_file
        if debug:
            print(f"Debug: Compressed and saved {compressed_file}.")

    for tmp_file in tmp_files:
        os.remove(tmp_file)

    if is_gzipped:
        os.remove(fragment_file[:-3])

    shutil.rmtree(temp_dir)

    if debug:
        print(f"Debug: Finished processing. Check the output files in {output_dir}.")
        processed_lines = 0
        for group, file_name in compressed_files.items():
            with gzip.open(file_name, 'rt') as f:
                for line in f:
                    processed_lines += 1
        print(f"Debug: Processed {processed_lines} lines out of {total_lines}.")
        if processed_lines != total_lines:
            print("Warning: The number of processed lines does not match the total expected lines.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Sort fragment file based on barcode groupings.')
    parser.add_argument('--fragment_file', type=str, help='Path to the fragment file (gzipped or not)')
    parser.add_argument('--barcode_group_file', type=str, help='Path to the barcode grouping file')
    parser.add_argument('--output_dir', type=str, default='.', help='Path to the output directory')
    parser.add_argument('--num_processes', type=int, default=1, help='Number of processes to use')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode for additional output verification')
    args = parser.parse_args()

    main(args.fragment_file, args.barcode_group_file, args.output_dir, args.num_processes, args.debug)
