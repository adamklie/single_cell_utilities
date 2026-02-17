import os
import sys
import glob
import logging
import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
_stderr = sys.stderr
null = open(os.devnull,'wb')


def main(args):
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    
    # Log arguments
    logging.info("Arguments: %s", args)

    # Dirs and paths
    logging.info("Setting up directories and paths")
    bed_dir = args.bed_dir
    bed_paths = sorted(glob.glob(os.path.join(bed_dir, '*.bed'))) + sorted(glob.glob(os.path.join(bed_dir, '*.bed.gz')))
    cellids = [os.path.basename(x).split('.')[0] for x in bed_paths]
    bed_path_dict = dict(zip(cellids, bed_paths))
    logging.info(f"Found {len(bed_paths)} bed files in {bed_dir}")
    logging.info(f"Dictionary of cellids and bed paths: {bed_path_dict}")
    tmp_dir = args.tmp_dir # default '/cellar/users/aklie/tmp/'
    output_dir = args.outdir_path
    if not os.path.exists(output_dir):
        logging.info(f"Creating output directory {output_dir}")
        os.makedirs(output_dir)

    # Get chromosome sizes
    import requests
    import pyranges as pr
    import pandas as pd
    logging.info(f"Fetching chromosome sizes from {args.chrom_sizes}")
    target_url = args.chrom_sizes  # default: 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
    chromsizes = pd.read_csv(target_url, sep='\t', header=None)
    chromsizes.columns = ['Chromosome', 'End']
    chromsizes['Start'] = [0]*chromsizes.shape[0]
    chromsizes = chromsizes.loc[:,['Chromosome', 'Start', 'End']]
    if args.cellranger_arc_annot:  # default is true
        # Exceptionally in this case, to agree with CellRangerARC annotations
        chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('v', '.') for x in range(len(chromsizes['Chromosome']))]
        chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].split('_')[1] if len(chromsizes['Chromosome'][x].split('_')) > 1 else chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))]
    chromsizes = pr.PyRanges(chromsizes)

    # Run peak calling
    from pycisTopic.pseudobulk_peak_calling import peak_calling
    macs_path = args.macs_path  # default: '/cellar/users/aklie/opt/miniconda3/envs/scenicplus/bin/macs2'
    logging.info(f"Running peak calling with macs_path={macs_path}")
    narrow_peaks_dict = peak_calling(
        macs_path=macs_path,
        bed_paths=bed_path_dict,
        outdir=output_dir,
        genome_size=args.genome_size,  # default: 'hs'
        n_cpu=args.num_cpus,  # default: 1
        input_format=args.input_format,  # default: 'BEDPE'
        shift=args.shift,  # default: 73
        ext_size=args.ext_size,  # default: 146
        keep_dup = args.keep_dup,  # default: 'all'
        q_value = args.q_value,  # default: 0.05
        _temp_dir = os.path.join(tmp_dir, 'ray_spill')
    )

    # Dump the return object to a pickle
    import pickle
    logging.info("Dumping narrow_peaks_dict to pickle")
    pickle.dump(
        narrow_peaks_dict,
        open(os.path.join(output_dir, 'narrow_peaks_dict.pkl'), 'wb')
    )

    # Get consensus peaks
    logging.info("Getting consensus peaks")
    from pycisTopic.iterative_peak_calling import get_consensus_peaks
    peak_half_width = args.peak_half_width  # default: 250
    path_to_blacklist= args.blacklist  # default: '/cellar/users/aklie/data/ref/blacklists/hg38/ENCFF356LFX.bed'
    consensus_peaks=get_consensus_peaks(
        narrow_peaks_dict=narrow_peaks_dict, 
        peak_half_width=peak_half_width, 
        chromsizes=chromsizes, 
        path_to_blacklist=path_to_blacklist
    )

    # Save to bed file
    logging.info(f"Saving consensus peaks to bed file {os.path.join(output_dir, 'consensus_regions.bed')}")
    consensus_peaks.to_bed(
        path = os.path.join(output_dir, 'consensus_regions.bed'),
        keep=True,
        compression='infer',
        chain=False
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--bed_dir', type=str, required=True, help='Directory containing fragment files in bed format.')
    parser.add_argument('--outdir_path', type=str, required=True, help='Path to output directory.')
    parser.add_argument('--tmp_dir', type=str, required=False, default='/cellar/users/aklie/tmp', help='Path to temporary directory.')
    parser.add_argument('--chrom_sizes', type=str, required=False, default='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes', help='Path to chromosome sizes file.')
    parser.add_argument('--macs_path', type=str, required=False, default='/cellar/users/aklie/opt/miniconda3/envs/scenicplus/bin/macs2', help='Path to macs2 executable.')
    parser.add_argument('--genome_size', type=str, required=False, default='hs', help='Genome size for macs2.')
    parser.add_argument('--num_cpus', type=int, required=False, default=1, help='Number of CPUs to use.')
    parser.add_argument('--input_format', type=str, required=False, default='BEDPE', help='Input format for macs2.')
    parser.add_argument('--shift', type=int, required=False, default=73, help='Shift for macs2.')
    parser.add_argument('--ext_size', type=int, required=False, default=146, help='Extension size for macs2.')
    parser.add_argument('--keep_dup', type=str, required=False, default='all', help='Keep duplicates for macs2.')
    parser.add_argument('--q_value', type=float, required=False, default=0.05, help='Q value for macs2.')
    parser.add_argument('--peak_half_width', type=int, required=False, default=250, help='Half width of peaks.')
    parser.add_argument('--blacklist', type=str, required=False, default='/cellar/users/aklie/data/ref/blacklists/hg38/ENCFF356LFX.bed', help='Path to blacklist file.')
    parser.add_argument('--cellranger_arc_annot', action='store_true', help='Whether to use CellRangerARC annotations for chromosomes.')
    args = parser.parse_args()
    main(args)