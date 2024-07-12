import sys
import os
import re
import subprocess
import signal
import pyranges as pr


def strip_ext_gz(f):
    return re.sub(r'\.gz$', '', str(f))


def run_shell_cmd(cmd):
    p = subprocess.Popen(
        ['/bin/bash', '-o', 'pipefail'],  # to catch error in pipe
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid)  # to make a new process with a new PGID
    pid = p.pid
    pgid = os.getpgid(pid)
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    err_str = 'PID={}, PGID={}, RC={}\nSTDERR={}\nSTDOUT={}'.format(
        pid, pgid, rc, stderr.strip(), stdout.strip())
    if rc:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            pass
        finally:
            raise Exception(err_str)
    return stdout.strip('\n')


def gunzip(f, suffix, out_dir):
    if not f.endswith('.gz'):
        raise Exception('Cannot gunzip a file without .gz extension.')
    gunzipped = os.path.join(out_dir,
                             os.path.basename(strip_ext_gz(f)))
    if suffix:
        gunzipped += '.{}'.format(suffix)
    # cmd = 'gzip -cd {} > {}'.format(f, gunzipped)
    cmd = 'zcat -f {} > {}'.format(f, gunzipped)
    run_shell_cmd(cmd)
    return gunzipped


def rm_f(files):
    if files:
        if type(files) == list:
            run_shell_cmd('rm -f {}'.format(' '.join(files)))
        else:
            run_shell_cmd('rm -f {}'.format(files))


def filter_peaks(
    peaks: pr.PyRanges,
    blacklist: pr.PyRanges,
    score_col: str = "Score",
    n_peaks: int = None,
    min_q_value: float = None,
) -> pr.PyRanges:
    """Filter peaks based on blacklist regions.

    Performs the following operations:
    1. Remove blacklist regions using PyRanges.overlap(invert=True)
    2. Keep only chr1-22, X, Y
    3. Remove duplicates based on Chromosome, Start, End, keeping row with the highest score_col value
    4. Optionally remove peaks with q_value > min_q_value
    5. Optionally take top n_peaks peaks

    Parameters
    ----------
    peaks : pr.PyRanges
        Peaks to filter as a PyRanges object.
    blacklist : pr.PyRanges
        Blacklist regions as a PyRanges object.
    score_col : str, optional
        Column to use for scoring peaks, by default "Score".
    n_peaks : int, optional
        Number of peaks to keep, by default None.
    
    Returns
    -------
    pr.PyRanges
        Filtered peaks as a PyRanges object.
    """
    # Remove blacklist regions
    peaks = peaks.overlap(blacklist, invert=True)

    # Keep only chr1-22, X, Y
    peaks = peaks[peaks.Chromosome.str.match("chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY")]

    # Remove duplicates based on Chromosome, Start, End, keeping row with the highest score_col value
    peaks_df = peaks.df.sort_values(score_col, ascending=False).drop_duplicates(["Chromosome", "Start", "End"])

    # Remove peaks with q_value > min_q_value. Need to convert min_q_value to 10^-min_q_value
    if min_q_value is not None:
        min_q_value = 10 ** -min_q_value
        peaks_df = peaks_df[peaks_df["q_value"] >= min_q_value]
    
    # Take top n_peaks peaks if specified
    if n_peaks is not None:
        peaks_df = peaks_df.head(n_peaks)


    # Convert back to pyranges
    peaks = pr.PyRanges(peaks_df)

    return peaks


def naive_overlap(
    peaks1: pr.PyRanges,
    peaks2: pr.PyRanges,
    frac: float = None,
    return_overlap: bool = False
) -> tuple:
    """Calculate the overlap between two peaksets.

    Calculates the overlap between two peaksets and returns 3 statistics: 
    1) the number of peaks in the overlap, defined as a any peak in peaks1 that overlaps with a peak in peaks2
    2) the percent of peaks in peaks1 that overlap with peaks2, defined as the number of peaks in the overlap divided by the number of peaks in peaks1
    3) the fraction of the peaks in peaks1 that overlap with peaks2 that are intersecting. Defined as the length of the intersecting part of the peaks 
       in peaks1 that overlap with peaks2, divided by the length of the peaks in peaks1 that overlap with peaks2

    If frac is specified, only peaks in peaks1 that overlap with peaks2 and have a fraction overlap greater than frac will be considered.
    

    Parameters
    ----------
    peaks1 : pr.PyRanges
        First peakset as a PyRanges object.
    peaks2 : pr.PyRanges
        Second peakset as a PyRanges object.

    Returns
    -------
    tuple
        Tuple of overlap number and overlap percent.
    """

    # Get the number of peaks in each peakset
    num_peaks = len(peaks1)

    overlap = peaks1.overlap(peaks2, how="first")  # Peaks in peaks1 that overlap with peaks2
    intersect = peaks1.intersect(peaks2, how="first")  # Part of peaks in peaks1 that overlap with peaks2
    overlap_lengths = overlap.lengths()  # Lengths of overlapping peaks in peaks1
    intersect_lengths = intersect.lengths()  # Lengths of intersecting peaks in peaks1
    frac_overlap = (intersect_lengths/overlap_lengths).values  # Fraction of overlapping peaks in peaks1 that are intersecting peaks

    # Filter by fraction overlap if applicable
    if frac:
        frac_mask = (frac_overlap >= frac)
        overlap = overlap[frac_mask]
        frac_overlap = frac_overlap[frac_mask]


    # Get the number of peaks in the overlap
    overlap_num = len(overlap)

    # Get the percent overlap
    overlap_percent = overlap_num / num_peaks
    
    # Return
    if return_overlap:
        return overlap_num, overlap_percent, frac_overlap, overlap
    else:
        return overlap_num, overlap_percent, frac_overlap
    