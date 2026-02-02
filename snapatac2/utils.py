#!/usr/bin/env python3
"""
utils.py

Utility functions for SnapATAC2 workflows including peak filtering and overlap analysis.

Author: Adam Klie
"""

import pyranges as pr
from typing import Optional, Tuple, Union


def filter_peaks(
    peaks: pr.PyRanges,
    blacklist: pr.PyRanges,
    score_col: str = "Score",
) -> pr.PyRanges:
    """Filter peaks based on blacklist regions.

    Performs the following operations:
    1. Remove blacklist regions using PyRanges.overlap(invert=True)
    2. Keep only chr1-22, X, Y
    3. Remove duplicates based on Chromosome, Start, End, keeping row with the highest score_col value

    Parameters
    ----------
    peaks : pr.PyRanges
        Peaks to filter as a PyRanges object.
    blacklist : pr.PyRanges
        Blacklist regions as a PyRanges object.
    score_col : str, optional
        Column name containing peak scores for deduplication. Default is "Score".

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
    peaks_df = peaks.df.sort_values(score_col, ascending=False).drop_duplicates(
        ["Chromosome", "Start", "End"]
    )

    # Convert back to pyranges
    peaks = pr.PyRanges(peaks_df)

    return peaks


def naive_overlap(
    peaks1: pr.PyRanges,
    peaks2: pr.PyRanges,
    frac: Optional[float] = None,
    return_overlap: bool = False,
) -> Union[Tuple[int, float, float], Tuple[int, float, float, pr.PyRanges]]:
    """Calculate the overlap between two peaksets.

    Calculates the overlap between two peaksets and returns 3 statistics:
    1) the number of peaks in the overlap, defined as any peak in peaks1 that overlaps with a peak in peaks2
    2) the percent of peaks in peaks1 that overlap with peaks2, defined as the number of peaks in the
       overlap divided by the number of peaks in peaks1
    3) the fraction of the peaks in peaks1 that overlap with peaks2 that are intersecting. Defined as
       the length of the intersecting part of the peaks in peaks1 that overlap with peaks2, divided by
       the length of the peaks in peaks1 that overlap with peaks2

    If frac is specified, only peaks in peaks1 that overlap with peaks2 and have a fraction overlap
    greater than frac will be considered.

    Parameters
    ----------
    peaks1 : pr.PyRanges
        First peakset as a PyRanges object.
    peaks2 : pr.PyRanges
        Second peakset as a PyRanges object.
    frac : float, optional
        Minimum fraction overlap threshold. If specified, only peaks with overlap >= frac are counted.
    return_overlap : bool, optional
        If True, also return the overlapping peaks PyRanges object. Default is False.

    Returns
    -------
    tuple
        If return_overlap is False: (overlap_num, overlap_percent, frac_overlap)
        If return_overlap is True: (overlap_num, overlap_percent, frac_overlap, overlap_pyranges)
    """
    # Get the number of peaks in each peakset
    num_peaks = len(peaks1)

    overlap = peaks1.overlap(peaks2, how="first")  # Peaks in peaks1 that overlap with peaks2
    intersect = peaks1.intersect(peaks2, how="first")  # Part of peaks in peaks1 that overlap with peaks2
    overlap_lengths = overlap.lengths()  # Lengths of overlapping peaks in peaks1
    intersect_lengths = intersect.lengths()  # Lengths of intersecting peaks in peaks1
    frac_overlap = (intersect_lengths / overlap_lengths).values  # Fraction of overlapping peaks that are intersecting

    # Filter by fraction overlap if applicable
    if frac:
        frac_mask = frac_overlap >= frac
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
