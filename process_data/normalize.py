import pandas as pd


def normalize_cpm(counts_matrix):
    """
    Normalize an RNA-seq counts matrix using Counts Per Million (CPM).
    
    Counts Per Million (CPM) is a widely used method in RNA-seq data analysis to normalize 
    gene expression levels. CPM aims to adjust for differences in sequencing depths across 
    samples and provide relative expression values on a comparable scale by scaling the raw 
    read counts of each gene by a sample-specific sequencing depth (total counts) and 
    multiplying by a scaling factor of one million (to obtain counts per million).

    Pros:
    - Relatively simple and intuitive to implement.
    - Allows for direct comparisons of gene expression levels between samples.
    
    Cons:
    - Does not account for gene length, potentially leading to biases in gene expression estimation.
    - May be sensitive to extreme values or outliers.
    
    Parameters:
    counts_matrix (pd.DataFrame): A pandas DataFrame containing raw read counts 
                                  with rows as genes and columns as samples.
    
    Returns:
    pd.DataFrame: A pandas DataFrame with CPM-normalized expression values.
    
    Raises:
    ValueError: If the input is not a pandas DataFrame.
    """
    
    if not isinstance(counts_matrix, pd.DataFrame):
        raise ValueError("Input must be a pandas DataFrame.")
    
    # Calculate total counts per sample
    total_counts = counts_matrix.sum(axis=0)
    
    # CPM normalization
    cpm_matrix = counts_matrix.div(total_counts, axis=1) * 1e6
    
    return cpm_matrix


def normalize_tpm(counts_matrix, gene_lengths):
    """
    Normalize an RNA-seq counts matrix using Transcripts Per Kilobase Million (TPM).
    
    Transcripts Per Kilobase Million (TPM) is an improvement over RPKM/FPKM. 
    TPM first normalizes for gene length, then for sequencing depth, making the sum 
    of all TPMs in each sample identical. This allows for a more accurate comparison 
    of gene expression between samples.
    
    Pros:
    - More accurate than RPKM/FPKM for comparing gene expression levels between samples.
    - Accounts for the total sum of normalized expression levels, allowing for a more balanced comparison.
    
    Cons:
    - Can be affected by highly expressed genes and depends on accurate estimates of gene lengths and accurate read mapping.
    - Cannot be used for differential expression analysis.
    
    Parameters:
    counts_matrix (pd.DataFrame): A pandas DataFrame containing raw read counts 
                                  with rows as genes and columns as samples.
    gene_lengths (pd.Series): A pandas Series containing the lengths of genes in kilobases.
                              The index should match the row index of the counts_matrix.
    
    Returns:
    pd.DataFrame: A pandas DataFrame with TPM-normalized expression values.
    
    Raises:
    ValueError: If the input counts_matrix is not a pandas DataFrame or if gene_lengths is not a pandas Series.
    """
    
    if not isinstance(counts_matrix, pd.DataFrame):
        raise ValueError("counts_matrix must be a pandas DataFrame.")
    
    if not isinstance(gene_lengths, pd.Series):
        raise ValueError("gene_lengths must be a pandas Series.")
    
    if not all(counts_matrix.index == gene_lengths.index):
        raise ValueError("The index of counts_matrix and gene_lengths must match.")
    
    # Normalize by gene length in kilobases
    length_normalized = counts_matrix.div(gene_lengths, axis=0)
    
    # Calculate per-sample scaling factors (total sum of length-normalized counts)
    scaling_factors = length_normalized.sum(axis=0)
    
    # Normalize by scaling factors and multiply by 1e6 to get TPM
    tpm_matrix = length_normalized.div(scaling_factors, axis=1) * 1e6
    
    return tpm_matrix


def normalize_rpkm(counts_matrix, gene_lengths):
    """
    Normalize an RNA-seq counts matrix using Reads Per Kilobase of transcript per Million mapped reads (RPKM).
    
    RPKM (Reads Per Kilobase of transcript per Million mapped reads) and FPKM (Fragments Per Kilobase 
    of transcript per Million mapped reads) normalize for both the length of the gene and the total 
    number of reads (i.e., the library size), making expression level comparisons between genes 
    in the same sample possible.
    
    Pros:
    - Widely used and established/incorporated in several software tools.
    - Makes it possible to compare gene expression levels within the same sample and between different samples.
    
    Cons:
    - Assumes that the total number of reads is the same across all samples, which isn't always accurate, 
      particularly when comparing different conditions or tissues.
    - Can be biased by highly expressed genes or transcripts, making it less reliable for datasets 
      with a high level of expression variation.
    - Like TPM, it cannot be used for differential expression analysis.
    
    Parameters:
    counts_matrix (pd.DataFrame): A pandas DataFrame containing raw read counts 
                                  with rows as genes and columns as samples.
    gene_lengths (pd.Series): A pandas Series containing the lengths of genes in kilobases.
                              The index should match the row index of the counts_matrix.
    
    Returns:
    pd.DataFrame: A pandas DataFrame with RPKM-normalized expression values.
    
    Raises:
    ValueError: If the input counts_matrix is not a pandas DataFrame or if gene_lengths is not a pandas Series.
    """
    
    if not isinstance(counts_matrix, pd.DataFrame):
        raise ValueError("counts_matrix must be a pandas DataFrame.")
    
    if not isinstance(gene_lengths, pd.Series):
        raise ValueError("gene_lengths must be a pandas Series.")
    
    if not all(counts_matrix.index == gene_lengths.index):
        raise ValueError("The index of counts_matrix and gene_lengths must match.")
    
    # Divide by both gene length and total counts, then multiply by 1e9 to get RPKM
    rpkm_matrix = counts_matrix.div(gene_lengths, axis=0).div(counts_matrix.sum(axis=0), axis=1) * 1e9
    
    return rpkm_matrix