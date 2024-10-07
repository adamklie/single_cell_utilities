import pandas as pd
from pybiomart import Server

def fetch_gene_info(assembly='hg38'):
    """
    Fetches gene information including gene name, start coordinate, 
    end coordinate, and gene length using the pybiomart package.

    Args:
    - assembly (str): Either 'hg19' or 'hg38' to specify the genome assembly (default is 'hg38').

    Returns:
    - pd.DataFrame: DataFrame containing the gene name, start coordinate, end coordinate, and gene length.
    """

    # Connect to the Ensembl BioMart server
    server = Server(host='http://www.ensembl.org')
    
    # Choose dataset based on assembly
    if assembly == 'hg19':
        dataset_name = 'hsapiens_gene_ensembl_grch37'
    elif assembly == 'hg38':
        dataset_name = 'hsapiens_gene_ensembl'
    else:
        raise ValueError("Invalid assembly version. Choose 'hg19' or 'hg38'.")
    
    # Fetch the dataset
    dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets[dataset_name]

    # Query BioMart for gene coordinates and gene names
    results = dataset.query(attributes=['external_gene_name', 'chromosome_name', 'start_position', 'end_position'])

    # Rename columns for clarity
    results.columns = ['gene_name', 'chromosome', 'start', 'end']
    
    # Calculate gene length
    results['length'] = results['end'] - results['start'] + 1

    # Return the final dataframe
    return results[['gene_name', 'chromosome', 'start', 'end', 'length']]