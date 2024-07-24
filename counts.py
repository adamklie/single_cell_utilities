import numpy as np

def sanitize_mtx(mtx: np.ndarray):
    cell_count_mask = mtx.sum(1) > 0  # count for each cell
    gene_count_mask = mtx.sum(0) > 0  # count for each gene

    genes_detected_mask = (mtx > 0).sum(1) > 0  # n genes per cell
    cells_detected_mask = (mtx > 0).sum(0) > 0  # n cells per gene
    row_mask = np.logical_and(cell_count_mask, genes_detected_mask)
    col_mask = np.logical_and(gene_count_mask, cells_detected_mask)

    return (row_mask, col_mask)


def is_counts(data):
    if np.all(data >= 0) and np.all(data.astype(int) == data):
        print("The matrix contains count data.")
    else:
        print("The matrix does not contain count data.")
        
        
def is_mostly_counts(data, percent=0.9):
    """Want to check if some percent of the data is counts

    Args:
        data (_type_): _description_
    """
    if np.all(data >= 0) and np.all(data.astype(int) == data):
        print("The matrix contains all count data.")
    elif np.sum(data >= 0) / data.size >= percent and np.sum(data.astype(int) == data) / data.size >= percent:
        greater_than_0 = (np.sum(data >= 0) / data.size)*100
        int_equals = (np.sum(data.astype(int) == data) / data.size)*100
        print(f"The matrix contains mostly count data. {greater_than_0}% of the data is greater than 0 and {int_equals}% of the data is equal to its integer value.")
    else:
        print("The matrix does not contain mostly count data.")