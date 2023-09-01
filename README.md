# single_cell_utilities
Work in progress libarary of utility functions and scripts for handling single cell data. Everything from acquiring it from databases to visualizing results of analysis.

# Configuring environments
Coming soon

# `data_acquisition`
Environment: `get_data`
Scripts and functions for acquiring data from databases and other sources. Relies heavily on the following packages
- `GEOparse`
- `pysradb`
- `parallel-fastq-dump`

# `data_wrangling`
Environments: `scverse-lite-py38`, `scverse-py38`, `scverse-R433`, `scenicplus`
Scripts and functions for wrangling data into a format that is ready for analysis. Relies heavily on the following packages orgainzed by data modality. Hope to get this to a less package dependent state in the future.

**scRNA-seq**
- `scanpy` -
- `Seurat` -

**scATAC-seq**
- `Signac` -
- `SnapATAC2` -
- `pycisTopic` -
- `ArchR` -

# `data_processing`
Environments: `scverse-lite-py38`, `scverse-py38`, `scverse-R433`, `scenicplus`
Scripts and functions to process the data prior to analysis. Functionality here is less package dependent and can be organized into the following categories
- cleaning
- quality control
- normalizations
- feature selection
- metacell aggregation
- imputation
- batch correction

# `data_analysis`
Environments: `scverse-lite-py38`, `scverse-py38`, `scverse-R433`, `scenicplus`
Scripts and functions to analyze the data that will often use the output of the `data_processing` functions.

# `data_visualization`
TODO