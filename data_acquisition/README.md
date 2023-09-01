# Downloading fastqs from SRA
I have found the greatest success (among numerous other attempts) at using pysradb to do the metadata downloading and the download of `.sra` files. 
I then use `parallel-fastq-dump` to convert the `.sra` files to `.fastq` files on a multi-core machine.
Each dataset is unique and usually requires some manual coding to make sure you are pulling the right metadata. I often do this in a Jupyter notebooks before updating the download script to reflect the changes I made to the metadata acquisition.

# Downloading processed data from GEO
I have found the greatest success (among numerous other attempts) at using GEOparse for downloading processed data from GEO. More details to come.

# Download validation


# File transferring
I have found the greatest success (among numerous other attempts) at using `rsync` to transfer files from one machine to another. More details to come. TODO: `/cellar/users/aklie/data/igvf/beta_cell_networks/download/igvf_sc-islet_10X-Multiome` make the code in this directory more generalizable and move it to `single_cell_utilities` repo.