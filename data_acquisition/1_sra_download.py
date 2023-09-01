import os
import sys
import glob
import subprocess
from tqdm.autonotebook import tqdm
from pysradb.sraweb import SRAweb

# Choose the current dataset we are working with
dataset_name = "Wang2023_islet_snATAC-seq"
srp_id = "SRP311849"

# Set-up directories
base_dir = "/cellar/users/aklie/data/igvf/beta_cell_networks"
cwd = os.path.join(base_dir, "download", dataset_name)
fastq_dir = os.path.join(base_dir, "fastq", dataset_name)

# Connect to SRA
db = SRAweb()

# Grab the metadata for the SRP
metadata = db.sra_metadata(srp_id, detailed=True)

# Save the metadata and the list of srr ids
metadata.to_csv(os.path.join(cwd, f"{srp_id}_metadata.tsv"), index=False, sep="\t")
metadata["run_accession"].to_csv(os.path.join(cwd, f"{srp_id}_srr_ids.txt"), index=False, header=False)

# Subset to conditions of interest
t2d_metadata = metadata[metadata["diagnosis"] == "T2D"]
nd_metadata = metadata[metadata["diagnosis"] == "Non-diabetic"]
pre_t2d_metadata = metadata[metadata["diagnosis"] == "Pre-T2D"]

# Save each subset along with the list of srr ids
t2d_metadata.to_csv(os.path.join(cwd, f"{srp_id}_t2d_metadata.tsv"), index=False, sep="\t")
t2d_metadata["run_accession"].to_csv(os.path.join(cwd, f"{srp_id}_t2d_srr_ids.txt"), index=False, header=False)

nd_metadata.to_csv(os.path.join(cwd, f"{srp_id}_nd_metadata.tsv"), index=False, sep="\t")
nd_metadata["run_accession"].to_csv(os.path.join(cwd, f"{srp_id}_nd_srr_ids.txt"), index=False, header=False)

pre_t2d_metadata.to_csv(os.path.join(cwd, f"{srp_id}_pre_t2d_metadata.tsv"), index=False, sep="\t")
pre_t2d_metadata["run_accession"].to_csv(os.path.join(cwd, f"{srp_id}_pre_t2d_srr_ids.txt"), index=False, header=False)


# # Download non-diabetic samples (`sra` files) to start
db.download(df=nd_metadata, out_dir=fastq_dir, skip_confirmation=True)

# Parameters for download
tmp_dir = "/cellar/users/aklie/tmp/fastq-dump"
gzip = True
split_files = True
threads = 32

# Loop through and fastqdump out each SRA download file within the subdirectories of the fastq_dir
for sra_file in glob.glob(os.path.join(fastq_dir, srp_id, "*", "*.sra")):
    sra_dir = os.path.dirname(sra_file)
    if gzip:
        cmd = f"parallel-fastq-dump --threads {threads} --outdir {sra_dir} --split-files --tmpdir {tmp_dir} --gzip -s {sra_file}"
    else:
        cmd = f"parallel-fastq-dump --threads {threads} --outdir {sra_dir} --split-files --tmpdir {tmp_dir} -s {sra_file}"
    print(cmd)
    if len(glob.glob(os.path.join(sra_dir, "*.fastq*"))) > 0:
        print(f"Files already downloaded for {sra_dir}")
    else:
        subprocess.run(cmd, shell=True)
