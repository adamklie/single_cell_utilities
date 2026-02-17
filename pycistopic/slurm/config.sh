#!/usr/bin/env bash
# config.sh — Central configuration for all pycisTopic SLURM scripts.
#
# Edit the variables below, then run:
#   sbatch slurm/01_build_cistopic_obj.sh
#   sbatch slurm/02_run_models_cgs.sh
#   ...
#
# All SLURM scripts source this file automatically.

# === ENVIRONMENT ===
# Option A: conda/mamba environment
CONDA_ENV="/path/to/your/conda/env"
# Option B: uv virtual environment (uncomment and set to use instead of conda)
# UV_VENV="/path/to/pycistopic/.venv"

# === PATHS ===
# Resolved automatically — points to the pycistopic/ project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Set after step 1 (output of build_cistopic_obj.py)
CISTOPIC_OBJ=""

# Set after step 2 or 3 (output directory of run_models.py)
MODELS_DIR=""

# === INPUT DATA ===
# MTX mode (set all three)
MTX_PATH=""
BARCODES_PATH=""
REGIONS_PATH=""

# H5AD mode (uncomment and set; comment out MTX paths above)
# H5AD_PATH=""

# Optional inputs
CELL_METADATA_PATH=""
BLACKLIST_PATH=""

# === PROJECT ===
PROJECT_NAME="my_project"
OUTPUT_BASE="/path/to/output"

# === SLURM DEFAULTS ===
SLURM_PARTITION="compute"
SLURM_ACCOUNT=""
SLURM_TIME_SHORT="4:00:00"
SLURM_TIME_LONG="2-00:00:00"
SLURM_MEM_BUILD="32G"
SLURM_MEM_TRAIN="64G"
SLURM_MEM_EVAL="32G"
SLURM_CPUS_TRAIN=8
LOG_DIR="${OUTPUT_BASE}/logs"

# === MODEL TRAINING ===
N_TOPICS="10 20 30 40 50 60 80"
N_ITER=500
ALPHA=50.0
ETA=0.1
SEED=555

# === MALLET (only needed for 03_run_models_mallet.sh) ===
MALLET_PATH="/path/to/mallet"
MALLET_MEMORY="100G"

# === EVALUATION ===
GROUPBY="cell_type"
