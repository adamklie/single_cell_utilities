#!/usr/bin/env bash
#SBATCH --job-name=cistopic_build
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# 01_build_cistopic_obj.sh — Build a cisTopic object from a peak matrix.
#
# Usage:
#   sbatch slurm/01_build_cistopic_obj.sh
#
# After completion, set CISTOPIC_OBJ in config.sh to the output path printed
# at the end of this script.

set -euo pipefail

# --- Source config ---
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CONFIG_DIR}/config.sh"

# --- Apply SLURM overrides from config ---
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --time=${SLURM_TIME_SHORT}
#SBATCH --mem=${SLURM_MEM_BUILD}
if [ -n "${SLURM_ACCOUNT}" ]; then
    #SBATCH --account=${SLURM_ACCOUNT}
    true
fi

# --- Create log and output directories ---
mkdir -p "${LOG_DIR}"
mkdir -p "${OUTPUT_BASE}/build"

# --- Activate environment ---
if [ -n "${UV_VENV:-}" ]; then
    echo "Activating uv venv: ${UV_VENV}"
    source "${UV_VENV}/bin/activate"
elif [ -n "${CONDA_ENV:-}" ]; then
    echo "Activating conda env: ${CONDA_ENV}"
    source activate "${CONDA_ENV}" 2>/dev/null || conda activate "${CONDA_ENV}"
fi

# --- Build command ---
OUTPUT_DIR="${OUTPUT_BASE}/build"
OUTPUT_PKL="${OUTPUT_DIR}/${PROJECT_NAME}.pkl"

echo "============================================"
echo " Building cisTopic Object"
echo "============================================"
echo "  Project    : ${PROJECT_NAME}"
echo "  Output     : ${OUTPUT_PKL}"

CMD="python ${SCRIPT_DIR}/build_cistopic_obj.py"
CMD="${CMD} --output_dir ${OUTPUT_DIR}"
CMD="${CMD} --project_name ${PROJECT_NAME}"

# MTX mode vs H5AD mode
if [ -n "${H5AD_PATH:-}" ]; then
    echo "  Mode       : H5AD"
    echo "  H5AD       : ${H5AD_PATH}"
    CMD="${CMD} --h5ad_path ${H5AD_PATH}"
else
    echo "  Mode       : MTX"
    echo "  Matrix     : ${MTX_PATH}"
    echo "  Barcodes   : ${BARCODES_PATH}"
    echo "  Regions    : ${REGIONS_PATH}"
    CMD="${CMD} --matrix_path ${MTX_PATH}"
    CMD="${CMD} --barcodes_path ${BARCODES_PATH}"
    CMD="${CMD} --regions_path ${REGIONS_PATH}"
fi

# Optional arguments
if [ -n "${BLACKLIST_PATH:-}" ]; then
    echo "  Blacklist  : ${BLACKLIST_PATH}"
    CMD="${CMD} --blacklist_path ${BLACKLIST_PATH}"
fi

if [ -n "${CELL_METADATA_PATH:-}" ]; then
    echo "  Metadata   : ${CELL_METADATA_PATH}"
    CMD="${CMD} --cell_metadata_path ${CELL_METADATA_PATH}"
fi

echo "============================================"
echo ""
echo "Running: ${CMD}"
echo ""

eval ${CMD}

echo ""
echo "Done. Set CISTOPIC_OBJ in config.sh to:"
echo "  CISTOPIC_OBJ=\"${OUTPUT_PKL}\""
