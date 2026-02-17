#!/usr/bin/env bash
#SBATCH --job-name=cistopic_eval
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# 04_evaluate_models.sh — Evaluate trained LDA models with metrics and optional UMAPs.
#
# Usage:
#   sbatch slurm/04_evaluate_models.sh
#
# Requires CISTOPIC_OBJ and MODELS_DIR to be set in config.sh.
# Works with output from either 02 (CGS) or 03 (MALLET).

set -euo pipefail

# --- Source config ---
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CONFIG_DIR}/config.sh"

# --- Apply SLURM overrides from config ---
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --time=${SLURM_TIME_SHORT}
#SBATCH --mem=${SLURM_MEM_EVAL}
if [ -n "${SLURM_ACCOUNT}" ]; then
    #SBATCH --account=${SLURM_ACCOUNT}
    true
fi

# --- Validate inputs ---
if [ -z "${CISTOPIC_OBJ}" ]; then
    echo "ERROR: CISTOPIC_OBJ is not set in config.sh."
    exit 1
fi
if [ ! -f "${CISTOPIC_OBJ}" ]; then
    echo "ERROR: CISTOPIC_OBJ not found: ${CISTOPIC_OBJ}"
    exit 1
fi
if [ -z "${MODELS_DIR}" ]; then
    echo "ERROR: MODELS_DIR is not set in config.sh. Run step 02 or 03 first."
    exit 1
fi

# Find models.pkl — either directly or inside the directory
MODELS_PKL=""
if [ -f "${MODELS_DIR}/models.pkl" ]; then
    MODELS_PKL="${MODELS_DIR}/models.pkl"
elif [ -f "${MODELS_DIR}" ]; then
    MODELS_PKL="${MODELS_DIR}"
else
    echo "ERROR: No models.pkl found in MODELS_DIR: ${MODELS_DIR}"
    exit 1
fi

# --- Create directories ---
mkdir -p "${LOG_DIR}"
EVAL_OUTPUT="${OUTPUT_BASE}/evaluation"
mkdir -p "${EVAL_OUTPUT}"

# --- Activate environment ---
if [ -n "${UV_VENV:-}" ]; then
    echo "Activating uv venv: ${UV_VENV}"
    source "${UV_VENV}/bin/activate"
elif [ -n "${CONDA_ENV:-}" ]; then
    echo "Activating conda env: ${CONDA_ENV}"
    source activate "${CONDA_ENV}" 2>/dev/null || conda activate "${CONDA_ENV}"
fi

# --- Summary ---
echo "============================================"
echo " Evaluating LDA Models"
echo "============================================"
echo "  cisTopic obj : ${CISTOPIC_OBJ}"
echo "  Models       : ${MODELS_PKL}"
echo "  Output dir   : ${EVAL_OUTPUT}"
echo "  Group by     : ${GROUPBY:-none}"
echo "============================================"
echo ""

CMD="python ${SCRIPT_DIR}/evaluate_models.py"
CMD="${CMD} --cistopic_obj ${CISTOPIC_OBJ}"
CMD="${CMD} --models ${MODELS_PKL}"
CMD="${CMD} --output_dir ${EVAL_OUTPUT}"

if [ -n "${GROUPBY:-}" ]; then
    CMD="${CMD} --groupby ${GROUPBY}"
fi

echo "Running: ${CMD}"
echo ""

eval ${CMD}

echo ""
echo "Done. Results in: ${EVAL_OUTPUT}/"
