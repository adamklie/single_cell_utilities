#!/usr/bin/env bash
#SBATCH --job-name=cistopic_mallet
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1

# 03_run_models_mallet.sh — Train LDA topic models using MALLET.
#
# Usage:
#   sbatch slurm/03_run_models_mallet.sh
#
# Requires CISTOPIC_OBJ and MALLET_PATH to be set in config.sh.

set -euo pipefail

# --- Source config ---
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CONFIG_DIR}/config.sh"

# --- Apply SLURM overrides from config ---
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --time=${SLURM_TIME_LONG}
#SBATCH --mem=${SLURM_MEM_TRAIN}
#SBATCH --cpus-per-task=${SLURM_CPUS_TRAIN}
if [ -n "${SLURM_ACCOUNT}" ]; then
    #SBATCH --account=${SLURM_ACCOUNT}
    true
fi

# --- Validate inputs ---
if [ -z "${CISTOPIC_OBJ}" ]; then
    echo "ERROR: CISTOPIC_OBJ is not set in config.sh. Run step 01 first."
    exit 1
fi
if [ ! -f "${CISTOPIC_OBJ}" ]; then
    echo "ERROR: CISTOPIC_OBJ not found: ${CISTOPIC_OBJ}"
    exit 1
fi
if [ -z "${MALLET_PATH}" ] || [ ! -x "${MALLET_PATH}" ]; then
    echo "ERROR: MALLET_PATH is not set or not executable: ${MALLET_PATH}"
    exit 1
fi

# --- Create directories ---
mkdir -p "${LOG_DIR}"
MODELS_OUTPUT="${OUTPUT_BASE}/models_mallet"
mkdir -p "${MODELS_OUTPUT}"

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
echo " Training LDA Models (MALLET)"
echo "============================================"
echo "  cisTopic obj  : ${CISTOPIC_OBJ}"
echo "  Output dir    : ${MODELS_OUTPUT}"
echo "  Topics        : ${N_TOPICS}"
echo "  Iterations    : ${N_ITER}"
echo "  Alpha         : ${ALPHA}"
echo "  Eta           : ${ETA}"
echo "  Seed          : ${SEED}"
echo "  CPUs          : ${SLURM_CPUS_TRAIN}"
echo "  MALLET binary : ${MALLET_PATH}"
echo "  MALLET memory : ${MALLET_MEMORY}"
echo "============================================"
echo ""

python "${SCRIPT_DIR}/run_models.py" \
    --input "${CISTOPIC_OBJ}" \
    --output "${MODELS_OUTPUT}" \
    --n_topics ${N_TOPICS} \
    --method mallet \
    --mallet_path "${MALLET_PATH}" \
    --mallet_memory "${MALLET_MEMORY}" \
    --n_cpu "${SLURM_CPUS_TRAIN}" \
    --n_iter "${N_ITER}" \
    --alpha "${ALPHA}" \
    --eta "${ETA}" \
    --seed "${SEED}"

echo ""
echo "Done. Set MODELS_DIR in config.sh to:"
echo "  MODELS_DIR=\"${MODELS_OUTPUT}\""
