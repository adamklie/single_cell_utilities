#!/usr/bin/env bash
# test_local.sh — Run the full pycisTopic pipeline on synthetic test data.
#
# Usage:
#   uv run bash tests/test_local.sh        # from the pycistopic/ directory
#   bash tests/test_local.sh               # if venv is already active
#
# All outputs go to tests/test_output/. Exit code 0 = all tests passed.

set -euo pipefail

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
TEST_DATA="${SCRIPT_DIR}/test_data"
TEST_OUTPUT="${SCRIPT_DIR}/test_output"

PASS=0
FAIL=0
SKIP=0

pass() { echo "  [PASS] $1"; PASS=$((PASS + 1)); }
fail() { echo "  [FAIL] $1"; FAIL=$((FAIL + 1)); }
skip() { echo "  [SKIP] $1"; SKIP=$((SKIP + 1)); }

echo "============================================"
echo " pycisTopic Local Test Suite"
echo "============================================"
echo "Project dir : ${PROJECT_DIR}"
echo "Test data   : ${TEST_DATA}"
echo "Test output : ${TEST_OUTPUT}"
echo ""

# Clean previous test output
rm -rf "${TEST_OUTPUT}"
mkdir -p "${TEST_OUTPUT}"

# ---------------------------------------------------------------------------
# Step 1: Generate test data
# ---------------------------------------------------------------------------
echo "--- Step 1: Generate test data ---"
if python "${SCRIPT_DIR}/generate_test_data.py" --output_dir "${TEST_DATA}"; then
    pass "generate_test_data.py"
else
    fail "generate_test_data.py"
    echo "Cannot continue without test data. Exiting."
    exit 1
fi
echo ""

# ---------------------------------------------------------------------------
# Step 2: Build cisTopic object (MTX mode)
# ---------------------------------------------------------------------------
echo "--- Step 2: Build cisTopic object (MTX mode) ---"
if python "${PROJECT_DIR}/build_cistopic_obj.py" \
    --matrix_path "${TEST_DATA}/counts.mtx" \
    --barcodes_path "${TEST_DATA}/barcodes.tsv" \
    --regions_path "${TEST_DATA}/regions.tsv" \
    --output_dir "${TEST_OUTPUT}/build_mtx" \
    --project_name test_mtx \
    --blacklist_path "${TEST_DATA}/blacklist.bed" \
    --cell_metadata_path "${TEST_DATA}/cell_metadata.tsv"; then
    if [ -f "${TEST_OUTPUT}/build_mtx/test_mtx.pkl" ]; then
        pass "build_cistopic_obj.py (MTX mode)"
    else
        fail "build_cistopic_obj.py (MTX mode) — output .pkl not found"
    fi
else
    fail "build_cistopic_obj.py (MTX mode)"
fi
echo ""

# ---------------------------------------------------------------------------
# Step 3: Build cisTopic object (H5AD mode)
# ---------------------------------------------------------------------------
echo "--- Step 3: Build cisTopic object (H5AD mode) ---"
if python "${PROJECT_DIR}/build_cistopic_obj.py" \
    --h5ad_path "${TEST_DATA}/peaks.h5ad" \
    --output_dir "${TEST_OUTPUT}/build_h5ad" \
    --project_name test_h5ad \
    --blacklist_path "${TEST_DATA}/blacklist.bed" \
    --cell_metadata_path "${TEST_DATA}/cell_metadata.tsv"; then
    if [ -f "${TEST_OUTPUT}/build_h5ad/test_h5ad.pkl" ]; then
        pass "build_cistopic_obj.py (H5AD mode)"
    else
        fail "build_cistopic_obj.py (H5AD mode) — output .pkl not found"
    fi
else
    fail "build_cistopic_obj.py (H5AD mode)"
fi
echo ""

# ---------------------------------------------------------------------------
# Step 4: Train models (CGS)
# ---------------------------------------------------------------------------
echo "--- Step 4: Train LDA models (CGS) ---"
CISTOPIC_OBJ="${TEST_OUTPUT}/build_mtx/test_mtx.pkl"
if [ ! -f "${CISTOPIC_OBJ}" ]; then
    skip "run_models.py (CGS) — no cisTopic object from step 2"
else
    if python "${PROJECT_DIR}/run_models.py" \
        --input "${CISTOPIC_OBJ}" \
        --output "${TEST_OUTPUT}/models_cgs" \
        --n_topics 2 3 5 \
        --method cgs \
        --n_iter 10 \
        --seed 42; then
        if [ -f "${TEST_OUTPUT}/models_cgs/models.pkl" ]; then
            pass "run_models.py (CGS)"
        else
            fail "run_models.py (CGS) — models.pkl not found"
        fi
    else
        fail "run_models.py (CGS)"
    fi
fi
echo ""

# ---------------------------------------------------------------------------
# Step 5: Evaluate models (CGS)
# ---------------------------------------------------------------------------
echo "--- Step 5: Evaluate models (CGS) ---"
if [ ! -f "${TEST_OUTPUT}/models_cgs/models.pkl" ]; then
    skip "evaluate_models.py (CGS) — no models from step 4"
else
    # Use only Arun_2010 + Cao_Juan_2009 + loglikelihood — Minmo_2011 (coherence)
    # produces NaN with very sparse test data and few topic counts.
    if python "${PROJECT_DIR}/evaluate_models.py" \
        --cistopic_obj "${CISTOPIC_OBJ}" \
        --models "${TEST_OUTPUT}/models_cgs/models.pkl" \
        --output_dir "${TEST_OUTPUT}/eval_cgs" \
        --groupby cell_type \
        --metrics Arun_2010 Cao_Juan_2009 loglikelihood; then
        pass "evaluate_models.py (CGS)"
    else
        fail "evaluate_models.py (CGS)"
    fi
fi
echo ""

# ---------------------------------------------------------------------------
# Step 6: Train + evaluate models (MALLET) — optional
# ---------------------------------------------------------------------------
echo "--- Step 6: Train LDA models (MALLET) ---"

# Find mallet binary
MALLET_BIN=""
if [ -n "${MALLET_PATH:-}" ]; then
    MALLET_BIN="${MALLET_PATH}"
elif command -v mallet &>/dev/null; then
    MALLET_BIN="$(command -v mallet)"
fi

if [ -n "${MALLET_BIN}" ]; then
    echo "  Using MALLET at: ${MALLET_BIN}"
    if [ ! -f "${CISTOPIC_OBJ}" ]; then
        skip "run_models.py (MALLET) — no cisTopic object"
    else
        if python "${PROJECT_DIR}/run_models.py" \
            --input "${CISTOPIC_OBJ}" \
            --output "${TEST_OUTPUT}/models_mallet" \
            --n_topics 2 3 \
            --method mallet \
            --mallet_path "${MALLET_BIN}" \
            --n_iter 10 \
            --seed 42; then
            if [ -f "${TEST_OUTPUT}/models_mallet/models.pkl" ]; then
                pass "run_models.py (MALLET)"
            else
                fail "run_models.py (MALLET) — models.pkl not found"
            fi
        else
            fail "run_models.py (MALLET)"
        fi

        # Evaluate MALLET models
        echo ""
        echo "--- Step 6b: Evaluate models (MALLET) ---"
        if [ -f "${TEST_OUTPUT}/models_mallet/models.pkl" ]; then
            if python "${PROJECT_DIR}/evaluate_models.py" \
                --cistopic_obj "${CISTOPIC_OBJ}" \
                --models "${TEST_OUTPUT}/models_mallet/models.pkl" \
                --output_dir "${TEST_OUTPUT}/eval_mallet" \
                --groupby cell_type \
                --metrics Arun_2010 Cao_Juan_2009 loglikelihood; then
                pass "evaluate_models.py (MALLET)"
            else
                fail "evaluate_models.py (MALLET)"
            fi
        else
            skip "evaluate_models.py (MALLET) — no models"
        fi
    fi
else
    echo "  MALLET not found (set MALLET_PATH or install mallet on PATH)."

    # Offer to download if Java is available
    if command -v java &>/dev/null; then
        echo "  Java is available. To download MALLET for testing:"
        echo "    mkdir -p ${SCRIPT_DIR}/mallet"
        echo "    curl -L https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz \\"
        echo "      | tar xz -C ${SCRIPT_DIR}/mallet --strip-components=1"
        echo "    export MALLET_PATH=${SCRIPT_DIR}/mallet/bin/mallet"
        echo "    bash ${BASH_SOURCE[0]}"
    else
        echo "  Java not found — install Java and MALLET to test MALLET mode."
    fi
    skip "run_models.py (MALLET) — MALLET not available"
    skip "evaluate_models.py (MALLET) — MALLET not available"
fi
echo ""

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "============================================"
echo " Test Summary"
echo "============================================"
echo "  PASS: ${PASS}"
echo "  FAIL: ${FAIL}"
echo "  SKIP: ${SKIP}"
echo ""

if [ "${FAIL}" -gt 0 ]; then
    echo "Some tests FAILED."
    exit 1
else
    echo "All tests passed!"
    exit 0
fi
