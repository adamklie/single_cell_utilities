# pycisTopic Pipeline

End-to-end cisTopic analysis: build a cisTopic object from any peak matrix, train LDA topic models, evaluate them, and interactively explore a selected model.

## Quick Start

```bash
# Step 1: Build cisTopic object from a peak matrix
python build_cistopic_obj.py \
    --matrix_path counts.mtx \
    --barcodes_path barcodes.tsv \
    --regions_path regions.tsv \
    --output_dir ./results \
    --project_name my_sample \
    --blacklist_path blacklist.bed

# Step 2: Train LDA models across multiple topic counts
python run_models.py \
    --input ./results/my_sample.pkl \
    --output ./results/models \
    --n_topics 10 20 30 40 50 60 \
    --method cgs \
    --n_cpu 8

# Step 3: Evaluate models (metrics + UMAPs)
python evaluate_models.py \
    --cistopic_obj ./results/my_sample.pkl \
    --models ./results/models/models.pkl \
    --output_dir ./results/evaluation \
    --groupby cell_type

# Step 4: Open the notebook to explore your chosen model
#   -> Fill in the config cell in notebooks/model_analysis.ipynb
#   -> Set selected_n_topics to your chosen topic count
jupyter notebook notebooks/model_analysis.ipynb
```

## Input Formats

### MTX Mode (3 files)

| File | Format | Description |
|------|--------|-------------|
| `counts.mtx` | Matrix Market (`.mtx`) | Sparse peak-by-cell or cell-by-peak matrix (auto-detected) |
| `barcodes.tsv` | One barcode per line, no header | Cell barcodes (e.g. `ACGTACGT-1`) |
| `regions.tsv` | One region per line, no header | Peak regions (e.g. `chr1:1000-2000`) |

The script auto-detects matrix orientation by matching dimensions against the barcodes and regions counts.

### H5AD Mode (1 file)

An AnnData `.h5ad` file with:
- **obs** (rows) = cells
- **var** (columns) = regions
- `.X` = cells-by-regions count matrix (transposed internally to regions-by-cells)

## Script Reference

### `build_cistopic_obj.py`

Build a cisTopic object from a peak matrix.

| Argument | Type | Required | Default | Description |
|----------|------|----------|---------|-------------|
| `--matrix_path` | str | MTX mode | — | Path to `.mtx` counts matrix |
| `--barcodes_path` | str | MTX mode | — | Path to barcodes file |
| `--regions_path` | str | MTX mode | — | Path to regions file |
| `--h5ad_path` | str | H5AD mode | — | Path to `.h5ad` AnnData file |
| `--output_dir` | str | Yes | — | Output directory |
| `--project_name` | str | No | `cistopic` | Project name (used as output filename) |
| `--blacklist_path` | str | No | `None` | Path to blacklist BED file |
| `--cell_metadata_path` | str | No | `None` | Path to cell metadata (`.tsv`, `.txt`, or `.csv`) |
| `--region_separator` | str | No | `:-` | Expected region name separator pattern |

**Output:** `{output_dir}/{project_name}.pkl`

### `run_models.py`

Train LDA topic models (CGS or MALLET).

| Argument | Type | Required | Default | Description |
|----------|------|----------|---------|-------------|
| `--input` | str | Yes | — | Path to cisTopic `.pkl` |
| `--output` | str | Yes | — | Output directory |
| `--n_topics` | int+ | Yes | — | Topic counts (space-separated) |
| `--method` | str | Yes | — | `cgs` or `mallet` |
| `--n_cpu` | int | No | `1` | Number of CPU cores |
| `--n_iter` | int | No | `150` | Number of iterations |
| `--alpha` | float | No | `50.0` | Alpha hyperparameter |
| `--alpha_by_topic` | flag | No | `True` | Divide alpha by n_topics |
| `--no_alpha_by_topic` | flag | No | — | Don't divide alpha by n_topics |
| `--eta` | float | No | `0.1` | Eta hyperparameter |
| `--eta_by_topic` | flag | No | `False` | Divide eta by n_topics |
| `--seed` | int | No | `555` | Random seed |
| `--save_path` | str | No | `None` | Save intermediate models |
| `--temp_dir` | str | No | `None` | Temporary directory |
| `--mallet_path` | str | MALLET | — | Path to MALLET binary |
| `--mallet_memory` | str | No | `100G` | Java heap memory for MALLET |
| `--reuse_corpus` | flag | No | `False` | Reuse existing MALLET corpus |
| `--top_topics_coh` | int | No | `5` | Top topics for coherence |

**Outputs:**
- `{output}/models.pkl` — list of trained models
- `{output}/models_params.yaml` — reproducibility parameters

### `evaluate_models.py`

Evaluate trained models with metrics and optional UMAPs.

| Argument | Type | Required | Default | Description |
|----------|------|----------|---------|-------------|
| `--cistopic_obj` | str | Yes | — | Path to cisTopic `.pkl` |
| `--models` | str | Yes | — | Path to models `.pkl` or directory |
| `--output_dir` | str | Yes | — | Output directory |
| `--n_topics` | int+ | No | all | Evaluate only these topic counts |
| `--groupby` | str | No | `None` | Metadata column for UMAP coloring |
| `--metrics` | str+ | No | all four | Metrics to compute |
| `--run_umap` / `--no_run_umap` | flag | No | `True` | Generate UMAPs per topic count |

**Outputs:**
- `{output_dir}/model_evaluation_metrics.pdf` — metrics plot
- `{output_dir}/model_evaluation_metrics.csv` — metrics table
- `{output_dir}/umap_per_topic_count.pdf` — UMAP panel (if enabled)
- `{output_dir}/evaluation_params.yaml` — parameters

## Setup (uv)

This pipeline uses [uv](https://docs.astral.sh/uv/) for environment management.

```bash
cd pycistopic
uv sync          # creates .venv/ and installs all dependencies
uv run python build_cistopic_obj.py --help   # verify it works
```

To activate the venv manually (e.g. for Jupyter):

```bash
source .venv/bin/activate
```

## Local Testing

Run the full pipeline on synthetic test data to verify everything works:

```bash
cd pycistopic
uv run bash tests/test_local.sh
```

This will:
1. Generate small test data (100 regions x 50 cells) in `tests/test_data/`
2. Build cisTopic objects from both MTX and H5AD inputs
3. Train CGS models with topic counts 2, 3, 5 (10 iterations)
4. Evaluate the trained models
5. Optionally test MALLET if available (set `MALLET_PATH` env var)

All outputs go to `tests/test_output/`. To test MALLET mode:

```bash
# Download MALLET (requires Java)
mkdir -p tests/mallet
curl -L https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz \
  | tar xz -C tests/mallet --strip-components=1
export MALLET_PATH=tests/mallet/bin/mallet
uv run bash tests/test_local.sh
```

## SLURM Submission

Pre-built SLURM scripts are in `slurm/`. All paths and settings are centralized in `slurm/config.sh`.

### Setup

1. Edit `slurm/config.sh` with your paths, SLURM partition, account, etc.
2. Submit jobs in order:

```bash
# Step 1: Build cisTopic object
sbatch slurm/01_build_cistopic_obj.sh
# -> Update CISTOPIC_OBJ in config.sh with the output path

# Step 2: Train models (CGS or MALLET — pick one or both)
sbatch slurm/02_run_models_cgs.sh
sbatch slurm/03_run_models_mallet.sh
# -> Update MODELS_DIR in config.sh with the output path

# Step 3: Evaluate models
sbatch slurm/04_evaluate_models.sh
```

Each script sources `config.sh` automatically. Supports both conda and uv venv activation (controlled in config).

## Dependencies

- Python >= 3.11
- [pycisTopic](https://github.com/aertslab/pycisTopic) (installed from git)
- numpy, pandas, scipy
- scanpy, anndata
- matplotlib, seaborn
- pyyaml, setuptools
- For MALLET: [MALLET](https://mimno.github.io/Mallet/) binary + Java

## Directory Structure

```
pycistopic/
    pyproject.toml                 # uv project config
    .python-version                # Python 3.11 pin
    build_cistopic_obj.py          # CLI: peak matrix -> cisTopic object
    run_models.py                  # CLI: cisTopic object -> trained LDA models
    evaluate_models.py             # CLI: evaluate trained models
    notebooks/
        model_analysis.ipynb       # Interactive: explore a selected model
    tests/
        generate_test_data.py      # Create synthetic test data
        test_local.sh              # Run full pipeline on test data
    slurm/
        config.sh                  # Central config for all SLURM scripts
        01_build_cistopic_obj.sh   # Build cisTopic object
        02_run_models_cgs.sh       # Train with CGS
        03_run_models_mallet.sh    # Train with MALLET
        04_evaluate_models.sh      # Evaluate models
    README.md                      # This file
    legacy/                        # All old files preserved here
```

## Troubleshooting

### "Matrix shape does not match barcodes or regions"
Your MTX file dimensions don't match the barcode/region counts. Verify:
```bash
wc -l barcodes.tsv  # should match one matrix dimension
wc -l regions.tsv   # should match the other dimension
```

### "Could not load pickle file"
The `.pkl` file may have been saved with a different Python/pycisTopic version. Ensure you're using the same environment that created the file.

### "MALLET binary not found"
Set `--mallet_path` to the full path of your MALLET binary (e.g. `/opt/Mallet/bin/mallet`).

### "Column not found in cell metadata"
The `--groupby` column name must exactly match a column in `cistopic_obj.cell_data`. Run `evaluate_models.py` without `--groupby` first — it will list available columns in its error message.

### Barcode index mismatch between cisTopic and AnnData
pycisTopic appends `___project_name` to barcodes internally. The notebook handles this by splitting on `___`. If you have `___` in your actual barcodes, set a project name that doesn't conflict.

### Out of memory during model training
- Reduce `--n_cpu` (each core needs its own copy of the model)
- For MALLET, increase `--mallet_memory` (e.g. `200G`)
- Train fewer topic counts per job and combine models later
