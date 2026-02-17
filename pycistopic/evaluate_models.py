#!/usr/bin/env python
"""Evaluate trained cisTopic LDA models: compute metrics and optionally generate UMAPs.

Outputs (all in --output_dir):
    model_evaluation_metrics.pdf   — multi-panel metrics plot
    model_evaluation_metrics.csv   — metrics values as a table
    umap_per_topic_count.pdf       — side-by-side UMAP panels (if --run_umap)
    evaluation_params.yaml         — parameters used
"""

import argparse
import glob
import logging
import os
import pickle
import sys

import yaml

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

ALL_METRICS = ["Arun_2010", "Cao_Juan_2009", "Minmo_2011", "loglikelihood"]


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------

def load_pickle(path, label):
    """Load a pickle file with error handling."""
    if not os.path.exists(path):
        logger.error("%s not found: %s", label, path)
        sys.exit(1)
    try:
        with open(path, "rb") as f:
            obj = pickle.load(f)
    except Exception as e:
        logger.error("Could not load %s: %s\n  %s", label, path, e)
        sys.exit(1)
    return obj


def load_models(models_path):
    """Load models from a single .pkl or a directory of .pkl files."""
    if os.path.isfile(models_path):
        logger.info("Loading models from file: %s", models_path)
        models = load_pickle(models_path, "models")
        if not isinstance(models, list):
            models = [models]
    elif os.path.isdir(models_path):
        logger.info("Loading models from directory: %s", models_path)
        pkl_files = sorted(glob.glob(os.path.join(models_path, "*.pkl")))
        if not pkl_files:
            logger.error("No .pkl files found in %s", models_path)
            sys.exit(1)
        models = []
        for pf in pkl_files:
            m = load_pickle(pf, f"model file {pf}")
            if isinstance(m, list):
                models.extend(m)
            else:
                models.append(m)
    else:
        logger.error("Models path not found: %s", models_path)
        sys.exit(1)

    logger.info("  Loaded %d model(s)", len(models))
    return models


def get_topic_counts(models):
    """Extract topic count from each model."""
    counts = []
    for m in models:
        counts.append(m.cell_topic.shape[0])
    return counts


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(args):
    # --- Load cisTopic object ---
    logger.info("Loading cisTopic object from %s", args.cistopic_obj)
    cistopic_obj = load_pickle(args.cistopic_obj, "cisTopic object")

    n_cells = cistopic_obj.fragment_matrix.shape[1]
    n_regions = cistopic_obj.fragment_matrix.shape[0]

    # --- Load models ---
    models = load_models(args.models)
    all_topic_counts = get_topic_counts(models)
    logger.info("  Available topic counts: %s", sorted(all_topic_counts))

    # --- Filter models by --n_topics if specified ---
    if args.n_topics is not None:
        requested = set(args.n_topics)
        available = set(all_topic_counts)
        missing = requested - available
        if missing:
            logger.error(
                "Requested topic counts not found in models: %s. "
                "Available: %s",
                sorted(missing),
                sorted(available),
            )
            sys.exit(1)
        models = [m for m in models if m.cell_topic.shape[0] in requested]
        all_topic_counts = get_topic_counts(models)
        logger.info("  Filtered to topic counts: %s", sorted(all_topic_counts))

    # --- Validate groupby ---
    if args.groupby is not None:
        cell_data = cistopic_obj.cell_data
        available_cols = list(cell_data.columns)
        if args.groupby not in available_cols:
            logger.error(
                "Column '%s' not found in cell metadata. "
                "Available columns: %s",
                args.groupby,
                available_cols,
            )
            sys.exit(1)

    # --- Validate metrics ---
    metrics = args.metrics
    for m in metrics:
        if m not in ALL_METRICS:
            logger.error(
                "Unknown metric '%s'. Available: %s", m, ALL_METRICS
            )
            sys.exit(1)

    # --- Summary ---
    logger.info("=" * 60)
    logger.info("Evaluation summary")
    logger.info("  n_models:     %d", len(models))
    logger.info("  Topic counts: %s", sorted(all_topic_counts))
    logger.info("  n_cells:      %d", n_cells)
    logger.info("  n_regions:    %d", n_regions)
    logger.info("  Metrics:      %s", metrics)
    logger.info("  Run UMAP:     %s", args.run_umap)
    if args.groupby:
        logger.info("  Group by:     %s", args.groupby)
    logger.info("  Output dir:   %s", args.output_dir)
    logger.info("=" * 60)

    # --- Prepare output ---
    os.makedirs(args.output_dir, exist_ok=True)

    # Save params
    params = {
        "cistopic_obj": os.path.abspath(args.cistopic_obj),
        "models": os.path.abspath(args.models),
        "output_dir": os.path.abspath(args.output_dir),
        "n_topics": args.n_topics,
        "metrics": metrics,
        "run_umap": args.run_umap,
        "groupby": args.groupby,
    }
    params_path = os.path.join(args.output_dir, "evaluation_params.yaml")
    with open(params_path, "w") as f:
        yaml.dump(params, f, default_flow_style=False, sort_keys=False)

    # --- Compute metrics ---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    from pycisTopic.lda_models import evaluate_models

    logger.info("Computing metrics...")

    # evaluate_models with plot=True generates the metrics plot
    metrics_fig_path = os.path.join(args.output_dir, "model_evaluation_metrics.pdf")
    model = evaluate_models(
        models,
        select_model=None,
        return_model=False,
        metrics=metrics,
        plot_metrics=True,
        plot=False,
    )
    plt.savefig(metrics_fig_path, bbox_inches="tight", dpi=150)
    plt.close()
    logger.info("  Metrics plot saved to %s", metrics_fig_path)

    # Build metrics CSV
    metrics_data = []
    for m in models:
        n_top = m.cell_topic.shape[0]
        row = {"n_topics": n_top}
        ll = getattr(m, "log_likelihood", None)
        if ll is not None:
            row["log_likelihood"] = ll
        metrics_data.append(row)
    metrics_df = pd.DataFrame(metrics_data).sort_values("n_topics")
    metrics_csv_path = os.path.join(args.output_dir, "model_evaluation_metrics.csv")
    metrics_df.to_csv(metrics_csv_path, index=False)
    logger.info("  Metrics CSV saved to %s", metrics_csv_path)

    # --- Generate UMAPs per topic count ---
    if args.run_umap:
        logger.info("Generating UMAPs per topic count...")
        from pycisTopic.clust_vis import run_umap, plot_metadata

        n_models = len(models)
        fig, axes = plt.subplots(
            1, n_models, figsize=(5 * n_models, 5), squeeze=False
        )

        for idx, m in enumerate(sorted(models, key=lambda x: x.cell_topic.shape[0])):
            n_top = m.cell_topic.shape[0]
            logger.info("  UMAP for n_topics=%d...", n_top)

            # Add model, run UMAP
            cistopic_obj.add_LDA_model(m)
            run_umap(cistopic_obj, target="cell", scale=False)

            umap_df = cistopic_obj.projections["cell"]["UMAP"]
            ax = axes[0, idx]
            ax.scatter(
                umap_df.iloc[:, 0],
                umap_df.iloc[:, 1],
                s=3,
                alpha=0.5,
                c="steelblue",
                rasterized=True,
            )
            ax.set_title(f"{n_top} topics")
            ax.set_xlabel("UMAP1")
            ax.set_ylabel("UMAP2")
            ax.set_xticks([])
            ax.set_yticks([])

            # If groupby is set, color by that variable
            if args.groupby is not None:
                cell_data = cistopic_obj.cell_data
                if args.groupby in cell_data.columns:
                    ax.clear()
                    categories = cell_data[args.groupby]
                    unique_cats = categories.unique()
                    colors = plt.cm.tab20.colors
                    color_map = {
                        cat: colors[i % len(colors)]
                        for i, cat in enumerate(unique_cats)
                    }
                    for cat in unique_cats:
                        mask = categories == cat
                        mask_idx = mask[mask].index
                        umap_subset = umap_df.loc[
                            umap_df.index.isin(mask_idx)
                        ]
                        ax.scatter(
                            umap_subset.iloc[:, 0],
                            umap_subset.iloc[:, 1],
                            s=3,
                            alpha=0.5,
                            label=cat,
                            c=[color_map[cat]],
                            rasterized=True,
                        )
                    ax.set_title(f"{n_top} topics")
                    ax.set_xlabel("UMAP1")
                    ax.set_ylabel("UMAP2")
                    ax.set_xticks([])
                    ax.set_yticks([])

        # Only add legend to last subplot
        if args.groupby is not None and n_models > 0:
            axes[0, -1].legend(
                bbox_to_anchor=(1.05, 1),
                loc="upper left",
                fontsize=6,
                markerscale=3,
            )

        plt.tight_layout()
        umap_path = os.path.join(args.output_dir, "umap_per_topic_count.pdf")
        plt.savefig(umap_path, bbox_inches="tight", dpi=150)
        plt.close()
        logger.info("  UMAPs saved to %s", umap_path)

    logger.info("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Evaluate trained cisTopic LDA models.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python evaluate_models.py \\
      --cistopic_obj cistopic.pkl \\
      --models ./models/models.pkl \\
      --output_dir ./evaluation \\
      --groupby cell_type

  # Evaluate only specific topic counts, skip UMAP
  python evaluate_models.py \\
      --cistopic_obj cistopic.pkl \\
      --models ./models/ \\
      --output_dir ./evaluation \\
      --n_topics 20 30 40 --no_run_umap
""",
    )

    parser.add_argument(
        "--cistopic_obj",
        type=str,
        required=True,
        help="Path to cisTopic object pickle file (.pkl)",
    )
    parser.add_argument(
        "--models",
        type=str,
        required=True,
        help="Path to models pickle file (.pkl) or directory of .pkl files",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory for evaluation results",
    )
    parser.add_argument(
        "--n_topics",
        type=int,
        nargs="+",
        default=None,
        help="Evaluate only these topic counts (default: all)",
    )
    parser.add_argument(
        "--groupby",
        type=str,
        default=None,
        help="Cell metadata column for coloring UMAPs (e.g. 'cell_type')",
    )
    parser.add_argument(
        "--metrics",
        type=str,
        nargs="+",
        default=ALL_METRICS,
        help="Metrics to compute (default: all four)",
    )
    parser.add_argument(
        "--run_umap",
        action="store_true",
        default=True,
        help="Generate UMAPs per topic count (default: True)",
    )
    parser.add_argument(
        "--no_run_umap",
        action="store_false",
        dest="run_umap",
        help="Skip UMAP generation (faster, metrics only)",
    )

    args = parser.parse_args()
    main(args)
