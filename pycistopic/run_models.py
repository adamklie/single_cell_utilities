#!/usr/bin/env python
"""Train LDA topic models on a cisTopic object (CGS or MALLET).

Outputs:
    {output}/models.pkl          — pickled List[CistopicLDAModel]
    {output}/models_params.yaml  — parameters used for reproducibility
"""

import argparse
import logging
import os
import pickle
import sys
import time

import yaml

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------

def load_cistopic_object(path):
    """Load and validate a CistopicObject from a pickle file."""
    if not os.path.isfile(path):
        logger.error("Input file not found: %s", path)
        sys.exit(1)
    try:
        with open(path, "rb") as f:
            obj = pickle.load(f)
    except Exception as e:
        logger.error("Could not load pickle file: %s\n  %s", path, e)
        sys.exit(1)

    # Verify it's a CistopicObject
    cls_name = type(obj).__name__
    if "CistopicObject" not in cls_name and "cistopic" not in cls_name.lower():
        logger.warning(
            "Loaded object is type '%s' — expected CistopicObject. "
            "Proceeding anyway.",
            cls_name,
        )
    return obj


def validate_mallet(args):
    """Validate MALLET-specific requirements."""
    if args.mallet_path is None:
        logger.error("--mallet_path is required when --method=mallet")
        sys.exit(1)
    if not os.path.isfile(args.mallet_path):
        logger.error("MALLET binary not found: %s", args.mallet_path)
        sys.exit(1)
    if not os.access(args.mallet_path, os.X_OK):
        logger.warning(
            "MALLET binary may not be executable: %s", args.mallet_path
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(args):
    # --- Load cisTopic object ---
    logger.info("Loading cisTopic object from %s", args.input)
    cistopic_obj = load_cistopic_object(args.input)

    n_cells = cistopic_obj.fragment_matrix.shape[1]
    n_regions = cistopic_obj.fragment_matrix.shape[0]
    logger.info("  cisTopic object: %d regions x %d cells", n_regions, n_cells)

    # --- Validate n_topics ---
    n_topics = args.n_topics
    for nt in n_topics:
        if nt <= 0:
            logger.error("n_topics must be positive integers, got %d", nt)
            sys.exit(1)
        if nt > n_cells:
            logger.warning(
                "n_topics=%d exceeds n_cells=%d — model may not converge", nt, n_cells
            )
        if nt > n_regions:
            logger.warning(
                "n_topics=%d exceeds n_regions=%d — model may not converge",
                nt,
                n_regions,
            )

    # --- Validate hyperparameters ---
    if args.alpha <= 0:
        logger.error("--alpha must be > 0, got %s", args.alpha)
        sys.exit(1)
    if args.eta <= 0:
        logger.error("--eta must be > 0, got %s", args.eta)
        sys.exit(1)

    # --- Validate method-specific args ---
    if args.method == "mallet":
        validate_mallet(args)

    # --- Prepare output directory ---
    os.makedirs(args.output, exist_ok=True)
    # Verify writability
    test_file = os.path.join(args.output, ".write_test")
    try:
        with open(test_file, "w") as f:
            f.write("test")
        os.remove(test_file)
    except OSError:
        logger.error("Output directory is not writable: %s", args.output)
        sys.exit(1)

    # --- Prepare temp directory ---
    if args.temp_dir is not None:
        os.makedirs(args.temp_dir, exist_ok=True)

    # --- Prepare save_path ---
    if args.save_path is not None:
        os.makedirs(args.save_path, exist_ok=True)

    # --- Pre-flight summary ---
    logger.info("=" * 60)
    logger.info("Pre-flight summary")
    logger.info("  Method:          %s", args.method)
    logger.info("  n_topics:        %s", n_topics)
    logger.info("  n_iter:          %d", args.n_iter)
    logger.info("  n_cpu:           %d", args.n_cpu)
    logger.info("  alpha:           %s", args.alpha)
    logger.info("  alpha_by_topic:  %s", args.alpha_by_topic)
    logger.info("  eta:             %s", args.eta)
    logger.info("  eta_by_topic:    %s", args.eta_by_topic)
    logger.info("  seed:            %d", args.seed)
    logger.info("  Input shape:     %d regions x %d cells", n_regions, n_cells)
    logger.info("  Output dir:      %s", args.output)
    if args.method == "mallet":
        logger.info("  MALLET binary:   %s", args.mallet_path)
        logger.info("  MALLET memory:   %s", args.mallet_memory)
        logger.info("  Reuse corpus:    %s", args.reuse_corpus)
        logger.info("  top_topics_coh:  %d", args.top_topics_coh)
    logger.info("=" * 60)

    # --- Save params ---
    params = {
        "method": args.method,
        "input": os.path.abspath(args.input),
        "output": os.path.abspath(args.output),
        "n_topics": n_topics,
        "n_iter": args.n_iter,
        "n_cpu": args.n_cpu,
        "alpha": args.alpha,
        "alpha_by_topic": args.alpha_by_topic,
        "eta": args.eta,
        "eta_by_topic": args.eta_by_topic,
        "seed": args.seed,
        "save_path": args.save_path,
        "temp_dir": args.temp_dir,
    }
    if args.method == "mallet":
        params["mallet_path"] = args.mallet_path
        params["mallet_memory"] = args.mallet_memory
        params["reuse_corpus"] = args.reuse_corpus
        params["top_topics_coh"] = args.top_topics_coh

    params_path = os.path.join(args.output, "models_params.yaml")
    with open(params_path, "w") as f:
        yaml.dump(params, f, default_flow_style=False, sort_keys=False)
    logger.info("Parameters saved to %s", params_path)

    # --- Run models ---
    from pycisTopic.lda_models import run_cgs_models, run_cgs_models_mallet

    start_time = time.time()

    if args.method == "cgs":
        logger.info("Running CGS models...")
        models = run_cgs_models(
            cistopic_obj,
            n_topics=n_topics,
            n_cpu=args.n_cpu,
            n_iter=args.n_iter,
            random_state=args.seed,
            alpha=args.alpha,
            alpha_by_topic=args.alpha_by_topic,
            eta=args.eta,
            eta_by_topic=args.eta_by_topic,
            save_path=args.save_path,
            _temp_dir=args.temp_dir,
        )
    else:  # mallet
        logger.info("Running MALLET models...")
        os.environ["MALLET_MEMORY"] = args.mallet_memory
        models = run_cgs_models_mallet(
            args.mallet_path,
            cistopic_obj,
            n_topics=n_topics,
            n_cpu=args.n_cpu,
            n_iter=args.n_iter,
            random_state=args.seed,
            alpha=args.alpha,
            alpha_by_topic=args.alpha_by_topic,
            eta=args.eta,
            eta_by_topic=args.eta_by_topic,
            save_path=args.save_path,
            top_topics_coh=args.top_topics_coh,
            tmp_path=args.temp_dir,
            reuse_corpus=args.reuse_corpus,
        )

    total_time = time.time() - start_time
    logger.info("Training completed in %.1f seconds", total_time)

    # --- Post-training summary ---
    logger.info("Model summary:")
    for model in models:
        n_top = model.cell_topic.shape[0]
        ll = getattr(model, "log_likelihood", None)
        ll_str = f"{ll:.2f}" if ll is not None else "N/A"
        logger.info("  n_topics=%d  log_likelihood=%s", n_top, ll_str)

    # --- Save models ---
    output_path = os.path.join(args.output, "models.pkl")
    logger.info("Saving models to %s", output_path)
    with open(output_path, "wb") as f:
        pickle.dump(models, f)
    logger.info("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Train LDA topic models on a cisTopic object.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # CGS
  python run_models.py \\
      --input cistopic.pkl --output ./models \\
      --n_topics 10 20 30 40 50 --method cgs --n_cpu 8

  # MALLET
  python run_models.py \\
      --input cistopic.pkl --output ./models \\
      --n_topics 10 20 30 40 50 --method mallet \\
      --mallet_path /path/to/mallet --n_cpu 8
""",
    )

    # --- Required ---
    parser.add_argument(
        "--input", "-i",
        type=str,
        required=True,
        help="Path to cisTopic object pickle file (.pkl)",
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        required=True,
        help="Output directory. Models saved as models.pkl inside this directory",
    )
    parser.add_argument(
        "--n_topics", "-nt",
        type=int,
        nargs="+",
        required=True,
        help="Topic counts to train (space-separated integers, e.g. 10 20 30)",
    )
    parser.add_argument(
        "--method", "-m",
        type=str,
        choices=["cgs", "mallet"],
        required=True,
        help="LDA inference method: 'cgs' (collapsed Gibbs sampling) or 'mallet'",
    )

    # --- Hyperparameters ---
    parser.add_argument(
        "--n_cpu", "-c",
        type=int,
        default=1,
        help="Number of CPU cores (default: 1)",
    )
    parser.add_argument(
        "--n_iter",
        type=int,
        default=150,
        help="Number of iterations (default: 150)",
    )
    parser.add_argument(
        "--alpha", "-a",
        type=float,
        default=50.0,
        help="Alpha hyperparameter (default: 50.0)",
    )
    parser.add_argument(
        "--alpha_by_topic",
        action="store_true",
        default=True,
        help="Divide alpha by number of topics (default: True)",
    )
    parser.add_argument(
        "--no_alpha_by_topic",
        action="store_false",
        dest="alpha_by_topic",
        help="Do NOT divide alpha by number of topics",
    )
    parser.add_argument(
        "--eta", "-e",
        type=float,
        default=0.1,
        help="Eta hyperparameter (default: 0.1)",
    )
    parser.add_argument(
        "--eta_by_topic",
        action="store_true",
        default=False,
        help="Divide eta by number of topics (default: False)",
    )
    parser.add_argument(
        "--seed", "-s",
        type=int,
        default=555,
        help="Random seed (default: 555)",
    )
    parser.add_argument(
        "--save_path",
        type=str,
        default=None,
        help="Path to save intermediate model files during training",
    )
    parser.add_argument(
        "--temp_dir",
        type=str,
        default=None,
        help="Temporary directory for intermediate files",
    )

    # --- MALLET-specific ---
    mallet_group = parser.add_argument_group("MALLET-specific options")
    mallet_group.add_argument(
        "--mallet_path",
        type=str,
        default=None,
        help="Path to MALLET binary (required if --method=mallet)",
    )
    mallet_group.add_argument(
        "--mallet_memory",
        type=str,
        default="100G",
        help="Java heap memory for MALLET (default: 100G)",
    )
    mallet_group.add_argument(
        "--reuse_corpus",
        action="store_true",
        default=False,
        help="Reuse existing MALLET corpus if available",
    )
    mallet_group.add_argument(
        "--top_topics_coh",
        type=int,
        default=5,
        help="Number of top topics for coherence (default: 5)",
    )

    args = parser.parse_args()
    main(args)
