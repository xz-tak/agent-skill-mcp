#!/usr/bin/env python3
"""CLI wrapper for Geneformer Classifier.

Usage:
    python classify_cells.py prepare --s3 --input-data tokenized.dataset --output-dir clf/ --output-prefix my_clf \
        --classifier-type cell --cell-state-dict '{"state_key":"disease","states":["healthy","diseased"]}'

    python classify_cells.py validate --s3 --prepared-data clf/my_clf_labeled.dataset \
        --id-class-dict clf/my_clf_id_class_dict.pkl --output-dir clf/ --output-prefix my_clf

    python classify_cells.py train-all --s3 --prepared-data clf/my_clf_labeled.dataset \
        --id-class-dict clf/my_clf_id_class_dict.pkl --output-dir clf/ --output-prefix my_clf

    python classify_cells.py evaluate --model-dir fine_tuned_model/ \
        --id-class-dict clf/my_clf_id_class_dict.pkl --test-data test.dataset \
        --output-dir clf/ --output-prefix eval
"""

import argparse
import json
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Cell/gene classification using Geneformer Classifier"
    )

    parser.add_argument(
        "subcommand",
        choices=["prepare", "validate", "train-all", "evaluate"],
        help="Pipeline step to run",
    )

    # Model source
    model_group = parser.add_mutually_exclusive_group()
    model_group.add_argument("--model-dir", help="Path to model directory")
    model_group.add_argument("--s3", action="store_true", help="Auto-cache model from S3")

    parser.add_argument(
        "--s3-tier",
        choices=["V2-104M", "V2-316M", "V2-104M-CLcancer"],
        default="V2-104M",
        help="S3 model tier (default: V2-104M)",
    )

    # Data inputs (subcommand-dependent)
    parser.add_argument("--input-data", help="Input .dataset path (for prepare)")
    parser.add_argument(
        "--prepared-data", help="Prepared labeled .dataset path (for validate/train-all)"
    )
    parser.add_argument(
        "--id-class-dict",
        help="Path to id_class_dict.pkl (for validate/train-all/evaluate)",
    )
    parser.add_argument("--test-data", help="Test .dataset path (for evaluate)")

    # Output
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--output-prefix", required=True, help="Output prefix")

    # Classifier config
    parser.add_argument(
        "--classifier-type",
        choices=["cell", "gene"],
        default="cell",
        help="Classifier type (default: cell)",
    )
    parser.add_argument(
        "--cell-state-dict",
        default=None,
        help='JSON string, e.g. \'{"state_key":"disease","states":["healthy","diseased"]}\'',
    )
    parser.add_argument(
        "--gene-class-dict",
        default=None,
        help='JSON string, e.g. \'{"Label_A":["ENSG..."],"Label_B":["ENSG..."]}\'',
    )
    parser.add_argument(
        "--filter-data",
        default=None,
        help='JSON string for cell filtering, e.g. \'{"cell_type":["Cardiomyocyte"]}\'',
    )
    parser.add_argument(
        "--freeze-layers",
        type=int,
        default=0,
        help="Number of transformer layers to freeze (default: 0)",
    )
    parser.add_argument(
        "--num-crossval-splits",
        type=int,
        choices=[0, 1, 5],
        default=1,
        help="Cross-validation splits: 0=no split, 1=train/eval, 5=5-fold (default: 1)",
    )
    parser.add_argument(
        "--forward-batch-size",
        type=int,
        default=200,
        help="Forward batch size (default: 200)",
    )
    parser.add_argument(
        "--nproc", type=int, default=4, help="Number of processes (default: 4)"
    )
    parser.add_argument(
        "--learning-rate", type=float, default=None, help="Learning rate"
    )
    parser.add_argument(
        "--num-epochs", type=int, default=None, help="Number of training epochs"
    )
    parser.add_argument(
        "--quantize",
        action="store_true",
        help="Enable 8-bit quantization",
    )

    return parser.parse_args()


def _parse_json_arg(value, arg_name):
    """Parse a JSON string argument, exiting on error."""
    if value is None:
        return None
    try:
        return json.loads(value)
    except json.JSONDecodeError as exc:
        print(f"Error parsing --{arg_name} JSON: {exc}", file=sys.stderr)
        sys.exit(1)


def main():
    args = parse_args()

    # Resolve model path (not needed for prepare)
    model_dir = None
    if args.subcommand != "prepare":
        if args.s3:
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            from s3_cache import get_model_path

            model_dir = get_model_path(tier=args.s3_tier)
        elif args.model_dir:
            model_dir = args.model_dir
        else:
            print(
                f"Error: --model-dir or --s3 is required for '{args.subcommand}'",
                file=sys.stderr,
            )
            sys.exit(1)

    # Parse JSON arguments
    cell_state_dict = _parse_json_arg(args.cell_state_dict, "cell-state-dict")
    gene_class_dict = _parse_json_arg(args.gene_class_dict, "gene-class-dict")
    filter_data = _parse_json_arg(args.filter_data, "filter-data")

    # Build training_args dict from CLI flags
    training_args = {}
    if args.learning_rate is not None:
        training_args["learning_rate"] = args.learning_rate
    if args.num_epochs is not None:
        training_args["num_train_epochs"] = args.num_epochs
    training_args_param = training_args if training_args else None

    os.makedirs(args.output_dir, exist_ok=True)

    from geneformer import Classifier

    classifier = Classifier(
        classifier=args.classifier_type,
        cell_state_dict=cell_state_dict,
        gene_class_dict=gene_class_dict,
        filter_data=filter_data,
        training_args=training_args_param,
        freeze_layers=args.freeze_layers,
        num_crossval_splits=args.num_crossval_splits,
        forward_batch_size=args.forward_batch_size,
        nproc=args.nproc,
        quantize=args.quantize,
    )

    if args.subcommand == "prepare":
        if not args.input_data:
            print("Error: --input-data is required for 'prepare'", file=sys.stderr)
            sys.exit(1)
        print(f"Preparing data from: {args.input_data}")
        classifier.prepare_data(args.input_data, args.output_dir, args.output_prefix)
        print(f"Data prepared in: {args.output_dir}")
        print(f"  Labeled dataset: {args.output_prefix}_labeled.dataset")
        print(f"  ID-class dict:   {args.output_prefix}_id_class_dict.pkl")

    elif args.subcommand == "validate":
        if not args.prepared_data or not args.id_class_dict:
            print(
                "Error: --prepared-data and --id-class-dict are required for 'validate'",
                file=sys.stderr,
            )
            sys.exit(1)
        print(f"Validating with model: {model_dir}")
        print(f"Prepared data: {args.prepared_data}")
        all_metrics = classifier.validate(
            model_dir,
            args.prepared_data,
            args.id_class_dict,
            args.output_dir,
            args.output_prefix,
            predict_eval=True,
        )
        print(f"\nValidation complete. Results in: {args.output_dir}")
        if "conf_matrix" in all_metrics:
            classifier.plot_conf_mat(
                {"Geneformer": all_metrics["conf_matrix"]},
                args.output_dir,
                args.output_prefix,
            )
            print(f"Confusion matrix saved.")

    elif args.subcommand == "train-all":
        if not args.prepared_data or not args.id_class_dict:
            print(
                "Error: --prepared-data and --id-class-dict are required for 'train-all'",
                file=sys.stderr,
            )
            sys.exit(1)
        print(f"Training on all data with model: {model_dir}")
        print(f"Prepared data: {args.prepared_data}")
        classifier.train_all_data(
            model_dir,
            args.prepared_data,
            args.id_class_dict,
            args.output_dir,
            args.output_prefix,
        )
        print(f"\nTraining complete. Model saved to: {args.output_dir}")

    elif args.subcommand == "evaluate":
        if not args.test_data or not args.id_class_dict:
            print(
                "Error: --test-data, --id-class-dict, and --model-dir are required for 'evaluate'",
                file=sys.stderr,
            )
            sys.exit(1)
        print(f"Evaluating model: {model_dir}")
        print(f"Test data: {args.test_data}")
        classifier.evaluate_saved_model(
            model_dir,
            args.id_class_dict,
            args.test_data,
            args.output_dir,
            args.output_prefix,
        )
        print(f"\nEvaluation complete. Results in: {args.output_dir}")


if __name__ == "__main__":
    main()
