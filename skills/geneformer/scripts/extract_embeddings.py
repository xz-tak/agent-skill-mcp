#!/usr/bin/env python3
"""CLI wrapper for Geneformer EmbExtractor.

Usage:
    python extract_embeddings.py --s3 --input-data tokenized.dataset --output-dir embs/ --output-prefix my_embs
    python extract_embeddings.py --model-dir /path/to/model --input-data tokenized.dataset --output-dir embs/ --output-prefix my_embs --emb-mode gene
"""

import argparse
import json
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract cell/gene embeddings using Geneformer EmbExtractor"
    )

    model_group = parser.add_mutually_exclusive_group(required=True)
    model_group.add_argument("--model-dir", help="Path to model directory")
    model_group.add_argument("--s3", action="store_true", help="Auto-cache model from S3")

    parser.add_argument(
        "--s3-tier",
        choices=["V2-104M", "V2-316M", "V2-104M-CLcancer"],
        default="V2-104M",
        help="S3 model tier (default: V2-104M)",
    )
    parser.add_argument("--input-data", required=True, help="Path to tokenized .dataset")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--output-prefix", required=True, help="Output prefix")
    parser.add_argument(
        "--model-type",
        choices=["Pretrained", "CellClassifier", "GeneClassifier"],
        default="Pretrained",
        help="Model type (default: Pretrained)",
    )
    parser.add_argument(
        "--num-classes", type=int, default=0, help="Number of classes (default: 0)"
    )
    parser.add_argument(
        "--emb-mode",
        choices=["cls", "cell", "gene"],
        default="cls",
        help="Embedding mode (default: cls)",
    )
    parser.add_argument(
        "--emb-layer",
        type=int,
        choices=[-1, 0],
        default=-1,
        help="Embedding layer: -1=2nd-to-last (general), 0=last (task-specific) (default: -1)",
    )
    parser.add_argument(
        "--filter-data",
        default=None,
        help='JSON string for cell filtering, e.g. \'{"cell_type":["T cell"]}\'',
    )
    parser.add_argument(
        "--max-ncells", type=int, default=None, help="Max cells to process (default: all)"
    )
    parser.add_argument(
        "--emb-label",
        default=None,
        help="Comma-separated column names to add as labels",
    )
    parser.add_argument(
        "--labels-to-plot",
        default=None,
        help="Comma-separated labels for coloring plots",
    )
    parser.add_argument(
        "--plot-style",
        choices=["heatmap", "umap", "none"],
        default="none",
        help="Plot style (default: none)",
    )
    parser.add_argument(
        "--summary-stat",
        choices=["none", "mean", "median", "exact_mean", "exact_median"],
        default="none",
        help="Summary statistic for embeddings (default: none)",
    )
    parser.add_argument(
        "--forward-batch-size",
        type=int,
        default=100,
        help="Forward batch size (default: 100)",
    )
    parser.add_argument(
        "--nproc", type=int, default=4, help="Number of processes (default: 4)"
    )
    parser.add_argument(
        "--output-torch-embs",
        action="store_true",
        help="Output embeddings as torch tensors",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Resolve model path
    if args.s3:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from s3_cache import get_model_path

        model_dir = get_model_path(tier=args.s3_tier)
    else:
        model_dir = args.model_dir

    # Parse optional JSON arguments
    filter_data = None
    if args.filter_data is not None:
        try:
            filter_data = json.loads(args.filter_data)
        except json.JSONDecodeError as exc:
            print(f"Error parsing --filter-data JSON: {exc}", file=sys.stderr)
            sys.exit(1)

    emb_label = args.emb_label.split(",") if args.emb_label else None
    labels_to_plot = args.labels_to_plot.split(",") if args.labels_to_plot else None
    summary_stat = None if args.summary_stat == "none" else args.summary_stat

    os.makedirs(args.output_dir, exist_ok=True)

    from geneformer import EmbExtractor

    embex = EmbExtractor(
        model_type=args.model_type,
        num_classes=args.num_classes,
        emb_mode=args.emb_mode,
        filter_data=filter_data,
        max_ncells=args.max_ncells,
        emb_layer=args.emb_layer,
        emb_label=emb_label,
        labels_to_plot=labels_to_plot,
        forward_batch_size=args.forward_batch_size,
        nproc=args.nproc,
        summary_stat=summary_stat,
    )

    print(f"Extracting {args.emb_mode} embeddings from: {args.input_data}")
    print(f"Model: {model_dir}")
    print(f"Output: {args.output_dir}/{args.output_prefix}")

    embs = embex.extract_embs(
        model_dir,
        args.input_data,
        args.output_dir,
        args.output_prefix,
        output_torch_embs=args.output_torch_embs,
    )

    print(f"\nEmbedding extraction complete.")
    if hasattr(embs, "shape"):
        print(f"  Shape: {embs.shape}")

    if args.plot_style != "none":
        print(f"Generating {args.plot_style} plot...")
        embex.plot_embs(
            embs,
            plot_style=args.plot_style,
            output_directory=args.output_dir,
            output_prefix=args.output_prefix,
        )
        print(f"Plot saved to {args.output_dir}")


if __name__ == "__main__":
    main()
