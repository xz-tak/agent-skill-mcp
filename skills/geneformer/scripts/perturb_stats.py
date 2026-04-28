#!/usr/bin/env python3
"""CLI wrapper for Geneformer InSilicoPerturberStats.

Computes statistical analysis of in silico perturbation results.
No model needed -- operates on raw perturbation pickle files.

Usage:
    conda activate geneformer

    # Goal state shift analysis (requires cell_states)
    python perturb_stats.py \
        --mode goal_state_shift \
        --input-dir isp_output/ \
        --output-dir stats/ --output-prefix my_stats \
        --cell-states '{"state_key":"disease","start_state":"dcm","goal_state":"nf","alt_states":["hcm"]}'

    # Mixture model (undirected impact)
    python perturb_stats.py \
        --mode mixture_model \
        --input-dir isp_output/ \
        --output-dir stats/ --output-prefix my_stats

    # Perturbation vs null distribution
    python perturb_stats.py \
        --mode vs_null \
        --input-dir isp_output/ \
        --null-dir null_output/ \
        --output-dir stats/ --output-prefix my_stats

    # Aggregate data for a single perturbation
    python perturb_stats.py \
        --mode aggregate_data \
        --input-dir isp_output/ \
        --output-dir stats/ --output-prefix my_stats

    # Aggregate gene-level shifts across perturbations
    python perturb_stats.py \
        --mode aggregate_gene_shifts \
        --input-dir isp_output/ \
        --output-dir stats/ --output-prefix my_stats
"""

import argparse
import json


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute statistics from Geneformer in silico perturbation results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Mode
    parser.add_argument(
        "--mode",
        choices=[
            "goal_state_shift",
            "vs_null",
            "mixture_model",
            "aggregate_data",
            "aggregate_gene_shifts",
        ],
        default="mixture_model",
        help="Statistics mode (default: mixture_model)",
    )

    # I/O directories
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing raw perturbation pickle files from InSilicoPerturber",
    )
    parser.add_argument(
        "--null-dir",
        type=str,
        default=None,
        help="Directory with null distribution data (required for vs_null mode)",
    )
    parser.add_argument("--output-dir", required=True, help="Output directory for statistics")
    parser.add_argument("--output-prefix", required=True, help="Output file prefix")

    # Gene settings
    parser.add_argument(
        "--genes-perturbed",
        type=str,
        default="all",
        help='Genes perturbed: "all" or comma-separated Ensembl IDs (default: all)',
    )
    parser.add_argument(
        "--combos",
        type=int,
        choices=[0, 1],
        default=0,
        help="0: individual genes, 1: pairwise combinations (default: 0)",
    )
    parser.add_argument("--anchor-gene", type=str, default=None, help="Ensembl ID for combination anchor")

    # Cell state modeling
    parser.add_argument(
        "--cell-states",
        type=str,
        default=None,
        help='JSON string for cell_states_to_model (required for goal_state_shift mode)',
    )

    # Pickle settings
    parser.add_argument(
        "--pickle-suffix",
        type=str,
        default="_raw.pickle",
        help='Suffix for pickle files (default: "_raw.pickle")',
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # Parse JSON arguments
    cell_states = json.loads(args.cell_states) if args.cell_states else None

    # Validate mode-specific requirements
    if args.mode == "goal_state_shift" and cell_states is None:
        raise ValueError("--cell-states is required for goal_state_shift mode")
    if args.mode == "vs_null" and args.null_dir is None:
        raise ValueError("--null-dir is required for vs_null mode")

    # Parse genes_perturbed
    if args.genes_perturbed == "all":
        genes_perturbed = "all"
    else:
        genes_perturbed = [g.strip() for g in args.genes_perturbed.split(",")]

    # Create and run stats
    from geneformer import InSilicoPerturberStats

    ispstats = InSilicoPerturberStats(
        mode=args.mode,
        genes_perturbed=genes_perturbed,
        combos=args.combos,
        anchor_gene=args.anchor_gene,
        cell_states_to_model=cell_states,
        pickle_suffix=args.pickle_suffix,
    )

    print(f"Computing {args.mode} statistics from: {args.input_dir}")
    ispstats.get_stats(args.input_dir, args.null_dir, args.output_dir, args.output_prefix)
    print(f"Statistics complete. Results saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
