#!/usr/bin/env python3
"""CLI wrapper for Geneformer InSilicoPerturber.

Runs in silico gene perturbation (delete, overexpress, inhibit, activate)
and measures embedding shifts for goal state analysis.

Usage:
    conda activate geneformer

    # Delete perturbation with pre-computed state embeddings
    python perturb_cells.py --s3 \
        --input-data tokenized.dataset \
        --output-dir isp/ --output-prefix my_isp \
        --perturb-type delete \
        --model-type CellClassifier --num-classes 3 \
        --cell-states '{"state_key":"disease","start_state":"dcm","goal_state":"nf","alt_states":["hcm"]}' \
        --state-embs-file state_embs.pkl

    # Compute state embeddings on the fly, then perturb
    python perturb_cells.py --model-dir /path/to/model \
        --input-data tokenized.dataset \
        --output-dir isp/ --output-prefix my_isp \
        --perturb-type delete \
        --model-type CellClassifier --num-classes 3 \
        --cell-states '{"state_key":"disease","start_state":"dcm","goal_state":"nf","alt_states":["hcm"]}' \
        --compute-state-embs

    # Overexpress specific genes
    python perturb_cells.py --s3 \
        --input-data tokenized.dataset \
        --output-dir isp/ --output-prefix my_isp \
        --perturb-type overexpress \
        --genes-to-perturb ENSG00000141510,ENSG00000171862

    # Pairwise combination perturbation with anchor gene
    python perturb_cells.py --s3 \
        --input-data tokenized.dataset \
        --output-dir isp/ --output-prefix combo_isp \
        --perturb-type delete --combos 1 \
        --anchor-gene ENSG00000141510
"""

import argparse
import json
import os
import pickle
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run Geneformer in silico perturbation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Model source (mutually exclusive)
    model_group = parser.add_mutually_exclusive_group(required=True)
    model_group.add_argument("--model-dir", type=str, help="Direct path to model directory")
    model_group.add_argument("--s3", action="store_true", help="Auto-cache model from S3")
    parser.add_argument(
        "--s3-tier",
        type=str,
        default="V2-104M",
        help="S3 model tier (default: V2-104M)",
    )

    # Required I/O
    parser.add_argument("--input-data", required=True, help="Path to tokenized .dataset")
    parser.add_argument("--output-dir", required=True, help="Output directory for perturbation results")
    parser.add_argument("--output-prefix", required=True, help="Output file prefix")

    # Perturbation settings
    parser.add_argument(
        "--perturb-type",
        choices=["delete", "overexpress", "inhibit", "activate"],
        default="delete",
        help="Perturbation type (default: delete)",
    )
    parser.add_argument(
        "--perturb-rank-shift",
        type=int,
        choices=[1, 2, 3],
        default=None,
        help="Rank shift magnitude for inhibit/activate",
    )
    parser.add_argument(
        "--genes-to-perturb",
        type=str,
        default="all",
        help='Genes to perturb: "all" or comma-separated Ensembl IDs (default: all)',
    )
    parser.add_argument(
        "--combos",
        type=int,
        choices=[0, 1],
        default=0,
        help="0: individual genes, 1: pairwise combinations (default: 0)",
    )
    parser.add_argument("--anchor-gene", type=str, default=None, help="Ensembl ID for combination anchor")

    # Model settings
    parser.add_argument(
        "--model-type",
        choices=[
            "Pretrained",
            "CellClassifier",
            "GeneClassifier",
            "MTLCellClassifier",
            "Pretrained-Quantized",
            "MTLCellClassifier-Quantized",
        ],
        default="Pretrained",
        help="Model type (default: Pretrained)",
    )
    parser.add_argument("--num-classes", type=int, default=0, help="Number of classes (default: 0)")
    parser.add_argument(
        "--emb-mode",
        choices=["cls", "cell", "cls_and_gene", "cell_and_gene"],
        default="cls",
        help="Embedding mode (default: cls)",
    )

    # Filtering and state modeling
    parser.add_argument("--filter-data", type=str, default=None, help="JSON string for cell filtering")
    parser.add_argument(
        "--cell-states",
        type=str,
        default=None,
        help='JSON string for cell_states_to_model, e.g. \'{"state_key":"disease","start_state":"dcm","goal_state":"nf","alt_states":["hcm"]}\'',
    )
    parser.add_argument(
        "--state-embs-file",
        type=str,
        default=None,
        help="Path to state_embs_dict pickle (from EmbExtractor.get_state_embs)",
    )
    parser.add_argument(
        "--compute-state-embs",
        action="store_true",
        help="Compute state embeddings on the fly (requires --cell-states)",
    )

    # Performance
    parser.add_argument("--max-ncells", type=int, default=None, help="Max number of cells to process")
    parser.add_argument(
        "--emb-layer",
        type=int,
        choices=[-1, 0],
        default=-1,
        help="Embedding layer: -1 (general) or 0 (task-specific) (default: -1)",
    )
    parser.add_argument("--forward-batch-size", type=int, default=100, help="Forward batch size (default: 100)")
    parser.add_argument("--nproc", type=int, default=4, help="Number of processes (default: 4)")

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

    # Parse JSON arguments
    filter_data = json.loads(args.filter_data) if args.filter_data else None
    cell_states = json.loads(args.cell_states) if args.cell_states else None

    # Parse genes_to_perturb
    if args.genes_to_perturb == "all":
        genes_to_perturb = "all"
    else:
        genes_to_perturb = [g.strip() for g in args.genes_to_perturb.split(",")]

    # Resolve state embeddings
    state_embs_dict = None

    if args.compute_state_embs:
        if cell_states is None:
            raise ValueError("--compute-state-embs requires --cell-states")
        from geneformer import EmbExtractor

        print("Computing state embeddings on the fly...")
        embex = EmbExtractor(
            model_type=args.model_type,
            num_classes=args.num_classes,
            filter_data=filter_data,
            max_ncells=args.max_ncells,
            emb_layer=args.emb_layer,
            summary_stat="exact_mean",
            forward_batch_size=args.forward_batch_size,
            nproc=args.nproc,
        )
        state_embs_dict = embex.get_state_embs(
            cell_states, model_dir, args.input_data, args.output_dir, args.output_prefix,
        )
        # Save for reuse
        embs_path = os.path.join(args.output_dir, f"{args.output_prefix}_state_embs.pkl")
        os.makedirs(args.output_dir, exist_ok=True)
        with open(embs_path, "wb") as f:
            pickle.dump(state_embs_dict, f)
        print(f"State embeddings saved to: {embs_path}")

    elif args.state_embs_file:
        with open(args.state_embs_file, "rb") as f:
            state_embs_dict = pickle.load(f)
        print(f"Loaded state embeddings from: {args.state_embs_file}")

    # Create and run perturbation
    from geneformer import InSilicoPerturber

    isp = InSilicoPerturber(
        perturb_type=args.perturb_type,
        perturb_rank_shift=args.perturb_rank_shift,
        genes_to_perturb=genes_to_perturb,
        combos=args.combos,
        anchor_gene=args.anchor_gene,
        model_type=args.model_type,
        num_classes=args.num_classes,
        emb_mode=args.emb_mode,
        filter_data=filter_data,
        cell_states_to_model=cell_states,
        state_embs_dict=state_embs_dict,
        max_ncells=args.max_ncells,
        emb_layer=args.emb_layer,
        forward_batch_size=args.forward_batch_size,
        nproc=args.nproc,
    )

    print(f"Running {args.perturb_type} perturbation on {args.input_data}...")
    isp.perturb_data(model_dir, args.input_data, args.output_dir, args.output_prefix)
    print(f"Perturbation complete. Results saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
