#!/usr/bin/env python3
"""Gene attribution analysis using scimilarity Interpreter.

Usage:
    python interpret_genes.py --input data.h5ad --model-dir /path/to/model_v1.1 \
        --anchor-col celltype --anchor-val macrophage \
        --negative-col celltype --negative-val fibroblast \
        --output-prefix results

    python interpret_genes.py --input data.h5ad --model-dir /path/to/model_v1.1 \
        --anchor-col condition --anchor-val disease \
        --negative-col condition --negative-val normal \
        --top-n 20 --output-prefix disease_vs_normal
"""

import argparse


def main():
    parser = argparse.ArgumentParser(description="Gene attribution analysis with scimilarity")
    parser.add_argument("--input", required=True, help="Input h5ad file path")
    parser.add_argument("--model-dir", default=None, help="Path to scimilarity model directory")
    parser.add_argument("--s3", action="store_true", help="Auto-cache model from S3 (embedding tier)")
    parser.add_argument("--anchor-col", required=True, help="obs column for selecting anchor cells")
    parser.add_argument("--anchor-val", required=True, help="Value in anchor-col to select anchors")
    parser.add_argument("--negative-col", required=True, help="obs column for selecting negative cells")
    parser.add_argument("--negative-val", required=True, help="Value in negative-col to select negatives")
    parser.add_argument("--output-prefix", required=True, help="Output prefix for results")
    parser.add_argument("--gpu", action="store_true", help="Use GPU")
    parser.add_argument("--top-n", type=int, default=15, help="Number of top genes to plot (default: 15)")
    parser.add_argument("--max-cells", type=int, default=500,
                        help="Max cells per group for attribution (default: 500)")
    args = parser.parse_args()

    import numpy as np
    import scanpy as sc
    from scimilarity import CellEmbedding, Interpreter
    from scimilarity.utils import align_dataset, lognorm_counts

    # Resolve model path
    model_dir = args.model_dir
    if args.s3:
        import sys, os
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from s3_cache import get_model_path
        model_dir = get_model_path(tier="embedding")
    if not model_dir:
        print("ERROR: provide --model-dir or --s3")
        return

    print(f"Loading data: {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Cells: {adata.n_obs}, Genes: {adata.n_vars}")

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # Initialize model
    print(f"Loading model: {model_dir}")
    ce = CellEmbedding(model_path=model_dir, use_gpu=args.gpu)

    # Prepare data
    print("Aligning and normalizing...")
    adata_aligned = align_dataset(adata, ce.gene_order)
    adata_aligned = lognorm_counts(adata_aligned)

    # Select anchor and negative cells
    anchor_mask = adata.obs[args.anchor_col] == args.anchor_val
    negative_mask = adata.obs[args.negative_col] == args.negative_val

    n_anchors = anchor_mask.sum()
    n_negatives = negative_mask.sum()
    print(f"  Anchors ({args.anchor_val}): {n_anchors} cells")
    print(f"  Negatives ({args.negative_val}): {n_negatives} cells")

    if n_anchors == 0 or n_negatives == 0:
        print("ERROR: No anchor or negative cells found")
        return

    # Subsample if too many cells
    n_use = min(n_anchors, n_negatives, args.max_cells)
    print(f"  Using {n_use} cells per group")

    anchor_idx = np.where(anchor_mask)[0]
    negative_idx = np.where(negative_mask)[0]

    rng = np.random.default_rng(42)
    if len(anchor_idx) > n_use:
        anchor_idx = rng.choice(anchor_idx, n_use, replace=False)
    if len(negative_idx) > n_use:
        negative_idx = rng.choice(negative_idx, n_use, replace=False)

    # Match counts (need same shape)
    n_matched = min(len(anchor_idx), len(negative_idx))
    anchor_idx = anchor_idx[:n_matched]
    negative_idx = negative_idx[:n_matched]

    anchor_X = adata_aligned[anchor_idx].X
    negative_X = adata_aligned[negative_idx].X

    # Compute attributions
    print("Computing attributions (Integrated Gradients)...")
    interp = Interpreter(encoder=ce.model, gene_order=ce.gene_order)
    attrs = interp.get_attributions(anchor_X, negative_X)

    # Rank genes
    ranked = interp.get_ranked_genes(attrs)
    print(f"\nTop {args.top_n} genes by attribution:")
    print(ranked.head(args.top_n).to_string(index=False))

    # Save results
    csv_path = f"{args.output_prefix}_ranked_genes.csv"
    ranked.to_csv(csv_path, index=False)
    print(f"\nSaved ranked genes: {csv_path}")

    # Plot
    plot_path = f"{args.output_prefix}_top_genes.pdf"
    interp.plot_ranked_genes(ranked, n_plot=args.top_n, filename=plot_path)
    print(f"Saved plot: {plot_path}")

    print("Done.")


if __name__ == "__main__":
    main()
