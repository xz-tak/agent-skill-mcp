#!/usr/bin/env python3
"""End-to-end cell type annotation pipeline using scimilarity.

Usage:
    python annotate_cells.py --input data.h5ad --model-dir /path/to/model_v1.1 --output annotated.h5ad
    python annotate_cells.py --input data.h5ad --model-dir /path/to/model_v1.1 --output annotated.h5ad --gpu --k 100
    python annotate_cells.py --input data.h5ad --model-dir /path/to/model_v1.1 --output annotated.h5ad --blocklist "T cell,B cell"
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Annotate cell types using scimilarity")
    parser.add_argument("--input", required=True, help="Input h5ad file path")
    parser.add_argument("--model-dir", default=None, help="Path to scimilarity model directory")
    parser.add_argument("--s3", action="store_true", help="Auto-cache model from S3 (annotation tier)")
    parser.add_argument("--output", required=True, help="Output h5ad file path")
    parser.add_argument("--gpu", action="store_true", help="Use GPU for embedding")
    parser.add_argument("--k", type=int, default=50, help="Number of nearest neighbors (default: 50)")
    parser.add_argument("--weighting", action="store_true", help="Use distance weighting for predictions")
    parser.add_argument("--blocklist", type=str, default=None, help="Comma-separated cell types to exclude")
    parser.add_argument("--safelist", type=str, default=None, help="Comma-separated cell types to include only")
    parser.add_argument("--filter-cells", action="store_true", help="Apply QC filtering before annotation")
    parser.add_argument("--min-genes", type=int, default=400, help="Min genes per cell for filtering")
    parser.add_argument("--mito-percent", type=float, default=30.0, help="Max mitochondrial percent")
    args = parser.parse_args()

    import scanpy as sc
    from scimilarity import CellAnnotation
    from scimilarity.utils import align_dataset, lognorm_counts, filter_cells

    # Resolve model path
    model_dir = args.model_dir
    if args.s3:
        import sys, os
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from s3_cache import get_model_path
        model_dir = get_model_path(tier="annotation")
    if not model_dir:
        print("ERROR: provide --model-dir or --s3")
        return

    print(f"Loading data: {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Cells: {adata.n_obs}, Genes: {adata.n_vars}")

    # Ensure raw counts are available
    if "counts" not in adata.layers:
        print("  Storing X as layers['counts']")
        adata.layers["counts"] = adata.X.copy()

    # Optional QC filtering
    if args.filter_cells:
        n_before = adata.n_obs
        adata = filter_cells(adata, min_genes=args.min_genes, mito_percent=args.mito_percent)
        print(f"  Filtered: {n_before} -> {adata.n_obs} cells")

    # Initialize annotator
    print(f"Loading model: {model_dir}")
    ca = CellAnnotation(model_path=model_dir, use_gpu=args.gpu)
    print(f"  Genes in model: {ca.n_genes}, Latent dim: {ca.latent_dim}")

    # Apply blocklist/safelist
    if args.blocklist:
        labels = [l.strip() for l in args.blocklist.split(",")]
        ca.blocklist_celltypes(labels)
        print(f"  Blocklisted: {labels}")
    elif args.safelist:
        labels = [l.strip() for l in args.safelist.split(",")]
        ca.safelist_celltypes(labels)
        print(f"  Safelisted: {labels}")

    # Align and normalize
    print("Aligning gene space...")
    adata_aligned = align_dataset(adata, ca.gene_order)
    print("Log normalizing...")
    adata_aligned = lognorm_counts(adata_aligned)

    # Compute embeddings
    print("Computing embeddings...")
    embeddings = ca.get_embeddings(adata_aligned.X)
    adata.obsm["X_scimilarity"] = embeddings
    print(f"  Embedding shape: {embeddings.shape}")

    # Get predictions
    print(f"Predicting cell types (k={args.k})...")
    predictions, nn_idxs, nn_dists, stats = ca.get_predictions_knn(
        embeddings, k=args.k, weighting=args.weighting
    )

    # Store results
    adata.obs["celltype_hint"] = predictions.values
    adata.obs["min_dist"] = stats["min_dist"].values
    adata.obs["celltype_hits"] = stats["hits"].values
    adata.obs["celltype_hits_weighted"] = stats["hits_weighted"].values
    adata.obs["celltype_hint_stat"] = stats["vsAll"].values
    adata.obs["celltype_hint_weighted_stat"] = stats["vsAll_weighted"].values

    # Summary
    print("\nAnnotation summary:")
    print(adata.obs["celltype_hint"].value_counts().head(20).to_string())
    print(f"\nMean confidence (vsAll): {stats['vsAll'].mean():.3f}")
    print(f"Mean min_dist: {stats['min_dist'].mean():.4f}")

    # Save
    print(f"\nSaving to: {args.output}")
    adata.write_h5ad(args.output)
    print("Done.")


if __name__ == "__main__":
    main()
