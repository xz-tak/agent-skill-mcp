#!/usr/bin/env python3
"""Extract cell embeddings using scimilarity.

Usage:
    python embed_cells.py --input data.h5ad --model-dir /path/to/model_v1.1 --output embedded.h5ad
    python embed_cells.py --input data.h5ad --model-dir /path/to/model_v1.1 --output embedded.h5ad --gpu
"""

import argparse


def main():
    parser = argparse.ArgumentParser(description="Compute scimilarity cell embeddings")
    parser.add_argument("--input", required=True, help="Input h5ad file path")
    parser.add_argument("--model-dir", default=None, help="Path to scimilarity model directory")
    parser.add_argument("--s3", action="store_true", help="Auto-cache model from S3 (embedding tier)")
    parser.add_argument("--output", required=True, help="Output h5ad file path")
    parser.add_argument("--gpu", action="store_true", help="Use GPU for embedding")
    parser.add_argument("--buffer-size", type=int, default=10000, help="Batch size for embedding")
    args = parser.parse_args()

    import scanpy as sc
    from scimilarity import CellEmbedding
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

    # Ensure raw counts
    if "counts" not in adata.layers:
        print("  Storing X as layers['counts']")
        adata.layers["counts"] = adata.X.copy()

    # Initialize
    print(f"Loading model: {model_dir}")
    ce = CellEmbedding(model_path=model_dir, use_gpu=args.gpu)
    print(f"  Genes: {ce.n_genes}, Latent dim: {ce.latent_dim}")

    # Prepare data
    print("Aligning gene space...")
    adata_aligned = align_dataset(adata, ce.gene_order)
    print("Log normalizing...")
    adata_aligned = lognorm_counts(adata_aligned)

    # Embed
    print(f"Computing embeddings (buffer_size={args.buffer_size})...")
    embeddings = ce.get_embeddings(adata_aligned.X, buffer_size=args.buffer_size)
    adata.obsm["X_scimilarity"] = embeddings
    print(f"  Embedding shape: {embeddings.shape}")

    # Save
    print(f"Saving to: {args.output}")
    adata.write_h5ad(args.output)
    print("Done.")


if __name__ == "__main__":
    main()
