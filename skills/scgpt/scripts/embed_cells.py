#!/usr/bin/env python3
"""Extract cell embeddings using scGPT.

Usage:
    python embed_cells.py --input data.h5ad --model-dir /path/to/scgpt_model --output embedded.h5ad
    python embed_cells.py --input data.h5ad --model-dir /path/to/scgpt_model --output embedded.h5ad --device cuda
"""

import argparse


def main():
    parser = argparse.ArgumentParser(description="Compute scGPT cell embeddings")
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--model-dir", default=None, help="scGPT model directory")
    parser.add_argument("--s3", action="store_true", help="Auto-cache model from S3 to /tmp/scgpt")
    parser.add_argument("--output", required=True, help="Output h5ad file")
    parser.add_argument("--device", default="cuda", help="Device (default: cuda)")
    parser.add_argument("--batch-size", type=int, default=64, help="Batch size")
    parser.add_argument("--max-length", type=int, default=1200, help="Max sequence length")
    parser.add_argument("--gene-col", default="feature_name", help="Gene column in var")
    args = parser.parse_args()

    from scgpt.tasks.cell_emb import embed_data

    model_dir = args.model_dir
    if args.s3:
        import sys, os
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from s3_cache import get_model_path
        model_dir = get_model_path()
    if not model_dir:
        print("ERROR: provide --model-dir or --s3")
        return

    print(f"Loading and embedding: {args.input}")
    print(f"Model: {model_dir}, Device: {args.device}")

    adata = embed_data(
        adata_or_file=args.input,
        model_dir=model_dir,
        gene_col=args.gene_col,
        max_length=args.max_length,
        batch_size=args.batch_size,
        device=args.device,
    )

    print(f"Embeddings shape: {adata.obsm['X_scGPT'].shape}")
    print(f"Saving to: {args.output}")
    adata.write_h5ad(args.output)
    print("Done.")


if __name__ == "__main__":
    main()
