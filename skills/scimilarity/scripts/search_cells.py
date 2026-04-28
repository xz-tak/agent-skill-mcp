#!/usr/bin/env python3
"""Cell search pipeline using scimilarity.

Usage:
    # Nearest neighbor search
    python search_cells.py --input data.h5ad --model-dir /path/to/model_v1.1 --mode nearest --k 10000 --output-prefix results

    # Centroid search
    python search_cells.py --input data.h5ad --model-dir /path/to/model_v1.1 --mode centroid --centroid-key query_cells --output-prefix results

    # Cluster search
    python search_cells.py --input data.h5ad --model-dir /path/to/model_v1.1 --mode cluster --cluster-key leiden --output-prefix results

    # Exhaustive search with max distance
    python search_cells.py --input data.h5ad --model-dir /path/to/model_v1.1 --mode exhaustive --max-dist 0.03 --output-prefix results
"""

import argparse


def main():
    parser = argparse.ArgumentParser(description="Search for similar cells using scimilarity")
    parser.add_argument("--input", required=True, help="Input h5ad file path")
    parser.add_argument("--model-dir", default=None, help="Path to scimilarity model directory")
    parser.add_argument("--s3", action="store_true", help="Auto-cache model from S3 (full tier)")
    parser.add_argument("--mode", required=True, choices=["nearest", "centroid", "cluster", "exhaustive"],
                        help="Search mode")
    parser.add_argument("--output-prefix", required=True, help="Output prefix for results files")
    parser.add_argument("--gpu", action="store_true", help="Use GPU for embedding")
    parser.add_argument("--k", type=int, default=10000, help="Number of nearest neighbors (default: 10000)")
    parser.add_argument("--max-dist", type=float, default=None, help="Maximum distance filter")
    parser.add_argument("--centroid-key", type=str, default=None, help="obs column with 0/1 for centroid cells")
    parser.add_argument("--cluster-key", type=str, default=None, help="obs column with cluster labels")
    parser.add_argument("--cluster-label", type=str, default=None, help="Specific cluster to search")
    parser.add_argument("--metadata-filter", type=str, default=None,
                        help="Metadata filter as key=value (e.g., tissue=lung)")
    parser.add_argument("--buffer-size", type=int, default=100000, help="Buffer size for exhaustive search")
    args = parser.parse_args()

    import numpy as np
    import scanpy as sc
    from scimilarity import CellQuery
    from scimilarity.utils import align_dataset, lognorm_counts

    # Resolve model path
    model_dir = args.model_dir
    if args.s3:
        import sys, os
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from s3_cache import get_model_path
        model_dir = get_model_path(tier="full")
    if not model_dir:
        print("ERROR: provide --model-dir or --s3")
        return

    print(f"Loading data: {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Cells: {adata.n_obs}, Genes: {adata.n_vars}")

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # Parse metadata filter
    metadata_filter = None
    if args.metadata_filter:
        key, value = args.metadata_filter.split("=")
        metadata_filter = {key.strip(): value.strip()}

    # Initialize
    load_knn = args.mode in ["nearest", "centroid", "cluster"]
    print(f"Loading model: {model_dir} (knn={load_knn})")
    cq = CellQuery(model_path=model_dir, use_gpu=args.gpu, load_knn=load_knn)

    if args.mode == "nearest":
        print("Aligning and normalizing...")
        adata_aligned = align_dataset(adata, cq.gene_order)
        adata_aligned = lognorm_counts(adata_aligned)
        embeddings = cq.get_embeddings(adata_aligned.X)

        print(f"Searching nearest (k={args.k}, max_dist={args.max_dist})...")
        nn_idxs, nn_dists, metadata = cq.search_nearest(
            embeddings, k=args.k, max_dist=args.max_dist
        )

    elif args.mode == "centroid":
        if not args.centroid_key:
            print("ERROR: --centroid-key required for centroid mode")
            return
        print(f"Searching centroid (key={args.centroid_key})...")
        centroid_emb, nn_idxs, nn_dists, metadata, qc_stats = cq.search_centroid_nearest(
            adata, centroid_key=args.centroid_key, k=args.k, max_dist=args.max_dist
        )
        print(f"  QC stats: {qc_stats}")
        np.save(f"{args.output_prefix}_centroid_embedding.npy", centroid_emb)

    elif args.mode == "cluster":
        if not args.cluster_key:
            print("ERROR: --cluster-key required for cluster mode")
            return
        print(f"Searching cluster centroids (key={args.cluster_key})...")
        centroid_embs, cluster_idx, nn_idxs, nn_dists, metadata = \
            cq.search_cluster_centroids_nearest(
                adata, cluster_key=args.cluster_key,
                cluster_label=args.cluster_label, k=args.k, max_dist=args.max_dist
            )
        print(f"  Clusters: {cluster_idx}")
        np.save(f"{args.output_prefix}_centroid_embeddings.npy", centroid_embs)

    elif args.mode == "exhaustive":
        print("Aligning and normalizing...")
        adata_aligned = align_dataset(adata, cq.gene_order)
        adata_aligned = lognorm_counts(adata_aligned)
        embeddings = cq.get_embeddings(adata_aligned.X)

        max_dist = args.max_dist if args.max_dist else 0.03
        print(f"Exhaustive search (max_dist={max_dist})...")
        nn_idxs, nn_dists, metadata = cq.search_exhaustive(
            embeddings, max_dist=max_dist,
            metadata_filter=metadata_filter, buffer_size=args.buffer_size
        )

    # Save results
    print(f"\nResults: {len(metadata)} hits")
    metadata_path = f"{args.output_prefix}_metadata.csv"
    metadata.to_csv(metadata_path, index=False)
    print(f"  Saved metadata: {metadata_path}")

    # Compile sample-level summary
    if args.mode in ["nearest", "exhaustive"] and len(nn_idxs) > 0:
        for i, (idxs, dists) in enumerate(zip(nn_idxs[:5], nn_dists[:5])):
            print(f"  Query {i}: {len(idxs)} hits, min_dist={min(dists):.4f}")

    print("Done.")


if __name__ == "__main__":
    main()
