#!/usr/bin/env python3
"""S3-to-local cache for scimilarity model files.

Syncs model files from S3 to a local cache directory on first use.
Subsequent calls return the cached path instantly.

Usage as module:
    from s3_cache import get_model_path
    model_path = get_model_path(tier="annotation")
    ca = CellAnnotation(model_path=model_path)

Usage as CLI:
    python s3_cache.py                          # default: annotation tier
    python s3_cache.py --tier full              # sync everything (~28GB)
    python s3_cache.py --tier embedding         # just encoder (~250MB)
    python s3_cache.py --clear                  # remove local cache
"""

import argparse
import os
import subprocess
import sys

S3_URI = "s3://tec-rnd-sci-dev-gi2/gi2-xz/models/scimilarity/model_v1.1/"
LOCAL_CACHE = "/tmp/scimilarity/model_v1.1"

# Files per tier (cumulative)
TIER_FILES = {
    "embedding": [
        "encoder.ckpt",
        "gene_order.tsv",
        "layer_sizes.json",
        "label_ints.csv",
    ],
    "annotation": [
        "encoder.ckpt",
        "gene_order.tsv",
        "layer_sizes.json",
        "label_ints.csv",
        "annotation/labelled_kNN.bin",
        "annotation/reference_labels.tsv",
    ],
    "full": None,  # sync everything
}

# Marker file used to check if a tier is already cached
TIER_MARKERS = {
    "embedding": "encoder.ckpt",
    "annotation": "annotation/labelled_kNN.bin",
    "full": "cellsearch/full_kNN.bin",
}


def is_cached(tier: str = "annotation") -> bool:
    """Check if the model tier is already cached locally."""
    marker = TIER_MARKERS.get(tier, "encoder.ckpt")
    return os.path.exists(os.path.join(LOCAL_CACHE, marker))


def sync_from_s3(tier: str = "annotation", verbose: bool = True) -> str:
    """Sync model files from S3 to local cache.

    Parameters
    ----------
    tier : str
        One of "embedding", "annotation", or "full".
    verbose : bool
        Print progress messages.

    Returns
    -------
    str
        Local path to the cached model directory.
    """
    os.makedirs(LOCAL_CACHE, exist_ok=True)

    files = TIER_FILES.get(tier)

    if files is None:
        # Full sync
        if verbose:
            print(f"Syncing full model from {S3_URI} -> {LOCAL_CACHE}")
        cmd = ["aws", "s3", "sync", S3_URI, LOCAL_CACHE]
        subprocess.run(cmd, check=True)
    else:
        # Selective sync — copy individual files
        for f in files:
            local_path = os.path.join(LOCAL_CACHE, f)
            if os.path.exists(local_path):
                if verbose:
                    print(f"  [cached] {f}")
                continue

            os.makedirs(os.path.dirname(local_path) or LOCAL_CACHE, exist_ok=True)
            s3_path = S3_URI + f
            if verbose:
                print(f"  [syncing] {f}")
            cmd = ["aws", "s3", "cp", s3_path, local_path]
            subprocess.run(cmd, check=True)

    if verbose:
        print(f"Model cached at: {LOCAL_CACHE}")
    return LOCAL_CACHE


def get_model_path(tier: str = "annotation", verbose: bool = True) -> str:
    """Get local model path, syncing from S3 if not cached.

    Parameters
    ----------
    tier : str
        "embedding" (~250MB), "annotation" (~9GB), or "full" (~28GB).
    verbose : bool
        Print progress messages.

    Returns
    -------
    str
        Local filesystem path to model directory.

    Examples
    --------
    >>> from s3_cache import get_model_path
    >>> model_path = get_model_path(tier="annotation")
    >>> from scimilarity import CellAnnotation
    >>> ca = CellAnnotation(model_path=model_path)
    """
    if is_cached(tier):
        if verbose:
            print(f"Using cached model: {LOCAL_CACHE} (tier={tier})")
        return LOCAL_CACHE

    return sync_from_s3(tier=tier, verbose=verbose)


def clear_cache(verbose: bool = True):
    """Remove the local cache directory."""
    import shutil

    if os.path.exists(LOCAL_CACHE):
        shutil.rmtree(LOCAL_CACHE)
        if verbose:
            print(f"Cleared cache: {LOCAL_CACHE}")
    else:
        if verbose:
            print("No cache to clear.")


def main():
    parser = argparse.ArgumentParser(description="S3 cache for scimilarity model")
    parser.add_argument(
        "--tier",
        choices=["embedding", "annotation", "full"],
        default="annotation",
        help="Model tier to sync (default: annotation)",
    )
    parser.add_argument("--clear", action="store_true", help="Clear local cache")
    parser.add_argument("--quiet", action="store_true", help="Suppress output")
    parser.add_argument("--path-only", action="store_true", help="Only print the path")
    args = parser.parse_args()

    verbose = not args.quiet and not args.path_only

    if args.clear:
        clear_cache(verbose=verbose)
        return

    path = get_model_path(tier=args.tier, verbose=verbose)

    if args.path_only:
        print(path)
    elif verbose:
        print(f"\nReady. Use with:")
        print(f'  model_path = "{path}"')


if __name__ == "__main__":
    main()
