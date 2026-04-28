#!/usr/bin/env python3
"""S3-to-local cache for scGPT model files.

Syncs model files from S3 to /tmp/scgpt/ on first use.
Subsequent calls return the cached path instantly.

Usage as module:
    from s3_cache import get_model_path
    model_dir = get_model_path()
    # model_dir = "/tmp/scgpt"

Usage as CLI:
    python s3_cache.py                  # sync model (~210MB)
    python s3_cache.py --clear          # remove local cache
    python s3_cache.py --path-only      # just print path
"""

import argparse
import os
import subprocess
import sys

S3_URI = "s3://tec-rnd-sci-dev-gi2/gi2-xz/models/scgpt/"
LOCAL_CACHE = "/tmp/scgpt"

MODEL_FILES = [
    "best_model.pt",
    "vocab.json",
    "args.json",
]

MARKER_FILE = "best_model.pt"


def is_cached() -> bool:
    """Check if the model is already cached locally."""
    return os.path.exists(os.path.join(LOCAL_CACHE, MARKER_FILE))


def sync_from_s3(verbose: bool = True) -> str:
    """Sync model files from S3 to local cache.

    Returns
    -------
    str
        Local path to the cached model directory.
    """
    os.makedirs(LOCAL_CACHE, exist_ok=True)

    for f in MODEL_FILES:
        local_path = os.path.join(LOCAL_CACHE, f)
        if os.path.exists(local_path):
            if verbose:
                print(f"  [cached] {f}")
            continue

        s3_path = S3_URI + f
        if verbose:
            print(f"  [syncing] {f}")
        cmd = ["aws", "s3", "cp", s3_path, local_path]
        subprocess.run(cmd, check=True)

    if verbose:
        print(f"Model cached at: {LOCAL_CACHE}")
    return LOCAL_CACHE


def get_model_path(verbose: bool = True) -> str:
    """Get local model path, syncing from S3 if not cached.

    Returns
    -------
    str
        Local filesystem path to model directory (/tmp/scgpt).

    Examples
    --------
    >>> from s3_cache import get_model_path
    >>> model_dir = get_model_path()
    >>> import torch
    >>> ckpt = torch.load(f"{model_dir}/best_model.pt", map_location="cpu")
    """
    if is_cached():
        if verbose:
            print(f"Using cached model: {LOCAL_CACHE}")
        return LOCAL_CACHE

    return sync_from_s3(verbose=verbose)


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
    parser = argparse.ArgumentParser(description="S3 cache for scGPT model")
    parser.add_argument("--clear", action="store_true", help="Clear local cache")
    parser.add_argument("--quiet", action="store_true", help="Suppress output")
    parser.add_argument("--path-only", action="store_true", help="Only print the path")
    args = parser.parse_args()

    verbose = not args.quiet and not args.path_only

    if args.clear:
        clear_cache(verbose=verbose)
        return

    path = get_model_path(verbose=verbose)

    if args.path_only:
        print(path)
    elif verbose:
        print(f'\nReady. Use with:\n  model_dir = "{path}"')


if __name__ == "__main__":
    main()
