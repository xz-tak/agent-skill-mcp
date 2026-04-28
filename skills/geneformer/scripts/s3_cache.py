#!/usr/bin/env python3
"""S3-to-local cache for Geneformer V2 model checkpoints.

Syncs model files from S3 to a local cache directory on first use.
Subsequent calls return the cached path instantly.

Usage as module:
    from s3_cache import get_model_path
    model_path = get_model_path(tier="V2-104M")  # default ~399MB

Usage as CLI:
    python s3_cache.py                          # default: V2-104M (~399MB)
    python s3_cache.py --tier V2-316M           # larger model (~1.2GB)
    python s3_cache.py --tier V2-104M-CLcancer  # cancer continual learning variant
    python s3_cache.py --clear                  # remove local cache
    python s3_cache.py --path-only              # just print the path
"""

import argparse
import os
import subprocess
import sys

S3_BASE = "s3://tec-rnd-sci-dev-gi2/gi2-xz/models/geneformer/"
LOCAL_CACHE = "/tmp/geneformer"

# V2 model tiers (V1 not supported)
TIER_VARIANTS = {
    "V2-104M": "Geneformer-V2-104M",                # ~399MB, DEFAULT
    "V2-316M": "Geneformer-V2-316M",                # ~1.2GB
    "V2-104M-CLcancer": "Geneformer-V2-104M_CLcancer",  # ~399MB, cancer CL variant
}

DEFAULT_TIER = "V2-104M"

# Marker file to detect cached models
MARKER_FILE = "model.safetensors"


def is_cached(tier: str = DEFAULT_TIER) -> bool:
    """Check if the model tier is already cached locally."""
    variant = TIER_VARIANTS.get(tier)
    if variant is None:
        raise ValueError(f"Unknown tier: {tier}. Choose from: {list(TIER_VARIANTS.keys())}")
    return os.path.exists(os.path.join(LOCAL_CACHE, variant, MARKER_FILE))


def sync_from_s3(tier: str = DEFAULT_TIER, verbose: bool = True) -> str:
    """Sync model files from S3 to local cache.

    Parameters
    ----------
    tier : str
        One of "V2-104M", "V2-316M", or "V2-104M-CLcancer".
    verbose : bool
        Print progress messages.

    Returns
    -------
    str
        Local path to the cached model directory.
    """
    variant = TIER_VARIANTS.get(tier)
    if variant is None:
        raise ValueError(f"Unknown tier: {tier}. Choose from: {list(TIER_VARIANTS.keys())}")

    local_dir = os.path.join(LOCAL_CACHE, variant)
    os.makedirs(local_dir, exist_ok=True)

    s3_uri = S3_BASE + variant + "/"
    if verbose:
        print(f"Syncing {tier} from {s3_uri} -> {local_dir}")

    cmd = ["aws", "s3", "sync", s3_uri, local_dir]
    subprocess.run(cmd, check=True)

    if verbose:
        print(f"Model cached at: {local_dir}")
    return local_dir


def get_model_path(tier: str = DEFAULT_TIER, verbose: bool = True) -> str:
    """Get local model path, syncing from S3 if not cached.

    Parameters
    ----------
    tier : str
        "V2-104M" (~399MB, default), "V2-316M" (~1.2GB), or "V2-104M-CLcancer" (~399MB).
    verbose : bool
        Print progress messages.

    Returns
    -------
    str
        Local filesystem path to model directory.

    Examples
    --------
    >>> from s3_cache import get_model_path
    >>> model_path = get_model_path(tier="V2-104M")
    >>> from geneformer import EmbExtractor
    >>> embex = EmbExtractor(model_type="Pretrained")
    >>> embex.extract_embs(model_path, ...)
    """
    if is_cached(tier):
        variant = TIER_VARIANTS[tier]
        local_dir = os.path.join(LOCAL_CACHE, variant)
        if verbose:
            print(f"Using cached model: {local_dir} (tier={tier})")
        return local_dir

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
    parser = argparse.ArgumentParser(description="S3 cache for Geneformer V2 models")
    parser.add_argument(
        "--tier",
        choices=list(TIER_VARIANTS.keys()),
        default=DEFAULT_TIER,
        help=f"Model tier to sync (default: {DEFAULT_TIER})",
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
