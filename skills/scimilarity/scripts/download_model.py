#!/usr/bin/env python3
"""Download scimilarity pretrained model from Zenodo.

Usage:
    python download_model.py --output-dir /path/to/destination
    python download_model.py --output-dir ./models --force
"""

import argparse
import os
import sys
import hashlib
from pathlib import Path

ZENODO_RECORD = "10685499"
ZENODO_URL = f"https://zenodo.org/records/{ZENODO_RECORD}/files/model_v1.1.tar.gz"
FILENAME = "model_v1.1.tar.gz"
EXPECTED_MD5 = "546251b7c435f3b1dbe38e2e420ad57f"
EXPECTED_SIZE_GB = 28.2


def download_with_progress(url: str, dest: str, force: bool = False) -> str:
    """Download file with progress bar, supports resume."""
    import requests
    from tqdm import tqdm

    filepath = os.path.join(dest, FILENAME)

    if os.path.exists(filepath) and not force:
        existing_size = os.path.getsize(filepath)
        print(f"File exists: {filepath} ({existing_size / 1e9:.1f} GB)")
        response = input("Resume/overwrite/skip? [r/o/s]: ").strip().lower()
        if response == "s":
            return filepath
        elif response == "o":
            os.remove(filepath)

    headers = {}
    mode = "wb"
    initial_size = 0
    if os.path.exists(filepath) and not force:
        initial_size = os.path.getsize(filepath)
        headers["Range"] = f"bytes={initial_size}-"
        mode = "ab"
        print(f"Resuming from {initial_size / 1e9:.1f} GB")

    print(f"Downloading from: {url}")
    print(f"Destination: {filepath}")

    response = requests.get(f"{url}?download=1", headers=headers, stream=True)
    total_size = int(response.headers.get("content-length", 0)) + initial_size

    with open(filepath, mode) as f, tqdm(
        total=total_size,
        initial=initial_size,
        unit="B",
        unit_scale=True,
        desc=FILENAME,
    ) as pbar:
        for chunk in response.iter_content(chunk_size=8192 * 1024):
            if chunk:
                f.write(chunk)
                pbar.update(len(chunk))

    print(f"Download complete: {filepath}")
    return filepath


def verify_md5(filepath: str, expected: str) -> bool:
    """Verify file MD5 checksum."""
    from tqdm import tqdm

    print("Verifying MD5 checksum...")
    md5 = hashlib.md5()
    size = os.path.getsize(filepath)
    with open(filepath, "rb") as f, tqdm(
        total=size, unit="B", unit_scale=True, desc="Verifying"
    ) as pbar:
        for chunk in iter(lambda: f.read(8192 * 1024), b""):
            md5.update(chunk)
            pbar.update(len(chunk))

    computed = md5.hexdigest()
    if computed == expected:
        print(f"MD5 verified: {computed}")
        return True
    else:
        print(f"MD5 MISMATCH: expected {expected}, got {computed}")
        return False


def extract_archive(filepath: str, dest: str) -> str:
    """Extract tar.gz archive."""
    import tarfile

    print(f"Extracting {filepath}...")
    with tarfile.open(filepath, "r:gz") as tar:
        tar.extractall(path=dest)

    model_dir = os.path.join(dest, "model_v1.1")
    print(f"Extracted to: {model_dir}")
    return model_dir


def verify_model_files(model_dir: str) -> bool:
    """Verify all expected model files exist."""
    required_files = [
        "encoder.ckpt",
        "gene_order.tsv",
        "layer_sizes.json",
        "label_ints.csv",
    ]
    optional_dirs = [
        "annotation",
        "cellsearch",
    ]

    all_ok = True
    for f in required_files:
        path = os.path.join(model_dir, f)
        if os.path.exists(path):
            size = os.path.getsize(path)
            print(f"  [OK] {f} ({size / 1e6:.1f} MB)")
        else:
            print(f"  [MISSING] {f}")
            all_ok = False

    for d in optional_dirs:
        path = os.path.join(model_dir, d)
        if os.path.isdir(path):
            contents = os.listdir(path)
            print(f"  [OK] {d}/ ({len(contents)} items)")
        else:
            print(f"  [MISSING] {d}/")

    return all_ok


def main():
    parser = argparse.ArgumentParser(description="Download scimilarity model from Zenodo")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory to download and extract model into",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-download even if file exists",
    )
    parser.add_argument(
        "--skip-extract",
        action="store_true",
        help="Skip extraction (download only)",
    )
    parser.add_argument(
        "--skip-verify",
        action="store_true",
        help="Skip MD5 verification",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Download
    filepath = download_with_progress(ZENODO_URL, args.output_dir, force=args.force)

    # Verify
    if not args.skip_verify:
        if not verify_md5(filepath, EXPECTED_MD5):
            print("WARNING: MD5 mismatch. File may be corrupted or partially downloaded.")
            sys.exit(1)

    # Extract
    if not args.skip_extract:
        model_dir = extract_archive(filepath, args.output_dir)

        print("\nVerifying model files:")
        if verify_model_files(model_dir):
            print(f"\nModel ready at: {model_dir}")
        else:
            print("\nWARNING: Some files are missing.")
            sys.exit(1)
    else:
        print(f"\nArchive saved at: {filepath}")


if __name__ == "__main__":
    main()
