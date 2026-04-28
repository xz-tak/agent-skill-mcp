#!/usr/bin/env python3
"""Validate the Geneformer conda environment.

Checks Python version, required packages, GPU availability, and Geneformer version.

Usage:
    python ensure_env.py
    python ensure_env.py --quiet   # exit code only: 0=ready, 1=not ready
"""

import argparse
import importlib
import sys


REQUIRED_PACKAGES = [
    "geneformer",
    "torch",
    "transformers",
    "datasets",
    "scanpy",
    "anndata",
    "numpy",
    "pandas",
]

MIN_PYTHON = (3, 10)


def check_python_version():
    """Return (ok, version_str)."""
    current = sys.version_info[:2]
    version_str = f"{current[0]}.{current[1]}"
    return current >= MIN_PYTHON, version_str


def check_package(name):
    """Return (ok, version_str)."""
    try:
        mod = importlib.import_module(name)
        version = getattr(mod, "__version__", "installed")
        return True, str(version)
    except ImportError:
        return False, "not found"


def check_cuda():
    """Return (ok, info_str)."""
    try:
        import torch
        available = torch.cuda.is_available()
        if available:
            device_name = torch.cuda.get_device_name(0)
            return True, f"yes ({device_name})"
        return False, "no"
    except Exception as exc:
        return False, f"error ({exc})"


def main():
    parser = argparse.ArgumentParser(description="Validate Geneformer environment")
    parser.add_argument(
        "--quiet", action="store_true", help="Exit code only: 0=ready, 1=not ready"
    )
    args = parser.parse_args()

    all_ok = True

    python_ok, python_ver = check_python_version()
    if not python_ok:
        all_ok = False

    results = []
    for pkg in REQUIRED_PACKAGES:
        ok, ver = check_package(pkg)
        results.append((pkg, ver, ok))
        if not ok:
            all_ok = False

    cuda_ok, cuda_info = check_cuda()
    if not cuda_ok:
        all_ok = False

    if args.quiet:
        sys.exit(0 if all_ok else 1)

    print(f"{'Package':<16} {'Version':<24} {'Status'}")
    print("-" * 52)
    print(f"{'python':<16} {python_ver:<24} {'OK' if python_ok else 'FAIL (need >=3.10)'}")
    for pkg, ver, ok in results:
        status = "OK" if ok else "MISSING"
        print(f"{pkg:<16} {ver:<24} {status}")
    print(f"{'CUDA':<16} {cuda_info:<24} {'OK' if cuda_ok else 'WARN'}")
    print("-" * 52)
    print(f"Environment: {'READY' if all_ok else 'NOT READY'}")

    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
