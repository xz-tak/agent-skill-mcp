#!/usr/bin/env python3
"""
Validate that a file is a valid AnnData h5ad object.

Usage:
    python validate_h5ad.py <file.h5ad>

Exit codes:
    0 - Valid AnnData object
    1 - File not found or not a valid h5ad file
    2 - Invalid arguments

Output:
    JSON with validation results and structure summary
"""

import sys
import json
import os


def validate_h5ad(file_path: str) -> dict:
    """
    Validate an h5ad file and return structure information.

    Returns:
        dict with:
        - valid: bool
        - error: str or None
        - structure: dict with AnnData structure info (if valid)
    """
    result = {
        "valid": False,
        "file_path": file_path,
        "error": None,
        "structure": None
    }

    # Check file exists
    if not os.path.exists(file_path):
        result["error"] = f"File not found: {file_path}"
        return result

    # Check file extension
    if not file_path.lower().endswith('.h5ad'):
        result["error"] = f"File does not have .h5ad extension: {file_path}"
        return result

    # Try to import anndata
    try:
        import anndata as ad
    except ImportError:
        result["error"] = "anndata package not installed. Install with: pip install anndata"
        return result

    # Try to load the h5ad file
    try:
        adata = ad.read_h5ad(file_path)
    except Exception as e:
        result["error"] = f"Failed to load h5ad file: {str(e)}"
        return result

    # File is valid - extract structure info
    result["valid"] = True

    structure = {
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
        "obs_columns": list(adata.obs.columns),
        "var_columns": list(adata.var.columns),
        "layers": list(adata.layers.keys()) if adata.layers else [],
        "obsm_keys": list(adata.obsm.keys()) if adata.obsm else [],
        "varm_keys": list(adata.varm.keys()) if adata.varm else [],
        "obsp_keys": list(adata.obsp.keys()) if adata.obsp else [],
        "uns_keys": list(adata.uns.keys()) if adata.uns else [],
        "has_raw": adata.raw is not None,
        "X_dtype": str(adata.X.dtype) if adata.X is not None else None,
        "X_shape": list(adata.X.shape) if adata.X is not None else None,
    }

    # Check if X appears to be normalized or raw counts
    if adata.X is not None:
        import numpy as np
        # Sample some values to check
        if hasattr(adata.X, 'data'):  # sparse
            sample = adata.X.data[:min(10000, len(adata.X.data))]
        else:  # dense
            sample = adata.X.flatten()[:10000]

        sample = np.array(sample)
        max_val = float(np.max(sample)) if len(sample) > 0 else 0
        is_integer = np.allclose(sample, np.round(sample), atol=1e-6) if len(sample) > 0 else False

        structure["X_max_value"] = max_val
        structure["X_appears_normalized"] = max_val < 20 and not is_integer
        structure["X_appears_counts"] = is_integer and max_val > 10

    result["structure"] = structure
    return result


def main():
    if len(sys.argv) < 2:
        print(json.dumps({
            "valid": False,
            "error": "No file path provided. Usage: python validate_h5ad.py <file.h5ad>"
        }))
        sys.exit(2)

    file_path = sys.argv[1]
    result = validate_h5ad(file_path)

    # Print JSON result
    print(json.dumps(result, indent=2))

    # Exit with appropriate code
    if result["valid"]:
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
