#!/usr/bin/env python3
"""
Export h5ad metadata (obs) as TSV for R import.
Fallback for when hdf5r can't parse complex h5ad categorical encodings.

Usage:
    python export_metadata.py <input.h5ad> <output_metadata.tsv>
"""
import sys
import anndata as ad

h5ad_path = sys.argv[1]
output_path = sys.argv[2]

print(f"Reading obs from {h5ad_path} (backed mode)...", flush=True)
adata = ad.read_h5ad(h5ad_path, backed='r')
meta = adata.obs.copy()
print(f"  Shape: {meta.shape}", flush=True)
meta.to_csv(output_path, sep="\t")
print(f"  Saved to {output_path}", flush=True)
