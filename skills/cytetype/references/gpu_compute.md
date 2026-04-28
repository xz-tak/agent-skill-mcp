# GPU-Aware Compute Patterns for CyteType Clustering

## Table of Contents
1. [GPU VRAM Detection](#1-gpu-vram-detection)
2. [3-Tier Backend Selection](#2-3-tier-backend-selection)
3. [Multi-GPU LocalCUDACluster Setup](#3-multi-gpu-localcudacluster-setup)
4. [h5ad ↔ zarr Conversion](#4-h5ad--zarr-conversion)
5. [Zarr-Based Multi-GPU Loading (Official rsc Pattern)](#5-zarr-based-multi-gpu-loading)
6. [X-Swap Pattern (Leiden Only)](#6-x-swap-pattern-leiden-only)
7. [Neighbors Auto-Compute](#7-neighbors-auto-compute)
8. [Leiden Resolution Iteration](#8-leiden-resolution-iteration)
9. [rank_genes_groups Method Selection](#9-rank_genes_groups-method-selection)
10. [Cleanup](#10-cleanup)
11. [Complete Workflow Example](#11-complete-workflow-example)

---

## 1. GPU VRAM Detection

Strict <10% VRAM usage threshold. If no GPU qualifies, fall back to scanpy CPU.

```python
import subprocess
import logging

def detect_available_gpus(max_vram_usage_pct=10):
    """Detect GPUs with VRAM usage below threshold.

    Returns list of dicts with keys: index, mem_used_mb, mem_total_mb, utilization_pct.
    Returns empty list if no GPUs qualify or nvidia-smi unavailable.
    """
    try:
        result = subprocess.run(
            ['nvidia-smi', '--query-gpu=index,memory.used,memory.total,utilization.gpu',
             '--format=csv,noheader,nounits'],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode != 0:
            logging.warning("nvidia-smi failed, no GPU available")
            return []

        idle_gpus = []
        for line in result.stdout.strip().split('\n'):
            parts = [x.strip() for x in line.split(',')]
            if len(parts) != 4:
                continue
            idx, mem_used, mem_total, util = parts
            usage_pct = float(mem_used) / float(mem_total) * 100
            if usage_pct < max_vram_usage_pct:
                idle_gpus.append({
                    'index': idx,
                    'mem_used_mb': int(float(mem_used)),
                    'mem_total_mb': int(float(mem_total)),
                    'utilization_pct': int(float(util)),
                })
                logging.info(f"GPU {idx}: {usage_pct:.1f}% VRAM used — eligible")
            else:
                logging.info(f"GPU {idx}: {usage_pct:.1f}% VRAM used — skipped (>{max_vram_usage_pct}%)")

        return idle_gpus

    except FileNotFoundError:
        logging.warning("nvidia-smi not found, no GPU available")
        return []
    except subprocess.TimeoutExpired:
        logging.warning("nvidia-smi timed out, no GPU available")
        return []
```

## 2. 3-Tier Backend Selection

```python
def select_compute_backend(idle_gpus):
    """Select compute backend based on available GPUs and packages.

    Returns: ('multi_gpu', idle_gpus) | ('single_gpu', [gpu]) | ('cpu', [])
    """
    rsc_available = False
    dask_cuda_available = False

    try:
        import rapids_singlecell as rsc
        rsc_available = True
    except ImportError:
        logging.warning("rapids_singlecell not available")

    try:
        import dask_cuda
        dask_cuda_available = True
    except ImportError:
        logging.info("dask_cuda not available, multi-GPU disabled")

    if len(idle_gpus) >= 2 and rsc_available and dask_cuda_available:
        logging.info(f"Backend: multi-GPU ({len(idle_gpus)} GPUs)")
        return ('multi_gpu', idle_gpus)

    if len(idle_gpus) >= 1 and rsc_available:
        logging.info(f"Backend: single-GPU (GPU {idle_gpus[0]['index']})")
        return ('single_gpu', [idle_gpus[0]])

    logging.info("Backend: scanpy CPU")
    return ('cpu', [])
```

**CPU fallback warning** — if dataset has >500K cells and backend is CPU:
```python
if backend == 'cpu' and adata.n_obs > 500_000:
    # Use AskUserQuestion to warn user:
    # "No GPU with <10% VRAM available. CPU on {adata.n_obs:,} cells will be very slow.
    #  (1) Proceed with CPU scanpy
    #  (2) Abort and wait for GPU availability"
```

## 3. Multi-GPU LocalCUDACluster Setup

```python
import rapids_singlecell as rsc
from dask_cuda import LocalCUDACluster
from dask.distributed import Client

def setup_multi_gpu(idle_gpus):
    """Initialize LocalCUDACluster with idle GPUs. Returns (client, cluster)."""
    gpu_ids = ",".join(g['index'] for g in idle_gpus)
    logging.info(f"Initializing LocalCUDACluster on GPUs: {gpu_ids}")

    cluster = LocalCUDACluster(
        CUDA_VISIBLE_DEVICES=gpu_ids,
        threads_per_worker=10,
        protocol="tcp",
        rmm_pool_size="10GB",
        rmm_maximum_pool_size="20GB",  # per-worker ceiling, ~87% of 23GB GPU
        rmm_allocator_external_lib_list="cupy",
    )
    client = Client(cluster)
    logging.info(f"Dask client ready: {client.dashboard_link}")
    return client, cluster
```

## 4. h5ad ↔ zarr Conversion

Zarr format enables dask-chunked loading, which is required for distributing large
expression matrices across multiple GPUs without staging through one GPU.

### h5ad → zarr (preserves full object: X, obs, var, obsm, obsp, uns, layers)

```python
import anndata as ad
import scanpy as sc
import logging
from pathlib import Path

def h5ad_to_zarr(h5ad_path, zarr_path=None, overwrite=False):
    """Convert h5ad to zarr format, preserving all AnnData components.

    Parameters
    ----------
    h5ad_path : str
        Path to input h5ad file.
    zarr_path : str, optional
        Path to output zarr store. Defaults to same name with .zarr extension.
    overwrite : bool
        If True, overwrite existing zarr store.

    Returns
    -------
    str : path to zarr store
    """
    h5ad_path = Path(h5ad_path)
    if zarr_path is None:
        zarr_path = h5ad_path.with_suffix('.zarr')
    zarr_path = Path(zarr_path)

    if zarr_path.exists() and not overwrite:
        logging.info(f"Zarr store already exists: {zarr_path}")
        return str(zarr_path)

    logging.info(f"Converting {h5ad_path} → {zarr_path}")
    adata = sc.read_h5ad(str(h5ad_path))
    logging.info(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # anndata.write_zarr preserves X, obs, var, obsm, obsp, uns, layers
    adata.write_zarr(str(zarr_path))
    logging.info(f"Zarr store written: {zarr_path}")
    return str(zarr_path)
```

### zarr → h5ad (preserves full object)

```python
def zarr_to_h5ad(zarr_path, h5ad_path=None, compression='gzip'):
    """Convert zarr store back to h5ad, preserving all AnnData components.

    Parameters
    ----------
    zarr_path : str
        Path to input zarr store.
    h5ad_path : str, optional
        Path to output h5ad file. Defaults to same name with .h5ad extension.
    compression : str
        Compression for h5ad. Default 'gzip'.

    Returns
    -------
    str : path to h5ad file
    """
    zarr_path = Path(zarr_path)
    if h5ad_path is None:
        h5ad_path = zarr_path.with_suffix('.h5ad')
    h5ad_path = Path(h5ad_path)

    logging.info(f"Converting {zarr_path} → {h5ad_path}")
    adata = ad.read_zarr(str(zarr_path))
    logging.info(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    adata.write_h5ad(str(h5ad_path), compression=compression)
    logging.info(f"h5ad written: {h5ad_path}")
    return str(h5ad_path)
```

### Load zarr with dask-chunked obs metadata (for inspection without full load)

```python
import zarr

def inspect_zarr(zarr_path):
    """Quick inspection of zarr store without loading X."""
    f = zarr.open(str(zarr_path), mode='r')
    info = {
        'X_shape': f['X'].attrs.get('shape', f['X'].shape),
        'obs_columns': list(f['obs'].keys()),
        'var_columns': list(f['var'].keys()),
        'obsm_keys': list(f['obsm'].keys()) if 'obsm' in f else [],
        'obsp_keys': list(f['obsp'].keys()) if 'obsp' in f else [],
        'uns_keys': list(f['uns'].keys()) if 'uns' in f else [],
        'layers_keys': list(f['layers'].keys()) if 'layers' in f else [],
    }
    return info
```

## 5. Zarr-Based Multi-GPU Loading

**This is the official rsc multi-GPU pattern.** It loads X as a dask array with
chunking so each GPU gets a portion of the data without staging through one GPU.

This is **required** when the expression matrix X exceeds single-GPU VRAM (e.g.,
2.98M cells x 38.8K genes = ~36 GB on 23 GB GPUs).

```python
import anndata as ad
import zarr
import numpy as np
from packaging.version import parse as parse_version

def load_zarr_for_multi_gpu(zarr_path, sparse_chunk_size=50_000):
    """Load zarr store as AnnData with dask-chunked X for multi-GPU distribution.

    This follows the official rsc multi-GPU pattern from:
    https://rapids-singlecell.readthedocs.io/en/latest/notebooks/06-multi_gpu_show.html

    Parameters
    ----------
    zarr_path : str
        Path to zarr store (created by adata.write_zarr()).
    sparse_chunk_size : int
        Number of rows per dask chunk. Default 50,000 (from official rsc docs).
        Adjust based on available GPU VRAM per worker.

    Returns
    -------
    AnnData with X as a lazy dask array, obs/var/obsm/obsp/uns loaded eagerly.
    """
    # Select correct dask reader based on anndata version
    if parse_version(ad.__version__) < parse_version("0.12.0rc1"):
        from anndata.experimental import read_elem_as_dask as read_dask
    else:
        from anndata.experimental import read_elem_lazy as read_dask

    f = zarr.open(str(zarr_path), mode='r')

    # Read X as dask array (lazy, chunked — never materializes on one GPU)
    X = f["X"]
    shape = X.attrs.get("shape", X.shape)
    X_dask = read_dask(X, (sparse_chunk_size, shape[1]))

    # Read metadata eagerly (small)
    obs = ad.io.read_elem(f["obs"])
    var = ad.io.read_elem(f["var"])

    # Build AnnData with dask X
    adata = ad.AnnData(X=X_dask, obs=obs, var=var)

    # Load obsm, obsp, uns, layers if present
    if "obsm" in f:
        for key in f["obsm"].keys():
            adata.obsm[key] = ad.io.read_elem(f["obsm"][key])
    if "obsp" in f:
        for key in f["obsp"].keys():
            adata.obsp[key] = ad.io.read_elem(f["obsp"][key])
    if "uns" in f:
        for key in f["uns"].keys():
            try:
                adata.uns[key] = ad.io.read_elem(f["uns"][key])
            except Exception:
                logging.warning(f"Could not read uns['{key}'], skipping")
    if "layers" in f:
        for key in f["layers"].keys():
            adata.layers[key] = read_dask(f["layers"][key], (sparse_chunk_size, shape[1]))

    logging.info(f"Loaded zarr: {adata.n_obs:,} x {adata.n_vars:,}, X is dask array")
    return adata


def transfer_zarr_to_multi_gpu(adata, n_gpus=8):
    """Transfer dask-chunked AnnData to multi-GPU, distribute across workers.

    Call AFTER setup_multi_gpu() and load_zarr_for_multi_gpu().
    """
    rsc.get.anndata_to_GPU(adata)
    adata.X = adata.X.persist()  # distribute chunks across GPU workers
    adata.X.compute_chunk_sizes()

    # Optional: rechunk evenly across workers for balanced compute
    n_rows = adata.shape[0]
    n_cols = adata.shape[1]
    rows_per_worker = (n_rows + n_gpus - 1) // n_gpus
    adata.X = adata.X.rechunk((rows_per_worker, n_cols)).persist()
    adata.X.compute_chunk_sizes()

    logging.info(f"X distributed across {n_gpus} GPUs, "
                 f"~{n_rows // n_gpus:,} rows per worker")
```

## 6. X-Swap Pattern (Leiden Only)

For leiden clustering, X is not needed — only obsp connectivities (~700 MB for 3M cells).
When loading from h5ad (not zarr), swap out X before GPU transfer to avoid OOM.

```python
import scipy.sparse as sp

def transfer_to_gpu_without_x(adata):
    """Transfer only obsp to GPU, swapping out X to avoid OOM.
    Returns (X_backup, layers_backup) for later restoration."""
    X_backup = adata.X
    layers_backup = dict(adata.layers)
    adata.X = sp.csr_matrix(adata.shape, dtype=np.float32)  # empty placeholder
    adata.layers.clear()
    logging.info(f"X swapped out. obsp nnz={adata.obsp['connectivities'].nnz:,}")

    rsc.get.anndata_to_GPU(adata)

    # Guard: verify obsp stayed on CPU (leiden's _create_graph_dask expects scipy)
    import scipy.sparse as sp_check
    assert sp_check.issparse(adata.obsp['connectivities']), \
        "obsp converted to GPU unexpectedly — X-swap pattern requires obsp on CPU"

    logging.info("Connectivities on GPU (X excluded)")
    return X_backup, layers_backup


def restore_x_after_gpu(adata, X_backup, layers_backup):
    """Restore X and layers after GPU computation."""
    rsc.get.anndata_to_CPU(adata)
    adata.X = X_backup
    for k, v in layers_backup.items():
        adata.layers[k] = v
    logging.info("X restored from backup")
```

## 7. Neighbors Auto-Compute

Check for connectivities. If absent, auto-detect PCA representation.

```python
def ensure_neighbors(adata, backend):
    """Compute neighbors if connectivities missing. Auto-detects PCA key."""
    if 'connectivities' in adata.obsp:
        logging.info("Neighbors already computed, skipping")
        return

    use_rep = None
    for key in ['X_pca_harmony', 'X_pca']:
        if key in adata.obsm:
            use_rep = key
            logging.info(f"Using {key} for neighbors computation")
            break

    if use_rep is None:
        raise ValueError("No suitable PCA representation found in adata.obsm. "
                         "Available keys: " + str(list(adata.obsm.keys())))

    if backend in ('multi_gpu', 'single_gpu'):
        rsc.pp.neighbors(adata, use_rep=use_rep)
    else:
        import scanpy as sc
        sc.pp.neighbors(adata, use_rep=use_rep)

    logging.info(f"Neighbors computed using {use_rep}")
```

## 8. Leiden Resolution Iteration

Start at 1.0, decrease by 0.1, floor at 0.2. Report summary stats each round.

```python
import numpy as np

def leiden_iteration(adata, backend, key_added='cytetype_leiden'):
    """Iterate leiden resolution until min cluster > 100 cells.

    Resolution: 1.0 → 0.9 → 0.8 → ... floor 0.2
    Reports: resolution, n_clusters, min/median/max cluster size each round.
    Returns: final resolution used.
    """
    resolution = 1.0
    floor = 0.2

    while resolution >= floor:
        if backend == 'multi_gpu':
            rsc.tl.leiden(adata, resolution=resolution, random_state=42,
                          key_added=key_added, use_dask=True)
        elif backend == 'single_gpu':
            rsc.tl.leiden(adata, resolution=resolution, random_state=42,
                          key_added=key_added)
        else:
            import scanpy as sc
            sc.tl.leiden(adata, resolution=resolution, random_state=42,
                         key_added=key_added)

        sizes = adata.obs[key_added].value_counts()
        n_clusters = len(sizes)
        min_size = sizes.min()
        median_size = int(np.median(sizes.values))
        max_size = sizes.max()
        min_cluster = sizes.idxmin()

        report = (
            f"Resolution {resolution:.2f}: {n_clusters} clusters | "
            f"min={min_size} (cluster {min_cluster}) | "
            f"median={median_size} | max={max_size}"
        )
        logging.info(report)
        print(report)

        if min_size > 100:
            return resolution

        resolution = round(resolution - 0.1, 1)

    return resolution
```

## 9. rank_genes_groups Method Selection

Critical: `wilcoxon` raises ValueError with Dask arrays. Must use `wilcoxon_binned` for multi-GPU.

**For multi-GPU: load from zarr** to get dask-chunked X, then use `wilcoxon_binned`.
**For single-GPU: use `wilcoxon`** (no dask arrays).
**For CPU: use `wilcoxon`**.

```python
def run_rank_genes_groups(adata, groupby, backend):
    """Run rank_genes_groups with correct method for backend.

    IMPORTANT: For multi-GPU, adata.X must be a dask array (loaded via zarr).
    If X is a scipy sparse matrix (loaded via h5ad), convert to zarr first.
    """
    if backend == 'multi_gpu':
        # wilcoxon_binned supports Dask arrays — required for multi-GPU
        rsc.tl.rank_genes_groups(adata, groupby=groupby,
                                 method='wilcoxon_binned', pts=True)
        logging.info("rank_genes_groups computed (rsc, wilcoxon_binned, multi-GPU)")

    elif backend == 'single_gpu':
        rsc.tl.rank_genes_groups(adata, groupby=groupby,
                                 method='wilcoxon', pts=True)
        logging.info("rank_genes_groups computed (rsc, wilcoxon, single-GPU)")

    else:
        import scanpy as sc
        sc.tl.rank_genes_groups(adata, groupby=groupby,
                                method='wilcoxon', pts=True)
        logging.info("rank_genes_groups computed (scanpy, wilcoxon, CPU)")
```

## 10. Cleanup

Always transfer data back to CPU and close dask resources.

```python
def cleanup_gpu(adata, backend, client=None, cluster=None):
    """Transfer data to CPU and close GPU resources."""
    if backend in ('multi_gpu', 'single_gpu'):
        try:
            rsc.get.anndata_to_CPU(adata)
            logging.info("AnnData transferred back to CPU")
        except Exception as e:
            logging.warning(f"anndata_to_CPU failed (may already be on CPU): {e}")

    if backend == 'multi_gpu':
        try:
            import cupy as cp
            cp.get_default_memory_pool().free_all_blocks()
            cp.get_default_pinned_memory_pool().free_all_blocks()
        except Exception:
            pass
        if client is not None:
            client.close()
            logging.info("Dask client closed")
        if cluster is not None:
            cluster.close()
            logging.info("LocalCUDACluster closed")
```

## 11. Complete Workflow Example

Full pipeline: h5ad → zarr → multi-GPU leiden + rank_genes_groups → h5ad.

```python
import logging
import rapids_singlecell as rsc
import scanpy as sc

# --- Step 1: Convert h5ad to zarr (one-time) ---
zarr_path = h5ad_to_zarr("integration_scanpy.h5ad")

# --- Step 2: Detect GPUs and setup cluster ---
idle_gpus = detect_available_gpus(max_vram_usage_pct=10)
backend, selected_gpus = select_compute_backend(idle_gpus)

client, dask_cluster = None, None
if backend == 'multi_gpu':
    client, dask_cluster = setup_multi_gpu(selected_gpus)

# --- Step 3: Load from zarr with dask chunking ---
if backend == 'multi_gpu':
    adata = load_zarr_for_multi_gpu(zarr_path, sparse_chunk_size=50_000)
    transfer_zarr_to_multi_gpu(adata, n_gpus=len(selected_gpus))
elif backend == 'single_gpu':
    adata = sc.read_h5ad("integration_scanpy.h5ad")
    rsc.get.anndata_to_GPU(adata)
else:
    adata = sc.read_h5ad("integration_scanpy.h5ad")

# --- Step 4: Ensure neighbors ---
ensure_neighbors(adata, backend)

# --- Step 5: Leiden iteration ---
final_resolution = leiden_iteration(adata, backend, key_added='cytetype_leiden')

# --- Step 6: rank_genes_groups (multi-GPU uses wilcoxon_binned on dask X) ---
run_rank_genes_groups(adata, groupby='cytetype_leiden', backend=backend)

# --- Step 7: Cleanup and save ---
cleanup_gpu(adata, backend, client, dask_cluster)
adata.write_h5ad("integration_cytetype.h5ad", compression="gzip")
```
