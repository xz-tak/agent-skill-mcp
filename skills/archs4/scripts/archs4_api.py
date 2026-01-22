"""
Simple Python client for the ARCHS4 /sigpy API.

Each function wraps a specific endpoint and returns the result
as one or more pandas.DataFrame objects where possible.

Endpoints covered (see ARCHS4 help/API docs):
- GET  /sigpy/status
- GET  /sigpy/meta/quicksearch
- GET  or POST /sigpy/meta/genes
- GET  or POST /sigpy/data/samples
- GET  /sigpy/data/samples/status/<task_id>
- GET  /sigpy/data/samples/download/<task_id>
- POST /sigpy/data/knn/signature
- POST /sigpy/data/correlation
- POST /sigpy/data/diffexp
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional

import io
import os
import zipfile

import pandas as pd
import requests

BASE_URL = "https://maayanlab.cloud"


# --------- low-level helpers -------------------------------------------------


def _get_json(path: str, params: Optional[Dict[str, Any]] = None, timeout: int = 60) -> Any:
    """Internal helper to issue a GET request and parse JSON."""
    url = f"{BASE_URL}{path}"
    resp = requests.get(url, params=params, timeout=timeout)
    resp.raise_for_status()
    return resp.json()


def _post_json(path: str, payload: Optional[Dict[str, Any]] = None, timeout: int = 60) -> Any:
    """Internal helper to issue a POST request and parse JSON."""
    url = f"{BASE_URL}{path}"
    resp = requests.post(url, json=payload, timeout=timeout)
    resp.raise_for_status()
    return resp.json()


# --------- Test & service status --------------------------------------------


def get_service_status() -> pd.DataFrame:
    """
    Wraps: GET /sigpy/status

    Returns
    -------
    pandas.DataFrame
        One-row DataFrame with the service status and last update time.
    """
    data = _get_json("/sigpy/status")
    # Expected to be a small dict like {"status": "...", "last_update": "..."}
    return pd.DataFrame([data])


# --------- Metadata endpoints ------------------------------------------------


def quicksearch_metadata(query: str, species: str = "human") -> pd.DataFrame:
    """
    Wraps: GET /sigpy/meta/quicksearch

    Parameters
    ----------
    query : str
        Search term (e.g. GEO accession ID, tissue, cell type).
    species : {"human", "mouse"}, optional
        Species (defaults to "human").

    Returns
    -------
    pandas.DataFrame
        DataFrame constructed from the JSON response. If the response has a
        top-level "samples" list, that list is returned as rows; otherwise
        the full JSON is normalized.
    """
    params = {"query": query, "species": species}
    data = _get_json("/sigpy/meta/quicksearch", params=params)

    if isinstance(data, list):
        return pd.DataFrame(data)

    if isinstance(data, dict):
        if "samples" in data and isinstance(data["samples"], list):
            return pd.DataFrame(data["samples"])
        # generic flattening fallback
        return pd.json_normalize(data)

    # Last resort, just wrap the object in a single-row DataFrame
    return pd.DataFrame([{"result": data}])


def list_genes(species: str = "human", method: str = "GET") -> pd.DataFrame:
    """
    Wraps: GET/POST /sigpy/meta/genes

    Parameters
    ----------
    species : {"human", "mouse"}, optional
        Species to list genes for.
    method : {"GET", "POST"}, optional
        Whether to use the GET or POST version of the endpoint.

    Returns
    -------
    pandas.DataFrame
        DataFrame with one column, 'gene', containing gene symbols.
    """
    method = method.upper()
    if method == "GET":
        data = _get_json("/sigpy/meta/genes", params={"species": species})
    elif method == "POST":
        data = _post_json("/sigpy/meta/genes", payload={"species": species})
    else:
        raise ValueError("method must be 'GET' or 'POST'")

    # According to the docs this is a flat list of gene symbols.
    if isinstance(data, dict) and "genes" in data:
        genes = data["genes"]
    else:
        genes = data

    return pd.DataFrame({"gene": list(genes)})


# --------- Data extraction: samples -----------------------------------------


def request_sample_data(
    gsm_ids: Iterable[str],
    species: str = "human",
    use_post: bool = True,
) -> pd.DataFrame:
    """
    Wraps: GET/POST /sigpy/data/samples

    Parameters
    ----------
    gsm_ids : iterable of str
        GSM sample accessions to request (max 10,000).
    species : {"human", "mouse"}, optional
        Species (defaults to "human").
    use_post : bool, optional
        If True, use POST with JSON body; otherwise use GET with
        comma-separated GSM IDs.

    Returns
    -------
    pandas.DataFrame
        One-row DataFrame with whatever JSON the endpoint returns, which
        typically includes a 'task_id' that you can poll with
        ``check_samples_task_status``.
    """
    gsm_ids = list(gsm_ids)
    if not gsm_ids:
        raise ValueError("gsm_ids must contain at least one GSM accession")

    if use_post:
        payload = {"gsm_ids": gsm_ids, "species": species}
        data = _post_json("/sigpy/data/samples", payload=payload)
    else:
        params = {"gsm_ids": ",".join(gsm_ids), "species": species}
        data = _get_json("/sigpy/data/samples", params=params)

    return pd.DataFrame([data])


def check_samples_task_status(task_id: str) -> pd.DataFrame:
    """
    Wraps: GET /sigpy/data/samples/status/<task_id>

    Parameters
    ----------
    task_id : str
        Task ID returned by request_sample_data.

    Returns
    -------
    pandas.DataFrame
        One-row DataFrame describing the task status.
    """
    path = f"/sigpy/data/samples/status/{task_id}"
    data = _get_json(path)
    return pd.DataFrame([data])


def download_samples_zip(task_id: str, out_file: str = "matrix.zip") -> str:
    """
    Wraps: GET /sigpy/data/samples/download/<task_id>

    This endpoint returns a ZIP archive, not JSON, so this helper just
    downloads it and saves it to disk.

    Parameters
    ----------
    task_id : str
        Task ID returned by request_sample_data.
    out_file : str, optional
        Path to save the ZIP file.

    Returns
    -------
    str
        Path to the downloaded ZIP file.
    """
    url = f"{BASE_URL}/sigpy/data/samples/download/{task_id}"
    resp = requests.get(url, timeout=300)
    resp.raise_for_status()
    with open(out_file, "wb") as f:
        f.write(resp.content)
    return out_file


def extract_first_matrix_from_zip(zip_bytes: bytes) -> pd.DataFrame:
    """
    Utility: given raw ZIP bytes from the samples download endpoint,
    return the first tab-delimited file inside as a DataFrame.

    This is optional; use it if you prefer to work in-memory instead
    of saving the ZIP to disk first.

    Parameters
    ----------
    zip_bytes : bytes
        Content returned by the /samples/download endpoint.

    Returns
    -------
    pandas.DataFrame
        DataFrame parsed from the first .txt/.tsv/.tab file found.

    Notes
    -----
    The exact file naming inside the ZIP may change; this function
    simply grabs the first tabular text file it finds.
    """
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
        # heuristic: pick the first text-like file
        candidates = [
            name
            for name in zf.namelist()
            if name.lower().endswith((".txt", ".tsv", ".tab"))
        ]
        if not candidates:
            raise ValueError("No tabular text file found inside ZIP archive")

        with zf.open(candidates[0]) as fh:
            return pd.read_csv(fh, sep="\t")


# --------- Data query: k-NN signature search --------------------------------


def knn_signature_from_gene_sets(
    up_genes: Iterable[str],
    down_genes: Optional[Iterable[str]] = None,
    species: str = "human",
    k: int = 10,
    signame: str = "gene-set similarity search",
) -> pd.DataFrame:
    """
    Wraps: POST /sigpy/data/knn/signature using up/down marker genes.

    Parameters
    ----------
    up_genes : iterable of str
        Genes characteristically up-regulated.
    down_genes : iterable of str, optional
        Genes characteristically down-regulated (can be empty / None).
    species : {"human", "mouse"}, optional
        Species (defaults to "human").
    k : int, optional
        Number of nearest neighbors to return.
    signame : str, optional
        A human-readable name for the query.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: 'sample', 'index', 'distance'.
    """
    up_genes = list(up_genes)
    down_genes_list = list(down_genes) if down_genes is not None else []

    payload: Dict[str, Any] = {
        "signatures": [{"up_genes": up_genes, "down_genes": down_genes_list}],
        "species": species,
        "k": k,
        "signame": signame,
    }

    data = _post_json("/sigpy/data/knn/signature", payload=payload)

    # Example JSON:
    # {
    #   "distances": [...],
    #   "indexes": [...],
    #   "samples": [...],
    #   "series_count": 130,
    #   "signame": "...",
    #   "species": "human"
    # }
    distances = data.get("distances", [])
    indexes = data.get("indexes", [])
    samples = data.get("samples", [])

    return pd.DataFrame(
        {
            "sample": samples,
            "index": indexes,
            "distance": distances,
        }
    )


def knn_signature_from_profile(
    genes: Iterable[str],
    values: Iterable[float],
    species: str = "human",
    k: int = 10,
    signame: str = "expression-profile similarity search",
) -> pd.DataFrame:
    """
    Wraps: POST /sigpy/data/knn/signature using a full gene expression profile.

    Parameters
    ----------
    genes : iterable of str
        Genes in the signature.
    values : iterable of float
        Expression values aligned with 'genes'.
    species : {"human", "mouse"}, optional
        Species (defaults to "human").
    k : int, optional
        Number of nearest neighbors to return.
    signame : str, optional
        A human-readable name for the query.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: 'sample', 'index', 'distance'.
    """
    genes = list(genes)
    values = list(values)
    if len(genes) != len(values):
        raise ValueError("genes and values must have the same length")

    payload: Dict[str, Any] = {
        "signatures": [{"genes": genes, "values": values}],
        "species": species,
        "k": k,
        "signame": signame,
    }

    data = _post_json("/sigpy/data/knn/signature", payload=payload)

    distances = data.get("distances", [])
    indexes = data.get("indexes", [])
    samples = data.get("samples", [])

    return pd.DataFrame(
        {
            "sample": samples,
            "index": indexes,
            "distance": distances,
        }
    )


# --------- Data query: gene correlation -------------------------------------


def gene_correlation(
    gene: str,
    meta: str,
    species: str = "human",
    k: int = 200,
    samples: Optional[List[str]] = None,
    filter_genes: Optional[List[str]] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Wraps: POST /sigpy/data/correlation

    Parameters
    ----------
    gene : str
        Gene symbol to examine.
    meta : str
        Metadata filter/search term (e.g. "keratinocyte").
    species : {"human", "mouse"}, optional
        Species (defaults to "human").
    k : int, optional
        Number of top correlated genes to return (default: 200).
        Ignored if filter_genes is specified.
    samples : list of str, optional
        Specific GSM sample IDs to analyze. If None, uses all samples
        matching the meta search term.
    filter_genes : list of str, optional
        Specific genes to filter for in results. When provided, requests
        comprehensive correlation data and filters to only these genes.
        Example: ["GENE1", "GENE2"] to find correlation between query gene
        and these specific genes. Ignores k parameter when set.

    Returns
    -------
    dict of str -> pandas.DataFrame
        Dictionary with keys:
        - 'summary'  : gene-level summary (1 row)
        - 'positive' : positively correlated genes (filtered if filter_genes set)
        - 'negative' : negatively correlated genes (filtered if filter_genes set)
        - 'samples'  : samples used in the calculation
        - 'all_correlations' : (only if filter_genes) all correlations combined

    Examples
    --------
    # Get top 200 correlations (default)
    results = gene_correlation("KRT14", "keratinocyte")

    # Get top 50 correlations
    results = gene_correlation("KRT14", "keratinocyte", k=50)

    # Find specific gene correlations (comprehensive)
    results = gene_correlation("KRT14", "keratinocyte", filter_genes=["KRT5", "KRT6A"])

    # Use specific samples
    results = gene_correlation("CD3D", "T cell", samples=["GSM1234567", "GSM1234568"])
    """
    # Build payload
    payload: Dict[str, Any] = {"gene": gene, "meta": meta, "species": species}

    # If filtering for specific genes, request comprehensive data (large k)
    if filter_genes is not None:
        # Request a large number to get comprehensive results
        payload["k"] = 10000
    else:
        payload["k"] = k

    # Add samples if specified
    if samples is not None:
        payload["samples"] = samples

    data = _post_json("/sigpy/data/correlation", payload=payload)

    summary_df = pd.DataFrame(
        [
            {
                "gene": data.get("gene"),
                "mean_log_expression": data.get("mean_log_expression"),
                "searchterm": data.get("searchterm"),
            }
        ]
    )

    pos_df = pd.DataFrame(data.get("positive_correlated_genes", []))
    neg_df = pd.DataFrame(data.get("negative_correlated_genes", []))
    samples_df = pd.DataFrame({"sample": data.get("samples", [])})

    result = {
        "summary": summary_df,
        "positive": pos_df,
        "negative": neg_df,
        "samples": samples_df,
    }

    # If filtering for specific genes, filter the results
    if filter_genes is not None and len(filter_genes) > 0:
        # Combine positive and negative correlations
        all_corr = pd.concat([pos_df, neg_df], ignore_index=True)

        # Filter for specific genes
        if "gene" in all_corr.columns:
            filtered = all_corr[all_corr["gene"].isin(filter_genes)].copy()

            # Sort by absolute correlation value
            if "correlation" in filtered.columns:
                filtered = filtered.sort_values(
                    by="correlation",
                    key=lambda x: abs(x),
                    ascending=False
                )

            # Split back into positive and negative
            if "correlation" in filtered.columns:
                result["positive"] = filtered[filtered["correlation"] >= 0].reset_index(drop=True)
                result["negative"] = filtered[filtered["correlation"] < 0].reset_index(drop=True)
                result["all_correlations"] = filtered.reset_index(drop=True)
            else:
                result["all_correlations"] = filtered.reset_index(drop=True)

    return result


def gene_correlation_pairwise(
    genes: List[str],
    meta: str,
    species: str = "human",
    samples: Optional[List[str]] = None,
    n_samples: Optional[int] = None,
    output_prefix: Optional[str] = None,
    output_dir: Optional[str] = None,
    generate_heatmap: bool = True,
) -> Dict[str, Any]:
    """
    Perform pairwise gene correlation analysis for multiple genes with bootstrap aggregation.

    For a list of genes [A, B, C], computes all pairwise correlations:
    A-B, A-C, B-C. When n_samples > 200, uses bootstrap aggregation by dividing
    samples into batches of 200 and aggregating results using Fisher's z-transformation.
    Returns correlation matrix and optionally generates an HTML heatmap.

    Parameters
    ----------
    genes : list of str
        List of gene symbols (minimum 2 genes required).
    meta : str
        Metadata filter/search term (e.g. "keratinocyte", "T cell").
    species : {"human", "mouse"}, optional
        Species (defaults to "human").
    samples : list of str, optional
        Specific GSM sample IDs to analyze. If None, uses all samples
        matching the meta search term.
    n_samples : int, optional
        Number of samples to use. If None, uses API default (~200 samples).
        If > 200, fetches all available samples, randomly selects n_samples,
        divides into batches of 200, and aggregates correlations across batches
        using Fisher's z-transformation. This overcomes the API's 200-sample
        limitation per query.
    output_prefix : str, optional
        Prefix for output files. If None, uses "gene_correlation_pairwise".
        Files generated: {prefix}_matrix.csv, {prefix}_pairwise.csv, {prefix}_heatmap.html
    output_dir : str, optional
        Directory where output files will be saved. If None, uses current working
        directory (where Python was launched). Can be absolute or relative path.
    generate_heatmap : bool, optional
        Whether to generate HTML heatmap visualization (default: True).

    Returns
    -------
    dict
        Dictionary containing:
        - 'correlation_matrix': DataFrame with pairwise correlations
        - 'pairwise_results': DataFrame with detailed pairwise results
        - 'p_value_matrix': DataFrame with p-values
        - 'heatmap_file': Path to HTML heatmap (if generate_heatmap=True)
        - 'matrix_file': Path to CSV matrix file (if output_prefix given)
        - 'n_batches': Number of bootstrap batches used (if > 1)
        - 'total_samples_used': Total number of unique samples analyzed

    Examples
    --------
    # Analyze with default samples (~200)
    genes = ["KRT14", "KRT5", "KRT6A", "KRT17"]
    results = gene_correlation_pairwise(genes, "keratinocyte")

    # Use bootstrap aggregation with 1600 samples (8 batches of 200)
    results = gene_correlation_pairwise(genes, "keratinocyte", n_samples=1600)
    print(f"Used {results['n_batches']} batches with {results['total_samples_used']} samples")

    # Access results
    print(results['correlation_matrix'])
    print(f"Heatmap saved to: {results['heatmap_file']}")

    Notes
    -----
    - Performs N*(N-1)/2 correlation queries for N genes
    - **Bootstrap Aggregation**: When n_samples > 200, samples are divided into
      non-overlapping batches of 200. Correlations are computed for each batch
      and aggregated using Fisher's z-transformation for robust estimates.
    - P-values are estimated using t-distribution based on effective sample size
    - Heatmap uses hierarchical clustering (average linkage)
    - Significance annotations: * (p<0.05), ** (p<0.01), *** (p<0.001), **** (p<0.0001)
    - Random sampling uses seed=42 for reproducibility
    """
    import numpy as np
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.stats import t as t_dist
    import plotly.graph_objects as go
    import plotly.figure_factory as ff

    # Validate input
    if len(genes) < 2:
        raise ValueError(f"At least 2 genes required, got {len(genes)}")

    # Convert to uppercase and remove duplicates
    genes = sorted(list(set([g.upper() for g in genes])))
    n_genes = len(genes)

    print(f"Performing pairwise correlation for {n_genes} genes: {genes}")
    print(f"Total pairwise comparisons: {n_genes * (n_genes - 1) // 2}")

    # Handle sample selection and batching
    sample_batches = []
    use_bootstrap = False

    if samples is None and n_samples is not None and n_samples > 200:
        # Fetch all available samples for this meta query
        print(f"\nFetching available samples for: {meta}")
        all_samples_df = quicksearch_metadata(meta, species=species)

        if 'sample' in all_samples_df.columns:
            available_samples = all_samples_df['sample'].tolist()
        elif len(all_samples_df.columns) > 0:
            # Take first column if 'sample' not found
            available_samples = all_samples_df.iloc[:, 0].tolist()
        else:
            available_samples = []

        print(f"Found {len(available_samples)} available samples")

        # Randomly sample if we have enough samples
        if len(available_samples) >= n_samples:
            np.random.seed(42)  # Set seed for reproducibility
            selected_samples = list(np.random.choice(available_samples, size=n_samples, replace=False))
            print(f"Randomly sampled {n_samples} samples (seed=42)")
        elif len(available_samples) > 0:
            selected_samples = available_samples
            print(f"Warning: Requested {n_samples} samples but only {len(available_samples)} available")
            print(f"Using all {len(available_samples)} available samples")
        else:
            print(f"Warning: No samples found for meta query '{meta}'")
            print(f"Will use API default sample selection")
            selected_samples = None

        # Divide into batches of 200
        if selected_samples and len(selected_samples) > 200:
            use_bootstrap = True
            batch_size = 200
            n_batches = (len(selected_samples) + batch_size - 1) // batch_size
            print(f"\nBootstrap aggregation: Dividing {len(selected_samples)} samples into {n_batches} batches of ~{batch_size}")

            for i in range(n_batches):
                start_idx = i * batch_size
                end_idx = min(start_idx + batch_size, len(selected_samples))
                batch = selected_samples[start_idx:end_idx]
                sample_batches.append(batch)
                print(f"  Batch {i+1}: {len(batch)} samples")
        elif selected_samples:
            sample_batches.append(selected_samples)
        else:
            sample_batches.append(None)
    elif samples is not None:
        # User provided specific samples
        if len(samples) > 200:
            use_bootstrap = True
            batch_size = 200
            n_batches = (len(samples) + batch_size - 1) // batch_size
            print(f"\nBootstrap aggregation: Dividing {len(samples)} provided samples into {n_batches} batches of ~{batch_size}")

            for i in range(n_batches):
                start_idx = i * batch_size
                end_idx = min(start_idx + batch_size, len(samples))
                batch = samples[start_idx:end_idx]
                sample_batches.append(batch)
                print(f"  Batch {i+1}: {len(batch)} samples")
        else:
            sample_batches.append(samples)
    else:
        # No sample specification, use API default
        sample_batches.append(None)

    # Initialize correlation matrix and p-value matrix
    corr_matrix = pd.DataFrame(np.eye(n_genes), index=genes, columns=genes)
    pval_matrix = pd.DataFrame(np.zeros((n_genes, n_genes)), index=genes, columns=genes)

    # Store detailed pairwise results
    pairwise_data = []

    # Perform pairwise correlations
    for i, gene1 in enumerate(genes):
        for j, gene2 in enumerate(genes):
            if i >= j:  # Skip diagonal and lower triangle (symmetric)
                continue

            print(f"  Computing correlation: {gene1} - {gene2}")

            # Collect correlations from all batches
            batch_correlations = []
            batch_n_samples = []

            for batch_idx, batch_samples in enumerate(sample_batches):
                if use_bootstrap:
                    print(f"    Batch {batch_idx + 1}/{len(sample_batches)}...", end=" ")

                try:
                    # Query correlation using gene1 as reference, filtering for gene2
                    result = gene_correlation(
                        gene=gene1,
                        meta=meta,
                        species=species,
                        samples=batch_samples,
                        filter_genes=[gene2],
                    )

                    # Extract correlation value
                    found_data = False
                    if 'all_correlations' in result and not result['all_correlations'].empty:
                        corr_value = result['all_correlations']['correlation'].iloc[0]
                        n_samp = len(result['samples']) if 'samples' in result else 200

                        batch_correlations.append(corr_value)
                        batch_n_samples.append(n_samp)
                        found_data = True

                        if use_bootstrap:
                            print(f"r={corr_value:.4f} (n={n_samp})")

                    # If forward direction failed, try reverse direction (handles API asymmetry)
                    if not found_data:
                        if use_bootstrap:
                            print("empty, trying reverse...", end=" ")
                        else:
                            print(f"    Forward query empty, trying reverse direction...")

                        result_reverse = gene_correlation(
                            gene=gene2,
                            meta=meta,
                            species=species,
                            samples=batch_samples,
                            filter_genes=[gene1],
                        )

                        if 'all_correlations' in result_reverse and not result_reverse['all_correlations'].empty:
                            corr_value = result_reverse['all_correlations']['correlation'].iloc[0]
                            n_samp = len(result_reverse['samples']) if 'samples' in result_reverse else 200

                            batch_correlations.append(corr_value)
                            batch_n_samples.append(n_samp)
                            found_data = True

                            if use_bootstrap:
                                print(f"r={corr_value:.4f} (n={n_samp}) [reversed]")
                            else:
                                print(f"    ✓ Reverse direction succeeded: r={corr_value:.4f}")
                        else:
                            if use_bootstrap:
                                print("No correlation found")
                            else:
                                print(f"    ✗ Both directions failed - no correlation data available")

                except Exception as e:
                    if use_bootstrap:
                        print(f"Error: {e}")
                    else:
                        print(f"    Error: {e}")

            # Aggregate correlations using Fisher's z-transformation
            if len(batch_correlations) > 0:
                if len(batch_correlations) == 1:
                    # Single batch, use directly
                    corr_value = batch_correlations[0]
                    total_n = batch_n_samples[0]
                else:
                    # Multiple batches, aggregate using Fisher's z-transformation
                    # z = arctanh(r), then average z values weighted by sample size
                    z_values = [np.arctanh(np.clip(r, -0.9999, 0.9999)) for r in batch_correlations]
                    weights = np.array(batch_n_samples)
                    total_n = int(np.sum(weights))

                    # Weighted average of z-scores
                    z_avg = np.average(z_values, weights=weights)

                    # Transform back to correlation
                    corr_value = np.tanh(z_avg)

                    print(f"    Aggregated: r={corr_value:.4f} (from {len(batch_correlations)} batches, total n={total_n})")

                # Store in matrix
                corr_matrix.loc[gene1, gene2] = corr_value
                corr_matrix.loc[gene2, gene1] = corr_value  # Symmetric

                # Estimate p-value from correlation using t-distribution
                # t = r * sqrt(n-2) / sqrt(1-r^2)
                if abs(corr_value) < 0.9999:  # Avoid division by zero
                    t_stat = corr_value * np.sqrt(total_n - 2) / np.sqrt(1 - corr_value**2)
                    p_value = 2 * (1 - t_dist.cdf(abs(t_stat), total_n - 2))
                else:
                    p_value = 0.0

                pval_matrix.loc[gene1, gene2] = p_value
                pval_matrix.loc[gene2, gene1] = p_value

                # Store detailed results
                pairwise_data.append({
                    'gene1': gene1,
                    'gene2': gene2,
                    'correlation': corr_value,
                    'p_value': p_value,
                    'n_samples': total_n,
                    'n_batches': len(batch_correlations) if use_bootstrap else 1,
                })

                if not use_bootstrap:
                    print(f"    Correlation: {corr_value:.4f}, p-value: {p_value:.2e}, n={total_n}")

            else:
                print(f"    Warning: No correlation found for {gene1}-{gene2}")
                corr_matrix.loc[gene1, gene2] = np.nan
                corr_matrix.loc[gene2, gene1] = np.nan

    # Create pairwise results DataFrame
    pairwise_df = pd.DataFrame(pairwise_data)

    # Determine output directory (default to current working directory)
    if output_dir is None:
        output_dir = os.getcwd()

    # Save matrix to CSV if output_prefix provided
    matrix_file = None
    if output_prefix:
        matrix_file = os.path.join(output_dir, f"{output_prefix}_matrix.csv")
        corr_matrix.to_csv(matrix_file)
        print(f"\nCorrelation matrix saved to: {matrix_file}")

        pairwise_file = os.path.join(output_dir, f"{output_prefix}_pairwise.csv")
        pairwise_df.to_csv(pairwise_file, index=False)
        print(f"Pairwise results saved to: {pairwise_file}")

    # Generate heatmap if requested
    heatmap_file = None
    if generate_heatmap:
        heatmap_file = os.path.join(output_dir, f"{output_prefix or 'gene_correlation_pairwise'}_heatmap.html")

        # Perform hierarchical clustering
        # Replace NaN with 0 for clustering
        corr_matrix_filled = corr_matrix.fillna(0)

        # Compute linkage
        linkage_matrix = linkage(corr_matrix_filled.values, method='average')

        # Get dendrogram order
        dend = dendrogram(linkage_matrix, no_plot=True)
        reorder_idx = dend['leaves']
        reordered_genes = [genes[i] for i in reorder_idx]

        # Reorder correlation and p-value matrices
        corr_reordered = corr_matrix.loc[reordered_genes, reordered_genes]
        pval_reordered = pval_matrix.loc[reordered_genes, reordered_genes]

        # Create significance annotations
        annot_text = []
        for i in range(n_genes):
            row_annot = []
            for j in range(n_genes):
                if i == j:
                    row_annot.append('')
                else:
                    p = pval_reordered.iloc[i, j]
                    if p < 0.0001:
                        row_annot.append('****')
                    elif p < 0.001:
                        row_annot.append('***')
                    elif p < 0.01:
                        row_annot.append('**')
                    elif p < 0.05:
                        row_annot.append('*')
                    else:
                        row_annot.append('')
            annot_text.append(row_annot)

        # Create heatmap with plotly
        fig = ff.create_annotated_heatmap(
            z=corr_reordered.values,
            x=reordered_genes,
            y=reordered_genes,
            annotation_text=annot_text,
            colorscale='RdBu_r',  # Blue-White-Red reversed (red for positive)
            zmid=0,  # Center at 0
            showscale=True,
            hovertemplate='%{x} - %{y}<br>Correlation: %{z:.3f}<extra></extra>',
        )

        # Update layout
        fig.update_layout(
            title=f'Gene Correlation Heatmap - {meta}<br>{"" if species == "human" else f"({species})"}',
            xaxis_title='',
            yaxis_title='',
            width=max(600, n_genes * 80),
            height=max(600, n_genes * 80),
            font=dict(size=10),
        )

        # Rotate x-axis labels
        fig.update_xaxes(tickangle=45)

        # Save to HTML
        fig.write_html(heatmap_file)
        print(f"Heatmap saved to: {heatmap_file}")

    # Calculate total samples used
    total_samples = sum([len(batch) for batch in sample_batches if batch is not None])

    return {
        'correlation_matrix': corr_matrix,
        'pairwise_results': pairwise_df,
        'p_value_matrix': pval_matrix,
        'heatmap_file': heatmap_file,
        'matrix_file': matrix_file,
        'n_batches': len(sample_batches),
        'total_samples_used': total_samples if total_samples > 0 else None,
    }


# --------- Data query: differential expression ------------------------------


def diffexp(
    gene: str,
    meta: str,
    species: str = "human",
    fdr_cutoff: float = 0.1,
    return_samples: bool = False,
) -> pd.DataFrame:
    """
    Wraps: POST /sigpy/data/diffexp

    Parameters
    ----------
    gene : str
        Gene to test for differential expression.
    meta : str
        Metadata/filter criteria (e.g. tissue or cell type term).
    species : {"human", "mouse"}, optional
        Species (defaults to "human").
    fdr_cutoff : float, optional
        False discovery rate threshold (default 0.1).
    return_samples : bool, optional
        If True, returns a dictionary with 'genes', 'control_samples', and
        'search_samples'. If False (default), returns only the genes DataFrame.

    Returns
    -------
    pandas.DataFrame or dict
        If return_samples=False (default):
            DataFrame with genes as rows and the following columns:
            - gene: gene symbol
            - fdr: false discovery rate (adjusted p-value)
            - t: t-statistic
            - mean_expression_control: mean log expression in control samples
            - mean_expression_search: mean log expression in search samples
            - log2_fold_change: log2(search/control) fold change
            - abs_log2_fold_change: absolute value of log2 fold change

        If return_samples=True:
            Dictionary with keys:
            - 'genes': DataFrame described above
            - 'control_samples': list of control sample IDs
            - 'search_samples': list of search sample IDs

    Notes
    -----
    The API returns differential expression results comparing samples matching
    the search term against control samples. The genes DataFrame is sorted by
    FDR (most significant first).

    Examples
    --------
    # Get differentially expressed genes
    df = diffexp("KRT14", "keratinocyte", fdr_cutoff=0.05)
    print(df.head())

    # Get results with sample information
    results = diffexp("KRT14", "keratinocyte", return_samples=True)
    print(f"Control samples: {len(results['control_samples'])}")
    print(results['genes'].head())
    """
    payload: Dict[str, Any] = {
        "gene": gene,
        "meta": meta,
        "species": species,
        "fdr_cutoff": fdr_cutoff,
    }
    data = _post_json("/sigpy/data/diffexp", payload=payload)

    # Extract genes list and convert to DataFrame
    if isinstance(data, dict) and "genes" in data:
        genes_list = data.get("genes", [])
        genes_df = pd.DataFrame(genes_list)

        # Calculate log2 fold change if we have the mean expression columns
        if "mean_expression_search" in genes_df.columns and "mean_expression_control" in genes_df.columns:
            # log2FC = log2(search) - log2(control) since values are already in log space
            genes_df["log2_fold_change"] = (
                genes_df["mean_expression_search"] - genes_df["mean_expression_control"]
            )
            genes_df["abs_log2_fold_change"] = genes_df["log2_fold_change"].abs()

        # Sort by FDR (most significant first)
        if "fdr" in genes_df.columns:
            genes_df = genes_df.sort_values("fdr", ascending=True).reset_index(drop=True)

        # Return based on return_samples flag
        if return_samples:
            return {
                "genes": genes_df,
                "control_samples": data.get("control_samples", []),
                "search_samples": data.get("search_samples", []),
            }
        else:
            return genes_df

    # Fallback for unexpected response format
    elif isinstance(data, list):
        genes_df = pd.DataFrame(data)
        if return_samples:
            return {"genes": genes_df, "control_samples": [], "search_samples": []}
        else:
            return genes_df
    else:
        # Last resort normalization
        genes_df = pd.json_normalize(data)
        if return_samples:
            return {"genes": genes_df, "control_samples": [], "search_samples": []}
        else:
            return genes_df


# --------- Data query: tissue expression atlas -----------------------------


def gene_expression(
    gene: str,
    species: str = "human",
) -> pd.DataFrame:
    """
    Get tissue expression atlas for a gene from ARCHS4.

    Wraps: POST /archs4/search/loadExpressionTissue.php

    Parameters
    ----------
    gene : str
        Gene symbol (e.g., "TP53", "BRCA1").
    species : {"human", "mouse"}, optional
        Species (defaults to "human").

    Returns
    -------
    pandas.DataFrame
        DataFrame with tissue expression data, sorted by median expression.
        Columns typically include:
        - tissue: tissue/organ name
        - mean: mean expression across samples
        - median: median expression across samples
        - min: minimum expression
        - max: maximum expression
        - std: standard deviation

    Examples
    --------
    # Get tissue expression for TP53
    df = gene_expression("TP53", species="human")
    print(df.head())

    # Get tissue expression for mouse gene
    df = gene_expression("Tp53", species="mouse")

    Notes
    -----
    - Gene names are case-insensitive (converted to uppercase internally)
    - Returns expression data across different tissues/organs
    - Data is sorted by median expression (highest first)
    """
    # Validate species
    if species not in ["human", "mouse"]:
        raise ValueError(f"species must be 'human' or 'mouse', got '{species}'")

    # Make gene uppercase
    gene = gene.upper()

    # Build query URL
    query = f"search={gene}&species={species}&type=tissue"
    url = f"{BASE_URL}/archs4/search/loadExpressionTissue.php?{query}"

    # Submit request
    response = requests.post(url, headers={"Content-Type": "application/json"}, timeout=60)
    response.raise_for_status()

    # Parse response as CSV
    df = pd.read_csv(io.StringIO(response.content.decode("utf-8")))

    # Check if any results were returned
    if len(df) < 2:
        raise ValueError(
            f"Gene '{gene}' did not return any tissue expression results. "
            "Please check the gene symbol and species."
        )

    # Clean up dataframe
    # Drop NaN rows
    df = df.dropna()

    # Drop color column if present
    if "color" in df.columns:
        df = df.drop(["color"], axis=1)

    # Rename 'id' column to 'tissue' if present
    if "id" in df.columns:
        df = df.rename(columns={"id": "tissue"})

    # Add 'mean' column if not present (calculate from q1, median, q3)
    if "mean" not in df.columns and "q1" in df.columns and "median" in df.columns and "q3" in df.columns:
        # Approximate mean as average of q1, median, q3
        df["mean"] = (df["q1"] + df["median"] + df["q3"]) / 3

    # Sort by median expression (highest first)
    if "median" in df.columns:
        df = df.sort_values("median", ascending=False)

    df = df.reset_index(drop=True)

    return df


def gene_expression_analysis(
    genes: List[str],
    species: str = "human",
    tissue_filter: Optional[str] = None,
    n_tissues: int = 10,
    output_prefix: Optional[str] = None,
    output_dir: Optional[str] = None,
    generate_plot: bool = True,
) -> Dict[str, Any]:
    """
    Analyze tissue expression for multiple genes with integrated visualization.

    Queries tissue expression for one or more genes, combines results into
    an integrated table, and generates an interactive boxplot showing expression
    patterns across tissues/cell types.

    Parameters
    ----------
    genes : list of str or str
        Gene symbol(s) to analyze. Can be a single gene or list of genes.
    species : {"human", "mouse"}, optional
        Species (defaults to "human").
    tissue_filter : str or list of str, optional
        Filter tissues by name (substring/fuzzy match). Can be single string
        or list of strings. If None, shows top N tissues by expression.
        Examples: "liver", ["liver", "hepatocyte"], ["T cell", "immune"].
    n_tissues : int, optional
        Number of top tissues to show (default: 10). For multiple genes,
        selects top N by mean expression across genes.
    output_prefix : str, optional
        Prefix for output files. If provided, saves:
        - {prefix}_expression_table.csv: Integrated expression data
        - {prefix}_boxplot.html: Interactive visualization
    output_dir : str, optional
        Directory where output files will be saved. If None, uses current working
        directory (where Python was launched). Can be absolute or relative path.
    generate_plot : bool, optional
        Whether to generate interactive boxplot (default: True).

    Returns
    -------
    dict
        Dictionary containing:
        - 'expression_data': DataFrame with expression for all genes/tissues
        - 'top_tissues': List of selected tissue names
        - 'plot_file': Path to HTML plot (if generate_plot=True)
        - 'table_file': Path to CSV table (if output_prefix given)

    Examples
    --------
    # Single gene analysis
    results = gene_expression_analysis("TP53")

    # Multiple genes with default top 10 tissues
    results = gene_expression_analysis(["KRT14", "KRT5", "KRT6A"])

    # Filter for specific tissue type
    results = gene_expression_analysis(
        ["ALB", "AFP", "TTR"],
        tissue_filter="liver",
        output_prefix="liver_markers"
    )

    # Custom number of tissues
    results = gene_expression_analysis(
        ["CD3D", "CD3E", "CD4", "CD8A"],
        n_tissues=15,
        tissue_filter="immune"
    )

    Notes
    -----
    - For multiple genes, top N tissues are selected by mean expression
    - Tissue filtering uses case-insensitive substring matching
    - Boxplot shows expression distribution with different color per gene
    - Interactive HTML plot allows zooming, panning, and data inspection
    """
    from difflib import get_close_matches
    import plotly.graph_objects as go

    # Handle single gene input
    if isinstance(genes, str):
        genes = [genes]

    # Validate input
    if len(genes) == 0:
        raise ValueError("At least one gene must be provided")

    genes = [g.upper() for g in genes]
    n_genes = len(genes)

    print(f"Analyzing expression for {n_genes} gene(s): {genes}")
    print(f"Species: {species}")

    # Query expression for each gene
    all_data = []
    for gene in genes:
        print(f"  Querying: {gene}")
        try:
            df = gene_expression(gene, species=species)
            df['gene'] = gene
            all_data.append(df)
        except Exception as e:
            print(f"  Warning: Failed to get expression for {gene}: {e}")

    if len(all_data) == 0:
        raise ValueError("No expression data retrieved for any gene")

    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)
    print(f"\nRetrieved expression data: {len(combined_df)} rows")

    # Filter tissues if requested
    if tissue_filter is not None:
        # Handle single string or list of strings
        if isinstance(tissue_filter, str):
            tissue_filters = [tissue_filter]
        else:
            tissue_filters = list(tissue_filter)

        print(f"\nFiltering tissues with: {tissue_filters}")

        # Build vocabulary from available tissues for suggestions
        unique_tissues = combined_df['tissue'].unique()
        # Extract meaningful terms from tissue names (system, organ, cell type)
        vocabulary = set()
        for tissue in unique_tissues:
            parts = tissue.lower().split('.')
            vocabulary.update(parts)
            # Add cell type as whole
            if len(parts) >= 3:
                vocabulary.add(parts[-1])

        # Add common filter terms
        common_terms = [
            'nervous', 'immune', 'digestive', 'respiratory', 'cardiovascular',
            'urogenital', 'integumentary', 'muscular', 'connective',
            'brain', 'liver', 'kidney', 'heart', 'lung', 'skin', 'blood',
            'neuron', 'hepatocyte', 'lymphocyte', 'epithelial', 'fibroblast'
        ]
        vocabulary.update(common_terms)

        # Try substring matching for each filter
        all_matches = set()
        for filt in tissue_filters:
            filt_lower = filt.lower()
            # Substring match (case-insensitive)
            mask = combined_df['tissue'].str.lower().str.contains(filt_lower, na=False)
            matches = combined_df[mask]['tissue'].unique()

            if len(matches) > 0:
                print(f"  '{filt}': {len(matches)} tissues matched")
                all_matches.update(matches)
            else:
                # Try fuzzy matching on full tissue names
                fuzzy_matches = get_close_matches(filt, unique_tissues, n=10, cutoff=0.4)
                if fuzzy_matches:
                    print(f"  '{filt}': {len(fuzzy_matches)} fuzzy matches")
                    all_matches.update(fuzzy_matches)
                else:
                    # No matches found - suggest alternatives
                    print(f"  '{filt}': No matches found")

                    # Find closest terms in vocabulary
                    suggestions = get_close_matches(filt_lower, vocabulary, n=5, cutoff=0.3)
                    if suggestions:
                        print(f"    💡 Did you mean: {', '.join(suggestions[:3])}?")
                        print(f"    💡 Try one of these filters: {', '.join(suggestions)}")
                    else:
                        # Suggest by category
                        print(f"    💡 Common filters by category:")
                        print(f"       Organs: liver, brain, kidney, heart, lung, skin")
                        print(f"       Cells: hepatocyte, neuron, lymphocyte, keratinocyte")
                        print(f"       Systems: immune, nervous, digestive, cardiovascular")

        if len(all_matches) > 0:
            filtered_df = combined_df[combined_df['tissue'].isin(all_matches)].copy()
            print(f"  Total: {len(all_matches)} unique tissues matched")
        else:
            print(f"  Warning: No tissues matching any filter. Using all tissues.")
            filtered_df = combined_df.copy()
    else:
        filtered_df = combined_df.copy()

    # Select top N tissues
    # Calculate mean expression per tissue (across all genes)
    tissue_mean = filtered_df.groupby('tissue')['median'].mean().sort_values(ascending=False)

    n_available = len(tissue_mean)
    n_select = min(n_tissues, n_available)

    if tissue_filter is not None:
        print(f"  Using {n_select} tissue(s) from filtered set")
    else:
        print(f"\nSelecting top {n_select} tissues by mean expression")

    top_tissues = tissue_mean.head(n_select).index.tolist()

    # Filter data to selected tissues
    plot_data = filtered_df[filtered_df['tissue'].isin(top_tissues)].copy()

    # Extract short tissue names (last component)
    plot_data['tissue_short'] = plot_data['tissue'].apply(
        lambda x: x.split('.')[-1] if '.' in x else x
    )

    # Sort tissues by mean expression
    tissue_order = plot_data.groupby('tissue_short')['median'].mean().sort_values(ascending=False).index.tolist()

    print(f"\nSelected tissues:")
    for i, tissue in enumerate(top_tissues, 1):
        short_name = tissue.split('.')[-1] if '.' in tissue else tissue
        mean_expr = tissue_mean[tissue]
        print(f"  {i:2d}. {short_name:30s} (mean={mean_expr:.2f})")

    # Determine output directory (default to current working directory)
    if output_dir is None:
        output_dir = os.getcwd()

    # Save integrated table
    table_file = None
    if output_prefix:
        table_file = os.path.join(output_dir, f"{output_prefix}_expression_table.csv")
        plot_data.to_csv(table_file, index=False)
        print(f"\n✓ Expression table saved: {table_file}")

    # Generate boxplot
    plot_file = None
    if generate_plot:
        plot_file = os.path.join(output_dir, f"{output_prefix or 'gene_expression'}_boxplot.html")

        fig = go.Figure()

        # Add trace for each gene
        colors = [
            '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
            '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
        ]

        for i, gene in enumerate(genes):
            gene_data = plot_data[plot_data['gene'] == gene]

            # Prepare data for boxplot (need values for each tissue)
            y_values = []
            x_labels = []
            for tissue in tissue_order:
                tissue_data = gene_data[gene_data['tissue_short'] == tissue]
                if not tissue_data.empty:
                    # Use median as representative value (boxplot needs list)
                    # For true distribution, we'd need raw samples, but we only have summary stats
                    # Create synthetic distribution from quartiles
                    row = tissue_data.iloc[0]
                    q1, median, q3 = row['q1'], row['median'], row['q3']
                    # Approximate distribution
                    values = [q1, q1, median, median, median, q3, q3]
                    y_values.extend(values)
                    x_labels.extend([tissue] * len(values))

            fig.add_trace(go.Box(
                y=y_values,
                x=x_labels,
                name=gene,
                marker_color=colors[i % len(colors)],
                boxmean='sd'
            ))

        # Update layout
        title_text = f"Gene Expression Across Tissues/Cell Types"
        if len(genes) == 1:
            title_text = f"{genes[0]} Expression Across Tissues"
        elif tissue_filter:
            # Handle both string and list cases
            if isinstance(tissue_filter, str):
                title_text = f"Gene Expression in {tissue_filter.title()} Tissues"
            else:
                # For list, use first term or generic label
                filter_label = tissue_filter[0].title() if len(tissue_filter) == 1 else "Filtered"
                title_text = f"Gene Expression in {filter_label} Tissues"

        fig.update_layout(
            title=title_text,
            xaxis_title="Tissue / Cell Type",
            yaxis_title="Expression Level",
            boxmode='group',  # Group boxes for each tissue
            xaxis={'categoryorder': 'array', 'categoryarray': tissue_order},
            hovermode='closest',
            width=max(800, n_select * 80),
            height=600,
            legend=dict(
                title="Gene",
                orientation="v",
                yanchor="top",
                y=1,
                xanchor="left",
                x=1.02
            )
        )

        # Rotate x-axis labels if many tissues
        if n_select > 5:
            fig.update_xaxes(tickangle=-45)

        fig.write_html(plot_file)
        print(f"✓ Interactive boxplot saved: {plot_file}")

    return {
        'expression_data': plot_data,
        'top_tissues': top_tissues,
        'plot_file': plot_file,
        'table_file': table_file,
        'n_genes': n_genes,
        'n_tissues': n_select,
    }
