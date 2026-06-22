#!/usr/bin/env python3
"""
SOMA Loader Module for Omicsoft DEG Analysis

Provides functions for loading TileDB-SOMA experiments from S3 URLs,
auto-detecting structure type (DEG vs Census), and translating filters.

Usage:
    from soma_loader import is_soma_uri, load_soma_deg, detect_soma_structure

Dependencies:
    - tiledbsoma (optional, for SOMA support)
    - anndata
    - pandas
"""

import os
from enum import Enum
from typing import Dict, List, Optional, Tuple, Any, Callable
import warnings

# Check if tiledbsoma is available
try:
    import tiledbsoma
    SOMA_AVAILABLE = True
except ImportError:
    SOMA_AVAILABLE = False

import pandas as pd
import numpy as np


class SOMAStructureType(Enum):
    """SOMA experiment structure types."""
    DEG = "deg"           # Has log2fc, padj layers
    EXPRESSION = "expression"  # Has tpm/fpkm/raw_counts layers (Omicsoft expression)
    CENSUS = "census"     # Has raw/normalized expression (CellxGene Census)
    UNKNOWN = "unknown"   # Unknown structure, requires user clarification


# Type alias for context
ContextType = Any


def _create_context(config: Dict[str, str]) -> ContextType:
    """Create context using SOMAContext (preferred) or SOMATileDBContext (fallback)."""
    if not SOMA_AVAILABLE:
        raise ImportError("tiledbsoma is not installed. Install with: pip install tiledbsoma")

    if hasattr(tiledbsoma, 'SOMAContext'):
        return tiledbsoma.SOMAContext(tiledb_config=config)
    else:
        # Fallback for older versions
        context = tiledbsoma.SOMATileDBContext()
        return context.replace(tiledb_config=config)


def get_soma_context(
    region: str = "us-west-2",
    access_key: Optional[str] = None,
    secret_key: Optional[str] = None,
    session_token: Optional[str] = None
) -> ContextType:
    """
    Build SOMA context for S3 access.

    Uses environment credentials if explicit credentials not provided.

    Args:
        region: AWS region (default: us-west-2)
        access_key: AWS access key ID (optional, uses env if not provided)
        secret_key: AWS secret access key (optional)
        session_token: AWS session token for STS (optional)

    Returns:
        Configured SOMA context
    """
    config = {"vfs.s3.region": region}

    if access_key:
        config["vfs.s3.aws_access_key_id"] = access_key
    if secret_key:
        config["vfs.s3.aws_secret_access_key"] = secret_key
    if session_token:
        config["vfs.s3.aws_session_token"] = session_token

    return _create_context(config)


def is_soma_uri(uri: str) -> bool:
    """
    Check if a URI is likely a SOMA experiment.

    SOMA URIs are S3 paths that don't end in .h5ad.

    Args:
        uri: URI to check

    Returns:
        True if likely a SOMA URI, False otherwise

    Examples:
        >>> is_soma_uri("s3://bucket/experiment_soma")
        True
        >>> is_soma_uri("s3://bucket/data.h5ad")
        False
        >>> is_soma_uri("/local/path/data.h5ad")
        False
    """
    if not isinstance(uri, str):
        return False

    # Must be an S3 URI
    if not uri.startswith("s3://"):
        return False

    # Should not end with .h5ad
    if uri.lower().endswith(".h5ad"):
        return False

    return True


def try_soma_open(uri: str, context: Optional[ContextType] = None) -> Tuple[bool, Optional[Any]]:
    """
    Try to open a URI as a SOMA experiment.

    Args:
        uri: URI to open
        context: Optional SOMA context for S3 access

    Returns:
        Tuple of (success: bool, experiment or None)
        Note: Caller is responsible for closing the experiment

    Example:
        success, exp = try_soma_open("s3://bucket/experiment")
        if success:
            print(f"Opened SOMA with {exp.obs.count} observations")
            exp.close()
    """
    if not SOMA_AVAILABLE:
        return False, None

    try:
        exp = tiledbsoma.open(uri, mode="r", context=context)
        # Verify it's a valid experiment by checking for obs
        _ = exp.obs.count
        return True, exp
    except Exception as e:
        return False, None


def detect_soma_structure(uri: str, context: Optional[ContextType] = None) -> SOMAStructureType:
    """
    Detect the structure type of a SOMA experiment.

    Args:
        uri: SOMA experiment URI
        context: Optional SOMA context for S3 access

    Returns:
        SOMAStructureType indicating DEG, CENSUS, or UNKNOWN

    Structure detection logic:
        - DEG: Has 'log2fc' and 'padj' in X layers
        - Census: Has 'raw' layer with Census-style obs columns
        - Unknown: Neither of the above
    """
    if not SOMA_AVAILABLE:
        raise ImportError("tiledbsoma is not installed. Install with: pip install tiledbsoma")

    with tiledbsoma.open(uri, mode="r", context=context) as exp:
        # Get X layer names
        x_layers = list(exp.ms["RNA"].X.keys())

        # DEG: has log2fc and padj
        if 'log2fc' in x_layers and 'padj' in x_layers:
            return SOMAStructureType.DEG

        # Omicsoft Expression: has tpm/fpkm/raw_counts layers
        expr_indicators = ['tpm', 'fpkm', 'raw_counts']
        if any(layer in x_layers for layer in expr_indicators):
            return SOMAStructureType.EXPRESSION

        # URI-based fallback for expression
        if '_expr' in uri.lower():
            return SOMAStructureType.EXPRESSION

        # Census: has raw + Census-style obs columns
        if 'raw' in x_layers:
            obs_cols = list(exp.obs.keys())
            census_indicators = ['cell_type', 'tissue_general', 'assay', 'cell_type_ontology_term_id']
            if any(col in obs_cols for col in census_indicators):
                return SOMAStructureType.CENSUS

        return SOMAStructureType.UNKNOWN


def prefetch_obs_values(exp, columns: List[str]) -> Dict[str, List[str]]:
    """
    Pre-fetch unique values from obs columns.

    Args:
        exp: Opened SOMA experiment
        columns: List of column names to fetch unique values for

    Returns:
        Dict mapping column name to list of unique values

    Example:
        unique_vals = prefetch_obs_values(exp, ['disease', 'tissue'])
        print(f"Diseases: {unique_vals['disease'][:10]}")
    """
    result = {}

    # Get available columns
    available_cols = list(exp.obs.keys())

    for col in columns:
        if col in available_cols:
            try:
                # Read only this column
                df = exp.obs.read(column_names=[col]).concat().to_pandas()
                # Get unique values, dropping NaN
                unique_values = df[col].dropna().unique().tolist()
                result[col] = unique_values
            except Exception as e:
                warnings.warn(f"Failed to fetch unique values for column '{col}': {e}")
                result[col] = []
        else:
            result[col] = []

    return result


def translate_substring_filter(
    available_values: List[str],
    search_patterns: List[str],
    case_insensitive: bool = True
) -> List[str]:
    """
    Translate substring search patterns to exact matching values.

    Args:
        available_values: List of available values in the column
        search_patterns: List of substring patterns to search for
        case_insensitive: Whether to match case-insensitively (default: True)

    Returns:
        List of exact values that match any of the patterns

    Example:
        >>> available = ["Crohn's disease (CD)", "Ulcerative colitis (UC)", "Normal"]
        >>> translate_substring_filter(available, ["crohn", "colitis"])
        ["Crohn's disease (CD)", "Ulcerative colitis (UC)"]
    """
    matches = []

    for value in available_values:
        value_str = str(value)
        value_check = value_str.lower() if case_insensitive else value_str

        for pattern in search_patterns:
            pattern_check = pattern.lower() if case_insensitive else pattern
            if pattern_check in value_check:
                matches.append(value_str)
                break  # No need to check more patterns for this value

    return matches


def _escape_soma_value(value: str) -> str:
    """Escape a value for use in SOMA filter strings."""
    # Escape single quotes by doubling them
    return value.replace("'", "\\'")


def build_soma_filter(filter_config: Dict[str, Dict[str, Any]]) -> Optional[str]:
    """
    Build a SOMA filter string from filter configuration.

    Args:
        filter_config: Dict with structure:
            {
                'column_name': {
                    'values': List[str],  # Exact values to match
                    'negate': bool,       # If True, exclude these values
                    'operator': str       # 'in' or 'not in'
                }
            }

    Returns:
        SOMA filter string or None if no filters

    Example:
        config = {
            'disease': {'values': ["Crohn's disease", "UC"], 'negate': False},
            'tissue': {'values': ["colon", "ileum"], 'negate': False}
        }
        filter_str = build_soma_filter(config)
        # Returns: "disease in ['Crohn\\'s disease', 'UC'] and tissue in ['colon', 'ileum']"
    """
    clauses = []

    for col, spec in filter_config.items():
        values = spec.get('values', [])
        negate = spec.get('negate', False)

        if not values:
            continue

        # Escape values for SOMA filter syntax
        escaped_values = [f"'{_escape_soma_value(v)}'" for v in values]
        values_str = f"[{', '.join(escaped_values)}]"

        if negate:
            clause = f"{col} not in {values_str}"
        else:
            clause = f"{col} in {values_str}"

        clauses.append(clause)

    if not clauses:
        return None

    return " and ".join(clauses)


def load_soma_deg(
    uri: str,
    filters: Optional[Dict] = None,
    context: Optional[ContextType] = None,
    ask_user_callback: Optional[Callable] = None,
    x_layer: str = "log2fc",
    include_layers: Optional[List[str]] = None
):
    """
    Load DEG data from a SOMA experiment as AnnData.

    Args:
        uri: SOMA experiment URI (local or S3)
        filters: Optional filter configuration dict with keys matching CLI args:
            - diseases: str (comma-separated patterns)
            - exclude_diseases: str (comma-separated patterns)
            - tissues: str (comma-separated patterns)
            - studies: str (comma-separated exact values)
            - comparison_category: str (comma-separated exact values)
            - case_treatment: str (comma-separated exact values)
            - comparison: str (comma-separated patterns)
        context: Optional SOMA context for S3 access
        ask_user_callback: Optional callback for user clarification (for UNKNOWN structure)
        x_layer: Primary X layer name (default: "log2fc")
        include_layers: List of additional layers to include (default: ["padj", "data", "sig_score"])

    Returns:
        AnnData object with DEG data

    Example:
        adata = load_soma_deg(
            "s3://bucket/deg_soma",
            filters={"diseases": "crohn,colitis", "tissues": "colon"}
        )
    """
    if not SOMA_AVAILABLE:
        raise ImportError("tiledbsoma is not installed. Install with: pip install tiledbsoma")

    import anndata

    # Default layers to include
    if include_layers is None:
        include_layers = ["padj", "data", "sig_score"]

    print(f"Opening SOMA experiment: {uri}")

    with tiledbsoma.open(uri, mode="r", context=context) as exp:
        # Detect structure type
        structure_type = detect_soma_structure(uri, context)
        print(f"Detected SOMA structure: {structure_type.value}")

        if structure_type == SOMAStructureType.UNKNOWN:
            if ask_user_callback:
                # Allow user to specify how to handle
                user_choice = ask_user_callback(
                    "SOMA experiment has unknown structure. How should we proceed?",
                    ["Treat as DEG data (log2fc/padj)", "Treat as expression data", "Cancel"]
                )
                if user_choice == 2:  # Cancel
                    raise ValueError("User cancelled loading of unknown SOMA structure")
                elif user_choice == 1:  # Expression
                    x_layer = "data" if "data" in exp.ms["RNA"].X.keys() else list(exp.ms["RNA"].X.keys())[0]
            else:
                warnings.warn(
                    f"SOMA experiment has unknown structure. Attempting to load as DEG data. "
                    f"Available X layers: {list(exp.ms['RNA'].X.keys())}"
                )

        # Build SOMA filter from user filters
        obs_filter = None
        if filters:
            filter_config = _translate_filters_to_soma(exp, filters)
            obs_filter = build_soma_filter(filter_config)

            if obs_filter:
                print(f"Applying SOMA filter: {obs_filter[:200]}..." if len(obs_filter) > 200 else f"Applying SOMA filter: {obs_filter}")

        # Build query
        obs_query = None
        if obs_filter:
            obs_query = tiledbsoma.AxisQuery(value_filter=obs_filter)

        query = exp.axis_query(
            measurement_name="RNA",
            obs_query=obs_query
        )

        # Convert to AnnData
        print(f"Converting to AnnData with X_name='{x_layer}'...")
        adata = query.to_anndata(X_name=x_layer)

        # Add additional layers
        x_layers = list(exp.ms["RNA"].X.keys())
        for layer_name in include_layers:
            if layer_name != x_layer and layer_name in x_layers:
                try:
                    print(f"  Adding layer: {layer_name}")
                    layer_adata = query.to_anndata(X_name=layer_name)
                    adata.layers[layer_name] = layer_adata.X
                except Exception as e:
                    warnings.warn(f"Failed to load layer '{layer_name}': {e}")

        query.close()

        # Ensure the primary layer is also named in layers dict
        if x_layer not in adata.layers:
            adata.layers[x_layer] = adata.X.copy()

        print(f"Loaded AnnData: {adata.shape[0]} observations x {adata.shape[1]} variables")
        print(f"  Layers: {list(adata.layers.keys())}")

        return adata


def _translate_filters_to_soma(exp, filters: Dict) -> Dict[str, Dict[str, Any]]:
    """
    Translate CLI-style filters to SOMA filter configuration.

    Performs substring matching for fuzzy filters and exact matching for
    exact filters, then builds a SOMA-compatible filter config.

    Args:
        exp: Opened SOMA experiment
        filters: Dict with filter names as keys and filter values as strings

    Returns:
        Filter configuration dict suitable for build_soma_filter()
    """
    # Column mapping: CLI arg -> SOMA column
    column_mapping = {
        'diseases': 'disease',
        'exclude_diseases': 'disease',
        'tissues': 'tissue',
        'studies': 'study',
        'comparison_category': 'comparison_category',
        'case_treatment': 'case_treatment',
        'comparison': 'comparison'
    }

    # Fuzzy (substring) vs exact matching
    fuzzy_filters = {'diseases', 'exclude_diseases', 'tissues', 'comparison'}

    # Pre-fetch unique values for columns we need
    needed_columns = set()
    for filter_key in filters.keys():
        if filter_key in column_mapping:
            needed_columns.add(column_mapping[filter_key])

    print(f"Pre-fetching unique values for columns: {list(needed_columns)}")
    unique_values = prefetch_obs_values(exp, list(needed_columns))

    # Build filter config
    filter_config = {}

    for filter_key, filter_value in filters.items():
        if filter_key not in column_mapping or not filter_value:
            continue

        column = column_mapping[filter_key]
        is_fuzzy = filter_key in fuzzy_filters
        is_negated = filter_key == 'exclude_diseases'

        # Parse comma-separated values
        search_terms = [v.strip() for v in filter_value.split(',') if v.strip()]

        if is_fuzzy:
            # Substring matching
            matched_values = translate_substring_filter(
                unique_values.get(column, []),
                search_terms,
                case_insensitive=True
            )

            if matched_values:
                print(f"  {filter_key}: {len(search_terms)} patterns -> {len(matched_values)} matches")

                # Handle the special case of both include and exclude diseases
                if filter_key == 'diseases':
                    if column not in filter_config:
                        filter_config[column] = {'values': matched_values, 'negate': False}
                    else:
                        # Append to existing
                        filter_config[column]['values'].extend(matched_values)
                        filter_config[column]['values'] = list(set(filter_config[column]['values']))
                elif filter_key == 'exclude_diseases':
                    # For exclude, we need to handle this differently
                    # We'll exclude these from the include list if it exists
                    if column in filter_config:
                        existing = set(filter_config[column]['values'])
                        to_exclude = set(matched_values)
                        filter_config[column]['values'] = list(existing - to_exclude)
                    else:
                        # No include filter, create a negate filter
                        filter_config[f"{column}_exclude"] = {'values': matched_values, 'negate': True}
                else:
                    filter_config[column] = {'values': matched_values, 'negate': is_negated}
            else:
                warnings.warn(f"No matches found for {filter_key}: {search_terms}")
        else:
            # Exact matching - verify values exist
            available = set(unique_values.get(column, []))
            valid_values = [v for v in search_terms if v in available]

            if valid_values:
                print(f"  {filter_key}: {len(valid_values)} exact matches")
                filter_config[column] = {'values': valid_values, 'negate': False}
            else:
                warnings.warn(f"No exact matches found for {filter_key}: {search_terms}")

    return filter_config


def validate_soma_access(uri: str, context: Optional[ContextType] = None) -> bool:
    """
    Validate access to SOMA experiment.

    Args:
        uri: SOMA experiment URI
        context: Optional SOMA context for S3 access

    Returns:
        True if accessible, False otherwise
    """
    if not SOMA_AVAILABLE:
        print("tiledbsoma is not installed")
        return False

    try:
        success, exp = try_soma_open(uri, context)
        if success:
            exp.close()
            print(f"SOMA access validated: {uri}")
            return True
        return False
    except Exception as e:
        print(f"SOMA access validation failed: {e}")
        return False


def get_soma_schema_info(uri: str, context: Optional[ContextType] = None) -> Dict[str, Any]:
    """
    Get schema information from a SOMA experiment.

    Args:
        uri: SOMA experiment URI
        context: Optional SOMA context for S3 access

    Returns:
        Dict with schema information including obs columns, var info, and X layers
    """
    if not SOMA_AVAILABLE:
        raise ImportError("tiledbsoma is not installed. Install with: pip install tiledbsoma")

    with tiledbsoma.open(uri, mode="r", context=context) as exp:
        structure_type = detect_soma_structure(uri, context)

        info = {
            'uri': uri,
            'structure_type': structure_type.value,
            'n_obs': exp.obs.count,
            'obs_columns': list(exp.obs.keys()),
            'measurements': {}
        }

        for ms_name in exp.ms.keys():
            ms = exp.ms[ms_name]
            ms_info = {
                'n_vars': ms.var.count,
                'var_columns': list(ms.var.keys()),
                'x_layers': list(ms.X.keys()),
            }
            info['measurements'][ms_name] = ms_info

        return info


if __name__ == "__main__":
    import sys

    print("SOMA Loader for Omicsoft DEG Analysis")
    print(f"tiledbsoma available: {SOMA_AVAILABLE}")

    if len(sys.argv) > 1:
        uri = sys.argv[1]

        print(f"\nChecking URI: {uri}")
        print(f"Is SOMA URI: {is_soma_uri(uri)}")

        if SOMA_AVAILABLE and is_soma_uri(uri):
            try:
                info = get_soma_schema_info(uri)
                print(f"\nSOMA Structure Type: {info['structure_type']}")
                print(f"Observations: {info['n_obs']}")
                print(f"Obs columns: {info['obs_columns'][:10]}...")

                for ms_name, ms_info in info['measurements'].items():
                    print(f"\nMeasurement: {ms_name}")
                    print(f"  Variables: {ms_info['n_vars']}")
                    print(f"  X layers: {ms_info['x_layers']}")
            except Exception as e:
                print(f"Error: {e}")
