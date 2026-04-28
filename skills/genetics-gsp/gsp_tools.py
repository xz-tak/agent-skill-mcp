"""
GSP Tools - Genetics Support Profiler Analysis Tools

Standalone tool functions for genetic evidence assessment using GSP data.
Supports two primary workflows:
1. Safety Risk Assessment - Identify potential safety concerns for a gene
2. Biological Risk Assessment - Evaluate target-disease associations

Data Access:
    Primary: REST API (configured in config.json)
    Testing: Excel files (for validation only, will be removed)

Usage:
    from gsp_tools import phenotypes_for_gene
    result = phenotypes_for_gene("NLRP3", "NLRP3")

Configuration:
    Edit config.json in this directory to set the API URL for your environment.
    Environment variable GSP_API_URL overrides config.json if set.
"""

from __future__ import annotations

import glob
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd


# ============================================================================
# CONFIGURATION
# ============================================================================

SCRIPT_DIR = Path(__file__).parent
CONFIG_PATH = SCRIPT_DIR / "config.json"


def _load_config() -> dict:
    """Load configuration from config.json."""
    if CONFIG_PATH.exists():
        try:
            with open(CONFIG_PATH) as f:
                return json.load(f)
        except Exception:
            pass
    return {}


_config = _load_config()

# API URL: environment variable overrides config.json
GSP_API_URL = os.environ.get("GSP_API_URL") or _config.get("api_url", "")
GSP_API_TIMEOUT = int(os.environ.get("GSP_API_TIMEOUT") or _config.get("api_timeout", 30))

# Try to import requests for API calls
try:
    import requests
    _HAS_REQUESTS = True
except ImportError:
    _HAS_REQUESTS = False


def _fetch_sheet_via_api(target: str, sheet_name: str) -> Optional[List[Dict[str, Any]]]:
    """Fetch raw sheet data via REST API.

    Args:
        target: Target name (e.g., "NLRP3")
        sheet_name: Sheet name (e.g., "GWAS Summary", "OMIM")

    Returns:
        List of row dictionaries if successful, None otherwise.
    """
    if not GSP_API_URL or not _HAS_REQUESTS:
        return None

    try:
        response = requests.get(
            f"{GSP_API_URL}/sheets/{target}/{sheet_name}",
            timeout=GSP_API_TIMEOUT
        )
        if response.ok:
            return response.json().get('data', [])
    except Exception:
        pass

    return None


def _fetch_all_sheets_via_api(sheet_name: str) -> Optional[List[Dict[str, Any]]]:
    """Fetch sheet data from ALL datasets via REST API.

    Used for cross-target queries like Study Info.

    Args:
        sheet_name: Sheet name (e.g., "Study Info")

    Returns:
        List of row dictionaries from all datasets if successful, None otherwise.
    """
    if not GSP_API_URL or not _HAS_REQUESTS:
        return None

    try:
        response = requests.get(
            f"{GSP_API_URL}/sheets-all/{sheet_name}",
            timeout=GSP_API_TIMEOUT
        )
        if response.ok:
            return response.json().get('data', [])
    except Exception:
        pass

    return None


# ============================================================================
# CONSTANTS
# ============================================================================

# Evidence strength scoring
EVIDENCE_SCORE_MAP = {
    'Very strong': 5.0,
    'Strong': 4.0,
    'Medium': 3.0,
    'Low': 2.0,
    'Very low': 1.0,
    'Unclear': 0.5,
    'None': 0.0
}

# Default data directory: ~/tmp (user home)
# Avoids polluting the skill installation directory with large data files
DATA_DIR = Path.home() / "tmp"
SCHEMA_PATH = SCRIPT_DIR / "gsp_schema.json"

REFERENCES_DIR = SCRIPT_DIR / "references"

# Module-level cache: data_dir_str -> Dict[sheet_name, DataFrame]
_GSP_REPORT_CACHE: Dict[str, Dict[str, pd.DataFrame]] = {}


# ============================================================================
# GSP REPORT LOADING (PICKLE / XLSX)
# ============================================================================

def _load_gsp_report(data_dir: Path) -> Optional[Dict[str, pd.DataFrame]]:
    """Load GSP report with caching and self-healing pickle.

    Priority:
    1. In-memory cache (instant)
    2. data_dir/output/GSP_v2.pkl.gz (compatible pickle, ~2-5s)
    3. data_dir/output/GSP.pkl.gz (original pickle, may fail)
    4. data_dir/output/GSP.xlsx (Excel fallback, ~30-60s)

    After xlsx load, writes GSP_v2.pkl.gz so future loads are fast.
    """
    cache_key = str(data_dir)
    if cache_key in _GSP_REPORT_CACHE:
        return _GSP_REPORT_CACHE[cache_key]

    output_dir = data_dir / "output"
    result = None

    # Try self-healed pickle first (our format, always compatible)
    pkl_v2 = output_dir / "GSP_v2.pkl.gz"
    if pkl_v2.exists():
        try:
            import gzip
            import pickle as _pickle
            with gzip.open(pkl_v2, 'rb') as f:
                result = _pickle.load(f)
            if not isinstance(result, dict):
                result = None
        except Exception:
            result = None

    # Try original pickle
    if result is None:
        pkl_path = output_dir / "GSP.pkl.gz"
        if pkl_path.exists():
            try:
                import gzip
                import pickle as _pickle
                with gzip.open(pkl_path, 'rb') as f:
                    result = _pickle.load(f)
                if not isinstance(result, dict):
                    result = None
            except Exception:
                result = None

    # Fallback: xlsx (slow but reliable)
    if result is None:
        xlsx_path = output_dir / "GSP.xlsx"
        if xlsx_path.exists():
            try:
                xls = pd.ExcelFile(xlsx_path, engine='openpyxl')
                result = {name: pd.read_excel(xls, sheet_name=name) for name in xls.sheet_names}
                # Self-heal: write compatible pickle for fast future loads
                try:
                    import gzip
                    import pickle as _pickle
                    import tempfile
                    tmp_path = pkl_v2.with_suffix('.tmp')
                    with gzip.open(tmp_path, 'wb') as f:
                        _pickle.dump(result, f, protocol=4)
                    tmp_path.rename(pkl_v2)
                except Exception:
                    pass
            except Exception:
                result = None

    if result is not None:
        _GSP_REPORT_CACHE[cache_key] = result
    return result


def clear_cache():
    """Clear the in-memory GSP report cache."""
    _GSP_REPORT_CACHE.clear()


def download_gsp_from_s3(
    s3_path: str,
    data_dir: Optional[str] = None,
    aws_profile: Optional[str] = None
) -> Path:
    """Download GSP report from S3.

    Downloads output/GSP.pkl.gz and output/GSP.xlsx from the given S3 path.

    Args:
        s3_path: S3 URI (e.g., s3://bucket/gsp/IBD_20260407/)
        data_dir: Local directory to store files (defaults to ~/tmp)
        aws_profile: AWS CLI profile name (e.g., 'cmp-dev')

    Returns:
        Path to the local data directory containing the downloaded files.
    """
    import subprocess

    if data_dir is None:
        data_dir = DATA_DIR
    data_dir = Path(data_dir)

    output_dir = data_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    stripped = s3_path.rstrip('/')
    if stripped.endswith('/output') or '/output/' in s3_path:
        base_s3 = stripped.rstrip('/output').rstrip('/') + '/output'
    else:
        base_s3 = stripped + '/output'

    for filename in ['GSP.pkl.gz', 'GSP.xlsx']:
        s3_uri = f"{base_s3}/{filename}"
        local_path = output_dir / filename
        if local_path.exists() and local_path.stat().st_size > 0:
            continue
        cmd = ['aws', 's3', 'cp', s3_uri, str(local_path)]
        if aws_profile:
            cmd.extend(['--profile', aws_profile])
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Warning: failed to download {s3_uri}: {e.stderr}")
            continue

    return data_dir


def load_dataset_registry() -> List[Dict[str, Any]]:
    """Load the GSP dataset registry from references/datasets.md.

    Returns:
        List of dataset dicts with keys: name, s3_pickle, s3_xlsx,
        aws_profile, local_dir, indications.
    """
    registry_path = REFERENCES_DIR / "datasets.md"
    if not registry_path.exists():
        return []

    try:
        import yaml
    except ImportError:
        # Parse YAML frontmatter manually if PyYAML not available
        text = registry_path.read_text()
        if not text.startswith('---'):
            return []
        end = text.index('---', 3)
        yaml_text = text[3:end].strip()
        # Minimal YAML parsing for our known structure
        import re
        datasets = []
        current = None
        current_indication = None
        for line in yaml_text.split('\n'):
            stripped = line.strip()
            if stripped.startswith('- name:'):
                # Flush pending indication
                if current_indication is not None and current is not None:
                    current['indications'].append(current_indication)
                    current_indication = None
                if current is not None:
                    datasets.append(current)
                current = {'name': stripped.split(':', 1)[1].strip(), 'indications': []}
            elif current is not None and stripped.startswith('- efo_id:'):
                # Flush pending indication before starting new one
                if current_indication is not None:
                    current['indications'].append(current_indication)
                current_indication = {'efo_id': stripped.split(':', 1)[1].strip()}
            elif current_indication is not None and stripped.startswith('name:'):
                current_indication['name'] = stripped.split(':', 1)[1].strip()
            elif current_indication is not None and stripped.startswith('aliases:'):
                current_indication['aliases'] = stripped.split(':', 1)[1].strip()
            elif current is not None and ':' in stripped and not stripped.startswith('-'):
                key, val = stripped.split(':', 1)
                current[key.strip()] = val.strip()
        # Flush remaining
        if current_indication is not None and current is not None:
            current['indications'].append(current_indication)
        if current is not None:
            datasets.append(current)
        return datasets

    text = registry_path.read_text()
    if not text.startswith('---'):
        return []
    end = text.index('---', 3)
    yaml_text = text[3:end].strip()
    data = yaml.safe_load(yaml_text)
    return data.get('datasets', [])


def resolve_dataset(user_input: str) -> Optional[Dict[str, Any]]:
    """Match user text to a dataset in the registry.

    Matches against dataset names and indication names (case-insensitive).

    Args:
        user_input: User text like "IBD", "Crohn's", "SSc", "IPF"

    Returns:
        Matched dataset dict, or None if no match.
    """
    registry = load_dataset_registry()
    query = user_input.strip().lower()

    # Exact dataset name match
    for ds in registry:
        if ds['name'].lower() == query:
            return ds

    # Partial/fuzzy indication name or alias match
    for ds in registry:
        for ind in ds.get('indications', []):
            if query in ind.get('name', '').lower():
                return ds
            if query == ind.get('efo_id', '').lower():
                return ds
            aliases = ind.get('aliases', '')
            if aliases:
                for alias in str(aliases).split(','):
                    if query == alias.strip().lower():
                        return ds

    return None


# ============================================================================
# DATA UTILITIES
# ============================================================================

def list_available_targets(
    data_dir: Optional[str] = None,
    verbose: bool = False
) -> Union[List[str], Dict[str, List[str]]]:
    """List all available GSP targets/datasets from both cloud and local sources.

    Args:
        data_dir: Optional data directory path for local xlsx files (defaults to ~/tmp)
        verbose: If True, return dict with sources separated; if False, return flat list

    Returns:
        If verbose=False (default): Sorted list of all available target names
        If verbose=True: Dict with 'cloud', 'local', and 'all' keys

    Example:
        >>> list_available_targets()
        ['ALS', 'NLRP3', 'TYK2']

        >>> list_available_targets(verbose=True)
        {'cloud': ['NLRP3', 'TYK2'], 'local': ['ALS', 'NLRP3'], 'all': ['ALS', 'NLRP3', 'TYK2']}
    """
    cloud_targets: set = set()
    local_targets: set = set()

    # Collect from REST API (cloud)
    if GSP_API_URL and _HAS_REQUESTS:
        try:
            response = requests.get(f"{GSP_API_URL}/datasets", timeout=GSP_API_TIMEOUT)
            if response.ok:
                cloud_targets = {d.get("gsp_prefix") for d in response.json() if d.get("gsp_prefix")}
        except Exception:
            pass

    # Collect from local xlsx files
    if data_dir is None:
        local_dir = DATA_DIR
    else:
        local_dir = Path(data_dir)

    if local_dir.exists():
        files = list(local_dir.glob("*_GSP.xlsx"))
        local_targets = {f.stem.replace("_GSP", "") for f in files}
        # New: directories containing output/GSP.pkl.gz or output/GSP.xlsx
        for subdir in local_dir.iterdir():
            if subdir.is_dir():
                output = subdir / "output"
                if (output / "GSP.pkl.gz").exists() or (output / "GSP.xlsx").exists():
                    local_targets.add(subdir.name)

    # Return based on verbose flag
    all_targets = cloud_targets | local_targets

    if verbose:
        return {
            "cloud": sorted(cloud_targets),
            "local": sorted(local_targets),
            "all": sorted(all_targets)
        }

    return sorted(all_targets)


def list_indications(
    data_dir: str
) -> List[Dict[str, Any]]:
    """List all disease indications available in a GSP dataset.

    Args:
        data_dir: Path to the dataset directory (e.g., data/IBD_20260407)

    Returns:
        List of dicts: [{efo_id, trait_name, study_count, gene_count}]
        sorted by study_count descending.
    """
    local_dir = Path(data_dir)
    report = _load_gsp_report(local_dir)
    if report is None:
        return []

    indications: Dict[str, Dict[str, Any]] = {}

    # From Study Info: study counts per disease
    study_df = report.get("Study Info", pd.DataFrame())
    if not study_df.empty and 'study_diseases' in study_df.columns:
        for efo_id, group in study_df.groupby('study_diseases'):
            efo_str = str(efo_id).strip()
            if efo_str and efo_str != 'nan':
                n_studies = int(group['study_id'].nunique()) if 'study_id' in group.columns else len(group)
                indications[efo_str] = {
                    'efo_id': efo_str,
                    'trait_name': '',
                    'study_count': n_studies,
                    'gene_count': 0
                }

    # From GWAS Summary: trait names and gene counts
    gwas_df = report.get("GWAS Summary", pd.DataFrame())
    if not gwas_df.empty and 'Trait EFO' in gwas_df.columns:
        for efo_id, group in gwas_df.groupby('Trait EFO'):
            efo_str = str(efo_id).strip()
            if efo_str and efo_str != 'nan':
                if efo_str not in indications:
                    indications[efo_str] = {
                        'efo_id': efo_str,
                        'trait_name': '',
                        'study_count': 0,
                        'gene_count': 0
                    }
                if 'Trait Reported' in group.columns:
                    names = group['Trait Reported'].dropna().unique()
                    if len(names) > 0:
                        indications[efo_str]['trait_name'] = str(names[0])
                if 'Gene Symbol' in group.columns:
                    indications[efo_str]['gene_count'] = int(group['Gene Symbol'].nunique())

    return sorted(indications.values(), key=lambda x: x['study_count'], reverse=True)


# ============================================================================
# TOOL 1: SAFETY RISK ASSESSMENT (phenotypes_for_gene)
# ============================================================================

def phenotypes_for_gene_xlsx(
    gsp_prefix: str,
    query_gene: str,
    cutoff: str = "Medium",
    data_dir: Optional[str] = None
) -> Dict[str, Any]:
    """Excel-based implementation - preserved for testing.

    Return all phenotypes associated with a gene for safety assessment.

    Loads the specified gsp_prefix GSP file and searches for phenotype
    associations from three data sources:
    - GWAS: Common variant associations
    - OMIM: Mendelian disease phenotypes with direction (LOF/GOF)
    - Gene Burden: Rare variant associations

    Args:
        gsp_prefix: Target protein name identifying which GSP file to use
                       (e.g., "NLRP3" loads NLRP3_GSP.xlsx)
        query_gene: Gene symbol to filter for (e.g., "NLRP3", "CASP1")
                   Can be any gene, including same as gsp_prefix
        cutoff: Minimum evidence strength threshold (default: "Medium")
                Options: "Strong", "Medium", "Low", "Very low", "Unclear", "None"
        data_dir: Optional data directory path (defaults to ~/tmp)

    Returns:
        Dictionary with:
        - gsp_prefix: The GSP file used
        - query_gene: The queried gene symbol
        - cutoff: The evidence cutoff used
        - gwas_results: List of GWAS associations meeting cutoff
        - omim_results: List of Mendelian phenotypes with direction
        - burden_results: List of rare variant burden associations

    Example:
        >>> result = phenotypes_for_gene_xlsx("NLRP3", "NLRP3", cutoff="Medium")
        >>> print(result['gwas_results'])
    """
    if data_dir is None:
        data_dir = DATA_DIR
    else:
        data_dir = Path(data_dir)

    # Normalize inputs
    gsp_prefix = gsp_prefix.strip().upper()
    query_gene = query_gene.strip().upper()

    if not query_gene:
        raise ValueError("query_gene is required and cannot be empty")

    if cutoff == "":
        cutoff = "Medium"

    if cutoff not in EVIDENCE_SCORE_MAP:
        raise ValueError(
            f"Invalid cutoff '{cutoff}'. Must be one of: {list(EVIDENCE_SCORE_MAP.keys())}"
        )

    # Load from the specific gsp_prefix file
    file_path = data_dir / f"{gsp_prefix}_GSP.xlsx"

    if not file_path.exists():
        available = list_available_targets(str(data_dir))
        raise FileNotFoundError(
            f"GSP file for '{gsp_prefix}' not found. "
            f"Available targets: {available}"
        )

    # Load and filter GWAS Summary
    gwas_df = pd.read_excel(file_path, sheet_name="GWAS Summary", engine="openpyxl")

    if "Final Evidence Strength" in gwas_df.columns:
        gwas_df["_evidence_num"] = (
            gwas_df["Final Evidence Strength"]
            .astype(str)
            .str.strip()
            .map(EVIDENCE_SCORE_MAP)
            .fillna(0.0)
        )
    else:
        gwas_df["_evidence_num"] = 0.0

    cutoff_score = EVIDENCE_SCORE_MAP[cutoff]
    filtered_gwas = gwas_df.query(
        '`Gene Symbol` == @query_gene & `_evidence_num` >= @cutoff_score'
    )
    gwas_results = (
        filtered_gwas[['Gene Symbol', 'Trait Reported', 'Final Evidence Strength']]
        .drop_duplicates()
        .to_dict('records')
    )

    # Load OMIM from the same file
    try:
        omim_df = pd.read_excel(file_path, sheet_name="OMIM", engine="openpyxl")
        omim_df = omim_df.query('`Approved Gene Symbol` == @query_gene')
        omim_results = (
            omim_df[['Approved Gene Symbol', 'Phenotype', 'Direction', 'Direction Source']]
            .rename(columns={'Approved Gene Symbol': 'Gene Symbol'})
            .drop_duplicates()
            .to_dict('records')
        )
    except Exception:
        omim_results = []

    # Load OT Gene Burden from the same file
    try:
        burden_df = pd.read_excel(file_path, sheet_name="OT Gene Burden", engine="openpyxl")
        burden_df = burden_df.query('`target.approvedSymbol` == @query_gene')
        burden_results = (
            burden_df[['target.approvedSymbol', 'diseaseFromSource', 'directionOnTarget', 'directionOnTrait']]
            .drop_duplicates()
            .to_dict('records')
        )
    except Exception:
        burden_results = []

    return {
        'gsp_prefix': gsp_prefix,
        'query_gene': query_gene,
        'cutoff': cutoff,
        'gwas_results': gwas_results,
        'omim_results': omim_results,
        'burden_results': burden_results,
    }


def phenotypes_for_gene(
    gsp_prefix: str,
    query_gene: str,
    cutoff: str = "Medium",
    data_dir: Optional[str] = None
) -> Dict[str, Any]:
    """Get all phenotypes associated with a gene for safety assessment.

    Data sources are tried in order:
    1. REST API (if GSP_API_URL environment variable is set)
    2. Django ORM (if running in Django context)
    3. Excel files (fallback)

    Args:
        gsp_prefix: Target protein name identifying which GSP dataset to use
                       (e.g., "NLRP3")
        query_gene: Gene symbol to filter for (e.g., "NLRP3", "CASP1")
                   Can be any gene, including same as gsp_prefix
        cutoff: Minimum evidence strength threshold (default: "Medium")
                Options: "Strong", "Medium", "Low", "Very low", "Unclear", "None"
        data_dir: Optional data directory path (defaults to ~/tmp)

    Returns:
        Dictionary with:
        - gsp_prefix: The GSP dataset used
        - query_gene: The queried gene symbol
        - cutoff: The evidence cutoff used
        - gwas_results: List of GWAS associations meeting cutoff
        - omim_results: List of Mendelian phenotypes with direction
        - burden_results: List of rare variant burden associations

    Example:
        >>> result = phenotypes_for_gene("NLRP3", "NLRP3", cutoff="Medium")
        >>> print(result['gwas_results'])
    """
    gsp_prefix = gsp_prefix.strip().upper()
    query_gene = query_gene.strip().upper()

    if not query_gene:
        raise ValueError("query_gene is required and cannot be empty")

    if cutoff == "":
        cutoff = "Medium"

    if cutoff not in EVIDENCE_SCORE_MAP:
        raise ValueError(
            f"Invalid cutoff '{cutoff}'. Must be one of: {list(EVIDENCE_SCORE_MAP.keys())}"
        )

    cutoff_score = EVIDENCE_SCORE_MAP.get(cutoff, 3.0)

    # Helper to process raw sheet data (same logic for API and ORM)
    def _process_phenotypes(gwas_data, omim_data, burden_data):
        # Process GWAS with deduplication
        gwas_seen = set()
        gwas_results = []
        for row in gwas_data:
            if (row.get('Gene Symbol') or '').upper() != query_gene:
                continue
            if EVIDENCE_SCORE_MAP.get(row.get('Final Evidence Strength', ''), 0) < cutoff_score:
                continue
            key = (row.get('Gene Symbol'), row.get('Trait Reported'), row.get('Final Evidence Strength'))
            if key not in gwas_seen:
                gwas_seen.add(key)
                gwas_results.append({
                    'Gene Symbol': row.get('Gene Symbol'),
                    'Trait Reported': row.get('Trait Reported'),
                    'Final Evidence Strength': row.get('Final Evidence Strength')
                })

        # Process OMIM with deduplication
        omim_seen = set()
        omim_results = []
        for row in omim_data:
            if (row.get('Approved Gene Symbol') or '').upper() != query_gene:
                continue
            key = (row.get('Approved Gene Symbol'), row.get('Phenotype'),
                   row.get('Direction'), row.get('Direction Source'))
            if key not in omim_seen:
                omim_seen.add(key)
                omim_results.append({
                    'Gene Symbol': row.get('Approved Gene Symbol'),
                    'Phenotype': row.get('Phenotype'),
                    'Direction': row.get('Direction'),
                    'Direction Source': row.get('Direction Source')
                })

        # Process Burden with deduplication
        burden_seen = set()
        burden_results = []
        for row in burden_data:
            if (row.get('target.approvedSymbol') or '').upper() != query_gene:
                continue
            key = (row.get('target.approvedSymbol'), row.get('diseaseFromSource'),
                   row.get('directionOnTarget'), row.get('directionOnTrait'))
            if key not in burden_seen:
                burden_seen.add(key)
                burden_results.append({
                    'target.approvedSymbol': row.get('target.approvedSymbol'),
                    'diseaseFromSource': row.get('diseaseFromSource'),
                    'directionOnTarget': row.get('directionOnTarget'),
                    'directionOnTrait': row.get('directionOnTrait')
                })

        return {
            'gsp_prefix': gsp_prefix,
            'query_gene': query_gene,
            'cutoff': cutoff,
            'gwas_results': gwas_results,
            'omim_results': omim_results,
            'burden_results': burden_results,
        }

    # Primary: REST API
    gwas_data = _fetch_sheet_via_api(gsp_prefix, "GWAS Summary")
    if gwas_data is not None:
        omim_data = _fetch_sheet_via_api(gsp_prefix, "OMIM") or []
        burden_data = _fetch_sheet_via_api(gsp_prefix, "OT Gene Burden") or []
        return _process_phenotypes(gwas_data, omim_data, burden_data)

    # Secondary: Pickle/XLSX report
    local_dir = Path(data_dir) if data_dir else DATA_DIR
    report = _load_gsp_report(local_dir)
    if report is not None:
        gwas_data = report.get("GWAS Summary", pd.DataFrame()).to_dict('records')
        omim_data = report.get("OMIM", pd.DataFrame()).to_dict('records')
        burden_data = report.get("OT Gene Burden", pd.DataFrame()).to_dict('records')
        return _process_phenotypes(gwas_data, omim_data, burden_data)

    # Fallback: Legacy xlsx files
    return phenotypes_for_gene_xlsx(gsp_prefix, query_gene, cutoff, data_dir)


def safety_assessment_summary(
    gsp_prefix: str,
    query_gene: str,
    cutoff: str = "Medium"
) -> str:
    """Generate a text summary of safety-relevant phenotypes for a gene.

    Args:
        gsp_prefix: Target protein name (GSP file to use)
        query_gene: Gene symbol to assess
        cutoff: Evidence strength cutoff (default: "Medium")

    Returns:
        Formatted text summary of safety findings
    """
    result = phenotypes_for_gene(gsp_prefix, query_gene, cutoff)

    lines = [
        f"Safety Assessment for {result['query_gene']} "
        f"(target: {result['gsp_prefix']}, cutoff: {result['cutoff']})"
    ]
    lines.append("=" * 70)

    lines.append(f"\nGWAS Associations ({len(result['gwas_results'])} found):")
    if result['gwas_results']:
        for r in result['gwas_results']:
            lines.append(f"  - {r['Trait Reported']}: {r['Final Evidence Strength']}")
    else:
        lines.append("  No associations meeting cutoff")

    lines.append(f"\nMendelian Phenotypes ({len(result['omim_results'])} found):")
    if result['omim_results']:
        for r in result['omim_results']:
            direction = r.get('Direction', 'Unknown')
            source = r.get('Direction Source', '')
            lines.append(f"  - {r['Phenotype']}: {direction} ({source})")
    else:
        lines.append("  No Mendelian phenotypes found")

    lines.append(f"\nGene Burden Evidence ({len(result['burden_results'])} found):")
    if result['burden_results']:
        for r in result['burden_results']:
            effect = r.get('directionOnTarget', 'Unknown')
            direction = r.get('directionOnTrait', 'Unknown')
            lines.append(f"  - {r['diseaseFromSource']}: {effect} - {direction}")
    else:
        lines.append("  No gene burden associations found")

    return "\n".join(lines)


# ============================================================================
# TOOL 2: BIOLOGICAL RISK ASSESSMENT (expand_gene_disease_pair)
# ============================================================================

def expand_gene_disease_pair_xlsx(
    gsp_prefix: str,
    query_gene: str,
    disease: str,
    data_dir: Optional[str] = None
) -> Dict[str, Any]:
    """Excel-based implementation - preserved for testing.

    Fetch all study-level associations for a gene-disease pair.

    Loads GWAS data from the specified gsp_prefix GSP file, and aggregates
    Study Info from ALL available GSP files (since studies may be recorded
    in multiple target files).

    Args:
        gsp_prefix: Target protein name identifying which GSP file to use
                       for GWAS data (e.g., "NLRP3" loads NLRP3_GSP.xlsx)
        query_gene: Gene symbol to filter for (e.g., "TYK2", "NLRP3")
        disease: Disease EFO/MONDO ID (e.g., "EFO_0000676" for Psoriasis)
        data_dir: Optional data directory path (defaults to ~/tmp)

    Returns:
        Dictionary with:
        - gsp_prefix: The GSP file used
        - gene: The queried gene symbol
        - disease_filter: The disease ID used
        - results: List of study-level associations with metadata

    Example:
        >>> result = expand_gene_disease_pair_xlsx("NLRP3", "NLRP3", "EFO_0004309")
        >>> print(f"Found {len(result['results'])} studies")
    """
    if data_dir is None:
        data_dir = DATA_DIR
    else:
        data_dir = Path(data_dir)

    # Normalize inputs
    gsp_prefix = gsp_prefix.strip().upper()
    query_gene = query_gene.strip().upper()
    disease = disease.strip()

    if not query_gene:
        raise ValueError("query_gene is required and cannot be empty")
    if not disease:
        raise ValueError("disease ID is required and cannot be empty")

    # Load from the specific gsp_prefix file
    file_path = data_dir / f"{gsp_prefix}_GSP.xlsx"

    if not file_path.exists():
        available = list_available_targets(str(data_dir))
        raise FileNotFoundError(
            f"GSP file for '{gsp_prefix}' not found. "
            f"Available targets: {available}"
        )

    # Load Study Info from ALL GSP files (studies may be in multiple files)
    gsp_files = glob.glob(str(data_dir / "*_GSP.xlsx"))
    study_dfs = []
    for fp in gsp_files:
        try:
            study_df = pd.read_excel(fp, sheet_name="Study Info", engine="openpyxl")
            study_df = study_df.query('study_diseases == @disease')
            if not study_df.empty:
                study_dfs.append(study_df[[
                    'study_id', 'study_diseases', 'study_publicationDate',
                    'study_hasSumstats', 'study_nSamples', 'study_nCases',
                    'study_cohorts', 'study_ldPopulationStructure'
                ]])
        except Exception:
            continue

    if study_dfs:
        study_info = pd.concat(study_dfs, join='inner', axis=0).drop_duplicates()
    else:
        study_info = pd.DataFrame()

    # Load and filter GWAS Summary from the target file
    gwas_df = pd.read_excel(file_path, sheet_name="GWAS Summary", engine="openpyxl")

    rows = (
        gwas_df[['Gene ID', 'Gene Symbol', 'Final Evidence Strength', 'Study ID', 'Trait Reported', 'Trait EFO']]
        .copy()
        .drop_duplicates()
        .query('`Gene Symbol` == @query_gene and `Trait EFO` == @disease')
    )

    if rows.empty:
        if study_info.empty:
            return {
                "gsp_prefix": gsp_prefix,
                "gene": query_gene,
                "disease_filter": disease,
                "results": [],
                "message": f"No evidence found for {query_gene} - {disease}"
            }
        # Return study info with "Not Sig." status
        study_info['Gene Symbol'] = query_gene
        study_info['Final Evidence Strength'] = 'Not Sig.'
        study_info['Trait Reported'] = ''
        study_info['Trait EFO'] = disease
        return {
            "gsp_prefix": gsp_prefix,
            "gene": query_gene,
            "disease_filter": disease,
            "results": study_info.to_dict('records')
        }

    # Merge with study info
    tdf = rows.merge(
        study_info,
        left_on='Study ID',
        right_on='study_id',
        how='outer'
    ).fillna('Not Sig.')

    tdf['study_cohorts'] = tdf['study_cohorts'].map(lambda x: '' if x == 'Not Sig.' else x)
    tdf['Gene ID'] = rows.iloc[0]['Gene ID']
    tdf['Gene Symbol'] = query_gene
    tdf['Trait Reported'] = rows.iloc[0]['Trait Reported']
    tdf['Trait EFO'] = disease

    return {
        "gsp_prefix": gsp_prefix,
        "gene": query_gene,
        "disease_filter": disease,
        "results": tdf.drop('Gene ID', axis=1).to_dict('records')
    }


def expand_gene_disease_pair(
    gsp_prefix: str,
    query_gene: str,
    disease: str,
    data_dir: Optional[str] = None
) -> Dict[str, Any]:
    """Fetch all study-level associations for a gene-disease pair.

    Data sources are tried in order:
    1. REST API (if GSP_API_URL environment variable is set)
    2. Django ORM (if running in Django context)
    3. Excel files (fallback)

    Args:
        gsp_prefix: Target protein name identifying which GSP dataset to use
                       for GWAS data (e.g., "NLRP3")
        query_gene: Gene symbol to filter for (e.g., "TYK2", "NLRP3")
        disease: Disease EFO/MONDO ID (e.g., "EFO_0000676" for Psoriasis)
        data_dir: Optional data directory path (defaults to ~/tmp)

    Returns:
        Dictionary with:
        - gsp_prefix: The GSP dataset used
        - gene: The queried gene symbol
        - disease_filter: The disease ID used
        - results: List of study-level associations with metadata

    Example:
        >>> result = expand_gene_disease_pair("NLRP3", "NLRP3", "EFO_0004309")
        >>> print(f"Found {len(result['results'])} studies")
    """
    gsp_prefix = gsp_prefix.strip().upper()
    query_gene = query_gene.strip().upper()
    disease = disease.strip()

    if not query_gene:
        raise ValueError("query_gene is required and cannot be empty")
    if not disease:
        raise ValueError("disease ID is required and cannot be empty")

    # Helper to process raw data (same logic for API and ORM)
    def _process_gene_disease(gwas_data, study_info_data):
        # Filter GWAS rows for gene and disease
        gwas_rows = [
            row for row in gwas_data
            if str(row.get('Gene Symbol', '')).strip().upper() == query_gene
            and str(row.get('Trait EFO', '')).strip() == disease
        ]

        # Filter and deduplicate study info
        study_info_rows = [
            {
                'study_id': row.get('study_id'),
                'study_diseases': row.get('study_diseases'),
                'study_publicationDate': row.get('study_publicationDate'),
                'study_hasSumstats': row.get('study_hasSumstats'),
                'study_nSamples': row.get('study_nSamples'),
                'study_nCases': row.get('study_nCases'),
                'study_cohorts': row.get('study_cohorts'),
                'study_ldPopulationStructure': row.get('study_ldPopulationStructure'),
            }
            for row in study_info_data
            if row.get('study_diseases') == disease
        ]

        seen_studies = set()
        unique_study_info = []
        for row in study_info_rows:
            if row['study_id'] not in seen_studies:
                seen_studies.add(row['study_id'])
                unique_study_info.append(row)

        if not gwas_rows:
            if not unique_study_info:
                return {
                    "gsp_prefix": gsp_prefix,
                    "gene": query_gene,
                    "disease_filter": disease,
                    "results": [],
                    "message": f"No evidence found for {query_gene} - {disease}"
                }
            results = []
            for study in unique_study_info:
                results.append({
                    'Gene Symbol': query_gene,
                    'Final Evidence Strength': 'Not Sig.',
                    'Trait Reported': '',
                    'Trait EFO': disease,
                    **study
                })
            return {
                "gsp_prefix": gsp_prefix,
                "gene": query_gene,
                "disease_filter": disease,
                "results": results
            }

        study_lookup = {s['study_id']: s for s in unique_study_info}
        results = []
        seen_study_ids = set()

        for row in gwas_rows:
            study_id = row.get('Study ID')
            seen_study_ids.add(study_id)

            result_row = {
                'Gene Symbol': row.get('Gene Symbol'),
                'Final Evidence Strength': row.get('Final Evidence Strength'),
                'Study ID': study_id,
                'Trait Reported': row.get('Trait Reported'),
                'Trait EFO': row.get('Trait EFO'),
            }

            if study_id in study_lookup:
                result_row.update(study_lookup[study_id])
            else:
                result_row.update({
                    'study_id': 'Not Sig.',
                    'study_diseases': 'Not Sig.',
                    'study_publicationDate': 'Not Sig.',
                    'study_hasSumstats': 'Not Sig.',
                    'study_nSamples': 'Not Sig.',
                    'study_nCases': 'Not Sig.',
                    'study_cohorts': '',
                    'study_ldPopulationStructure': 'Not Sig.',
                })

            results.append(result_row)

        trait_reported = gwas_rows[0].get('Trait Reported', '') if gwas_rows else ''
        for study_id, study in study_lookup.items():
            if study_id not in seen_study_ids:
                results.append({
                    'Gene Symbol': query_gene,
                    'Final Evidence Strength': 'Not Sig.',
                    'Study ID': study_id,
                    'Trait Reported': trait_reported,
                    'Trait EFO': disease,
                    **study
                })

        return {
            "gsp_prefix": gsp_prefix,
            "gene": query_gene,
            "disease_filter": disease,
            "results": results
        }

    # Primary: REST API
    gwas_data = _fetch_sheet_via_api(gsp_prefix, "GWAS Summary")
    if gwas_data is not None:
        study_info_data = _fetch_all_sheets_via_api("Study Info") or []
        return _process_gene_disease(gwas_data, study_info_data)

    # Secondary: Pickle/XLSX report
    local_dir = Path(data_dir) if data_dir else DATA_DIR
    report = _load_gsp_report(local_dir)
    if report is not None:
        gwas_data = report.get("GWAS Summary", pd.DataFrame()).to_dict('records')
        study_info_data = report.get("Study Info", pd.DataFrame()).to_dict('records')
        return _process_gene_disease(gwas_data, study_info_data)

    # Fallback: Legacy xlsx files
    return expand_gene_disease_pair_xlsx(gsp_prefix, query_gene, disease, data_dir)


def biological_assessment_summary(
    gsp_prefix: str,
    query_gene: str,
    disease: str,
    data_dir: Optional[str] = None
) -> str:
    """Generate a text summary of biological evidence for a gene-disease pair.

    Args:
        gsp_prefix: Target protein name (GSP file to use)
        query_gene: Gene symbol to assess
        disease: Disease EFO/MONDO ID
        data_dir: Optional data directory path

    Returns:
        Formatted text summary of biological evidence
    """
    result = expand_gene_disease_pair(gsp_prefix, query_gene, disease, data_dir)

    lines = [
        f"Biological Assessment for {result['gene']} - {result['disease_filter']} "
        f"(target: {result['gsp_prefix']})"
    ]
    lines.append("=" * 70)

    results = result['results']

    if not results:
        lines.append("\nNo studies found for this gene-disease pair.")
        return "\n".join(lines)

    evidence_counts = {}
    for r in results:
        strength = r.get('Final Evidence Strength', 'Unknown')
        evidence_counts[strength] = evidence_counts.get(strength, 0) + 1

    lines.append(f"\nTotal studies: {len(results)}")
    lines.append("Evidence distribution:")
    for strength in ['Strong', 'Medium', 'Low', 'Very low', 'Not Sig.']:
        if strength in evidence_counts:
            lines.append(f"  - {strength}: {evidence_counts[strength]}")

    lines.append("\nStudy Details:")
    lines.append("-" * 70)

    for r in results:
        study_id = r.get('Study ID', r.get('study_id', 'Unknown'))
        strength = r.get('Final Evidence Strength', 'Unknown')
        n_samples = r.get('study_nSamples', 'N/A')
        n_cases = r.get('study_nCases', 'N/A')
        pub_date = r.get('study_publicationDate', 'N/A')
        cohorts = r.get('study_cohorts', '')

        if isinstance(cohorts, list):
            cohorts = ', '.join(str(c) for c in cohorts)

        lines.append(f"\n  Study: {study_id}")
        lines.append(f"    Evidence: {strength}")
        lines.append(f"    Samples: {n_samples} (cases: {n_cases})")
        lines.append(f"    Published: {pub_date}")
        if cohorts:
            lines.append(f"    Cohorts: {cohorts}")

    return "\n".join(lines)


def biological_risk_score(
    results: List[Dict[str, Any]]
) -> Dict[str, Any]:
    """Apply the 7-rule AgenticBoost scoring algorithm to study results.

    Args:
        results: Study result dicts from expand_gene_disease_pair().

    Returns:
        Dict with risk_level, risk_label, cautionary_notes,
        high_quality_studies, total_studies, max_sample_size,
        evidence_distribution, ancestry_summary.
    """
    if not results:
        return {
            'risk_level': 'High risk',
            'risk_label': 'High risk: No studies found',
            'cautionary_notes': [],
            'high_quality_studies': [],
            'total_studies': 0,
            'max_sample_size': 0,
            'evidence_distribution': {},
            'ancestry_summary': {}
        }

    STRONG_LEVELS = {'Strong', 'Very strong'}
    WEAK_LEVELS = {'Not Sig.', 'Low', 'Very low'}

    def _get_sample_size(r):
        nc = r.get('study_nCases', 0)
        ns = r.get('study_nSamples', 0)
        try:
            nc = int(float(nc)) if str(nc) not in ('Not Sig.', 'nan', '', 'None') else 0
        except (ValueError, TypeError):
            nc = 0
        try:
            ns = int(float(ns)) if str(ns) not in ('Not Sig.', 'nan', '', 'None') else 0
        except (ValueError, TypeError):
            ns = 0
        return nc if nc > 0 else ns

    sample_sizes = [_get_sample_size(r) for r in results]
    max_ss = max(sample_sizes) if sample_sizes else 0
    threshold = 0.1 * max_ss

    hq_studies = [
        {**r, '_sample_size': ss}
        for r, ss in zip(results, sample_sizes)
        if ss >= threshold
    ]

    hq_strong = [s for s in hq_studies if s.get('Final Evidence Strength') in STRONG_LEVELS]
    hq_medium = [s for s in hq_studies if s.get('Final Evidence Strength') == 'Medium']

    evidence_dist = {}
    for r in results:
        strength = str(r.get('Final Evidence Strength', 'Unknown'))
        evidence_dist[strength] = evidence_dist.get(strength, 0) + 1

    ancestry_summary = {}
    for s in hq_studies:
        ld = s.get('study_ldPopulationStructure', '')
        if ld and str(ld) not in ('Not Sig.', 'nan', ''):
            sid = s.get('Study ID', s.get('study_id', 'unknown'))
            ancestry_summary[sid] = str(ld)

    cautionary_notes = []
    n_strong = len(hq_strong)
    n_medium = len(hq_medium)

    if n_strong == 0 and n_medium == 0:
        risk_level = 'High risk'
        risk_label = 'High risk: Low-to-no genetic evidence for this gene in this trait'
    elif n_strong == 0 and n_medium >= 2:
        risk_level = 'Medium risk'
        risk_label = 'Medium risk: preliminary associations'
    elif n_strong >= 1 and (n_strong + n_medium) >= 2:
        risk_level = 'Low risk'
        risk_label = 'Low risk: consistent genetic evidence across multiple high-quality studies'
    elif n_strong == 0 and n_medium == 1:
        risk_level = 'High risk'
        risk_label = 'High risk: single study with moderate evidence'
        cautionary_notes.append(
            'Only one high-quality study shows Medium evidence. '
            'Requires further expert review before drawing conclusions.'
        )
    elif n_strong >= 1 and n_medium == 0 and n_strong == 1:
        risk_level = 'Low risk'
        risk_label = 'Low risk: strong evidence from single high-quality study'
        cautionary_notes.append(
            'Only one high-quality study provides Strong/Very strong evidence. '
            'Further expert confirmation recommended.'
        )
    else:
        risk_level = 'Medium risk'
        risk_label = 'Medium risk: evidence pattern requires expert interpretation'
        cautionary_notes.append('Evidence pattern does not match standard rules. Expert review needed.')

    if (n_strong + n_medium) == 1 and len(ancestry_summary) > 1:
        cautionary_notes.append(
            'Discrepancies among studies may be attributable to '
            'differences in cohort ancestries.'
        )

    return {
        'risk_level': risk_level,
        'risk_label': risk_label,
        'cautionary_notes': cautionary_notes,
        'high_quality_studies': hq_studies,
        'total_studies': len(results),
        'max_sample_size': max_ss,
        'evidence_distribution': evidence_dist,
        'ancestry_summary': ancestry_summary
    }


_REPORT_DISCLAIMER = (
    "> **DISCLAIMER:** The risk assessment and evidence presented in this report "
    "are based on raw GSP (Genetics Support Profiler) scores derived from "
    "automated computational analysis. These results should be interpreted "
    "with caution and are not a substitute for expert review. The scoring "
    "reflects statistical associations only and does not establish causality."
)


def generate_biological_report(
    data_dir: str,
    query_gene: str,
    diseases: Optional[List[str]] = None,
    output_dir: Optional[str] = None,
    dataset_name: Optional[str] = None
) -> str:
    """Generate a Markdown biological risk report for a gene.

    Args:
        data_dir: Path to dataset directory (e.g., data/IBD_20260407)
        query_gene: Gene symbol to assess
        diseases: List of EFO IDs, or None to auto-discover all
        output_dir: Output directory (defaults to {cwd}/gsp_{dataset}/)
        dataset_name: Dataset label for filenames (e.g., "IBD")

    Returns:
        Path to generated Markdown report file.
    """
    import datetime

    local_dir = Path(data_dir)
    query_gene = query_gene.strip().upper()

    if dataset_name is None:
        dataset_name = local_dir.name

    if output_dir is None:
        out_path = Path.cwd() / f"gsp_{dataset_name}"
    else:
        out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    if diseases is None:
        indications = list_indications(data_dir)
        diseases = [ind['efo_id'] for ind in indications]

    sections = []
    sections.append(f"# Biological Risk Report: {query_gene}\n")
    sections.append(_REPORT_DISCLAIMER)
    sections.append(f"\n**Generated:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}")
    sections.append(f"**Dataset:** {dataset_name}")
    sections.append(f"**Indications assessed:** {len(diseases)}\n")

    sections.append("## Summary\n")
    sections.append("| Indication | EFO ID | Risk Level | HQ Studies | Max N |")
    sections.append("|---|---|---|---|---|")

    detail_sections = []

    for efo_id in diseases:
        result = expand_gene_disease_pair(
            gsp_prefix=dataset_name,
            query_gene=query_gene,
            disease=efo_id,
            data_dir=data_dir
        )

        studies = result.get('results', [])
        score = biological_risk_score(studies)

        trait_name = efo_id
        for s in studies:
            t = s.get('Trait Reported', '')
            if t and str(t) not in ('Not Sig.', 'nan', ''):
                trait_name = str(t)
                break

        sections.append(
            f"| {trait_name} | {efo_id} | {score['risk_level']} | "
            f"{len(score['high_quality_studies'])} | {score['max_sample_size']:,} |"
        )

        detail = [f"\n---\n\n## {trait_name} ({efo_id})\n"]
        detail.append(f"**Risk Level:** {score['risk_label']}")
        detail.append(f"**Total studies:** {score['total_studies']}, "
                      f"High-quality: {len(score['high_quality_studies'])}")
        detail.append(f"**Max sample size:** {score['max_sample_size']:,}")

        if score['evidence_distribution']:
            detail.append("\n**Evidence distribution:**")
            for level in ['Very strong', 'Strong', 'Medium', 'Low', 'Very low', 'Not Sig.']:
                count = score['evidence_distribution'].get(level, 0)
                if count > 0:
                    detail.append(f"- {level}: {count}")

        if score['cautionary_notes']:
            detail.append("\n**Cautionary notes:**")
            for note in score['cautionary_notes']:
                detail.append(f"- {note}")

        if score['ancestry_summary']:
            detail.append("\n**Ancestry coverage (HQ studies):**")
            for sid, anc in score['ancestry_summary'].items():
                detail.append(f"- {sid}: {anc}")

        if studies:
            detail.append("\n### Study Details\n")
            detail.append("| Study ID | Evidence | nSamples | nCases | Published | LD Population |")
            detail.append("|---|---|---|---|---|---|")
            for s in studies:
                sid = s.get('Study ID', s.get('study_id', ''))
                ev = s.get('Final Evidence Strength', '')
                ns = s.get('study_nSamples', '')
                nc = s.get('study_nCases', '')
                pub = s.get('study_publicationDate', '')
                ld = s.get('study_ldPopulationStructure', '')
                detail.append(f"| {sid} | {ev} | {ns} | {nc} | {pub} | {ld} |")

        detail_sections.extend(detail)

    sections.extend(detail_sections)

    report_file = out_path / f"{query_gene}_{dataset_name}_gsp_risk.md"
    report_file.write_text("\n".join(sections))

    return str(report_file)


# ============================================================================
# COMBO / MULTI-GENE SCORING
# ============================================================================

# Score mapping: biological_risk_score() risk_level -> (user_label, numeric)
# "Low risk" = strong genetic evidence = high confidence in target
RISK_TO_USER_SCORE = {
    'Low risk': ('High', 100),
    'Medium risk': ('Medium', 50),
    'High risk': ('Low', 25),
}


def discover_dataset_genes(data_dir: str) -> set:
    """Find which gene symbols exist in the GWAS Summary sheet.

    Args:
        data_dir: Path to dataset directory (e.g., data/IBD_20260407)

    Returns:
        Set of uppercase gene symbols found in GWAS Summary.
    """
    report = _load_gsp_report(Path(data_dir))
    gwas_df = report.get('GWAS Summary')
    if gwas_df is None or gwas_df.empty:
        return set()
    return set(gwas_df['Gene Symbol'].dropna().str.strip().str.upper().unique())


def load_combos_from_file(file_path: str) -> tuple:
    """Load gene lists and combo definitions from a TSV or JSON file.

    TSV format: columns combo_id, gene (one row per gene per combo)
    JSON format: list of lists, e.g. [["TYK2","JAK1"],["OSMR","IL23A"]]

    Args:
        file_path: Path to input file (.tsv, .csv, or .json)

    Returns:
        Tuple of (unique_genes: list[str], combos: list[list[str]])
    """
    import csv as _csv

    fp = Path(file_path)
    ext = fp.suffix.lower()

    if ext == '.json':
        with open(fp) as f:
            combos = json.load(f)
        if not isinstance(combos, list) or not all(isinstance(c, list) for c in combos):
            raise ValueError("JSON must be a list of lists, e.g. [[\"A\",\"B\"],[\"C\",\"D\"]]")
        combos = [[g.strip().upper() for g in c] for c in combos]
    elif ext in ('.tsv', '.csv'):
        delimiter = '\t' if ext == '.tsv' else ','
        combo_map = {}
        with open(fp, newline='') as f:
            reader = _csv.DictReader(f, delimiter=delimiter)
            for row in reader:
                cid = row.get('combo_id', '').strip()
                gene = row.get('gene', '').strip().upper()
                if cid and gene:
                    combo_map.setdefault(cid, []).append(gene)
        combos = list(combo_map.values())
    else:
        raise ValueError(f"Unsupported file format: {ext}. Use .json, .tsv, or .csv")

    unique_genes = sorted(set(g for combo in combos for g in combo))
    return unique_genes, combos


def score_gene_across_phenotypes(
    gene: str,
    phenotypes: List[tuple],
    data_dir: str,
    dataset_name: str,
    dataset_genes: set,
    aggregation: str = 'average',
) -> Dict[str, Any]:
    """Score a gene across multiple phenotypes using biological_risk_score.

    Args:
        gene: Gene symbol (e.g., "TYK2")
        phenotypes: List of (efo_id, short_name, full_name) tuples
        data_dir: Path to dataset directory
        dataset_name: Dataset label (e.g., "IBD")
        dataset_genes: Set of genes present in GWAS Summary
        aggregation: 'average' or 'best' across phenotypes

    Returns:
        Dict with gene, in_dataset, final_score, final_label, per_phenotype details.
    """
    gene_upper = gene.strip().upper()
    per_phenotype = {}

    if gene_upper not in dataset_genes:
        for efo_id, short_name, full_name in phenotypes:
            per_phenotype[efo_id] = {
                'phenotype_short': short_name,
                'phenotype_full': full_name,
                'risk_level': 'N/A',
                'risk_label': 'Gene not found in dataset',
                'label': 'None',
                'numeric': 0,
                'total_studies': 0,
                'hq_studies': 0,
                'evidence_distribution': {},
                'cautionary_notes': [],
                'ancestry_summary': {},
                'max_sample_size': 0,
                'hq_study_ids': [],
            }
        return {
            'gene': gene_upper,
            'in_dataset': False,
            'final_score': 0.0,
            'final_label': 'None',
            'per_phenotype': per_phenotype,
        }

    for efo_id, short_name, full_name in phenotypes:
        try:
            result = expand_gene_disease_pair(
                gsp_prefix=dataset_name,
                query_gene=gene_upper,
                disease=efo_id,
                data_dir=data_dir,
            )
            studies = result.get('results', [])
            score = biological_risk_score(studies)
            label, numeric = RISK_TO_USER_SCORE.get(score['risk_level'], ('Low', 25))

            hq_ids = [
                s.get('Study ID', s.get('study_id', ''))
                for s in score['high_quality_studies']
            ]

            per_phenotype[efo_id] = {
                'phenotype_short': short_name,
                'phenotype_full': full_name,
                'risk_level': score['risk_level'],
                'risk_label': score['risk_label'],
                'label': label,
                'numeric': numeric,
                'total_studies': score['total_studies'],
                'hq_studies': len(score['high_quality_studies']),
                'evidence_distribution': score['evidence_distribution'],
                'cautionary_notes': score['cautionary_notes'],
                'ancestry_summary': score['ancestry_summary'],
                'max_sample_size': score['max_sample_size'],
                'hq_study_ids': hq_ids,
            }
        except Exception as e:
            per_phenotype[efo_id] = {
                'phenotype_short': short_name,
                'phenotype_full': full_name,
                'risk_level': 'Error',
                'risk_label': str(e),
                'label': 'Unclear',
                'numeric': 0,
                'total_studies': 0,
                'hq_studies': 0,
                'evidence_distribution': {},
                'cautionary_notes': [],
                'ancestry_summary': {},
                'max_sample_size': 0,
                'hq_study_ids': [],
            }

    scores = [p['numeric'] for p in per_phenotype.values()]

    if aggregation == 'best':
        final_score = float(max(scores)) if scores else 0.0
    else:
        final_score = round(sum(scores) / len(scores), 1) if scores else 0.0

    if final_score >= 75:
        final_label = 'High'
    elif final_score >= 37.5:
        final_label = 'Medium'
    elif final_score > 0:
        final_label = 'Low'
    else:
        final_label = 'None'

    return {
        'gene': gene_upper,
        'in_dataset': True,
        'final_score': final_score,
        'final_label': final_label,
        'per_phenotype': per_phenotype,
    }


def compute_combo_scores(
    gene_results: Dict[str, Dict],
    combos: List[List[str]],
    phenotypes: List[tuple],
) -> List[Dict[str, Any]]:
    """Compute combo scores as mean of member gene final_scores.

    Args:
        gene_results: Dict mapping gene symbol -> score_gene_across_phenotypes result
        combos: List of gene-symbol lists (e.g., [["TYK2","JAK1"], ...])
        phenotypes: List of (efo_id, short_name, full_name) tuples

    Returns:
        List of combo result dicts, sorted by combo_score descending.
    """
    combo_results = []
    for combo in combos:
        members = []
        for gene in combo:
            g = gene.strip().upper()
            gr = gene_results.get(g, {
                'final_score': 0, 'final_label': 'None',
                'per_phenotype': {
                    efo: {'numeric': 0, 'label': 'None', 'phenotype_short': sn}
                    for efo, sn, _ in phenotypes
                },
            })
            members.append({
                'gene': g,
                'final_score': gr['final_score'],
                'final_label': gr['final_label'],
                'per_phenotype': {
                    efo: {
                        'score': gr['per_phenotype'].get(efo, {}).get('numeric', 0),
                        'label': gr['per_phenotype'].get(efo, {}).get('label', 'None'),
                    }
                    for efo, _, _ in phenotypes
                },
            })

        combo_score = round(
            sum(m['final_score'] for m in members) / len(members), 1
        )

        per_phenotype_combo = {}
        for efo_id, short_name, _ in phenotypes:
            pheno_scores = [m['per_phenotype'][efo_id]['score'] for m in members]
            pheno_avg = round(sum(pheno_scores) / len(pheno_scores), 1)
            per_phenotype_combo[efo_id] = {
                'phenotype_short': short_name,
                'combo_pheno_score': pheno_avg,
                'member_scores': [
                    {'gene': m['gene'], 'score': m['per_phenotype'][efo_id]['score'],
                     'label': m['per_phenotype'][efo_id]['label']}
                    for m in members
                ],
            }

        combo_results.append({
            'combo': combo,
            'combo_label': ' + '.join(combo),
            'size': len(combo),
            'combo_score': combo_score,
            'members': members,
            'per_phenotype': per_phenotype_combo,
        })

    return sorted(combo_results, key=lambda x: x['combo_score'], reverse=True)


def _fmt_evidence_dist(dist: dict) -> str:
    """Format evidence distribution dict as a compact string."""
    if not dist:
        return '--'
    order = ['Very strong', 'Strong', 'Medium', 'Low', 'Very low',
             'Unclear', 'Not Sig.', 'nan']
    parts = []
    for k in order:
        if k in dist:
            parts.append(f"{k}:{dist[k]}")
    for k in sorted(dist.keys()):
        if k not in order:
            parts.append(f"{k}:{dist[k]}")
    return ', '.join(parts) if parts else '--'


def export_individual_tsv(
    gene_results: Dict[str, Dict],
    phenotypes: List[tuple],
    output_path: str,
) -> str:
    """Export per-gene per-phenotype scores to TSV.

    Args:
        gene_results: Dict mapping gene symbol -> score_gene_across_phenotypes result
        phenotypes: List of (efo_id, short_name, full_name) tuples
        output_path: Path for the output TSV file

    Returns:
        Path to the written TSV file.
    """
    import csv as _csv

    pheno_ids = [efo for efo, _, _ in phenotypes]
    pheno_shorts = [sn for _, sn, _ in phenotypes]

    with open(output_path, 'w', newline='') as f:
        writer = _csv.writer(f, delimiter='\t')
        header = ['Rank', 'Gene', 'Final_Score', 'Final_Label', 'In_Dataset']
        for sn in pheno_shorts:
            header.extend([
                f'{sn}_Score', f'{sn}_Label', f'{sn}_RiskLevel',
                f'{sn}_RiskLabel', f'{sn}_TotalStudies',
                f'{sn}_HQStudies', f'{sn}_MaxSampleSize',
                f'{sn}_EvidenceDistribution', f'{sn}_CautionaryNotes',
                f'{sn}_Ancestry', f'{sn}_HQStudyIDs',
            ])
        header.append('Phenotypes')
        writer.writerow(header)

        sorted_genes = sorted(
            gene_results.values(),
            key=lambda g: g['final_score'],
            reverse=True,
        )

        phenotype_list_str = ', '.join(pheno_shorts)
        for rank, gr in enumerate(sorted_genes, 1):
            row = [rank, gr['gene'], gr['final_score'], gr['final_label'],
                   'Yes' if gr['in_dataset'] else 'No']
            for efo_id in pheno_ids:
                pp = gr['per_phenotype'].get(efo_id, {})
                anc = pp.get('ancestry_summary', {})
                anc_str = '; '.join(f"{k}={v}" for k, v in anc.items()) if anc else '--'
                hq_ids = pp.get('hq_study_ids', [])
                hq_str = '; '.join(str(s) for s in hq_ids) if hq_ids else '--'
                row.extend([
                    pp.get('numeric', 0),
                    pp.get('label', 'None'),
                    pp.get('risk_level', 'N/A'),
                    pp.get('risk_label', '--'),
                    pp.get('total_studies', 0),
                    pp.get('hq_studies', 0),
                    pp.get('max_sample_size', 0),
                    _fmt_evidence_dist(pp.get('evidence_distribution', {})),
                    '; '.join(pp.get('cautionary_notes', [])) or '--',
                    anc_str,
                    hq_str,
                ])
            row.append(phenotype_list_str)
            writer.writerow(row)

    return str(output_path)


def export_combo_tsv(
    combo_results: List[Dict],
    gene_results: Dict[str, Dict],
    phenotypes: List[tuple],
    output_path: str,
) -> str:
    """Export per-combo per-phenotype scores to TSV.

    Args:
        combo_results: Output from compute_combo_scores()
        gene_results: Dict mapping gene symbol -> score_gene_across_phenotypes result
        phenotypes: List of (efo_id, short_name, full_name) tuples
        output_path: Path for the output TSV file

    Returns:
        Path to the written TSV file.
    """
    import csv as _csv

    pheno_shorts = [sn for _, sn, _ in phenotypes]

    with open(output_path, 'w', newline='') as f:
        writer = _csv.writer(f, delimiter='\t')
        header = ['Rank', 'Combination', 'Size', 'Combo_Score', 'Member_Scores']
        for _, sn, _ in phenotypes:
            header.extend([
                f'{sn}_ComboScore', f'{sn}_MemberScores',
                f'{sn}_MemberEvidence',
            ])
        header.append('Phenotypes')
        writer.writerow(header)

        phenotype_list_str = ', '.join(pheno_shorts)
        for rank, cr in enumerate(combo_results, 1):
            member_scores_str = '; '.join(
                f"{m['gene']}:{m['final_score']}({m['final_label']})"
                for m in cr['members']
            )
            row = [rank, cr['combo_label'], cr['size'],
                   cr['combo_score'], member_scores_str]
            for efo_id, sn, _ in phenotypes:
                pp = cr['per_phenotype'][efo_id]
                pheno_member_str = '; '.join(
                    f"{ms['gene']}:{ms['score']}({ms['label']})"
                    for ms in pp['member_scores']
                )
                evidence_parts = []
                for ms in pp['member_scores']:
                    gr = gene_results.get(ms['gene'], {})
                    gpp = gr.get('per_phenotype', {}).get(efo_id, {})
                    ev = _fmt_evidence_dist(gpp.get('evidence_distribution', {}))
                    evidence_parts.append(f"{ms['gene']}=[{ev}]")
                row.extend([pp['combo_pheno_score'], pheno_member_str,
                            '; '.join(evidence_parts)])
            row.append(phenotype_list_str)
            writer.writerow(row)

    return str(output_path)


def generate_combo_report_template(
    gene_results: Dict[str, Dict],
    combo_results: List[Dict],
    phenotypes: List[tuple],
    dataset_name: str,
    output_dir: str,
    aggregation: str = 'average',
    include_disclaimer: bool = True,
) -> Dict[str, str]:
    """Generate structured markdown report templates for combo scoring.

    Produces two files:
    - gsp_{dataset}_individual_scores.md (gene-level)
    - gsp_{dataset}_combo_scores.md (combo-level)

    The reports contain structured tables and data sections. Claude should
    enhance them with per-gene biomedical interpretation and strategic analysis.

    Args:
        gene_results: Dict mapping gene -> score_gene_across_phenotypes result
        combo_results: Output from compute_combo_scores()
        phenotypes: List of (efo_id, short_name, full_name) tuples
        dataset_name: Dataset label (e.g., "IBD")
        output_dir: Output directory path
        aggregation: 'average' or 'best'
        include_disclaimer: Whether to include raw-scores disclaimer

    Returns:
        Dict with 'individual' and 'combo' keys pointing to file paths.
    """
    import datetime

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    pheno_shorts = [sn for _, sn, _ in phenotypes]
    pheno_list_str = ', '.join(pheno_shorts)
    agg_desc = 'Average' if aggregation == 'average' else 'Best'

    disclaimer = ""
    if include_disclaimer:
        disclaimer = (
            "\n> **RAW SCORES -- MANUAL REVIEW REQUIRED**\n>\n"
            "> All scores in this report are derived from **raw, automated GSP "
            "(Genetics Support Profiler) analysis** and have **NOT been manually "
            "curated**. The automated HQ-study filter and gene-to-locus mapping "
            "can significantly under- or over-estimate true genetic evidence. "
            "**These scores must be reviewed and adjusted by a genetics expert "
            "before use in target prioritization decisions.**\n"
        )

    now_str = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')

    header = (
        f"**Generated:** {now_str}\n"
        f"**Dataset:** {dataset_name}\n"
        f"**Phenotypes included:** {pheno_list_str}\n"
        f"**Aggregation:** {agg_desc} across {len(phenotypes)} phenotype(s)\n"
    )

    methodology = (
        "## Scoring Methodology\n\n"
        "| GSP Risk Level | Evidence Interpretation | Label | Score |\n"
        "|---|---|---|---|\n"
        "| Low risk | Strong, consistent genetic evidence across multiple HQ studies | High | 100 |\n"
        "| Medium risk | Preliminary associations; multiple medium-quality hits | Medium | 50 |\n"
        "| High risk | Weak or no genetic evidence in high-quality studies | Low | 25 |\n"
        "| Gene absent | Gene not found in GSP GWAS dataset | None | 0 |\n"
        "| Error/ambiguous | Cannot determine | Unclear | 0 |\n"
    )

    # --- Individual report ---
    ind_lines = [f"# {dataset_name} Individual Gene GSP Scoring Report\n"]
    if disclaimer:
        ind_lines.append(disclaimer)
    ind_lines.append(header)
    ind_lines.append(f"**Genes assessed:** {len(gene_results)}\n")
    ind_lines.append("---\n")
    ind_lines.append(methodology)
    ind_lines.append("---\n")

    ind_lines.append("## Individual Gene Scores -- Final Rankings\n")
    ind_header = "| Rank | Gene | Score | Label | In Dataset |"
    ind_sep = "|---|---|---|---|---|"
    for sn in pheno_shorts:
        ind_header += f" {sn} |"
        ind_sep += "---|"
    ind_lines.append(ind_header)
    ind_lines.append(ind_sep)

    sorted_genes = sorted(gene_results.values(), key=lambda g: g['final_score'], reverse=True)
    for rank, gr in enumerate(sorted_genes, 1):
        row = f"| {rank} | **{gr['gene']}** | {gr['final_score']} | {gr['final_label']} | {'Yes' if gr['in_dataset'] else 'No'} |"
        for efo_id, _, _ in phenotypes:
            pp = gr['per_phenotype'].get(efo_id, {})
            row += f" {pp.get('numeric', 0)} ({pp.get('label', 'None')}) |"
        ind_lines.append(row)

    ind_lines.append("\n---\n")
    ind_lines.append("## Per-Gene Detail\n")
    ind_lines.append("<!-- Claude: Add biomedical interpretation for each gene below -->\n")

    for gr in sorted_genes:
        if not gr['in_dataset']:
            continue
        ind_lines.append(f"### {gr['gene']} (Score: {gr['final_score']} -- {gr['final_label']})\n")
        ind_lines.append("| Phenotype | Score | Risk Level | Total Studies | HQ Studies | Evidence Distribution |")
        ind_lines.append("|---|---|---|---|---|---|")
        for efo_id, sn, _ in phenotypes:
            pp = gr['per_phenotype'].get(efo_id, {})
            ev = _fmt_evidence_dist(pp.get('evidence_distribution', {}))
            ind_lines.append(
                f"| **{sn}** | {pp.get('numeric', 0)} ({pp.get('label', 'None')}) "
                f"| {pp.get('risk_level', 'N/A')} | {pp.get('total_studies', 0)} "
                f"| {pp.get('hq_studies', 0)} | {ev} |"
            )
        ind_lines.append("")

    ind_lines.append("---\n")
    ind_lines.append("## Supplementary: Per-Phenotype Tables\n")
    for efo_id, sn, fn in phenotypes:
        ind_lines.append(f"### {fn} ({efo_id})\n")
        ind_lines.append("| Gene | Score | Label | Risk Level | Total Studies | HQ Studies | Evidence |")
        ind_lines.append("|---|---|---|---|---|---|---|")
        for gr in sorted_genes:
            pp = gr['per_phenotype'].get(efo_id, {})
            ev = _fmt_evidence_dist(pp.get('evidence_distribution', {}))
            ind_lines.append(
                f"| {gr['gene']} | {pp.get('numeric', 0)} | {pp.get('label', 'None')} "
                f"| {pp.get('risk_level', 'N/A')} | {pp.get('total_studies', 0)} "
                f"| {pp.get('hq_studies', 0)} | {ev} |"
            )
        ind_lines.append("")

    ind_lines.append("\n---\n")
    ind_lines.append("## Interpretation\n")
    ind_lines.append("<!-- Claude: Add comprehensive interpretation narrative here -->\n")

    ind_path = out / f"gsp_{dataset_name}_individual_scores.md"
    ind_path.write_text("\n".join(ind_lines))

    # --- Combo report ---
    combo_lines = [f"# {dataset_name} Gene Combination GSP Scoring Report\n"]
    if disclaimer:
        combo_lines.append(disclaimer)
    combo_lines.append(header)
    combo_lines.append(f"**Combinations assessed:** {len(combo_results)}\n")
    combo_lines.append("---\n")
    combo_lines.append(methodology)
    combo_lines.append("---\n")

    combo_lines.append("## Combination Rankings -- Final\n")
    c_header = "| Rank | Combination | Size | Combo Score | Member Scores |"
    c_sep = "|---|---|---|---|---|"
    for sn in pheno_shorts:
        c_header += f" {sn} |"
        c_sep += "---|"
    combo_lines.append(c_header)
    combo_lines.append(c_sep)

    for rank, cr in enumerate(combo_results, 1):
        ms_str = ', '.join(f"{m['gene']}:{m['final_score']}" for m in cr['members'])
        row = f"| {rank} | **{cr['combo_label']}** | {cr['size']} | **{cr['combo_score']}** | {ms_str} |"
        for efo_id, _, _ in phenotypes:
            pp = cr['per_phenotype'][efo_id]
            row += f" {pp['combo_pheno_score']} |"
        combo_lines.append(row)

    combo_lines.append("\n---\n")
    combo_lines.append("## Supplementary: Per-Phenotype Combo Scores\n")
    for efo_id, sn, fn in phenotypes:
        combo_lines.append(f"### {fn} ({efo_id})\n")
        combo_lines.append("| Rank | Combination | Score | Member Scores | Evidence |")
        combo_lines.append("|---|---|---|---|---|")
        for rank, cr in enumerate(combo_results, 1):
            pp = cr['per_phenotype'][efo_id]
            ms_str = ', '.join(f"{ms['gene']}:{ms['score']}({ms['label']})" for ms in pp['member_scores'])
            ev_parts = []
            for ms in pp['member_scores']:
                gr = gene_results.get(ms['gene'], {})
                gpp = gr.get('per_phenotype', {}).get(efo_id, {})
                ev_parts.append(f"{ms['gene']}=[{_fmt_evidence_dist(gpp.get('evidence_distribution', {}))}]")
            combo_lines.append(
                f"| {rank} | {cr['combo_label']} | {pp['combo_pheno_score']} | {ms_str} | {'; '.join(ev_parts)} |"
            )
        combo_lines.append("")

    combo_lines.append("\n---\n")
    combo_lines.append("## Interpretation\n")
    combo_lines.append("<!-- Claude: Add combo-level strategic interpretation here -->\n")

    combo_path = out / f"gsp_{dataset_name}_combo_scores.md"
    combo_path.write_text("\n".join(combo_lines))

    return {
        'individual': str(ind_path),
        'combo': str(combo_path),
    }


# ============================================================================
# CLI INTERFACE
# ============================================================================

def main():
    """CLI entry point for workflow-gsp."""
    import sys

    if len(sys.argv) < 2:
        print("GSP Tools - Genetics Support Profiler Analysis")
        print("=" * 50)
        print("\nUsage:")
        print("  workflow-gsp safety <gsp_prefix> <query_gene> [cutoff]")
        print("  workflow-gsp biological <gsp_prefix> <query_gene> <disease_id>")
        print("\nExamples:")
        print("  workflow-gsp safety NLRP3 NLRP3 Medium")
        print("  workflow-gsp biological NLRP3 NLRP3 EFO_0004309")
        print("\nAvailable targets:", list_available_targets())
        sys.exit(1)

    command = sys.argv[1].lower()

    if command == "safety":
        if len(sys.argv) < 4:
            print("Usage: workflow-gsp safety <gsp_prefix> <query_gene> [cutoff]")
            print("\nExample: workflow-gsp safety NLRP3 NLRP3 Medium")
            sys.exit(1)
        gsp_prefix = sys.argv[2]
        query_gene = sys.argv[3]
        cutoff = sys.argv[4] if len(sys.argv) > 4 else "Medium"
        print(safety_assessment_summary(gsp_prefix, query_gene, cutoff))

    elif command == "biological":
        if len(sys.argv) < 5:
            print("Usage: workflow-gsp biological <gsp_prefix> <query_gene> <disease_id>")
            print("\nExample: workflow-gsp biological NLRP3 NLRP3 EFO_0004309")
            print("\nCommon disease IDs:")
            print("  EFO_0000676  - Psoriasis")
            print("  EFO_0000685  - Rheumatoid arthritis")
            print("  EFO_0000384  - Crohn's disease")
            sys.exit(1)
        gsp_prefix = sys.argv[2]
        query_gene = sys.argv[3]
        disease = sys.argv[4]
        print(biological_assessment_summary(gsp_prefix, query_gene, disease))

    elif command == "indications":
        if len(sys.argv) < 3:
            print("Usage: workflow-gsp indications <data_dir>")
            sys.exit(1)
        data_dir = sys.argv[2]
        for ind in list_indications(data_dir):
            print(f"  {ind['efo_id']:20s} {ind['trait_name']:50s} ({ind['study_count']} studies, {ind['gene_count']} genes)")

    elif command == "report":
        if len(sys.argv) < 4:
            print("Usage: workflow-gsp report <data_dir> <query_gene> [--diseases EFO1,EFO2] [--output-dir DIR] [--dataset-name NAME]")
            sys.exit(1)
        data_dir = sys.argv[2]
        query_gene = sys.argv[3]
        diseases = None
        output_dir = None
        dataset_name = None
        i = 4
        while i < len(sys.argv):
            if sys.argv[i] == '--diseases' and i + 1 < len(sys.argv):
                diseases = [d.strip() for d in sys.argv[i + 1].split(',')]
                i += 2
            elif sys.argv[i] == '--output-dir' and i + 1 < len(sys.argv):
                output_dir = sys.argv[i + 1]
                i += 2
            elif sys.argv[i] == '--dataset-name' and i + 1 < len(sys.argv):
                dataset_name = sys.argv[i + 1]
                i += 2
            else:
                i += 1
        path = generate_biological_report(data_dir, query_gene, diseases, output_dir, dataset_name)
        print(f"Report generated: {path}")

    else:
        print(f"Unknown command: {command}")
        print("Use 'safety' or 'biological'")
        sys.exit(1)


if __name__ == "__main__":
    main()
