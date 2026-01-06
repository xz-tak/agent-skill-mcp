#!/usr/bin/env python3
"""
Census-based Gene Specificity Analysis

Query CELLxGENE Census and compute cell type-specific marker scores
using scanpy's rank_genes_groups with custom specificity metrics.

This script provides an alternative to the CellGuide-based marker extraction,
computing specificity scores directly from Census expression data.

Key Features:
- Query Census with tissue/disease/species filters (includes all cell types)
- Specify target cell type(s) for specificity calculation via name, ID, or substring
- Automatic cell type matching with confirmation
- Custom specificity scoring: lfc_logp, mean_specificity_ratio, normalized score

Usage examples:
    # Analyze specificity for fibroblasts in intestinal tissue
    python compute_specificity.py \
        --genes "IL11,GREM1,TYK2,JAK1" \
        --tissue intestine \
        --target-cell-type "fibroblast" \
        --output fibroblast_specificity

    # Use substring matching for plasma cells
    python compute_specificity.py \
        --genes "IL11,GREM1" \
        --tissue intestine \
        --target-cell-type "plasma" \
        --match-mode substring \
        --output plasma_specificity

    # Analyze with Cell Ontology ID
    python compute_specificity.py \
        --genes "APOE,APOC1" \
        --tissue liver \
        --target-cell-type "CL:0000057" \
        --match-mode id \
        --output hepatocyte_specificity

    # Analyze all cell types (no target filter)
    python compute_specificity.py \
        --genes "IL11,GREM1" \
        --tissue intestine \
        --output all_cells_specificity
"""

import argparse
import json
import logging
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Dict, Any, Union, Set

import numpy as np
import pandas as pd

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def match_cell_types(
    available_cell_types: List[str],
    target_pattern: str,
    match_mode: str = 'substring',
    cell_type_ids: Optional[Dict[str, str]] = None
) -> List[str]:
    """
    Match cell types based on pattern and mode.

    Parameters
    ----------
    available_cell_types : list
        List of cell type names available in the data
    target_pattern : str
        Pattern to match (name, substring, or Cell Ontology ID)
    match_mode : str
        Matching mode: 'exact', 'substring', 'regex', or 'id'
    cell_type_ids : dict, optional
        Mapping of cell type names to Cell Ontology IDs (for 'id' mode)

    Returns
    -------
    list
        List of matched cell type names
    """
    matched = []

    if match_mode == 'exact':
        # Exact match (case-insensitive)
        target_lower = target_pattern.lower()
        matched = [ct for ct in available_cell_types if ct.lower() == target_lower]

    elif match_mode == 'substring':
        # Substring match (case-insensitive)
        target_lower = target_pattern.lower()
        matched = [ct for ct in available_cell_types if target_lower in ct.lower()]

    elif match_mode == 'regex':
        # Regex pattern match
        try:
            pattern = re.compile(target_pattern, re.IGNORECASE)
            matched = [ct for ct in available_cell_types if pattern.search(ct)]
        except re.error as e:
            logger.error(f"Invalid regex pattern: {e}")
            return []

    elif match_mode == 'id':
        # Match by Cell Ontology ID
        if cell_type_ids:
            # Find cell types with matching ontology ID
            target_id = target_pattern.replace('_', ':')  # Normalize format
            matched = [ct for ct, ct_id in cell_type_ids.items()
                      if ct_id and target_id.lower() in ct_id.lower()]
        else:
            logger.warning("Cell type IDs not available for 'id' match mode")

    return matched


def confirm_cell_type_matches(
    matched_types: List[str],
    target_pattern: str,
    interactive: bool = True,
    group_cell_types: Optional[bool] = None
) -> tuple:
    """
    Display matched cell types and optionally get user confirmation.

    Parameters
    ----------
    matched_types : list
        List of matched cell type names
    target_pattern : str
        Original target pattern
    interactive : bool
        Whether to prompt for confirmation
    group_cell_types : bool, optional
        If True, group all matched cell types as one.
        If None and interactive, ask user.
        If False, treat each cell type separately.

    Returns
    -------
    tuple
        (confirmed_cell_types: list, should_group: bool)
    """
    if not matched_types:
        logger.warning(f"No cell types matched pattern: '{target_pattern}'")
        return [], False

    logger.info(f"\nMatched cell types for pattern '{target_pattern}':")
    for i, ct in enumerate(matched_types, 1):
        logger.info(f"  {i}. {ct}")

    confirmed_types = matched_types
    should_group = group_cell_types if group_cell_types is not None else False

    if interactive:
        print(f"\nFound {len(matched_types)} matching cell types.")
        print("Proceed with these cell types? (yes/no/select): ", end='')
        response = input().strip().lower()

        if response in ['no', 'n']:
            logger.info("Cell type selection cancelled")
            return [], False
        elif response == 'select':
            print("Enter numbers to select (comma-separated) or 'all': ", end='')
            selection = input().strip()
            if selection.lower() == 'all':
                confirmed_types = matched_types
            else:
                try:
                    indices = [int(x.strip()) - 1 for x in selection.split(',')]
                    confirmed_types = [matched_types[i] for i in indices if 0 <= i < len(matched_types)]
                except (ValueError, IndexError):
                    logger.error("Invalid selection")
                    return [], False

        # Ask about grouping if multiple cell types and not already specified
        if len(confirmed_types) > 1 and group_cell_types is None:
            print(f"\nGroup these {len(confirmed_types)} cell types as ONE for comparison? (yes/no): ", end='')
            group_response = input().strip().lower()
            should_group = group_response in ['yes', 'y']

    return confirmed_types, should_group


class SpecificityAnalyzer:
    """Compute cell type specificity from Census expression data."""

    def __init__(
        self,
        groupby: str = 'cell_type',
        target_cell_types: Optional[List[str]] = None,
        group_cell_types: bool = False,
        grouped_name: Optional[str] = None
    ):
        """
        Initialize the analyzer.

        Parameters
        ----------
        groupby : str
            Column to use for grouping cells (default: 'cell_type')
        target_cell_types : list, optional
            List of target cell types for specificity calculation.
            If None, calculate for all cell types.
        group_cell_types : bool
            If True, treat all target_cell_types as a single group for comparison.
        grouped_name : str, optional
            Name for the grouped cell type (default: 'target_group')
        """
        self.groupby = groupby
        self.target_cell_types = target_cell_types
        self.group_cell_types = group_cell_types
        self.grouped_name = grouped_name or 'target_group'

    def run_marker_analysis(
        self,
        adata,
        gene_subset: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """
        Run rank_genes_groups and extract marker DataFrame.

        Parameters
        ----------
        adata : AnnData
            AnnData object with expression data
        gene_subset : list, optional
            Subset of genes to analyze (if None, analyze all)

        Returns
        -------
        pd.DataFrame
            Marker gene results with statistics (unfiltered)
        """
        import scanpy as sc

        # Handle Census data where var_names may be numeric indices
        # Gene symbols are in 'feature_name' column
        if 'feature_name' in adata.var.columns:
            # Check if var_names are not gene symbols (e.g., numeric indices)
            sample_var_name = str(adata.var_names[0]) if len(adata.var_names) > 0 else ""
            if sample_var_name.isdigit() or sample_var_name.startswith('ENSG'):
                logger.info("Detected Census format - setting var_names to feature_name")
                adata.var_names = adata.var['feature_name'].values
                adata.var_names_make_unique()

        # Subset to genes of interest if specified
        if gene_subset:
            # Find genes that exist in the data (case-insensitive)
            var_names_upper = {str(g).upper(): g for g in adata.var_names}
            genes_found = []
            genes_missing = []
            for g in gene_subset:
                g_upper = str(g).upper()
                if g_upper in var_names_upper:
                    genes_found.append(var_names_upper[g_upper])
                else:
                    genes_missing.append(g)

            if genes_missing:
                logger.warning(f"Genes not found in data: {genes_missing}")

            if not genes_found:
                raise ValueError("None of the specified genes were found in the data")

            logger.info(f"Analyzing {len(genes_found)} of {len(gene_subset)} specified genes")
            adata = adata[:, genes_found].copy()

        # Ensure groupby column exists
        if self.groupby not in adata.obs.columns:
            available_cols = list(adata.obs.columns)
            raise ValueError(
                f"Groupby column '{self.groupby}' not found. "
                f"Available columns: {available_cols[:20]}..."
            )

        # Handle grouped cell type analysis
        analysis_groupby = self.groupby
        if self.group_cell_types and self.target_cell_types:
            logger.info(f"Grouping {len(self.target_cell_types)} cell types as '{self.grouped_name}'")

            # Create a new column for grouped analysis
            def assign_group(cell_type):
                """Assign cell type to grouped target or 'other'."""
                for target in self.target_cell_types:
                    if target.lower() in str(cell_type).lower():
                        return self.grouped_name
                return 'other'

            adata.obs['_grouped_cell_type'] = adata.obs[self.groupby].apply(assign_group)
            analysis_groupby = '_grouped_cell_type'

            # Log the grouping
            group_counts = adata.obs['_grouped_cell_type'].value_counts()
            logger.info(f"Grouped cell counts:")
            logger.info(f"  {self.grouped_name}: {group_counts.get(self.grouped_name, 0):,} cells")
            logger.info(f"  other: {group_counts.get('other', 0):,} cells")

        # Check number of groups
        n_groups = adata.obs[analysis_groupby].nunique()
        logger.info(f"Found {n_groups} unique groups in '{analysis_groupby}'")

        if n_groups < 2:
            raise ValueError(f"Need at least 2 groups for comparison, found {n_groups}")

        # Run rank_genes_groups with wilcoxon test
        logger.info(f"Running rank_genes_groups (method=wilcoxon, groupby={analysis_groupby})...")
        sc.tl.rank_genes_groups(
            adata,
            groupby=analysis_groupby,
            method='wilcoxon',
            pts=True,  # Include percentage of cells expressing
            use_raw=False
        )

        # Extract results as DataFrame - NO FILTERING
        logger.info("Extracting marker results (no filtering)...")
        marker_df = sc.get.rank_genes_groups_df(
            adata,
            group=None,  # Get all groups
            pval_cutoff=1.0,  # No p-value filter
            log2fc_min=None   # No fold change filter
        )

        logger.info(f"Found {len(marker_df)} marker entries (all genes × all cell types)")

        return marker_df

    def compute_lfc_logp(self, marker_df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate lfc_logp composite score and filter to target cell types.

        Scoring formula:
        - lfc_logp = logfoldchanges * -log(pvals_adj + 1e-323) *
                     (pct_nz_group - pct_nz_reference) * (1 - pct_nz_reference)

        Parameters
        ----------
        marker_df : pd.DataFrame
            Output from rank_genes_groups_df

        Returns
        -------
        pd.DataFrame
            DataFrame with lfc_logp scores, filtered to target cell types if specified
        """
        if len(marker_df) == 0:
            logger.warning("Empty marker DataFrame, returning empty result")
            return marker_df

        df = marker_df.copy()

        # Calculate lfc_logp composite score
        logger.info("Computing lfc_logp scores...")
        df['lfc_logp'] = (
            df['logfoldchanges']
            * -np.log(df['pvals_adj'] + 1e-323)
            * (df['pct_nz_group'] - df['pct_nz_reference'])
            * (1 - df['pct_nz_reference'])
        )

        # Filter to target cell types if specified
        if self.target_cell_types:
            logger.info(f"Filtering to target cell types: {self.target_cell_types}")

            # Create mask for target cell types
            def is_target_cell_type(cell_type):
                """Check if cell type matches any target pattern."""
                for target in self.target_cell_types:
                    if target.lower() in cell_type.lower():
                        return True
                return False

            # Mark and filter to target cell types
            df['is_target'] = df['group'].apply(is_target_cell_type)
            target_df = df[df['is_target']].copy()

            if len(target_df) == 0:
                logger.warning("No results found for target cell types")
                return pd.DataFrame()

            df = target_df.drop(columns=['is_target'])

        # Sort by lfc_logp descending
        df = df.sort_values('lfc_logp', ascending=False).reset_index(drop=True)

        # Reorder columns for clarity
        col_order = ['names', 'group', 'logfoldchanges', 'pvals', 'pvals_adj',
                     'pct_nz_group', 'pct_nz_reference', 'lfc_logp']
        # Only include columns that exist
        col_order = [c for c in col_order if c in df.columns]
        df = df[col_order]

        return df


def parse_genes(gene_input: str) -> List[str]:
    """
    Parse gene input from comma-separated string or file.

    Parameters
    ----------
    gene_input : str
        Comma-separated gene list or path to file with one gene per line

    Returns
    -------
    list
        List of gene symbols
    """
    gene_path = Path(gene_input)

    if gene_path.exists():
        # Read from file
        with open(gene_path, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
        logger.info(f"Loaded {len(genes)} genes from {gene_path}")
    else:
        # Parse comma-separated
        genes = [g.strip() for g in gene_input.split(',') if g.strip()]
        logger.info(f"Parsed {len(genes)} genes from input")

    return genes


def normalize_input(value: Union[str, List, Set, None]) -> Optional[Set[str]]:
    """Normalize input to a set of lowercase strings."""
    if value is None:
        return None
    if isinstance(value, str):
        return {value.lower()}
    if isinstance(value, (list, set)):
        return {str(v).lower() for v in value}
    return {str(value).lower()}


def query_census_data(
    species: str = 'human',
    tissue: Optional[List[str]] = None,
    disease: Optional[List[str]] = None,
    development_stage: Optional[List[str]] = None,
    gene_set: Optional[List[str]] = None
):
    """
    Query CELLxGENE Census for expression data.

    Parameters
    ----------
    species : str
        Species (human or mouse)
    tissue : list, optional
        Tissue filter(s)
    disease : list, optional
        Disease filter(s)
    development_stage : list, optional
        Development stage filter(s)
    gene_set : list, optional
        Gene symbols to query

    Returns
    -------
    AnnData
        Expression data with metadata
    """
    import cellxgene_census
    import tiledbsoma as soma

    # Map species to Census organism key
    species_map = {
        'human': 'homo_sapiens',
        'homo sapiens': 'homo_sapiens',
        'mouse': 'mus_musculus',
        'mus musculus': 'mus_musculus',
    }
    organism_key = species_map.get(species.lower(), 'homo_sapiens')

    logger.info(f"Opening Census connection (organism: {organism_key})...")

    with cellxgene_census.open_soma(census_version="stable") as census:
        experiment = census["census_data"][organism_key]

        # Load observation metadata to build filters
        logger.info("Loading observation metadata...")
        obs_df = experiment.obs.read().concat().to_pandas()
        logger.info(f"Total cells in Census: {len(obs_df):,}")

        # Build filter conditions
        filter_conditions = []

        # Development stage filter (default to adult)
        if development_stage is None:
            development_stage = ['adult']

        dev_stages = normalize_input(development_stage)
        if dev_stages:
            matching_stages = set()
            for stage in obs_df['development_stage'].unique():
                for ds in dev_stages:
                    if ds.lower() in str(stage).lower():
                        matching_stages.add(stage)

            if matching_stages:
                stage_filter = ' or '.join([f'development_stage == "{s}"' for s in matching_stages])
                filter_conditions.append(f'({stage_filter})')
                logger.info(f"Development stage filter: {len(matching_stages)} stages matched")

        # Tissue filter
        if tissue:
            tissue_set = normalize_input(tissue)
            matching_tissues = set()
            for t in obs_df['tissue_general'].unique():
                for ts in tissue_set:
                    if ts.lower() in str(t).lower():
                        matching_tissues.add(t)

            if matching_tissues:
                tissue_filter = ' or '.join([f'tissue_general == "{t}"' for t in matching_tissues])
                filter_conditions.append(f'({tissue_filter})')
                logger.info(f"Tissue filter: {matching_tissues}")
            else:
                logger.warning(f"No matching tissues found for: {tissue}")

        # Disease filter
        if disease:
            disease_set = normalize_input(disease)
            matching_diseases = set()
            for d in obs_df['disease'].unique():
                for ds in disease_set:
                    if ds.lower() in str(d).lower():
                        matching_diseases.add(d)

            if matching_diseases:
                disease_filter = ' or '.join([f'disease == "{d}"' for d in matching_diseases])
                filter_conditions.append(f'({disease_filter})')
                logger.info(f"Disease filter: {matching_diseases}")
            else:
                logger.warning(f"No matching diseases found for: {disease}")

        # Build value filter string
        value_filter = ' and '.join(filter_conditions) if filter_conditions else None

        # Build gene filter
        var_filter = None
        if gene_set:
            gene_symbols = [f'feature_name == "{gene}"' for gene in gene_set]
            var_filter = ' or '.join(gene_symbols)
            logger.info(f"Gene filter: {len(gene_set)} genes")

        logger.info(f"Query filter: {value_filter or 'None'}")

        # Execute query
        logger.info("Fetching data from Census...")
        with experiment.axis_query(
            measurement_name="RNA",
            obs_query=soma.AxisQuery(value_filter=value_filter) if value_filter else None,
            var_query=soma.AxisQuery(value_filter=var_filter) if var_filter else None,
        ) as query:
            adata = query.to_anndata(X_name="raw")

        logger.info(f"Retrieved {adata.n_obs:,} cells x {adata.n_vars:,} genes")

        return adata


def save_results(
    results_df: pd.DataFrame,
    summary: Dict[str, Any],
    output_prefix: str,
    output_dir: Path
):
    """Save analysis results to files."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save specificity results CSV
    csv_path = output_dir / f"{output_prefix}_specificity.csv"
    results_df.to_csv(csv_path, index=False)
    logger.info(f"Saved results to {csv_path}")

    # Save summary JSON
    summary_path = output_dir / f"{output_prefix}_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    logger.info(f"Saved summary to {summary_path}")

    # Generate and save markdown summary
    md_path = output_dir / f"{output_prefix}_summary.md"
    md_content = generate_markdown_summary(results_df, summary)
    with open(md_path, 'w') as f:
        f.write(md_content)
    logger.info(f"Saved markdown summary to {md_path}")


def generate_markdown_summary(results_df: pd.DataFrame, summary: Dict[str, Any]) -> str:
    """Generate a markdown summary report."""
    params = summary.get('parameters', {})
    results = summary.get('results', {})
    gene_summary = summary.get('gene_summary', {})

    # Check if grouped analysis
    is_grouped = params.get('grouped', False)
    grouped_name = params.get('grouped_name')

    lines = [
        "# Cell Type Specificity Analysis Report",
        "",
        f"**Generated:** {summary.get('timestamp', 'N/A')}",
        "",
        "## Query Parameters",
        "",
        f"- **Genes:** {', '.join(params.get('genes', []))}",
        f"- **Tissue:** {params.get('tissue', 'All')}",
        f"- **Disease:** {params.get('disease', 'All')}",
        f"- **Species:** {params.get('species', 'human')}",
        f"- **Target cell type:** {params.get('target_cell_type', 'All')}",
        f"- **Match mode:** {params.get('match_mode', 'substring')}",
        f"- **Grouped analysis:** {'Yes' if is_grouped else 'No'}",
    ]

    if is_grouped and grouped_name:
        lines.append(f"- **Grouped as:** {grouped_name}")

    lines.extend([
        "",
        "## Results Overview",
        "",
        f"- **Cells analyzed:** {results.get('n_cells', 0):,}",
        f"- **Genes queried:** {results.get('n_genes_queried', 0)}",
        f"- **Genes found:** {results.get('n_genes_found', 0)}",
        f"- **Cell type groups:** {results.get('n_groups', 0)}",
        f"- **Target groups:** {results.get('n_target_groups', 0)}",
        "",
    ])

    # Matched cell types
    if params.get('matched_cell_types'):
        lines.append("### Matched Cell Types")
        lines.append("")
        for ct in params.get('matched_cell_types', []):
            lines.append(f"- {ct}")
        lines.append("")

    # Per-gene summary
    if gene_summary:
        lines.append("## Per-Gene Summary")
        lines.append("")

        for gene, gdata in gene_summary.items():
            lines.append(f"### {gene}")
            lines.append("")
            if 'status' in gdata:
                lines.append(f"*{gdata['status']}*")
            elif 'grouped' in gdata or 'individual' in gdata:
                # Grouped + individual format
                if 'grouped' in gdata:
                    g = gdata['grouped']
                    lines.append(f"**GROUPED ({g.get('cell_type', 'N/A')}):**")
                    lines.append("")
                    lines.append(f"| Metric | Value |")
                    lines.append(f"|--------|-------|")
                    lines.append(f"| Log2FC | {g.get('log2fc', 0):.2f} |")
                    lines.append(f"| P-adj | {g.get('pval_adj', 1):.2e} |")
                    lines.append(f"| % Expressing | {g.get('pct_expressing', 0)*100:.1f}% |")
                    lines.append(f"| % Reference | {g.get('pct_reference', 0)*100:.2f}% |")
                    lines.append(f"| lfc_logp | {g.get('lfc_logp', 0):.2f} |")
                    lines.append("")

                if 'individual' in gdata:
                    ind = gdata['individual']
                    lines.append(f"**INDIVIDUAL BREAKDOWN** (best: {ind.get('best_cell_type', 'N/A')}):")
                    lines.append("")
                    breakdown = ind.get('breakdown', [])
                    if breakdown:
                        lines.append("| Cell Type | Log2FC | P-adj | % Expr | lfc_logp |")
                        lines.append("|-----------|--------|-------|--------|----------|")
                        for ct in breakdown:
                            pct = ct.get('pct_nz_group', 0) * 100
                            lines.append(
                                f"| {ct.get('group', 'N/A')} | "
                                f"{ct.get('logfoldchanges', 0):.2f} | "
                                f"{ct.get('pvals_adj', 1):.2e} | "
                                f"{pct:.1f}% | "
                                f"{ct.get('lfc_logp', 0):.2f} |"
                            )
                    lines.append("")
            else:
                # Standard format (no grouping)
                lines.append(f"- **Best cell type:** {gdata.get('best_cell_type', 'N/A')}")
                lines.append(f"- **Log2 fold change:** {gdata.get('best_log2fc', 0):.2f}")
                lines.append(f"- **Adjusted p-value:** {gdata.get('best_pval_adj', 1):.2e}")
                lines.append(f"- **% cells expressing:** {gdata.get('best_pct_expressing', 0)*100:.1f}%")
                lines.append(f"- **lfc_logp score:** {gdata.get('best_lfc_logp', 0):.2f}")
                lines.append(f"- **Detected in:** {gdata.get('n_cell_types_detected', 0)} cell type(s)")
                lines.append("")

                # Cell types table
                cell_types = gdata.get('cell_types', [])
                if cell_types:
                    lines.append("| Cell Type | Log2FC | P-adj | % Expr | lfc_logp |")
                    lines.append("|-----------|--------|-------|--------|----------|")
                    for ct in cell_types:
                        pct = ct.get('pct_nz_group', 0) * 100
                        lines.append(
                            f"| {ct.get('group', 'N/A')} | "
                            f"{ct.get('logfoldchanges', 0):.2f} | "
                            f"{ct.get('pvals_adj', 1):.2e} | "
                            f"{pct:.1f}% | "
                            f"{ct.get('lfc_logp', 0):.2f} |"
                        )
            lines.append("")

    # Top results table
    has_analysis_type = 'analysis_type' in results_df.columns if len(results_df) > 0 else False

    if len(results_df) > 0:
        lines.append("## All Results (by lfc_logp)")
        lines.append("")

        if has_analysis_type:
            lines.append("| Gene | Cell Type | Type | Log2FC | P-adj | % Expr | lfc_logp |")
            lines.append("|------|-----------|------|--------|-------|--------|----------|")
        else:
            lines.append("| Gene | Cell Type | Log2FC | P-adj | % Expr | lfc_logp |")
            lines.append("|------|-----------|--------|-------|--------|----------|")

        for _, row in results_df.head(25).iterrows():
            pct = row.get('pct_nz_group', 0) * 100
            if has_analysis_type:
                lines.append(
                    f"| {row.get('names', 'N/A')} | "
                    f"{row.get('group', 'N/A')} | "
                    f"{row.get('analysis_type', 'N/A')} | "
                    f"{row.get('logfoldchanges', 0):.2f} | "
                    f"{row.get('pvals_adj', 1):.2e} | "
                    f"{pct:.1f}% | "
                    f"{row.get('lfc_logp', 0):.2f} |"
                )
            else:
                lines.append(
                    f"| {row.get('names', 'N/A')} | "
                    f"{row.get('group', 'N/A')} | "
                    f"{row.get('logfoldchanges', 0):.2f} | "
                    f"{row.get('pvals_adj', 1):.2e} | "
                    f"{pct:.1f}% | "
                    f"{row.get('lfc_logp', 0):.2f} |"
                )
        lines.append("")

    # Interpretation notes
    lines.extend([
        "## Interpretation Guide",
        "",
        "- **Log2FC**: Log2 fold change comparing target cell type vs all other cell types",
        "- **P-adj**: Benjamini-Hochberg adjusted p-value from Wilcoxon rank-sum test",
        "- **% Expr**: Percentage of cells in the target cell type expressing the gene",
        "- **lfc_logp**: Composite score = log2FC × -log(p_adj) × (pct_target - pct_ref) × (1 - pct_ref)",
        "",
        "Higher lfc_logp indicates stronger, more specific expression in the target cell type.",
        ""
    ])

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description='Compute cell type specificity scores from CELLxGENE Census',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze specificity for fibroblasts in intestinal tissue
  python compute_specificity.py \\
    --genes "IL11,GREM1,TYK2,JAK1" \\
    --tissue intestine \\
    --target-cell-type "fibroblast" \\
    --output fibroblast_specificity

  # Use substring matching for plasma cells
  python compute_specificity.py \\
    --genes "IL11,GREM1" \\
    --tissue intestine \\
    --target-cell-type "plasma" \\
    --match-mode substring \\
    --output plasma_specificity

  # Analyze with Cell Ontology ID
  python compute_specificity.py \\
    --genes "APOE,APOC1" \\
    --tissue liver \\
    --target-cell-type "CL:0000057" \\
    --match-mode id \\
    --output hepatocyte_specificity

  # Use regex pattern for multiple cell types
  python compute_specificity.py \\
    --genes "IL11,GREM1" \\
    --tissue intestine \\
    --target-cell-type "fibroblast|myofibroblast|smooth muscle" \\
    --match-mode regex \\
    --output stromal_specificity

  # Analyze all cell types (no target filter)
  python compute_specificity.py \\
    --genes "IL11,GREM1" \\
    --tissue intestine \\
    --output all_cells_specificity
        """
    )

    # Required arguments
    parser.add_argument('--genes', type=str, required=True,
                        help='Gene list (comma-separated) or path to file with one gene per line')
    parser.add_argument('--output', type=str, required=True,
                        help='Output prefix for result files')

    # Filter arguments
    parser.add_argument('--tissue', type=str, default=None,
                        help='Tissue filter (comma-separated for multiple)')
    parser.add_argument('--disease', type=str, default=None,
                        help='Disease filter (comma-separated for multiple)')
    parser.add_argument('--species', type=str, default='human',
                        choices=['human', 'mouse'],
                        help='Species (default: human)')
    parser.add_argument('--development-stage', type=str, default=None,
                        help='Development stage filter (default: adult)')

    # Target cell type arguments
    parser.add_argument('--target-cell-type', type=str, default=None,
                        help='Target cell type(s) for specificity calculation. '
                             'Can be name, substring, regex, or Cell Ontology ID.')
    parser.add_argument('--match-mode', type=str, default='substring',
                        choices=['exact', 'substring', 'regex', 'id'],
                        help='Cell type matching mode (default: substring)')
    parser.add_argument('--group-cell-type', action='store_true',
                        help='Group all matched cell types as ONE for comparison against other cells. '
                             'If not set and multiple cell types match, each is analyzed separately.')
    parser.add_argument('--grouped-name', type=str, default=None,
                        help='Name for the grouped cell type (default: based on target pattern)')

    # Analysis arguments
    parser.add_argument('--group-by', type=str, default='cell_type',
                        help='Column for grouping cells (default: cell_type)')

    # Output arguments
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Output directory (default: current directory)')
    parser.add_argument('--no-interactive', action='store_true',
                        help='Skip interactive confirmation for cell type matching')

    args = parser.parse_args()

    # Parse inputs
    def parse_list_arg(arg):
        if arg is None:
            return None
        return [x.strip() for x in arg.split(',')]

    genes = parse_genes(args.genes)
    tissue = parse_list_arg(args.tissue)
    disease = parse_list_arg(args.disease)
    development_stage = parse_list_arg(args.development_stage)

    logger.info("=" * 80)
    logger.info("CELLxGENE Census Specificity Analysis")
    logger.info("=" * 80)
    logger.info(f"Genes: {genes}")
    logger.info(f"Tissue: {tissue}")
    logger.info(f"Disease: {disease}")
    logger.info(f"Species: {args.species}")
    logger.info(f"Development stage: {development_stage or ['adult']}")
    logger.info(f"Group by: {args.group_by}")
    logger.info(f"Target cell type: {args.target_cell_type or 'All'}")
    logger.info(f"Match mode: {args.match_mode}")
    logger.info("=" * 80)

    try:
        # Query Census data
        adata = query_census_data(
            species=args.species,
            tissue=tissue,
            disease=disease,
            development_stage=development_stage,
            gene_set=genes
        )

        if adata.n_obs == 0:
            logger.error("No cells found matching the filters")
            sys.exit(1)

        # Get available cell types from the data
        available_cell_types = adata.obs[args.group_by].unique().tolist()
        logger.info(f"Found {len(available_cell_types)} unique cell types in the data")

        # Match and confirm target cell types if specified
        target_cell_types = None
        should_group = False
        grouped_name = args.grouped_name

        if args.target_cell_type:
            # Build cell type ID mapping for 'id' mode
            cell_type_ids = None
            if args.match_mode == 'id' and 'cell_type_ontology_term_id' in adata.obs.columns:
                cell_type_ids = dict(zip(
                    adata.obs['cell_type'].unique(),
                    adata.obs.groupby('cell_type')['cell_type_ontology_term_id'].first()
                ))

            # Match cell types
            matched_types = match_cell_types(
                available_cell_types=available_cell_types,
                target_pattern=args.target_cell_type,
                match_mode=args.match_mode,
                cell_type_ids=cell_type_ids
            )

            # Confirm matches (and ask about grouping if interactive)
            # If --group-cell-type flag is set, force grouping
            group_cell_types_arg = True if args.group_cell_type else None
            target_cell_types, should_group = confirm_cell_type_matches(
                matched_types=matched_types,
                target_pattern=args.target_cell_type,
                interactive=not args.no_interactive,
                group_cell_types=group_cell_types_arg
            )

            if not target_cell_types:
                logger.error("No cell types selected. Exiting.")
                sys.exit(1)

            # Set default grouped name if grouping and not specified
            if should_group and not grouped_name:
                grouped_name = f"{args.target_cell_type}_grouped"

            if should_group:
                logger.info(f"GROUPING {len(target_cell_types)} cell types as '{grouped_name}'")
            else:
                logger.info(f"Selected {len(target_cell_types)} cell types for specificity analysis (separate)")

        # Run specificity analysis
        if should_group:
            # Run BOTH grouped and individual analyses
            logger.info("Running grouped + individual breakdown analysis...")

            # 1. Grouped analysis (combined cell types vs others)
            grouped_analyzer = SpecificityAnalyzer(
                groupby=args.group_by,
                target_cell_types=target_cell_types,
                group_cell_types=True,
                grouped_name=grouped_name
            )
            grouped_marker_df = grouped_analyzer.run_marker_analysis(
                adata.copy(),  # Use copy to avoid modifying original
                gene_subset=genes
            )
            if len(grouped_marker_df) > 0:
                grouped_results = grouped_analyzer.compute_lfc_logp(grouped_marker_df)
                grouped_results['analysis_type'] = 'grouped'
            else:
                grouped_results = pd.DataFrame()

            # 2. Individual analysis (each cell type separately)
            logger.info("Running individual cell type breakdown...")
            individual_analyzer = SpecificityAnalyzer(
                groupby=args.group_by,
                target_cell_types=target_cell_types,
                group_cell_types=False
            )
            individual_marker_df = individual_analyzer.run_marker_analysis(
                adata,
                gene_subset=genes
            )
            if len(individual_marker_df) > 0:
                individual_results = individual_analyzer.compute_lfc_logp(individual_marker_df)
                individual_results['analysis_type'] = 'individual'
            else:
                individual_results = pd.DataFrame()

            # Combine results: grouped first, then individual breakdown
            results_dfs = []
            if len(grouped_results) > 0:
                results_dfs.append(grouped_results)
            if len(individual_results) > 0:
                results_dfs.append(individual_results)

            if results_dfs:
                results_df = pd.concat(results_dfs, ignore_index=True)
                # Reorder columns to put analysis_type near the front
                cols = results_df.columns.tolist()
                cols.remove('analysis_type')
                cols.insert(2, 'analysis_type')
                results_df = results_df[cols]
            else:
                results_df = pd.DataFrame()
        else:
            # Standard individual analysis only
            analyzer = SpecificityAnalyzer(
                groupby=args.group_by,
                target_cell_types=target_cell_types,
                group_cell_types=False
            )
            marker_df = analyzer.run_marker_analysis(
                adata,
                gene_subset=genes
            )
            if len(marker_df) == 0:
                logger.warning("No marker results found")
                results_df = pd.DataFrame()
            else:
                results_df = analyzer.compute_lfc_logp(marker_df)

        # Generate per-gene summary for queried genes
        gene_summary = {}
        has_analysis_type = 'analysis_type' in results_df.columns if len(results_df) > 0 else False

        if len(results_df) > 0:
            for gene in genes:
                gene_data = results_df[results_df['names'].str.upper() == gene.upper()]
                if len(gene_data) > 0:
                    gene_entry = {}

                    if has_analysis_type and should_group:
                        # Separate grouped and individual results
                        grouped_data = gene_data[gene_data['analysis_type'] == 'grouped']
                        individual_data = gene_data[gene_data['analysis_type'] == 'individual']

                        # Grouped summary
                        if len(grouped_data) > 0:
                            grouped_row = grouped_data.iloc[0]
                            gene_entry['grouped'] = {
                                'cell_type': grouped_row['group'],
                                'log2fc': float(grouped_row['logfoldchanges']),
                                'pval_adj': float(grouped_row['pvals_adj']),
                                'pct_expressing': float(grouped_row['pct_nz_group']),
                                'pct_reference': float(grouped_row['pct_nz_reference']),
                                'lfc_logp': float(grouped_row['lfc_logp']),
                            }

                        # Individual breakdown
                        if len(individual_data) > 0:
                            best_indiv = individual_data.loc[individual_data['lfc_logp'].idxmax()]
                            gene_entry['individual'] = {
                                'best_cell_type': best_indiv['group'],
                                'best_log2fc': float(best_indiv['logfoldchanges']),
                                'best_pval_adj': float(best_indiv['pvals_adj']),
                                'best_pct_expressing': float(best_indiv['pct_nz_group']),
                                'best_lfc_logp': float(best_indiv['lfc_logp']),
                                'n_cell_types_detected': len(individual_data[individual_data['pct_nz_group'] > 0]),
                                'breakdown': individual_data[['group', 'logfoldchanges', 'pvals_adj',
                                                              'pct_nz_group', 'lfc_logp']].to_dict('records')
                            }

                        gene_summary[gene] = gene_entry
                    else:
                        # Standard summary (no grouping)
                        best_row = gene_data.loc[gene_data['lfc_logp'].idxmax()]
                        gene_summary[gene] = {
                            'best_cell_type': best_row['group'],
                            'best_log2fc': float(best_row['logfoldchanges']),
                            'best_pval_adj': float(best_row['pvals_adj']),
                            'best_pct_expressing': float(best_row['pct_nz_group']),
                            'best_lfc_logp': float(best_row['lfc_logp']),
                            'n_cell_types_detected': len(gene_data[gene_data['pct_nz_group'] > 0]),
                            'cell_types': gene_data[['group', 'logfoldchanges', 'pvals_adj',
                                                      'pct_nz_group', 'lfc_logp']].to_dict('records')
                        }
                else:
                    gene_summary[gene] = {'status': 'not_found_in_target_cell_types'}

        # Generate summary
        summary = {
            'timestamp': datetime.now().isoformat(),
            'parameters': {
                'genes': genes,
                'tissue': tissue,
                'disease': disease,
                'species': args.species,
                'development_stage': development_stage or ['adult'],
                'target_cell_type': args.target_cell_type,
                'match_mode': args.match_mode,
                'matched_cell_types': target_cell_types,
                'grouped': should_group,
                'grouped_name': grouped_name if should_group else None,
                'group_by': args.group_by,
            },
            'results': {
                'n_cells': int(adata.n_obs),
                'n_genes_queried': len(genes),
                'n_genes_found': int(adata.n_vars),
                'n_groups': int(adata.obs[args.group_by].nunique()),
                'n_target_groups': len(target_cell_types) if target_cell_types else int(adata.obs[args.group_by].nunique()),
                'n_result_entries': len(results_df),
                'unique_genes_in_results': int(results_df['names'].nunique()) if len(results_df) > 0 else 0,
                'unique_groups_in_results': int(results_df['group'].nunique()) if len(results_df) > 0 else 0,
            },
            'gene_summary': gene_summary,
            'cell_type_distribution': adata.obs[args.group_by].value_counts().head(30).to_dict(),
        }

        # Save results
        save_results(
            results_df=results_df,
            summary=summary,
            output_prefix=args.output,
            output_dir=Path(args.output_dir)
        )

        # Print summary
        logger.info("=" * 80)
        logger.info("ANALYSIS COMPLETE")
        logger.info("=" * 80)
        logger.info(f"Cells analyzed: {summary['results']['n_cells']:,}")
        logger.info(f"Genes found: {summary['results']['n_genes_found']}/{summary['results']['n_genes_queried']}")
        logger.info(f"Cell type groups: {summary['results']['n_groups']}")
        if target_cell_types:
            logger.info(f"Target cell types: {summary['results']['n_target_groups']}")
        logger.info(f"Result entries: {summary['results']['n_result_entries']}")

        # Print per-gene summary
        if gene_summary:
            logger.info("\n" + "=" * 80)
            logger.info("PER-GENE SUMMARY (in target cell types)")
            logger.info("=" * 80)
            for gene, gdata in gene_summary.items():
                if 'status' in gdata:
                    logger.info(f"\n{gene}: {gdata['status']}")
                elif 'grouped' in gdata or 'individual' in gdata:
                    # Grouped + individual format
                    logger.info(f"\n{gene}:")
                    if 'grouped' in gdata:
                        g = gdata['grouped']
                        logger.info(f"  [GROUPED] {g['cell_type']}:")
                        logger.info(f"    Log2FC: {g['log2fc']:.2f}, % expr: {g['pct_expressing']*100:.1f}%, "
                                   f"p-adj: {g['pval_adj']:.2e}, lfc_logp: {g['lfc_logp']:.2f}")
                    if 'individual' in gdata:
                        ind = gdata['individual']
                        logger.info(f"  [INDIVIDUAL BREAKDOWN] Best: {ind['best_cell_type']}")
                        logger.info(f"    Log2FC: {ind['best_log2fc']:.2f}, % expr: {ind['best_pct_expressing']*100:.1f}%, "
                                   f"lfc_logp: {ind['best_lfc_logp']:.2f}")
                        logger.info(f"    Detected in {ind['n_cell_types_detected']} cell type(s)")
                else:
                    # Standard format (no grouping)
                    logger.info(f"\n{gene}:")
                    logger.info(f"  Best cell type: {gdata['best_cell_type']}")
                    logger.info(f"  Log2FC: {gdata['best_log2fc']:.2f}")
                    logger.info(f"  Adj. p-value: {gdata['best_pval_adj']:.2e}")
                    logger.info(f"  % expressing: {gdata['best_pct_expressing']*100:.1f}%")
                    logger.info(f"  lfc_logp: {gdata['best_lfc_logp']:.2f}")
                    logger.info(f"  Detected in {gdata['n_cell_types_detected']} cell type(s)")

        # Print top results table
        if len(results_df) > 0:
            logger.info("\n" + "=" * 80)
            logger.info("ALL RESULTS (sorted by lfc_logp)")
            logger.info("=" * 80)

            # Include analysis_type if present
            if 'analysis_type' in results_df.columns:
                display_cols = ['names', 'group', 'analysis_type', 'logfoldchanges', 'pvals_adj', 'pct_nz_group', 'lfc_logp']
            else:
                display_cols = ['names', 'group', 'logfoldchanges', 'pvals_adj', 'pct_nz_group', 'lfc_logp']
            display_cols = [c for c in display_cols if c in results_df.columns]
            top_results = results_df.head(20)[display_cols].copy()

            # Format for display
            top_results['logfoldchanges'] = top_results['logfoldchanges'].apply(lambda x: f"{x:.2f}")
            top_results['pvals_adj'] = top_results['pvals_adj'].apply(lambda x: f"{x:.2e}")
            top_results['pct_nz_group'] = top_results['pct_nz_group'].apply(lambda x: f"{x*100:.1f}%")
            top_results['lfc_logp'] = top_results['lfc_logp'].apply(lambda x: f"{x:.2f}")

            if 'analysis_type' in top_results.columns:
                top_results.columns = ['Gene', 'Cell Type', 'Type', 'Log2FC', 'P-adj', '% Expr', 'lfc_logp']
            else:
                top_results.columns = ['Gene', 'Cell Type', 'Log2FC', 'P-adj', '% Expr', 'lfc_logp']
            print(top_results.to_string(index=False))

        logger.info("\n" + "=" * 80)
        logger.info(f"Output files:")
        logger.info(f"  {args.output_dir}/{args.output}_specificity.csv")
        logger.info(f"  {args.output_dir}/{args.output}_summary.json")
        logger.info(f"  {args.output_dir}/{args.output}_summary.md")
        logger.info("=" * 80)

    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
