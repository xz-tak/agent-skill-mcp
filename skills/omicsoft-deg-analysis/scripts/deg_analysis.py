#!/usr/bin/env python3
"""
Omicsoft Bulk DEG Analysis Script
Customizable differential gene expression analysis for h5ad files

Usage:
    python deg_analysis.py --file <path> --target <gene> --diseases <disease1,disease2> [options]
"""

import os
import sys
import argparse
from datetime import datetime
import tempfile
import subprocess

import anndata
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

import gseapy as gp
from tqdm import tqdm

import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
import plotly.figure_factory as ff
from plotly.colors import qualitative

import random

# Set random seed for reproducibility
random.seed(42)


def is_s3_uri(path):
    """Check if path is an S3 URI"""
    return isinstance(path, str) and path.startswith('s3://')


def validate_s3_access(s3_uri):
    """
    Validate access to S3 object

    Parameters
    ----------
    s3_uri : str
        S3 URI in format s3://bucket/key

    Returns
    -------
    bool
        True if accessible, False otherwise
    """
    try:
        import boto3
        from botocore.exceptions import ClientError, NoCredentialsError

        # Parse S3 URI
        if not s3_uri.startswith('s3://'):
            return False

        parts = s3_uri[5:].split('/', 1)
        if len(parts) != 2:
            print(f"Invalid S3 URI format: {s3_uri}")
            return False

        bucket, key = parts

        # Check access
        s3_client = boto3.client('s3')
        s3_client.head_object(Bucket=bucket, Key=key)
        print(f"✓ S3 access validated: {s3_uri}")
        return True

    except NoCredentialsError:
        print(f"✗ No AWS credentials found. Cannot access S3 URI: {s3_uri}")
        return False
    except ClientError as e:
        error_code = e.response['Error']['Code']
        if error_code == '404':
            print(f"✗ S3 object not found: {s3_uri}")
        elif error_code == '403':
            print(f"✗ Access denied to S3 object: {s3_uri}")
        else:
            print(f"✗ Error accessing S3 object: {e}")
        return False
    except Exception as e:
        print(f"✗ Unexpected error validating S3 access: {e}")
        return False


def load_h5ad(file_path):
    """
    Load h5ad file from local path or S3 URI and clean up expression data from uns

    Parameters
    ----------
    file_path : str
        Local file path or S3 URI (s3://bucket/key)

    Returns
    -------
    anndata.AnnData
        Loaded AnnData object with expression data removed from uns

    Notes
    -----
    If the file is named with pattern *_deg.h5ad and contains a corresponding
    *_expr key in uns, that key will be deleted to reduce memory usage and
    prevent confusion during DEG analysis.
    """
    if is_s3_uri(file_path):
        print(f"Detected S3 URI: {file_path}")

        # Validate S3 access
        if not validate_s3_access(file_path):
            raise ValueError(f"Cannot access S3 object: {file_path}")

        # Load directly from S3 using s3fs
        print("Loading h5ad file directly from S3...")
        try:
            import s3fs
            fs = s3fs.S3FileSystem()
            with fs.open(file_path, 'rb') as f:
                adata = sc.read_h5ad(f)
            print("✓ Successfully loaded from S3")
        except ImportError:
            raise ImportError(
                "s3fs package is required to read from S3. "
                "Install it with: pip install s3fs"
            )
    else:
        print(f"Detected local path: {file_path}")

        # Validate local file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Local file not found: {file_path}")

        print("Loading h5ad file from local path...")
        adata = sc.read_h5ad(file_path)
        print("✓ Successfully loaded from local path")

    # Clean up expression data from uns if present
    # Extract base name from file path (works for both local and S3)
    if is_s3_uri(file_path):
        # For S3 URIs, get the last part after the final /
        base_name = file_path.split('/')[-1]
    else:
        # For local paths, get the filename
        base_name = os.path.basename(file_path)

    # Remove .h5ad extension
    base_name = base_name.replace('.h5ad', '')

    # Check if file ends with _deg and construct potential expr key
    if base_name.endswith('_deg'):
        potential_expr_key = base_name.replace('_deg', '_expr')

        # Check if this key exists in uns
        if potential_expr_key in adata.uns:
            print(f"\nℹ Found expression data in uns['{potential_expr_key}']")
            print(f"  Removing to reduce memory usage and focus on DEG analysis...")
            del adata.uns[potential_expr_key]
            print(f"✓ Removed uns['{potential_expr_key}']")
        else:
            print(f"\nℹ No expression data found in uns (checked for '{potential_expr_key}')")

    return adata


def ad2long_df(adata, var_name='Gene', padj_col='padj'):
    """Convert AnnData to long-format DataFrame"""
    df_list = []
    adata.obs.index.name = 'index'
    id_vars = [adata.obs.index.name]

    for k, v in adata.layers.items():
        df = pd.DataFrame(
            index=adata.obs_names,
            columns=adata.var_names,
            data=v.toarray() if isinstance(v, sparse._csr.csr_matrix) else v,
        ).reset_index()

        df_long = pd.melt(
            df,
            id_vars=id_vars,
            value_vars=[col for col in df.columns if col not in id_vars],
            var_name=var_name,
            value_name=k
        ).set_index(id_vars[0])

        df_list.append(df_long)

    concat_df = pd.concat(df_list, axis=1)
    concat_df = concat_df.loc[:, ~concat_df.columns.duplicated()]

    # Replace zero with smallest nonzero value
    concat_df[padj_col] = concat_df[padj_col].replace(
        0,
        min(concat_df[padj_col][concat_df[padj_col] != 0]) / 10
    )

    concat_df['-logpadj'] = -np.log(concat_df[padj_col])
    concat_df['trans_logpadj'] = concat_df['-logpadj'].apply(lambda x: np.log2(x + 1))

    meta_df = pd.merge(concat_df, adata.obs, left_index=True, right_index=True, how='outer')

    return meta_df


def lfc_heatmap_traces(
    adata,
    geneset,
    fig_title,
    lfc="log2fc",
    padj="padj",
    lfc_limit=0,
    padj_limit=0.05,
):
    """
    Create a clustered annotated heatmap for the given anndata object and geneset,
    and return its traces, its layout annotations, and the title.
    No plotting or saving is done here; this is meant to be used in a combined
    figure (e.g., with a dropdown).
    """

    # Use the lfc layer as X and subset to genes in geneset present in adata
    adata = adata.copy()  # avoid mutating original
    adata.X = adata.layers[lfc].copy()
    adata = adata[:, [gene for gene in sorted(geneset) if gene in adata.var_names]]

    # Subset samples with any lfc beyond lfc_limit
    X_dense = adata.X.toarray()
    include_comp = adata.obs_names[(np.abs(X_dense) > lfc_limit).any(axis=1)]
    include_gene = adata.var_names[(np.abs(X_dense) > lfc_limit).any(axis=0)]
    adata = adata[include_comp, include_gene]

    # If nothing survives filtering, just return empty
    if adata.n_obs == 0 or adata.n_vars == 0:
        return [], [], fig_title

    lfc_df = adata.to_df()
    padj_matrix = adata.layers[padj].toarray()
    annotations = np.where(padj_matrix < padj_limit, "*", "")

    # Custom y-axis labels: {study}_{disease}_{comparison}
    y_labels = []
    for idx in lfc_df.index:
        study = adata.obs.loc[idx, "study"]
        disease = adata.obs.loc[idx, "disease_category"]
        comparison = adata.obs.loc[idx, "comparison"]
        y_labels.append(f"{study}_{disease}_{comparison}")

    n_rows, n_cols = lfc_df.shape

    # --- Hierarchical clustering: rows (comparisons) ---
    if n_rows > 1:
        row_linkage = linkage(
            pdist(lfc_df.values, metric="euclidean"),
            method="average",
        )
        row_order = leaves_list(row_linkage)
    else:
        row_order = np.arange(n_rows)

    # --- Hierarchical clustering: columns (genes) ---
    if n_cols > 1:
        col_linkage = linkage(
            pdist(lfc_df.values.T, metric="euclidean"),
            method="average",
        )
        col_order = leaves_list(col_linkage)
    else:
        col_order = np.arange(n_cols)

    # Reorder matrix, annotations, and labels according to clustering
    lfc_df = lfc_df.iloc[row_order, :]
    lfc_df = lfc_df.iloc[:, col_order]

    annotations = annotations[row_order, :][:, col_order]
    y_labels = [y_labels[i] for i in row_order]
    x_labels = list(lfc_df.columns)

    num_x_elements = len(x_labels)
    num_y_elements = len(y_labels)
    plot_width = 50 * num_x_elements
    plot_height = 200 + 30 * num_y_elements

    fig_heatmap = ff.create_annotated_heatmap(
        z=lfc_df.values,
        x=x_labels,
        y=y_labels,
        annotation_text=annotations,
        colorscale="RdBu_r",
        zmid=0,
        hoverinfo="z",
        showscale=True,
    )

    fig_heatmap.update_layout(
        title=fig_title,
        xaxis_title="Genes",
        yaxis_title="Samples",
        xaxis=dict(tickmode="linear"),
        yaxis=dict(tickmode="linear"),
        width=plot_width,
        height=plot_height,
    )

    # Make "*" green
    for annot in fig_heatmap.layout.annotations:
        if annot.text == "*":
            annot.font.color = "green"
            annot.font.size = 20

    # Return traces + annotations + title (for combined dropdown figure)
    return fig_heatmap.data, fig_heatmap.layout.annotations, fig_heatmap.layout.title.text


def create_enhanced_gsea_summary(res_deg_score, signature_set, target_genes=None, fdr_threshold=0.05, output_dir='.', output_filename='gsea_enhanced_summary.csv'):
    """
    Create enhanced GSEA summary table with significance counts, NES grouping, and leading edge annotation

    Parameters
    ----------
    res_deg_score : pd.DataFrame
        GSEA results dataframe with columns: Term, NES, FDR q-val, Lead_genes, study
    signature_set : dict
        Dictionary of signature names to gene lists
    target_genes : list, optional
        List of target genes to track in leading edge (if None, uses all signature genes)
    fdr_threshold : float
        FDR q-value threshold for significance (default: 0.05)
    output_dir : str
        Output directory path
    output_filename : str
        Name of output CSV file

    Returns
    -------
    pd.DataFrame
        Enhanced summary table
    """
    print(f"\nCreating enhanced GSEA summary (FDR < {fdr_threshold})...")

    # Count total comparisons (unique studies in results)
    total_comparisons = res_deg_score['study'].nunique()
    print(f"Total comparisons analyzed: {total_comparisons}")

    # Get genes to track in leading edge
    if target_genes is not None:
        # Use only target genes for leading edge tracking
        genes_to_track = set(target_genes)
        print(f"Tracking {len(genes_to_track)} target genes in leading edge")
    else:
        # Use all signature genes for leading edge tracking
        genes_to_track = set()
        for genes in signature_set.values():
            genes_to_track.update(genes)
        print(f"Tracking {len(genes_to_track)} signature genes in leading edge")

    # Process each term and comparison_category combination
    summary_rows = []

    # Check if comparison_category is available
    has_comparison_category = 'comparison_category' in res_deg_score.columns

    if has_comparison_category:
        # Group by both Term and comparison_category
        groupby_keys = ['Term', 'comparison_category']
    else:
        # Group by Term only
        groupby_keys = ['Term']

    for group_key, group_df in res_deg_score.groupby(groupby_keys):
        # Extract term and category
        if has_comparison_category:
            term, comparison_category = group_key
        else:
            term = group_key
            comparison_category = 'All'

        # Filter for significant results
        sig_df = group_df[group_df['FDR q-val'] < fdr_threshold]

        if len(sig_df) == 0:
            continue

        # Count significant comparisons for this term-category combination
        n_significant = len(sig_df)

        # Separate by NES direction
        pos_nes = sig_df[sig_df['NES'] > 0]
        neg_nes = sig_df[sig_df['NES'] < 0]

        # Get studies and NES values for each direction
        pos_studies = []
        neg_studies = []

        if len(pos_nes) > 0:
            for _, row in pos_nes.iterrows():
                pos_studies.append(f"{row['study']} (NES={row['NES']:.2f})")

        if len(neg_nes) > 0:
            for _, row in neg_nes.iterrows():
                neg_studies.append(f"{row['study']} (NES={row['NES']:.2f})")

        # Annotate leading edge genes
        # Check which target genes appear in leading edges
        leading_edge_targets = set()

        for _, row in sig_df.iterrows():
            if 'Lead_genes' in row and pd.notna(row['Lead_genes']):
                # Parse leading edge genes (format may vary)
                lead_genes_str = str(row['Lead_genes'])
                lead_genes = [g.strip() for g in lead_genes_str.split(';') if g.strip()]

                # Check intersection with target genes
                for gene in lead_genes:
                    if gene in genes_to_track:
                        leading_edge_targets.add(gene)

        # Create summary row
        summary_row = {
            'Term': term,
            'Comparison_Category': comparison_category,
            'N_Significant_Comparisons': n_significant,
            'Total_Comparisons': total_comparisons,
            'Percent_Significant': f"{(n_significant/total_comparisons*100):.1f}%",
            'N_Positive_NES': len(pos_nes),
            'N_Negative_NES': len(neg_nes),
            'Positive_NES_Studies': '; '.join(pos_studies) if pos_studies else 'None',
            'Negative_NES_Studies': '; '.join(neg_studies) if neg_studies else 'None',
            'Input_Targets_In_Leading_Edge': ', '.join(sorted(leading_edge_targets)) if leading_edge_targets else 'None',
            'N_Input_Targets_In_LE': len(leading_edge_targets)
        }

        summary_rows.append(summary_row)

    # Create summary dataframe
    summary_df = pd.DataFrame(summary_rows)

    if len(summary_df) > 0:
        # Sort by term and comparison category
        if has_comparison_category:
            summary_df = summary_df.sort_values(['N_Significant_Comparisons', 'Term', 'Comparison_Category'], ascending=[False, True, True])
        else:
            summary_df = summary_df.sort_values('N_Significant_Comparisons', ascending=False)

        # Save to file
        output_path = f"{output_dir}/{output_filename}"
        summary_df.to_csv(output_path, index=False)
        print(f"Saved enhanced GSEA summary to {output_path}")

        # Print summary statistics
        print(f"\nGSEA Summary Statistics:")
        if has_comparison_category:
            print(f"  Unique terms with FDR < {fdr_threshold}: {summary_df['Term'].nunique()}")
            print(f"  Term-Category combinations: {len(summary_df)}")
        else:
            print(f"  Terms with FDR < {fdr_threshold}: {len(summary_df)}")
        print(f"  Total comparisons: {total_comparisons}")
        print(f"\nTop 10 enriched term-category combinations by number of significant comparisons:")

        # Select columns to display
        display_cols = ['Term', 'Comparison_Category', 'N_Significant_Comparisons', 'N_Positive_NES', 'N_Negative_NES', 'N_Input_Targets_In_LE']

        print(summary_df[display_cols].head(10).to_string(index=False))
    else:
        print(f"No terms with FDR < {fdr_threshold} found")

    return summary_df


def run_analysis(args):
    """Main analysis function"""

    dt = datetime.now().strftime('%Y%m%d')
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 80)
    print(f"OMICSOFT DEG ANALYSIS - {args.target_name}")
    print("=" * 80)
    print(f"\nParameters:")
    print(f"  Target: {args.target_name}")
    print(f"  File: {args.file}")
    print(f"  Diseases: {args.diseases if args.diseases else 'None (all diseases)'}")
    print(f"  Studies: {args.studies if args.studies else 'None (all studies)'}")
    print(f"  Tissues: {args.tissues if args.tissues else 'None (all tissues)'}")
    print(f"  Comparison Category: {args.comparison_category if args.comparison_category else 'None (all categories)'}")
    print(f"  Case Treatment: {args.case_treatment if args.case_treatment else 'None (all treatments)'}")
    print(f"  Comparison: {args.comparison if args.comparison else 'None (all comparisons)'}")
    print(f"  Log2FC threshold: {args.lfc_threshold}")
    print(f"  Adjusted p-value threshold: {args.padj_threshold}")
    print(f"  Output directory: {output_dir}")

    # Load data
    print(f"\n[1/5] Loading h5ad file...")
    adata = load_h5ad(args.file)
    print(f"Loaded AnnData object: {adata.shape} (obs x vars)")
    print(f"  Observations: {adata.n_obs:,}")
    print(f"  Variables: {adata.n_vars:,}")

    # Filter data
    print(f"\n[2/5] Filtering data...")

    # Start with all observations (no filter)
    base_filter = pd.Series([True] * adata.n_obs, index=adata.obs.index)

    # Apply disease filter if provided
    if args.diseases:
        print(f"\nTotal diseases in dataset: {len(adata.obs.disease.unique())}")
        disease_pattern = '|'.join([d.strip() for d in args.diseases.split(',')])
        print(f"Filtering for diseases matching: {disease_pattern}")

        matching_diseases = [d for d in adata.obs.disease.unique()
                            if any(term.lower() in d.lower() for term in args.diseases.split(','))]

        print(f"Found {len(matching_diseases)} matching disease(s):")
        for disease in sorted(matching_diseases):
            count = (adata.obs.disease == disease).sum()
            print(f"  - {disease}: {count} observations")

        disease_filter = adata.obs.disease.str.contains(disease_pattern, case=False, na=False)
        base_filter = base_filter & disease_filter
        print(f"  Observations after disease filter: {base_filter.sum():,}")

    # Apply disease exclusion filter if provided
    if args.exclude_diseases:
        exclude_pattern = '|'.join([d.strip() for d in args.exclude_diseases.split(',')])
        print(f"\nExcluding diseases matching: {exclude_pattern}")

        # Get currently filtered observations
        temp_adata = adata[base_filter]
        excluded_diseases = [d for d in temp_adata.obs.disease.unique()
                            if any(term.lower() in d.lower() for term in args.exclude_diseases.split(','))]

        if excluded_diseases:
            print(f"Found {len(excluded_diseases)} disease(s) to exclude:")
            for disease in sorted(excluded_diseases):
                count = (temp_adata.obs.disease == disease).sum()
                print(f"  - {disease}: {count} observations")

        exclude_filter = ~adata.obs.disease.str.contains(exclude_pattern, case=False, na=False)
        base_filter = base_filter & exclude_filter
        print(f"  Observations after disease exclusion: {base_filter.sum():,}")

    # Apply study filter if provided
    if args.studies:
        study_list = [s.strip() for s in args.studies.split(',')]
        print(f"\nFiltering for studies: {study_list}")
        study_filter = adata.obs.study.isin(study_list)
        base_filter = base_filter & study_filter
        print(f"  Observations after study filter: {base_filter.sum():,}")

    # Apply tissue filter if provided
    if args.tissues:
        print(f"\nTotal tissues in dataset: {len(adata.obs.tissue.unique())}")
        tissue_pattern = '|'.join([t.strip() for t in args.tissues.split(',')])
        print(f"Filtering for tissues matching: {tissue_pattern}")

        matching_tissues = [t for t in adata.obs.tissue.unique()
                           if any(term.lower() in t.lower() for term in args.tissues.split(','))]

        print(f"Found {len(matching_tissues)} matching tissue(s):")
        for tissue in sorted(matching_tissues):
            count = (adata.obs.tissue == tissue).sum()
            print(f"  - {tissue}: {count} observations")

        tissue_filter = adata.obs.tissue.str.contains(tissue_pattern, case=False, na=False)
        base_filter = base_filter & tissue_filter
        print(f"  Observations after tissue filter: {base_filter.sum():,}")

    # Apply comparison category filter if provided
    if args.comparison_category:
        category_list = [c.strip() for c in args.comparison_category.split(',')]
        print(f"\nFiltering for comparison categories: {category_list}")
        category_filter = adata.obs.comparison_category.isin(category_list)
        base_filter = base_filter & category_filter
        print(f"  Observations after comparison category filter: {base_filter.sum():,}")

    # Apply case treatment filter if provided
    if args.case_treatment:
        treatment_list = [t.strip() for t in args.case_treatment.split(',')]
        print(f"\nFiltering for case treatments: {treatment_list}")
        treatment_filter = adata.obs.case_treatment.isin(treatment_list)
        base_filter = base_filter & treatment_filter
        print(f"  Observations after case treatment filter: {base_filter.sum():,}")

    # Apply comparison filter if provided
    if args.comparison:
        print(f"\nTotal comparisons in dataset: {len(adata.obs.comparison.unique())}")
        comparison_pattern = '|'.join([c.strip() for c in args.comparison.split(',')])
        print(f"Filtering for comparisons matching: {comparison_pattern}")

        matching_comparisons = [c for c in adata.obs.comparison.unique()
                               if any(term.lower() in c.lower() for term in args.comparison.split(','))]

        print(f"Found {len(matching_comparisons)} matching comparison(s):")
        for comparison in sorted(matching_comparisons):
            count = (adata.obs.comparison == comparison).sum()
            print(f"  - {comparison}: {count} observations")

        comparison_filter = adata.obs.comparison.str.contains(comparison_pattern, case=False, na=False)
        base_filter = base_filter & comparison_filter
        print(f"  Observations after comparison filter: {base_filter.sum():,}")

    # Remove duplicates
    base_filter = base_filter & ~adata.obs.index.duplicated()

    filtered_deg_adata = adata[base_filter]

    print(f"\nFiltered dataset: {filtered_deg_adata.shape}")
    print(f"  Observations: {filtered_deg_adata.n_obs:,}")
    print(f"  Variables: {filtered_deg_adata.n_vars:,}")

    if filtered_deg_adata.n_obs == 0:
        print("\nNo observations remain after filtering! Exiting.")
        return

    # Process disease categories
    print(f"\n[3/5] Processing disease categories...")
    filtered_deg_adata.obs['study'] = filtered_deg_adata.obs['study'].apply(lambda x: x.replace('/', '-'))
    filtered_deg_adata.obs['disease_category'] = filtered_deg_adata.obs.disease

    print(f"\nDisease categories: {len(filtered_deg_adata.obs.disease_category.unique())}")

    # Create study summary
    study_df = filtered_deg_adata.obs.reset_index(drop=True).loc[:, ['study', 'tissue', 'disease', 'disease_category']].drop_duplicates()
    print(f"Number of unique studies: {len(study_df)}")
    study_df.to_csv(f'{output_dir}/study_summary.csv', index=False)
    print(f"Saved study summary to {output_dir}/study_summary.csv")

    # Convert to long format
    print(f"\n[4/5] Converting to long format and creating visualizations...")
    long_df_comp = ad2long_df(filtered_deg_adata)
    print(f"Long format DataFrame shape: {long_df_comp.shape}")

    # Parse signature sets from user input
    # Format: "Sig1:Gene1,Gene2,Gene3;Sig2:GeneA,GeneB,GeneC"
    signature_set = {}
    if args.signatures:
        for sig_entry in args.signatures.split(';'):
            sig_entry = sig_entry.strip()
            if ':' not in sig_entry:
                print(f"Warning: Skipping invalid signature format: {sig_entry}")
                continue
            sig_name, genes = sig_entry.split(':', 1)
            gene_list = [g.strip() for g in genes.split(',') if g.strip()]
            if gene_list:
                signature_set[sig_name.strip()] = gene_list
        print(f"\nParsed {len(signature_set)} signature(s):")
        for sig_name, genes in signature_set.items():
            print(f"  - {sig_name}: {len(genes)} genes")

    # Parse target genes if provided
    target_genes = None
    if args.targets:
        target_genes = [g.strip() for g in args.targets.split(',') if g.strip()]
        print(f"\nParsed {len(target_genes)} target gene(s) for leading edge tracking:")
        print(f"  Target genes: {', '.join(target_genes)}")
        print(f"  Note: Targets will be plotted but NOT included in GSEA enrichment")

    # Create signature_set_for_plotting that includes targets
    signature_set_for_plotting = signature_set.copy()
    if target_genes:
        signature_set_for_plotting['targets'] = target_genes

    # Keep original signature_set for GSEA (without targets)
    signature_set_for_gsea = signature_set.copy()

    if not signature_set_for_plotting:
        print("Warning: No signatures or targets provided, skipping signature analysis")
    else:
        # Create scatter plots
        print(f"\nCreating scatter plots for {len(signature_set_for_plotting)} item(s)...")
        fig_dict = {}

        for sig, geneset in tqdm(signature_set_for_plotting.items(), desc="Creating scatter plots"):
            filtered_data = long_df_comp[
                (long_df_comp['Gene'].isin(sorted(geneset))) &
                (long_df_comp['padj'] < args.padj_threshold) &
                (abs(long_df_comp['log2fc']) > args.lfc_threshold)
            ]

            if len(filtered_data) == 0:
                print(f"No significant genes found for {sig}")
                continue

            fig_d = px.scatter(
                filtered_data,
                x="Gene", y="log2fc",
                size='-logpadj', color="disease_category", color_discrete_sequence=qualitative.Set3,
                hover_name="study",
                hover_data=['log2fc', 'padj', 'tissue', "disease_category", 'comparison', "case_treatment_status", "case_treatment_status"],
                size_max=25, title=f"{sig} log2fc", width=1200, height=600
            )
            fig_d.update_xaxes(categoryorder='category ascending')
            fig_d.add_hline(y=-args.lfc_threshold, line_dash="dot", line_color="green")
            fig_d.add_hline(y=args.lfc_threshold, line_dash="dot", line_color="green")
            fig_dict[sig] = fig_d

        # Create combined figure with dropdown
        if fig_dict:
            combined_fig = go.Figure()
            combined_fig.update_xaxes(tickfont=dict(size=20))
            combined_fig.update_yaxes(tickfont=dict(size=20))

            for i, (sig, fig) in enumerate(fig_dict.items()):
                for trace in fig.data:
                    combined_fig.add_trace(trace)
                    combined_fig.data[-1].visible = (i == 0)
                combined_fig.add_hline(y=-args.lfc_threshold, line_dash="dot", line_color="green", visible=(i == 0))
                combined_fig.add_hline(y=args.lfc_threshold, line_dash="dot", line_color="green", visible=(i == 0))

            # Create dropdown buttons
            dropdown_buttons = []
            for i, sig in enumerate(fig_dict.keys()):
                visibility = [False] * len(combined_fig.data)
                start_idx = sum(len(fig_dict[key].data) for key in list(fig_dict.keys())[:i])
                end_idx = start_idx + len(fig_dict[sig].data)
                for j in range(start_idx, end_idx):
                    visibility[j] = True

                dropdown_buttons.append(dict(
                    label=sig,
                    method="update",
                    args=[{"visible": visibility}, {"title": f"{sig} log2fc"}]
                ))

            combined_fig.update_layout(
                updatemenus=[dict(
                    active=0, buttons=dropdown_buttons, direction="down",
                    showactive=True, x=1, y=1
                )],
                legend=dict(itemsizing='constant', itemwidth=100),
                title="Select a signature from Dropdown"
            )

            combined_fig.write_html(f"{output_dir}/target_signature_{dt}.html")
            print(f"Saved interactive scatter plot to {output_dir}/target_signature_{dt}.html")

        # Create heatmaps for each signature
        print(f"\nCreating heatmaps...")
        heatmap_dict = {}

        # For each signature, generate heatmap traces
        for sig, geneset in tqdm(signature_set_for_plotting.items(), desc="Creating heatmaps"):
            title = f"{sig} log2fc"
            try:
                traces, annots, stored_title = lfc_heatmap_traces(
                    filtered_deg_adata, geneset, title,
                    lfc_limit=args.lfc_threshold,
                    padj_limit=args.padj_threshold
                )
                heatmap_dict[sig] = (traces, annots, stored_title)
            except Exception as e:
                print(f"Skipping {sig} due to error: {e}")
                continue

        # Create a combined figure
        if heatmap_dict:
            combined_heatmap = go.Figure()
            combined_fig.update_xaxes(tickfont=dict(size=20))
            combined_fig.update_yaxes(tickfont=dict(size=20))

            # Record the number of traces per signature for index calculation
            trace_counts = [len(heatmap_dict[key][0]) for key in heatmap_dict.keys()]
            total_traces = sum(trace_counts)

            # Add all traces to the combined figure with proper initial visibility
            for i, key in enumerate(heatmap_dict.keys()):
                for trace in heatmap_dict[key][0]:
                    combined_heatmap.add_trace(trace)
                    # Only show first signature initially
                    combined_heatmap.data[-1].visible = (i == 0)

            # Set initial layout annotations from the first signature
            first_key = list(heatmap_dict.keys())[0]
            combined_heatmap.update_layout(
                title=heatmap_dict[first_key][2],
                annotations=heatmap_dict[first_key][1],
                xaxis_title="Genes",
                yaxis_title="Samples"
            )

            # Create dropdown buttons to update visibility and layout annotations
            dropdown_buttons = []
            keys = list(heatmap_dict.keys())
            for i, sig in enumerate(keys):
                visibility = [False] * total_traces
                start_idx = sum(trace_counts[:i])
                end_idx = start_idx + trace_counts[i]
                for j in range(start_idx, end_idx):
                    visibility[j] = True
                dropdown_buttons.append(
                    dict(
                        label=sig,
                        method="update",
                        args=[
                            {"visible": visibility},
                            {"title": heatmap_dict[sig][2],
                             "annotations": heatmap_dict[sig][1]}
                        ]
                    )
                )

            combined_heatmap.update_layout(
                updatemenus=[dict(
                    active=0,
                    buttons=dropdown_buttons,
                    direction="down",
                    showactive=True,
                    x=1, y=1.05
                )],
                legend=dict(itemsizing='constant', itemwidth=100),
                title="Select a signature from Dropdown"
            )

            combined_heatmap.write_html(f"{output_dir}/comparison_signature_{dt}.html")
            print(f"Saved interactive heatmap to {output_dir}/comparison_signature_{dt}.html")
        else:
            print("No heatmaps generated")

        # Create summary table for all genes (signatures + targets)
        print(f"\nCreating {args.target_name} signature summary...")
        all_genes = list(set(g for dic in signature_set_for_plotting.values() for g in dic))

        sub_df = (
            long_df_comp[long_df_comp['Gene'].isin(all_genes)]
            .assign(sig=lambda x: np.where(
                (x['log2fc'] > args.lfc_threshold) & (x['padj'] < args.padj_threshold), 'up',
                np.where((x['log2fc'] < -args.lfc_threshold) & (x['padj'] < args.padj_threshold), 'dn', 'na')
            ))
            .groupby(['Gene', 'tissue', 'disease_category', 'comparison_category', 'sig'])
            .agg(count=('Gene', 'size'), study=('study', list))
            .reset_index()
            .assign(total_count=lambda df: df.groupby(['Gene', 'tissue', 'disease_category', 'comparison_category'])['count'].transform('sum'))
            .assign(threshold=f'|lfc|>{args.lfc_threshold};padj<{args.padj_threshold}')
            .query('count != 0 and sig != "na"')
        )

        if len(sub_df) > 0:
            sub_df.to_csv(f'{output_dir}/signature_summary.csv', index=False)
            print(f"Saved signature summary to {output_dir}/signature_summary.csv")
            print(f"\nSignature hits across diseases and comparison categories:")
            print(sub_df.groupby(['disease_category', 'comparison_category', 'tissue', 'sig'])['Gene'].count())
        else:
            print("No significant genes found in signature analysis")

        # Export detailed long table for signature/target genes
        print(f"\nExporting detailed long table for signature/target genes...")
        detailed_long_df = long_df_comp[long_df_comp['Gene'].isin(all_genes)].copy()
        detailed_long_df.to_csv(f'{output_dir}/detailed_gene_table.csv', index=False)
        print(f"Saved detailed gene table to {output_dir}/detailed_gene_table.csv")
        print(f"  Shape: {detailed_long_df.shape} (rows x columns)")
        print(f"  Genes included: {detailed_long_df['Gene'].nunique()}")

    # GSEA Analysis
    if args.run_gsea:
        print(f"\n[5/5] Running GSEA analysis...")
        long_df_comp['padj'] = long_df_comp['padj'].replace(
            0,
            long_df_comp.loc[long_df_comp['padj'] > 0, 'padj'].min() / 10
        )
        enrich_df = long_df_comp.reset_index(drop=True).drop_duplicates().set_index('Gene')
        enrich_df['sig_score'] = enrich_df['log2fc'] * enrich_df['-logpadj']

        res_deg_list = []

        # MSigDB_Hallmark_2020 is always included by default
        gsea_gene_sets = ['MSigDB_Hallmark_2020']

        print(f"GSEA will run with MSigDB_Hallmark_2020 (default) and {len(signature_set_for_gsea)} custom signature(s)")

        for c, df in tqdm(enrich_df.groupby(['study','comparison','comparison_category']), desc="Running GSEA"):
            study_comp = ':'.join(c[:2])  # study:comparison
            comp_category = c[2]  # comparison_category
            try:
                # Run GSEA on custom signatures if provided (excluding targets)
                if signature_set_for_gsea:
                    res_deg_enr = gp.prerank(
                        rnk=df['sig_score'][df['sig_score'] != 0].sort_values(ascending=False),
                        gene_sets=signature_set_for_gsea,
                        permutation_num=1000,
                        min_size=1,
                        max_size=8000,
                        threads=30,
                        seed=42,
                        outdir=None
                    )
                    res_deg_enr_df = res_deg_enr.res2d.assign(study=study_comp, comparison_category=comp_category)
                    res_deg_list.append(res_deg_enr_df)

                # Always run GSEA on MSigDB_Hallmark_2020
                res_deg = gp.prerank(
                    rnk=df['sig_score'][df['sig_score'] != 0].sort_values(ascending=False),
                    gene_sets=gsea_gene_sets,
                    permutation_num=1000,
                    min_size=1,
                    max_size=8000,
                    threads=30,
                    seed=42,
                    outdir=None
                )
                res_deg_df = res_deg.res2d.assign(study=study_comp, comparison_category=comp_category)
                res_deg_list.append(res_deg_df)
            except:
                print(f'skip {study_comp}')

        if res_deg_list:
            res_deg_score = pd.concat(res_deg_list)
            res_deg_score.to_csv(f'{output_dir}/gsea_all_results.csv', index=False)
            print(f"Saved GSEA results to {output_dir}/gsea_all_results.csv")

            # Create enhanced GSEA summary with leading edge annotation
            enhanced_summary_df = create_enhanced_gsea_summary(
                res_deg_score,
                signature_set_for_gsea,
                target_genes=target_genes,  # Pass target genes for leading edge tracking
                fdr_threshold=args.padj_threshold,
                output_dir=output_dir,
                output_filename='gsea_enhanced_summary.csv'
            )

            # Create GSEA scatter plot
            filtered_gsea = res_deg_score[(res_deg_score['NOM p-val'] < 0.05)]
            filtered_gsea['-logpadj'] = filtered_gsea['FDR q-val'].apply(lambda x: -np.log(x) if x != 0 else 0)

            fig_gsea = px.scatter(
                filtered_gsea,
                x="Term", y="NES",
                size='-logpadj', color="study", color_discrete_sequence=qualitative.Set3,
                hover_name="study",
                hover_data=['Term', 'NES', 'NOM p-val', "FDR q-val", "FWER p-val"],
                size_max=25, title=f"Pathway Enrichment Score", width=1600, height=800
            )
            fig_gsea.update_xaxes(categoryorder='category ascending')
            fig_gsea.add_hline(y=-1, line_dash="dot", line_color="green")
            fig_gsea.add_hline(y=1, line_dash="dot", line_color="green")
            fig_gsea.write_html(f"{output_dir}/gsea_score_pathway_{dt}.html")
            print(f"Saved GSEA scatter plot to {output_dir}/gsea_score_pathway_{dt}.html")

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE!")
    print("=" * 80)
    print(f"\nOutput files saved in {output_dir}/")


def main():
    parser = argparse.ArgumentParser(description='Omicsoft Bulk DEG Analysis')
    parser.add_argument('--file', required=True, help='Path to h5ad file')
    parser.add_argument('--target-name', required=True, help='Name of target gene/pathway for analysis')
    parser.add_argument('--signatures', required=True, help='Signatures in format "Name:Gene1,Gene2;Name2:GeneA,GeneB" (see references/signature_format.md)')

    # Optional target genes parameter
    parser.add_argument('--targets', default=None, help='Optional: Comma-separated list of target genes for tracking in leading edge and plotting (e.g., "GENE1,GENE2,GENE3"). These will be plotted but NOT included in GSEA enrichment.')

    # Optional filtering parameters
    parser.add_argument('--diseases', default=None, help='Optional: Comma-separated list of disease keywords to filter (e.g., "scleroderma,sclerosis")')
    parser.add_argument('--exclude-diseases', default=None, help='Optional: Comma-separated list of disease keywords to exclude (e.g., "ALS,amyotrophic lateral sclerosis")')
    parser.add_argument('--studies', default=None, help='Optional: Comma-separated list of study names to filter (e.g., "GSE12345,GSE67890")')
    parser.add_argument('--tissues', default=None, help='Optional: Comma-separated list of tissue keywords to filter (e.g., "skin,blood,lung")')
    parser.add_argument('--comparison-category', default=None, help='Optional: Comma-separated list of comparison categories to filter (e.g., "Disease vs. Normal")')
    parser.add_argument('--case-treatment', default=None, help='Optional: Comma-separated list of treatment keywords to filter (e.g., "none,NA")')
    parser.add_argument('--comparison', default=None, help='Optional: Comma-separated list of comparison keywords for fuzzy/partial matching (e.g., "response vs no response")')

    parser.add_argument('--lfc-threshold', type=float, default=0.0, help='Absolute log2 fold change threshold (default: 0.0)')
    parser.add_argument('--padj-threshold', type=float, default=0.05, help='Adjusted p-value threshold (default: 0.05)')
    parser.add_argument('--output-dir', default='deg_results', help='Output directory (default: deg_results)')
    parser.add_argument('--run-gsea', action='store_true', help='Run GSEA analysis (includes MSigDB_Hallmark_2020 by default)')

    args = parser.parse_args()
    run_analysis(args)


if __name__ == "__main__":
    main()
