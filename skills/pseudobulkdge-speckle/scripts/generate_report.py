#!/usr/bin/env python3
"""
Report Data Preparation Script for Pseudobulk DEG Analysis

================================================================================
IMPORTANT: This script ONLY prepares filtered data files.
The biological interpretation is done by Claude Code at runtime using:
- Ultrathink: For deep analysis of gene functions, pathway connections, etc.
- WebFetch: For literature support and validation

After running this script, Claude Code will:
1. Read each {comparison}_top_deg.tsv, {comparison}_top_pathways.tsv, etc.
2. Use ultrathink to generate professional biologist interpretation
3. Use WebFetch to query relevant literature
4. Append the interpretation directly to the markdown report
================================================================================

Usage:
    python generate_report.py \
        --output_dir <output_dir> \
        --fdr_cutoff 0.05 \
        --logfc_cutoff 0.5 \
        --nes_cutoff 0 \
        --top_k 30
"""

import argparse
import os
import glob
import pandas as pd
from pathlib import Path


def filter_top_deg(df: pd.DataFrame, top_k: int = 30, fdr_cutoff: float = 0.05,
                   logfc_cutoff: float = 0.5) -> pd.DataFrame:
    """
    Filter DEG results to top significant genes.

    Args:
        df: DataFrame with DEG results (must have FDR, logFC columns)
        top_k: Total number of top results to return
        fdr_cutoff: FDR threshold for significance
        logfc_cutoff: Absolute log2FC threshold

    Returns:
        DataFrame with top up and down regulated genes
    """
    df = df.copy()
    if 'FDR' not in df.columns and 'padj' in df.columns:
        df['FDR'] = df['padj']
    if 'logFC' not in df.columns and 'log2FoldChange' in df.columns:
        df['logFC'] = df['log2FoldChange']

    # Filter significant genes
    sig = df[(df['FDR'] < fdr_cutoff) & (df['logFC'].abs() > logfc_cutoff)].copy()

    if len(sig) == 0:
        print(f"  Warning: No genes pass filters (FDR < {fdr_cutoff}, |logFC| > {logfc_cutoff})")
        sig = df.nsmallest(top_k, 'FDR')
        return sig

    # Split into up and down regulated
    up = sig[sig['logFC'] > 0].nlargest(top_k // 2, 'logFC')
    down = sig[sig['logFC'] < 0].nsmallest(top_k // 2, 'logFC')

    return pd.concat([up, down])


def filter_top_pathways(df: pd.DataFrame, top_k: int = 30, fdr_cutoff: float = 0.05,
                        nes_cutoff: float = 0) -> pd.DataFrame:
    """
    Filter pathway enrichment results to top significant pathways.

    Args:
        df: DataFrame with pathway results (must have padj, NES columns)
        top_k: Total number of top results to return
        fdr_cutoff: FDR threshold for significance
        nes_cutoff: Absolute NES threshold

    Returns:
        DataFrame with top enriched and depleted pathways
    """
    df = df.copy()

    if 'padj' not in df.columns and 'FDR' in df.columns:
        df['padj'] = df['FDR']

    # Filter significant pathways
    sig = df[(df['padj'] < fdr_cutoff) & (df['NES'].abs() > nes_cutoff)].copy()

    if len(sig) == 0:
        print(f"  Warning: No pathways pass filters (padj < {fdr_cutoff}, |NES| > {nes_cutoff})")
        sig = df.nsmallest(top_k, 'padj')
        return sig

    # Split into positive and negative NES
    up = sig[sig['NES'] > 0].nlargest(top_k // 2, 'NES')
    down = sig[sig['NES'] < 0].nsmallest(top_k // 2, 'NES')

    return pd.concat([up, down])


def filter_speckle_results(df: pd.DataFrame, fdr_cutoff: float = 0.05) -> pd.DataFrame:
    """
    Filter speckle cell composition results to significant cell types.

    Args:
        df: DataFrame with speckle results
        fdr_cutoff: FDR threshold for significance

    Returns:
        DataFrame with significant cell type changes
    """
    df = df.copy()

    fdr_col = None
    for col in ['FDR', 'padj', 'P.Value', 'pval']:
        if col in df.columns:
            fdr_col = col
            break

    if fdr_col is None:
        print("  Warning: No FDR column found in speckle results")
        return df

    sig = df[df[fdr_col] < fdr_cutoff].copy()

    if len(sig) == 0:
        print(f"  Warning: No cell types pass FDR < {fdr_cutoff}, returning all results")
        return df

    return sig.sort_values(fdr_col)


def main():
    parser = argparse.ArgumentParser(
        description='Prepare filtered results for Claude Code ultrathink interpretation',
        epilog='''
================================================================================
NOTE: This script ONLY prepares data files.
Claude Code will read these files and use ultrathink + WebFetch for interpretation.
================================================================================
        '''
    )
    parser.add_argument('--output_dir', required=True, help='Pipeline output directory')
    parser.add_argument('--fdr_cutoff', type=float, default=0.05, help='FDR cutoff')
    parser.add_argument('--logfc_cutoff', type=float, default=0.5, help='|log2FC| cutoff')
    parser.add_argument('--nes_cutoff', type=float, default=0, help='|NES| cutoff')
    parser.add_argument('--top_k', type=int, default=30, help='Top k results per comparison')

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    report_dir = output_dir / 'report_data'
    report_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("PREPARING DATA FOR CLAUDE CODE ULTRATHINK INTERPRETATION")
    print("=" * 70)
    print()
    print("This script prepares filtered data files.")
    print("After this completes, Claude Code will:")
    print("  1. Read each filtered file")
    print("  2. Use ULTRATHINK for deep biological interpretation")
    print("  3. Use WEBFETCH for literature validation")
    print("  4. Generate professional biologist-level analysis")
    print()
    print(f"Output directory: {output_dir}")
    print(f"Report data directory: {report_dir}")
    print(f"\nFilter parameters (user-aligned):")
    print(f"  FDR cutoff: {args.fdr_cutoff}")
    print(f"  |log2FC| cutoff: {args.logfc_cutoff}")
    print(f"  |NES| cutoff: {args.nes_cutoff}")
    print(f"  Top k: {args.top_k}")

    # Find DEG results
    print("\n" + "-" * 60)
    print("Processing DEG results...")
    print("-" * 60)

    deg_files = list(output_dir.glob('pseudobulk_*_DE.tsv'))
    comparisons_found = set()

    if deg_files:
        deg_file = deg_files[0]
        print(f"  Reading: {deg_file}")
        deg_df = pd.read_csv(deg_file, sep='\t')

        comparisons = deg_df['comparison'].unique()
        print(f"  Found {len(comparisons)} comparisons")

        for comparison in comparisons:
            comparisons_found.add(comparison)
            print(f"\n  Processing comparison: {comparison}")
            comp_df = deg_df[deg_df['comparison'] == comparison].copy()

            clusters = comp_df['cluster'].unique()

            all_top_deg = []
            for cluster in clusters:
                cluster_df = comp_df[comp_df['cluster'] == cluster]
                top_df = filter_top_deg(cluster_df, args.top_k, args.fdr_cutoff, args.logfc_cutoff)
                all_top_deg.append(top_df)

            if all_top_deg:
                combined_top = pd.concat(all_top_deg)
                out_file = report_dir / f'{comparison}_top_deg.tsv'
                combined_top.to_csv(out_file, sep='\t', index=False)
                print(f"    -> {out_file.name} ({len(combined_top)} genes)")
    else:
        print("  Warning: No DEG files found")

    # Find pathway enrichment results
    print("\n" + "-" * 60)
    print("Processing pathway enrichment results...")
    print("-" * 60)

    fgsea_dir = output_dir / 'fgsea_groups_clusters'
    if fgsea_dir.exists():
        pathway_files = list(fgsea_dir.glob('*.tsv'))

        for pathway_file in pathway_files:
            print(f"\n  Reading: {pathway_file.name}")
            pathway_df = pd.read_csv(pathway_file, sep='\t')

            if 'comparison' in pathway_df.columns:
                comparisons = pathway_df['comparison'].unique()

                for comparison in comparisons:
                    comp_df = pathway_df[pathway_df['comparison'] == comparison].copy()
                    top_df = filter_top_pathways(comp_df, args.top_k, args.fdr_cutoff, args.nes_cutoff)

                    db_name = pathway_file.stem.split('_')[-1] if '_' in pathway_file.stem else 'pathway'
                    out_file = report_dir / f'{comparison}_top_pathways_{db_name}.tsv'
                    top_df.to_csv(out_file, sep='\t', index=False)
                    print(f"    -> {out_file.name} ({len(top_df)} pathways)")
            else:
                top_df = filter_top_pathways(pathway_df, args.top_k, args.fdr_cutoff, args.nes_cutoff)
                out_file = report_dir / f'top_pathways_{pathway_file.stem}.tsv'
                top_df.to_csv(out_file, sep='\t', index=False)
                print(f"    -> {out_file.name} ({len(top_df)} pathways)")
    else:
        print("  Warning: No fgsea_groups_clusters directory found")

    # Find custom signature results
    custom_sig_file = output_dir / 'fgsea_custom_signatures.tsv'
    if custom_sig_file.exists():
        print("\n  Processing custom signatures...")
        custom_df = pd.read_csv(custom_sig_file, sep='\t')

        if 'comparison' in custom_df.columns:
            comparisons = custom_df['comparison'].unique()
            for comparison in comparisons:
                comp_df = custom_df[custom_df['comparison'] == comparison].copy()
                top_df = filter_top_pathways(comp_df, args.top_k, args.fdr_cutoff, args.nes_cutoff)
                out_file = report_dir / f'{comparison}_top_custom_signatures.tsv'
                top_df.to_csv(out_file, sep='\t', index=False)
                print(f"    -> {out_file.name} ({len(top_df)} signatures)")

    # Find speckle results
    print("\n" + "-" * 60)
    print("Processing speckle cell composition results...")
    print("-" * 60)

    speckle_dir = output_dir / 'speckle_diffprop'
    if speckle_dir.exists():
        speckle_files = list(speckle_dir.glob('*.tsv'))

        for speckle_file in speckle_files:
            print(f"\n  Reading: {speckle_file.name}")
            speckle_df = pd.read_csv(speckle_file, sep='\t', index_col=0)

            comparison = speckle_file.stem.replace('comb_condition_', '')

            top_df = filter_speckle_results(speckle_df, args.fdr_cutoff)
            out_file = report_dir / f'{comparison}_speckle.tsv'
            top_df.to_csv(out_file, sep='\t')
            print(f"    -> {out_file.name} ({len(top_df)} cell types)")
    else:
        print("  Warning: No speckle_diffprop directory found")

    # Generate instruction file for Claude Code
    print("\n" + "-" * 60)
    print("Generating Claude Code instruction file...")
    print("-" * 60)

    instruction_file = report_dir / 'CLAUDE_CODE_INSTRUCTIONS.md'
    with open(instruction_file, 'w') as f:
        f.write("# Claude Code Ultrathink Interpretation Instructions\n\n")
        f.write("## Overview\n\n")
        f.write("The filtered data files have been prepared. Claude Code MUST now:\n\n")
        f.write("1. Read each comparison's filtered files\n")
        f.write("2. Use **ultrathink** to analyze and interpret the biological significance\n")
        f.write("3. Use **WebFetch** to query relevant literature for validation\n")
        f.write("4. Append the interpretation to the final markdown report\n\n")

        f.write("## Files to Process\n\n")
        for fpath in sorted(report_dir.glob('*.tsv')):
            f.write(f"- `{fpath.name}`\n")

        f.write("\n## For Each Comparison\n\n")
        f.write("Claude Code should analyze:\n\n")
        f.write("### DEG Interpretation\n")
        f.write("- Top upregulated genes: What are their functions? Disease associations?\n")
        f.write("- Top downregulated genes: What biological processes are reduced?\n")
        f.write("- Cross-cluster patterns: Any genes consistently DE across cell types?\n\n")

        f.write("### Pathway Interpretation\n")
        f.write("- What biological processes are enriched?\n")
        f.write("- How do the pathways connect mechanistically?\n")
        f.write("- What therapeutic implications emerge?\n\n")

        f.write("### Cell Composition Interpretation\n")
        f.write("- Which cell types show significant proportion changes?\n")
        f.write("- What does this suggest about the disease mechanism?\n")
        f.write("- How do composition changes relate to the DEG findings?\n\n")

        f.write("### Literature Validation (WebFetch)\n")
        f.write("- Query top genes in context of the disease/condition\n")
        f.write("- Find supporting evidence from recent publications\n")
        f.write("- Identify potential drug targets mentioned in literature\n\n")

        f.write("## Filter Parameters Used\n\n")
        f.write(f"- FDR cutoff: {args.fdr_cutoff}\n")
        f.write(f"- |log2FC| cutoff: {args.logfc_cutoff}\n")
        f.write(f"- |NES| cutoff: {args.nes_cutoff}\n")
        f.write(f"- Top k results: {args.top_k}\n")

    print(f"  -> {instruction_file.name}")

    # Generate summary
    summary_file = report_dir / 'report_summary.txt'
    files_generated = list(report_dir.glob('*.tsv'))
    with open(summary_file, 'w') as f:
        f.write("PSEUDOBULK DEG ANALYSIS - FILTERED RESULTS SUMMARY\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Output directory: {output_dir}\n")
        f.write(f"Report data directory: {report_dir}\n\n")
        f.write("Filter parameters:\n")
        f.write(f"  FDR cutoff: {args.fdr_cutoff}\n")
        f.write(f"  |log2FC| cutoff: {args.logfc_cutoff}\n")
        f.write(f"  |NES| cutoff: {args.nes_cutoff}\n")
        f.write(f"  Top k: {args.top_k}\n\n")
        f.write("Files generated:\n")
        for fname in sorted([f.name for f in files_generated]):
            f.write(f"  - {fname}\n")
        f.write("\nNEXT STEP: Claude Code reads these files and uses ultrathink/WebFetch\n")
        f.write("to generate professional biological interpretation.\n")

    print(f"  -> {summary_file.name}")

    print("\n" + "=" * 70)
    print("DATA PREPARATION COMPLETE")
    print("=" * 70)
    print(f"\nFiltered results saved to: {report_dir}")
    print(f"Total files generated: {len(files_generated)}")
    print()
    print("NEXT STEP: Claude Code will now read these files and use")
    print("ULTRATHINK + WEBFETCH to generate biological interpretation.")
    print()


if __name__ == '__main__':
    main()
