#!/usr/bin/env python3
"""
Analyze OFF-X safety data from Excel files
Calculate average safety scores per gene and generate comprehensive report
"""

import pandas as pd
import os
import glob
import json
from pathlib import Path
from collections import defaultdict
import numpy as np

# Score mapping
SCORE_MAP = {
    'Very high': -10,
    'High': -8,
    'Medium': -6,
    'Low': -2,
    'Very low': -1,
    'Not associated': 2,
    '': 0,
    'NA': 0,
    'None': 0,
    np.nan: 0
}

def map_score(value):
    """Map score label to numerical value"""
    if pd.isna(value) or value == '':
        return 0
    value_str = str(value).strip()
    return SCORE_MAP.get(value_str, 0)

def is_medium_high(value):
    """Check if score is Medium, High, or Very high"""
    if pd.isna(value) or value == '':
        return False
    value_str = str(value).strip()
    return value_str in ['Medium', 'High', 'Very high']

def is_na(value):
    """Check if value is NA/empty"""
    return pd.isna(value) or str(value).strip() in ['', 'NA', 'None']

def extract_gene_from_filepath(filepath):
    """Extract gene name from filepath by reading metadata.json"""
    # Look for metadata.json in the same directory as the Excel file
    excel_dir = Path(filepath).parent
    metadata_path = excel_dir / "metadata.json"

    if metadata_path.exists():
        try:
            with open(metadata_path, 'r') as f:
                metadata = json.load(f)
                entity = metadata.get('entity', '')
                if entity:
                    return entity.upper()
        except Exception as e:
            print(f"Warning: Could not read metadata for {filepath}: {e}")

    # Fallback to directory name parsing (e.g., "il17a_targets_2026-01-02_17-36-47")
    dir_name = excel_dir.name
    # Extract first part before underscore
    parts = dir_name.split('_')
    if parts:
        return parts[0].upper()

    return "UNKNOWN"

def analyze_excel_file(filepath):
    """Analyze a single Excel file for safety scores"""
    try:
        # Try to read the "Data" sheet
        df = pd.read_excel(filepath, sheet_name='Data')

        # Find the score column - look for column with "Score" and "label" in name
        score_col = None
        for col in df.columns:
            if 'score' in col.lower() and 'label' in col.lower():
                score_col = col
                break

        if score_col is None:
            print(f"Warning: No score column found in {filepath}")
            return None

        # Calculate statistics
        scores = df[score_col].apply(map_score)

        result = {
            'filepath': filepath,
            'filename': Path(filepath).name,
            'total_records': len(df),
            'scores': scores.tolist(),
            'average_score': scores.mean(),
            'medium_high_count': df[score_col].apply(is_medium_high).sum(),
            'medium_high_pct': (df[score_col].apply(is_medium_high).sum() / len(df) * 100) if len(df) > 0 else 0,
            'na_count': df[score_col].apply(is_na).sum(),
            'na_pct': (df[score_col].apply(is_na).sum() / len(df) * 100) if len(df) > 0 else 0,
            'score_distribution': df[score_col].value_counts().to_dict()
        }

        return result

    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return None

def group_by_gene(file_results):
    """Group results by gene"""
    gene_data = defaultdict(list)

    for result in file_results:
        if result is None:
            continue
        gene = extract_gene_from_filepath(result['filepath'])
        gene_data[gene].append(result)

    return gene_data

def calculate_gene_statistics(gene_results):
    """Calculate aggregated statistics for a gene across all files"""
    all_scores = []
    total_records = 0
    total_medium_high = 0
    total_na = 0

    for result in gene_results:
        all_scores.extend(result['scores'])
        total_records += result['total_records']
        total_medium_high += result['medium_high_count']
        total_na += result['na_count']

    avg_score = np.mean(all_scores) if all_scores else 0
    medium_high_pct = (total_medium_high / total_records * 100) if total_records > 0 else 0
    na_pct = (total_na / total_records * 100) if total_records > 0 else 0

    return {
        'average_score': avg_score,
        'total_records': total_records,
        'total_adverse_events': total_records,
        'medium_high_pct': medium_high_pct,
        'na_pct': na_pct,
        'file_count': len(gene_results)
    }

def interpret_score(score, na_pct):
    """Interpret safety score"""
    if na_pct > 90:
        return "Insufficient data (>90% NA)"
    elif score < -1.0:
        return "High concern - significant safety liabilities"
    elif score < -0.5:
        return "Moderate concern - notable adverse events"
    elif score < 0:
        return "Low concern - manageable safety profile"
    else:
        return "Excellent - minimal safety concerns"

def generate_markdown_report(gene_stats, lists):
    """Generate comprehensive markdown report"""

    # Sort genes by score
    sorted_genes = sorted(gene_stats.items(), key=lambda x: x[1]['average_score'])

    report = []
    report.append("# OFF-X Safety Analysis Report")
    report.append(f"\n## Executive Summary\n")
    report.append(f"Analyzed {len(gene_stats)} genes across {sum(g['file_count'] for g in gene_stats.values())} OFF-X Master view files.\n")

    # Overall rankings
    report.append("\n## Gene Safety Rankings (Best to Worst)\n")
    report.append("| Rank | Gene | Avg Score | Med-High % | NA % | Total AEs | Interpretation |")
    report.append("|------|------|-----------|------------|------|-----------|----------------|")

    for rank, (gene, stats) in enumerate(reversed(sorted_genes), 1):
        interpretation = interpret_score(stats['average_score'], stats['na_pct'])
        report.append(f"| {rank} | **{gene}** | {stats['average_score']:.2f} | {stats['medium_high_pct']:.1f}% | {stats['na_pct']:.1f}% | {stats['total_adverse_events']} | {interpretation} |")

    # Detailed gene analysis
    report.append("\n## Detailed Gene Analysis\n")

    for gene, stats in sorted_genes:
        report.append(f"\n### {gene}")
        report.append(f"- **Average Safety Score**: {stats['average_score']:.2f}")
        report.append(f"- **Total Adverse Events**: {stats['total_adverse_events']}")
        report.append(f"- **Medium-High Severity**: {stats['medium_high_pct']:.1f}% ({int(stats['medium_high_pct'] * stats['total_adverse_events'] / 100)} events)")
        report.append(f"- **NA/Unknown**: {stats['na_pct']:.1f}% ({int(stats['na_pct'] * stats['total_adverse_events'] / 100)} events)")
        report.append(f"- **Files Analyzed**: {stats['file_count']}")
        report.append(f"- **Interpretation**: {interpret_score(stats['average_score'], stats['na_pct'])}")

    # List-specific analysis
    report.append("\n## Combination Analysis\n")

    for list_name, genes in lists.items():
        report.append(f"\n### {list_name}: {', '.join(genes)}")

        # Calculate combo statistics
        list_genes_stats = [gene_stats[g] for g in genes if g in gene_stats]
        if not list_genes_stats:
            report.append("*No data available for this combination*\n")
            continue

        avg_combo_score = np.mean([g['average_score'] for g in list_genes_stats])
        avg_medium_high = np.mean([g['medium_high_pct'] for g in list_genes_stats])
        avg_na = np.mean([g['na_pct'] for g in list_genes_stats])
        total_aes = sum([g['total_adverse_events'] for g in list_genes_stats])

        report.append(f"- **Combined Average Score**: {avg_combo_score:.2f}")
        report.append(f"- **Average Medium-High %**: {avg_medium_high:.1f}%")
        report.append(f"- **Average NA %**: {avg_na:.1f}%")
        report.append(f"- **Total Adverse Events (all targets)**: {total_aes}")
        report.append(f"- **Clinical Validation Status**: {'Well-characterized' if avg_na < 30 else 'Limited data'}")
        report.append(f"- **Competitive Intelligence**: {interpret_combo(avg_combo_score, avg_medium_high)}")

        # Individual gene breakdown for combo
        report.append(f"\n**Individual Genes in {list_name}:**")
        for gene in genes:
            if gene in gene_stats:
                s = gene_stats[gene]
                report.append(f"- {gene}: Score {s['average_score']:.2f}, Med-High {s['medium_high_pct']:.1f}%, NA {s['na_pct']:.1f}%")
            else:
                report.append(f"- {gene}: No data")

    return "\n".join(report)

def interpret_combo(score, medium_high_pct):
    """Interpret combination for competitive intelligence"""
    if score < -1.0 and medium_high_pct > 10:
        return "High safety risk - may face development challenges. Competitors likely targeting due to mechanism, but safety concerns present significant hurdle."
    elif score < -0.5:
        return "Moderate safety risk - manageable with careful patient selection and monitoring. Competitive space with known safety profile."
    elif score < 0:
        return "Favorable safety profile - good target for development. Competitive advantages if efficacy demonstrated."
    else:
        return "Excellent safety profile - strong competitive position. High probability of clinical success if efficacy proven."

def main():
    # Find all Excel files
    excel_files = glob.glob("offx_playwright_result/**/*.xlsx", recursive=True)
    print(f"Found {len(excel_files)} Excel files to analyze\n")

    # Analyze each file
    file_results = []
    for filepath in excel_files:
        print(f"Analyzing: {Path(filepath).name}")
        result = analyze_excel_file(filepath)
        if result:
            file_results.append(result)

    print(f"\nSuccessfully analyzed {len(file_results)} files\n")

    # Group by gene
    gene_data = group_by_gene(file_results)
    print(f"Found data for {len(gene_data)} genes\n")

    # Calculate gene statistics
    gene_stats = {}
    for gene, results in gene_data.items():
        gene_stats[gene] = calculate_gene_statistics(results)
        print(f"{gene}: {len(results)} files, avg score {gene_stats[gene]['average_score']:.2f}")

    # Define gene lists
    lists = {
        "List 1": ["IL17A", "IL17F", "TNFSF13B", "TNFSF13"],
        "List 2": ["IL17A", "IL17F", "TNFSF13B"],
        "List 3": ["CXCR5", "CD19", "CD3D", "CD3E", "CD3G"]
    }

    # Normalize gene names in stats to match lists
    # (extract the core gene name, e.g., IL17A from IL-17A)
    normalized_stats = {}
    for gene, stats in gene_stats.items():
        # Try to match to list genes
        matched = False
        for list_genes in lists.values():
            for list_gene in list_genes:
                if list_gene in gene or gene in list_gene:
                    normalized_stats[list_gene] = stats
                    matched = True
                    break
            if matched:
                break
        if not matched:
            normalized_stats[gene] = stats

    # Generate report
    report = generate_markdown_report(normalized_stats, lists)

    # Save report
    report_path = "offx_safety_analysis_report.md"
    with open(report_path, 'w') as f:
        f.write(report)

    print(f"\n✅ Report saved to: {report_path}")

    # Save JSON
    json_path = "offx_safety_analysis_data.json"
    with open(json_path, 'w') as f:
        json.dump({
            'gene_stats': normalized_stats,
            'lists': lists
        }, f, indent=2, default=str)

    print(f"✅ Data saved to: {json_path}")

if __name__ == "__main__":
    main()
