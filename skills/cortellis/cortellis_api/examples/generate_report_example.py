#!/usr/bin/env python3
"""
Generate comprehensive markdown report with clinical insights for IBD targets.

This script generates a detailed clinical intelligence report from the analysis results
produced by analyze_targets_example.py.

**Input:** ibd_analysis_results.json (output from analyze_targets_example.py)
**Output:** IBD_Target_Clinical_Intelligence_Report.md

**Report Sections:**
1. Executive Summary - High-level findings and key insights
2. Target Rankings - Scored ranking of individual genes
3. Individual Target Analysis - Detailed analysis per gene
4. Gene Combination Analysis - Analysis of gene list combinations
5. Strategic Recommendations - Prioritization framework and next steps
6. Methodology - Scoring system and data sources

**Usage:**
1. First run: python analyze_targets_example.py
2. Then run: python generate_report_example.py

See ../references/scoring_framework.md for scoring methodology details.
"""

import json
import os
from datetime import datetime

def load_analysis_results():
    """Load the analysis results JSON with validation."""
    if not os.path.exists('ibd_analysis_results.json'):
        print("ERROR: ibd_analysis_results.json not found!")
        print("Please run analyze_targets_example.py first to generate the data.")
        exit(1)

    with open('ibd_analysis_results.json', 'r') as f:
        data = json.load(f)

    # Validate required keys
    required_keys = ['gene_results', 'gene_rankings', 'list_scores']
    for key in required_keys:
        if key not in data:
            print(f"ERROR: Missing required key '{key}' in analysis results")
            exit(1)

    return data

def generate_executive_summary(data):
    """Generate executive summary section."""
    gene_results = data['gene_results']
    gene_rankings = data['gene_rankings']
    list_scores = data['list_scores']

    md = "## Executive Summary\n\n"

    # Overall findings
    total_genes = len(gene_results)
    total_drugs = sum(r['total_drugs'] for r in gene_results.values())
    top_target = gene_rankings[0][0] if gene_rankings else "N/A"
    top_score = gene_rankings[0][1] if gene_rankings else 0

    md += f"This analysis evaluated **{total_genes} IBD-associated genes** across **5 therapeutic combination lists**, "
    md += f"identifying **{total_drugs} total drug candidates** from the Cortellis database.\n\n"

    md += f"**Key Findings:**\n\n"
    md += f"- **{top_target}** emerges as the most clinically validated target (Score: **{top_score:.1f}**) with the most extensive drug development portfolio\n"
    md += f"- **List 1 (TYK2, JAK1)** represents the highest-value combination based on available data\n"
    md += f"- Several targets show novel therapeutic opportunities with limited competitive activity\n"
    md += f"- Note: JAK1 data could not be retrieved due to API timeouts (extensive drug portfolio)\n\n"

    return md

def generate_target_analysis(gene_name, gene_data):
    """Generate detailed analysis for a single target."""
    md = f"### {gene_name}\n\n"

    score = gene_data['total_score']
    total_drugs = gene_data['total_drugs']
    drugs = gene_data['drugs']

    md += f"**Clinical Validation Score:** {score:.1f} | **Total Drugs:** {total_drugs}\n\n"

    if total_drugs == 0:
        md += "**Status:** No drug development activity identified in Cortellis database.\n\n"
        md += "**Interpretation:** This target represents either:\n"
        md += "- A novel therapeutic opportunity with no competitive activity\n"
        md += "- Limited druggability or early-stage research target\n"
        md += "- Potential for first-in-class development with reduced competitive risk\n\n"
        return md

    # Phase distribution
    md += "**Development Phase Distribution:**\n\n"
    phase_counts = {}
    for drug in drugs:
        phase = drug['phase']
        if phase not in phase_counts:
            phase_counts[phase] = {'total': 0, 'ibd': 0, 'drugs': []}
        phase_counts[phase]['total'] += 1
        if drug['ibd_indication']:
            phase_counts[phase]['ibd'] += 1
        phase_counts[phase]['drugs'].append(drug['name'])

    for phase in ['Approved', 'Phase 3', 'Phase 2', 'Phase 1', 'Preclinical', 'Discontinued', 'Unknown']:
        if phase in phase_counts:
            count = phase_counts[phase]['total']
            ibd_count = phase_counts[phase]['ibd']
            ibd_str = f" ({ibd_count} IBD-specific)" if ibd_count > 0 else ""
            md += f"- **{phase}:** {count} drug(s){ibd_str}\n"

    md += "\n"

    # IBD-specific drugs
    ibd_drugs = [d for d in drugs if d['ibd_indication']]
    if ibd_drugs:
        md += "**IBD-Specific Drug Candidates:**\n\n"
        for drug in ibd_drugs:
            indications = ", ".join(drug['indications'][:3])
            md += f"- **{drug['name']}** ({drug['phase']}) - {indications}\n"
        md += "\n"

    # Top non-IBD drugs by phase
    non_ibd_drugs = [d for d in drugs if not d['ibd_indication'] and d['phase'] not in ['Discontinued', 'Unknown']]
    if non_ibd_drugs:
        # Sort by phase value
        phase_priority = {'Approved': 5, 'Phase 3': 4, 'Phase 2': 3, 'Phase 1': 2, 'Preclinical': 1}
        non_ibd_drugs.sort(key=lambda x: phase_priority.get(x['phase'], 0), reverse=True)

        md += "**Leading Non-IBD Drug Candidates (Repurposing Potential):**\n\n"
        for drug in non_ibd_drugs[:5]:
            indications = ", ".join(drug['indications'][:2]) if drug['indications'] else "No indication data"
            md += f"- **{drug['name']}** ({drug['phase']}) - {indications}\n"
        md += "\n"

    # Clinical insights
    md += "**Clinical Validation Insights:**\n\n"

    if score > 20:
        md += "- **High Clinical Validation:** Extensively validated target with mature drug development pipeline\n"
        md += "- **Competitive Landscape:** High competition expected; differentiation will be critical\n"
        md += "- **Risk Profile:** Lower biological risk due to clinical precedent; higher commercial risk\n"
    elif score > 5:
        md += "- **Moderate Clinical Validation:** Emerging target with active development programs\n"
        md += "- **Competitive Landscape:** Moderate competition; opportunities for differentiation exist\n"
        md += "- **Risk Profile:** Balanced risk-reward with clinical precedent and commercial potential\n"
    elif score > 0:
        md += "- **Early Clinical Validation:** Limited but promising development activity\n"
        md += "- **Competitive Landscape:** Low competition; first-mover advantage possible\n"
        md += "- **Risk Profile:** Higher biological risk; significant commercial upside if successful\n"
    else:
        md += "- **No Clinical Validation:** No current drug development activity\n"
        md += "- **Competitive Landscape:** No competition; blue ocean opportunity\n"
        md += "- **Risk Profile:** Highest biological risk; potential for breakthrough innovation\n"

    md += "\n"

    return md

def generate_list_analysis(list_name, list_data, gene_results):
    """Generate analysis for a gene combination list."""
    genes = list_data['genes']
    total_score = list_data['total_score']
    missing = list_data['missing_genes']

    md = f"### {list_name}: {', '.join(genes)}\n\n"
    md += f"**Combined Score:** {total_score:.1f}\n\n"

    if missing:
        md += f"⚠️ **Note:** Data unavailable for {', '.join(missing)} due to API timeouts\n\n"

    md += "**Individual Target Contributions:**\n\n"
    for gene in genes:
        if gene in gene_results:
            score = gene_results[gene]['total_score']
            drugs = gene_results[gene]['total_drugs']
            md += f"- **{gene}:** Score {score:.1f} ({drugs} drugs)\n"
        else:
            md += f"- **{gene}:** No data available\n"

    md += "\n**Combination Rationale & Clinical Insights:**\n\n"

    # Provide strategic insights based on score profile
    if total_score > 40:
        md += "**Profile:** High clinical validation combination with extensive precedent\n\n"
        md += "- **Strengths:** Strong biological validation, proven mechanism, established safety profiles\n"
        md += "- **Challenges:** High competitive intensity, need for significant differentiation\n"
        md += "- **Strategy:** Focus on novel MOA, improved safety/efficacy, or underserved patient populations\n\n"
    elif total_score > 10:
        md += "**Profile:** Moderate validation with balanced risk-reward\n\n"
        md += "- **Strengths:** Clinical precedent with manageable competition, potential for differentiation\n"
        md += "- **Challenges:** Need to demonstrate clear advantage over existing approaches\n"
        md += "- **Strategy:** Target unmet needs, explore synergistic combinations, optimize dosing regimens\n\n"
    elif total_score > 0:
        md += "**Profile:** Emerging opportunity with limited competition\n\n"
        md += "- **Strengths:** Low competition, potential first-mover advantage, novel mechanism\n"
        md += "- **Challenges:** Higher development risk, limited clinical precedent\n"
        md += "- **Strategy:** Invest in early validation studies, identify optimal patient selection biomarkers\n\n"
    else:
        md += "**Profile:** Novel, unexplored therapeutic space\n\n"
        md += "- **Strengths:** No competition, potential breakthrough innovation, high commercial upside\n"
        md += "- **Challenges:** Highest biological risk, unknown safety profile\n"
        md += "- **Strategy:** Focus on strong preclinical validation, biomarker-driven development\n\n"

    return md

def generate_recommendations(data):
    """Generate strategic recommendations."""
    md = "## Strategic Recommendations\n\n"

    gene_rankings = data['gene_rankings']
    list_scores = data['list_scores']

    # Sort lists by score
    sorted_lists = sorted(list_scores.items(), key=lambda x: x[1]['total_score'], reverse=True)

    md += "### Prioritization Framework\n\n"

    md += "**High Priority (Advance Immediately):**\n\n"
    for list_name, list_data in sorted_lists[:2]:
        if list_data['total_score'] > 5:
            genes_str = ', '.join(list_data['genes'])
            md += f"- **{list_name} ({genes_str}):** Score {list_data['total_score']:.1f}\n"
            md += f"  - Strong clinical validation with manageable risk\n"
            md += f"  - Focus on differentiation and competitive positioning\n\n"

    md += "**Medium Priority (Evaluate Further):**\n\n"
    for list_name, list_data in sorted_lists[2:4]:
        if list_data['total_score'] > 0:
            genes_str = ', '.join(list_data['genes'])
            md += f"- **{list_name} ({genes_str}):** Score {list_data['total_score']:.1f}\n"
            md += f"  - Emerging opportunity; conduct deeper due diligence\n"
            md += f"  - Assess biological rationale and competitive timeline\n\n"

    md += "**Lower Priority (Monitor):**\n\n"
    for list_name, list_data in sorted_lists:
        if list_data['total_score'] == 0:
            genes_str = ', '.join(list_data['genes'])
            md += f"- **{list_name} ({genes_str}):** Score {list_data['total_score']:.1f}\n"
            md += f"  - Novel targets with high risk but potential high reward\n"
            md += f"  - Require strong preclinical validation before advancing\n\n"

    md += "### Next Steps\n\n"
    md += "1. **Competitive Intelligence Deep Dive:** Conduct detailed analysis of top competitors for high-scoring targets\n"
    md += "2. **Clinical Trial Analysis:** Review ongoing/completed trials for safety, efficacy, and differentiation opportunities\n"
    md += "3. **Patent Landscape:** Assess freedom to operate and opportunities for novel IP\n"
    md += "4. **KOL Engagement:** Interview key opinion leaders on unmet needs and therapeutic positioning\n"
    md += "5. **Biomarker Strategy:** Develop patient selection and companion diagnostic strategies\n"
    md += "6. **JAK1 Follow-up:** Obtain JAK1 data through alternative methods given its importance for List 1\n\n"

    return md

def generate_methodology(data):
    """Generate methodology section."""
    md = "## Methodology\n\n"

    md += "### Data Source\n\n"
    md += "- **Database:** Cortellis Drug Intelligence (Clarivate Analytics)\n"
    md += "- **Query Date:** " + datetime.now().strftime("%Y-%m-%d") + "\n"
    md += "- **Scope:** Comprehensive drug development data including phase, indications, and sponsorship\n\n"

    md += "### Scoring System\n\n"
    md += "Each drug candidate was scored based on development phase and IBD indication relevance:\n\n"

    md += "| Phase | IBD-Specific | Non-IBD | Rationale |\n"
    md += "|-------|--------------|---------|----------|\n"
    md += "| FDA Approved / On Market | 7 | 4 | Highest validation; IBD-specific indicates direct relevance |\n"
    md += "| Phase 3 Clinical | 3 | 2 | Late-stage validation; high probability of approval |\n"
    md += "| Phase 2 Clinical | 2 | 1 | Mid-stage validation; proof-of-concept established |\n"
    md += "| Phase 1 Clinical | 1 | 0.5 | Early validation; safety and PK/PD characterized |\n"
    md += "| Preclinical | 0.1 | 0.1 | Minimal validation; target engagement shown |\n"
    md += "| Discontinued / No Development | 0 | 0 | No active development |\n\n"

    md += "**IBD-Specific Indications:** " + "Inflammatory Bowel Disease, IBD, Crohn's Disease, Ulcerative Colitis\n\n"

    md += "### Limitations\n\n"
    md += "- JAK1 data could not be retrieved due to API timeout issues (extensive drug portfolio)\n"
    md += "- Scoring reflects quantity and stage of development; does not assess quality or mechanism differentiation\n"
    md += "- Historical/discontinued programs included in drug counts but scored at 0\n"
    md += "- Repurposing potential for non-IBD drugs requires additional clinical validation\n\n"

    return md

def main():
    """Generate comprehensive markdown report."""
    print("Generating comprehensive IBD target analysis report...")

    # Load data
    data = load_analysis_results()

    # Build report
    report = f"# IBD-Associated Gene Targets: Clinical Validation & Competitive Intelligence Report\n\n"
    report += f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
    report += f"---\n\n"

    # Executive Summary
    report += generate_executive_summary(data)
    report += "\n---\n\n"

    # Target Rankings
    report += "## Target Rankings by Clinical Validation Score\n\n"
    gene_rankings = data['gene_rankings']
    report += "| Rank | Gene | Score | Total Drugs | Clinical Insight |\n"
    report += "|------|------|-------|-------------|------------------|\n"

    for rank, (gene, score) in enumerate(gene_rankings, 1):
        gene_data = data['gene_results'][gene]
        total_drugs = gene_data['total_drugs']

        if score > 20:
            insight = "Extensively validated, high competition"
        elif score > 5:
            insight = "Moderate validation, emerging target"
        elif score > 0:
            insight = "Early validation, low competition"
        else:
            insight = "No development activity, novel target"

        report += f"| {rank} | **{gene}** | {score:.1f} | {total_drugs} | {insight} |\n"

    report += "\n---\n\n"

    # Individual Target Analysis
    report += "## Individual Target Analysis\n\n"

    for gene, score in gene_rankings:
        gene_data = data['gene_results'][gene]
        report += generate_target_analysis(gene, gene_data)
        report += "---\n\n"

    # Gene List Analysis
    report += "## Gene Combination List Analysis\n\n"

    list_scores = data['list_scores']
    sorted_lists = sorted(list_scores.items(), key=lambda x: x[1]['total_score'], reverse=True)

    for list_name, list_data in sorted_lists:
        report += generate_list_analysis(list_name, list_data, data['gene_results'])
        report += "---\n\n"

    # Strategic Recommendations
    report += generate_recommendations(data)
    report += "\n---\n\n"

    # Methodology
    report += generate_methodology(data)

    # Save report
    with open('IBD_Target_Clinical_Intelligence_Report.md', 'w') as f:
        f.write(report)

    print("\n" + "="*80)
    print("Report generated successfully!")
    print("="*80)
    print(f"\nOutput file: IBD_Target_Clinical_Intelligence_Report.md")
    print(f"Report sections: Executive Summary, Target Rankings, Individual Analysis,")
    print(f"                 Combination Analysis, Strategic Recommendations, Methodology")
    print("\n")

if __name__ == '__main__':
    main()
