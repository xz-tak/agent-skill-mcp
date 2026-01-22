#!/usr/bin/env python3
"""
Analyze Cortellis data for IBD-associated genes and calculate drug development scores.

This script demonstrates how to analyze Cortellis API v2.0 JSON data for target scoring.

**Schema Compatibility:**
- Works with Cortellis API v2.0 JSON structure
- Compatible with both API-generated and Excel-converted JSON files
- Handles both IndicationsPrimary (API v2.0) and IndicationsSecondary (legacy)
- Supports single object and array indication formats

**Key Schema Features:**
- DrugRecord is a dictionary keyed by drug ID (not drug name)
- Phase format: {"@id": "C3", "$": "Phase 3 Clinical"}
- Indications can be single object or array (script handles both)
- Values extracted from "$" field in nested objects

See ../references/json_schema.md for complete schema documentation.
"""

import json
import os
from pathlib import Path
from collections import defaultdict
import pandas as pd

# IBD-related indications (for targeted scoring)
IBD_INDICATIONS = [
    'inflammatory bowel disease', 'ibd', "crohn's disease", 'crohn',
    'ulcerative colitis', 'uc', 'colitis'
]

# Gene lists as specified
GENE_LISTS = {
    'List 1': ['TYK2', 'JAK1'],
    'List 2': ['TNFRSF25', 'GREM1'],
    'List 3': ['TNFRSF25', 'PCOLCE'],
    'List 4': ['CDKN2D', 'ITGA4', 'ITGB7'],
    'List 5': ['CDKN2D', 'PCOLCE']
}

# Helper functions for schema compatibility (from json_schema.md)
def extract_value(obj, default=''):
    """Extract value from dict with '$' key or return the object itself."""
    if isinstance(obj, dict):
        return obj.get('$', default)
    return obj if obj is not None else default

def as_list(x):
    """Convert to list if not already."""
    if x is None:
        return []
    return x if isinstance(x, list) else [x]

def is_ibd_indication(indication_text):
    """Check if indication text mentions IBD-related terms."""
    if not indication_text:
        return False
    indication_lower = str(indication_text).lower()
    return any(ibd_term in indication_lower for ibd_term in IBD_INDICATIONS)

def parse_phase(phase_str):
    """Parse phase string and return standardized phase."""
    if not phase_str:
        return 'Unknown'

    phase_lower = str(phase_str).lower()

    # Approved/On-market
    if any(term in phase_lower for term in ['approved', 'launched', 'marketed', 'on market']):
        return 'Approved'

    # Clinical phases
    if 'phase 3' in phase_lower or 'phase iii' in phase_lower:
        return 'Phase 3'
    if 'phase 2' in phase_lower or 'phase ii' in phase_lower:
        return 'Phase 2'
    if 'phase 1' in phase_lower or 'phase i' in phase_lower:
        return 'Phase 1'

    # Preclinical
    if 'preclinical' in phase_lower or 'discovery' in phase_lower or 'research' in phase_lower:
        return 'Preclinical'

    # Discontinued/No development
    if 'discontinued' in phase_lower or 'no development' in phase_lower or 'suspended' in phase_lower:
        return 'Discontinued'

    return 'Unknown'

def calculate_drug_score(phase, has_ibd_indication):
    """Calculate score for a single drug based on phase and indication."""
    if phase == 'Discontinued' or phase == 'Unknown':
        return 0

    # Scoring matrix
    if phase == 'Approved':
        return 7 if has_ibd_indication else 4
    elif phase == 'Phase 3':
        return 3 if has_ibd_indication else 2
    elif phase == 'Phase 2':
        return 2 if has_ibd_indication else 1
    elif phase == 'Phase 1':
        return 1 if has_ibd_indication else 0.5
    elif phase == 'Preclinical':
        return 0.1

    return 0

def analyze_gene(gene_name):
    """Analyze a single gene's Cortellis data."""
    json_file = f"{gene_name}_cortellis_data.json"

    if not os.path.exists(json_file):
        return None

    with open(json_file, 'r') as f:
        data = json.load(f)

    results = {
        'gene': gene_name,
        'total_drugs': 0,
        'drugs': [],
        'score_breakdown': defaultdict(int),
        'total_score': 0
    }

    # Analyze DrugRecord (comprehensive data from Investigational Drugs API)
    # IMPORTANT: DrugRecord is nested under data['drug']['DrugRecord']
    # IMPORTANT: DrugRecord is a dictionary keyed by drug ID (not drug name!)
    drug_data = data.get('drug', {})
    drug_records = drug_data.get('DrugRecord', {})

    if not isinstance(drug_records, dict):
        print(f"  Warning: DrugRecord is not a dict for {gene_name}, skipping")
        return results

    # Iterate through drug records (drug_key is the drug ID)
    for drug_key, drug in drug_records.items():
        drug_name = drug.get('DrugName', drug_key)

        # Get phase from PhaseHighest (use helper function)
        phase_data = drug.get('PhaseHighest', {})
        phase_str = extract_value(phase_data)
        phase = parse_phase(phase_str)

        # Check indications
        indications = []
        has_ibd = False

        # Check IndicationsPrimary (API v2.0) or IndicationsSecondary (legacy)
        ind_data = drug.get('IndicationsPrimary') or drug.get('IndicationsSecondary', {})

        # Use helper function to handle single/array pattern
        ind_list = as_list(ind_data.get('Indication', []))

        for ind in ind_list:
            # Extract indication text using helper function
            indication = extract_value(ind)
            if indication:
                indications.append(indication)
                if is_ibd_indication(indication):
                    has_ibd = True

        # Calculate score
        score = calculate_drug_score(phase, has_ibd)

        drug_info = {
            'name': drug_name,
            'phase': phase,
            'indications': indications,
            'ibd_indication': has_ibd,
            'score': score
        }

        results['drugs'].append(drug_info)
        results['total_score'] += score
        results['score_breakdown'][f"{phase}{'_IBD' if has_ibd else ''}"] += 1

    results['total_drugs'] = len(drug_records)

    return results

def main():
    print("=" * 80)
    print("IBD Target Drug Development Analysis")
    print("=" * 80)
    print()

    # Analyze all genes
    all_results = {}
    genes_to_analyze = set()
    for gene_list in GENE_LISTS.values():
        genes_to_analyze.update(gene_list)

    for gene in sorted(genes_to_analyze):
        print(f"Analyzing {gene}...")
        result = analyze_gene(gene)
        if result:
            all_results[gene] = result
            print(f"  ✓ Found {result['total_drugs']} drugs, Total Score: {result['total_score']:.1f}")
        else:
            print(f"  ✗ No data available")

    print()
    print("=" * 80)
    print("Gene Rankings by Score")
    print("=" * 80)
    print()

    # Sort genes by score
    sorted_genes = sorted(
        all_results.items(),
        key=lambda x: x[1]['total_score'],
        reverse=True
    )

    print(f"{'Rank':<6} {'Gene':<12} {'Total Score':<12} {'Total Drugs':<12} {'Key Phases'}")
    print("-" * 80)

    for rank, (gene, data) in enumerate(sorted_genes, 1):
        key_phases = ', '.join([f"{k}: {v}" for k, v in sorted(data['score_breakdown'].items(), key=lambda x: x[1], reverse=True)[:3]])
        print(f"{rank:<6} {gene:<12} {data['total_score']:<12.1f} {data['total_drugs']:<12} {key_phases}")

    print()
    print("=" * 80)
    print("Gene List Scores")
    print("=" * 80)
    print()

    # Calculate scores for each list
    list_scores = {}
    for list_name, genes in GENE_LISTS.items():
        total_score = sum(all_results.get(gene, {}).get('total_score', 0) for gene in genes)
        available_genes = [g for g in genes if g in all_results]
        missing_genes = [g for g in genes if g not in all_results]

        list_scores[list_name] = {
            'genes': genes,
            'available_genes': available_genes,
            'missing_genes': missing_genes,
            'total_score': total_score
        }

    # Sort lists by score
    sorted_lists = sorted(list_scores.items(), key=lambda x: x[1]['total_score'], reverse=True)

    print(f"{'Rank':<6} {'List':<12} {'Genes':<30} {'Total Score':<12} {'Status'}")
    print("-" * 100)

    for rank, (list_name, data) in enumerate(sorted_lists, 1):
        genes_str = ', '.join(data['genes'])
        status = "Complete" if not data['missing_genes'] else f"Missing: {', '.join(data['missing_genes'])}"
        print(f"{rank:<6} {list_name:<12} {genes_str:<30} {data['total_score']:<12.1f} {status}")

    print()

    # Save detailed results
    output_data = {
        'gene_results': all_results,
        'gene_rankings': [(gene, data['total_score']) for gene, data in sorted_genes],
        'list_scores': list_scores
    }

    with open('ibd_analysis_results.json', 'w') as f:
        json.dump(output_data, f, indent=2)

    print("Detailed results saved to: ibd_analysis_results.json")
    print()

    return output_data

if __name__ == '__main__':
    main()
