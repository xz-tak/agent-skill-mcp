#!/usr/bin/env python3
"""
Filter Validation Script for Omicsoft DEG Analysis
Tests filters step-by-step before running full analysis

Usage:
    python validate_filters.py --file <path> --diseases <terms> --tissues <terms> [options]
"""

import os
import sys
import argparse
import pandas as pd
import scanpy as sc

# Import from deg_analysis.py
sys.path.insert(0, os.path.dirname(__file__))
from deg_analysis import load_h5ad, is_s3_uri

# Check if SOMA support is available
try:
    from soma_loader import (
        SOMA_AVAILABLE,
        is_soma_uri,
        try_soma_open,
        detect_soma_structure,
        get_soma_context,
        get_soma_schema_info,
        prefetch_obs_values,
        translate_substring_filter,
        build_soma_filter,
        SOMAStructureType
    )
except ImportError:
    SOMA_AVAILABLE = False


def validate_filters_soma(args):
    """
    Validate filters for SOMA experiment using server-side queries.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments with filter parameters

    Returns
    -------
    dict
        Validation results with filter status and suggestions
    """
    import tiledbsoma

    print("=" * 80)
    print("FILTER VALIDATION - SOMA Experiment (Step-by-Step)")
    print("=" * 80)

    validation_results = {
        'total_obs': 0,
        'filters_applied': [],
        'final_obs': None,
        'success': True,
        'failed_filter': None,
        'suggestions': [],
        'data_source': 'soma'
    }

    # Get SOMA context (use AWS_DEFAULT_REGION env var or default to us-east-1)
    region = os.environ.get("AWS_DEFAULT_REGION", "us-east-1")
    context = get_soma_context(region=region)

    print(f"\n[STEP 0] Opening SOMA experiment...")
    print(f"  URI: {args.file}")

    try:
        # Get schema info
        schema_info = get_soma_schema_info(args.file, context)
        structure_type = schema_info['structure_type']

        print(f"  Structure type: {structure_type}")
        print(f"  Total observations: {schema_info['n_obs']:,}")

        validation_results['total_obs'] = schema_info['n_obs']
        validation_results['structure_type'] = structure_type

        if 'RNA' in schema_info['measurements']:
            print(f"  Total variables: {schema_info['measurements']['RNA']['n_vars']:,}")
            print(f"  X layers: {schema_info['measurements']['RNA']['x_layers']}")

    except Exception as e:
        print(f"  Error opening SOMA: {e}")
        validation_results['success'] = False
        validation_results['suggestions'].append(f"Failed to open SOMA experiment: {e}")
        return validation_results

    # Define filter columns and their match types
    filter_specs = [
        ('diseases', 'disease', 'substring', False),
        ('exclude_diseases', 'disease', 'substring', True),
        ('tissues', 'tissue', 'substring', False),
        ('studies', 'study', 'exact', False),
        ('comparison_category', 'comparison_category', 'exact', False),
        ('case_treatment', 'case_treatment', 'exact', False),
        ('comparison', 'comparison', 'substring', False),
    ]

    # Pre-fetch all needed columns
    needed_columns = ['disease', 'tissue', 'study', 'comparison_category', 'case_treatment', 'comparison']
    print(f"\n[STEP 1] Pre-fetching unique values from SOMA obs columns...")

    with tiledbsoma.open(args.file, mode="r", context=context) as exp:
        unique_values = prefetch_obs_values(exp, needed_columns)

    for col, values in unique_values.items():
        print(f"  {col}: {len(values)} unique values")

    # Build cumulative filter and validate step by step
    filter_config = {}
    step_num = 2

    for arg_name, col_name, match_type, is_negated in filter_specs:
        filter_value = getattr(args, arg_name.replace('-', '_'), None)

        if not filter_value:
            continue

        print(f"\n[STEP {step_num}] {arg_name.replace('_', ' ').title()} Filter")
        print(f"  Query: {filter_value}")
        print(f"  Column: {col_name}")
        print(f"  Match type: {match_type}")

        search_terms = [v.strip() for v in filter_value.split(',') if v.strip()]
        available = unique_values.get(col_name, [])

        if match_type == 'substring':
            # Substring matching
            matched_values = translate_substring_filter(available, search_terms, case_insensitive=True)
            print(f"  Search patterns: {search_terms}")
        else:
            # Exact matching
            matched_values = [v for v in search_terms if v in available]
            print(f"  Exact values: {search_terms}")

        if matched_values:
            print(f"\n  Found {len(matched_values)} matching value(s):")
            for val in sorted(set(matched_values))[:15]:
                print(f"    - {val}")
            if len(set(matched_values)) > 15:
                print(f"    ... and {len(set(matched_values)) - 15} more")

            # Handle exclude vs include for disease
            if arg_name == 'exclude_diseases':
                # Remove excluded from included
                if 'disease' in filter_config:
                    existing = set(filter_config['disease']['values'])
                    filter_config['disease']['values'] = list(existing - set(matched_values))
                    print(f"\n  Excluding {len(matched_values)} disease(s) from filter")
                else:
                    print(f"\n  No include disease filter to exclude from - creating negated filter")
                    filter_config[f'{col_name}_exclude'] = {'values': matched_values, 'negate': True}
            else:
                filter_config[col_name] = {'values': matched_values, 'negate': is_negated}

            validation_results['filters_applied'].append({
                'step': step_num,
                'filter': arg_name,
                'column': col_name,
                'query': filter_value,
                'matches': len(set(matched_values)),
                'status': 'success'
            })
        else:
            print(f"\n  No matching values found!")
            print(f"  Available {col_name} values (first 20):")
            for val in sorted(available)[:20]:
                print(f"    - {val}")

            validation_results['success'] = False
            validation_results['failed_filter'] = arg_name
            validation_results['suggestions'].append(
                f"No {col_name} values match '{filter_value}'. "
                f"Check available values in schema explorer."
            )
            return validation_results

        step_num += 1

    # Build and test the final filter
    print(f"\n[STEP {step_num}] Testing complete filter...")

    soma_filter = build_soma_filter(filter_config)

    if soma_filter:
        print(f"  SOMA filter: {soma_filter[:200]}..." if len(soma_filter) > 200 else f"  SOMA filter: {soma_filter}")

        # Count matching observations
        try:
            with tiledbsoma.open(args.file, mode="r", context=context) as exp:
                query = exp.axis_query(
                    measurement_name="RNA",
                    obs_query=tiledbsoma.AxisQuery(value_filter=soma_filter)
                )
                # Get obs to count
                obs_df = query.obs().concat().to_pandas()
                final_obs = len(obs_df)
                query.close()

            print(f"\n  Final observation count: {final_obs:,}")
            validation_results['final_obs'] = final_obs

            if final_obs == 0:
                validation_results['success'] = False
                validation_results['failed_filter'] = 'combined_filter'
                validation_results['suggestions'].append(
                    "Combined filters resulted in 0 observations. "
                    "Try relaxing some filter constraints."
                )
        except Exception as e:
            print(f"  Error testing filter: {e}")
            validation_results['success'] = False
            validation_results['suggestions'].append(f"Filter test failed: {e}")
    else:
        print("  No filters to apply - all data will be loaded")
        validation_results['final_obs'] = validation_results['total_obs']

    # Final summary
    if validation_results['success']:
        print(f"\n" + "=" * 80)
        print(f"FILTER VALIDATION SUMMARY - SOMA")
        print(f"=" * 80)
        print(f"\n All filters validated successfully!")
        print(f"\n  Initial observations: {validation_results['total_obs']:,}")
        print(f"  Final observations: {validation_results['final_obs']:,}")

        if validation_results['total_obs'] > 0:
            reduction = (validation_results['total_obs'] - validation_results['final_obs']) / validation_results['total_obs'] * 100
            print(f"  Reduction: {reduction:.1f}%")

        print(f"\n  Filters applied: {len(validation_results['filters_applied'])}")
        for filter_info in validation_results['filters_applied']:
            print(f"    {filter_info['step']}. {filter_info['filter']}: {filter_info['matches']} matches")

        print(f"\n Ready to proceed with SOMA analysis")

    return validation_results


def validate_filters_stepwise(args):
    """
    Validate filters step-by-step, showing intermediate results.

    Auto-detects SOMA vs h5ad and dispatches to appropriate validator.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments with filter parameters

    Returns
    -------
    dict
        Validation results with filter status and suggestions
    """
    # Check if this is a SOMA URI
    if SOMA_AVAILABLE and is_s3_uri(args.file) and is_soma_uri(args.file):
        print("Detected SOMA URI - using SOMA validation")
        return validate_filters_soma(args)

    # Standard h5ad validation
    print("=" * 80)
    print("FILTER VALIDATION - Step-by-Step Query")
    print("=" * 80)

    # Load data
    print(f"\n[STEP 0] Loading h5ad file...")
    print(f"  File: {args.file}")
    adata = load_h5ad(args.file)
    print(f"  ✓ Loaded: {adata.shape} (obs x vars)")
    print(f"  Total observations: {adata.n_obs:,}")
    print(f"  Total variables: {adata.n_vars:,}")

    # Initialize filter tracking
    current_obs = adata.n_obs
    base_filter = pd.Series([True] * adata.n_obs, index=adata.obs.index)
    validation_results = {
        'total_obs': adata.n_obs,
        'filters_applied': [],
        'final_obs': None,
        'success': True,
        'failed_filter': None,
        'suggestions': []
    }

    # Step 1: Disease filter
    if args.diseases:
        print(f"\n[STEP 1] Disease Filter")
        print(f"  Query: {args.diseases}")

        disease_terms = [d.strip() for d in args.diseases.split(',')]
        print(f"  Filter type: Substring match (case-insensitive)")
        print(f"  Search terms: {disease_terms}")

        # Show available diseases
        unique_diseases = adata.obs.disease.unique()
        print(f"\n  Available diseases in dataset: {len(unique_diseases)}")

        # Find matches
        matching_diseases = []
        for disease in unique_diseases:
            for term in disease_terms:
                if term.lower() in disease.lower():
                    matching_diseases.append(disease)
                    break

        if matching_diseases:
            print(f"\n  ✓ Found {len(matching_diseases)} matching disease(s):")
            for disease in sorted(set(matching_diseases)):
                count = (adata.obs.disease == disease).sum()
                print(f"    - {disease}: {count:,} obs")

            disease_pattern = '|'.join([d.strip() for d in args.diseases.split(',')])
            disease_filter = adata.obs.disease.str.contains(disease_pattern, case=False, na=False)
            base_filter = base_filter & disease_filter
            current_obs = base_filter.sum()

            print(f"\n  → Observations after disease filter: {current_obs:,}")
            validation_results['filters_applied'].append({
                'step': 1,
                'filter': 'disease',
                'query': args.diseases,
                'matches': len(set(matching_diseases)),
                'obs_remaining': current_obs,
                'status': 'success'
            })
        else:
            print(f"\n  ✗ No matching diseases found!")
            print(f"\n  Available disease keywords (showing first 20):")
            for i, disease in enumerate(sorted(unique_diseases)[:20]):
                print(f"    - {disease}")

            validation_results['success'] = False
            validation_results['failed_filter'] = 'disease'
            validation_results['suggestions'].append(
                f"No diseases match '{args.diseases}'. Please check available diseases in schema viewer or try broader terms."
            )
            validation_results['final_obs'] = 0
            return validation_results

    # Step 1.5: Disease exclusion filter
    if args.exclude_diseases:
        print(f"\n[STEP 1.5] Disease Exclusion Filter")
        print(f"  Exclude query: {args.exclude_diseases}")

        exclude_terms = [d.strip() for d in args.exclude_diseases.split(',')]
        print(f"  Filter type: Substring match (case-insensitive)")
        print(f"  Exclusion terms: {exclude_terms}")

        # Show diseases that will be excluded
        filtered_adata = adata[base_filter]
        unique_diseases = filtered_adata.obs.disease.unique()

        # Find matches to exclude
        excluded_diseases = []
        for disease in unique_diseases:
            for term in exclude_terms:
                if term.lower() in disease.lower():
                    excluded_diseases.append(disease)
                    break

        if excluded_diseases:
            print(f"\n  ✓ Found {len(excluded_diseases)} disease(s) to exclude:")
            for disease in sorted(set(excluded_diseases)):
                count = (filtered_adata.obs.disease == disease).sum()
                print(f"    - {disease}: {count:,} obs")

            # Apply exclusion filter
            exclude_pattern = '|'.join([d.strip() for d in args.exclude_diseases.split(',')])
            exclude_filter = ~adata.obs.disease.str.contains(exclude_pattern, case=False, na=False)
            base_filter = base_filter & exclude_filter
            current_obs = base_filter.sum()

            print(f"\n  → Observations after disease exclusion: {current_obs:,}")
            validation_results['filters_applied'].append({
                'step': 1.5,
                'filter': 'disease_exclusion',
                'query': args.exclude_diseases,
                'matches': len(set(excluded_diseases)),
                'obs_remaining': current_obs,
                'status': 'success'
            })
        else:
            print(f"\n  ℹ No diseases matched exclusion criteria (this is normal)")

    # Step 2: Tissue filter
    if args.tissues:
        print(f"\n[STEP 2] Tissue Filter")
        print(f"  Query: {args.tissues}")

        tissue_terms = [t.strip() for t in args.tissues.split(',')]
        print(f"  Filter type: Substring match (case-insensitive)")
        print(f"  Search terms: {tissue_terms}")

        # Show available tissues in current filtered data
        filtered_adata = adata[base_filter]
        unique_tissues = filtered_adata.obs.tissue.unique()
        print(f"\n  Available tissues in current dataset: {len(unique_tissues)}")

        # Find matches
        matching_tissues = []
        for tissue in unique_tissues:
            for term in tissue_terms:
                if term.lower() in tissue.lower():
                    matching_tissues.append(tissue)
                    break

        if matching_tissues:
            print(f"\n  ✓ Found {len(matching_tissues)} matching tissue(s):")
            for tissue in sorted(set(matching_tissues)):
                count = (filtered_adata.obs.tissue == tissue).sum()
                print(f"    - {tissue}: {count:,} obs")

            tissue_pattern = '|'.join([t.strip() for t in args.tissues.split(',')])
            tissue_filter = adata.obs.tissue.str.contains(tissue_pattern, case=False, na=False)
            base_filter = base_filter & tissue_filter
            current_obs = base_filter.sum()

            print(f"\n  → Observations after tissue filter: {current_obs:,}")
            validation_results['filters_applied'].append({
                'step': 2,
                'filter': 'tissue',
                'query': args.tissues,
                'matches': len(set(matching_tissues)),
                'obs_remaining': current_obs,
                'status': 'success'
            })

            if current_obs == 0:
                print(f"\n  ✗ Filter resulted in 0 observations!")
                print(f"  Tissue filter may be too restrictive for selected diseases.")
                validation_results['success'] = False
                validation_results['failed_filter'] = 'tissue'
                validation_results['suggestions'].append(
                    f"Tissue filter '{args.tissues}' combined with disease filter resulted in 0 observations. "
                    f"Try broader tissue terms or check schema viewer for tissue-disease combinations."
                )
                validation_results['final_obs'] = 0
                return validation_results
        else:
            print(f"\n  ✗ No matching tissues found in current filtered dataset!")
            print(f"\n  Available tissue keywords (showing first 20):")
            for i, tissue in enumerate(sorted(unique_tissues)[:20]):
                print(f"    - {tissue}")

            validation_results['success'] = False
            validation_results['failed_filter'] = 'tissue'
            validation_results['suggestions'].append(
                f"No tissues match '{args.tissues}' in the disease-filtered dataset. "
                f"Please check available tissues in schema viewer."
            )
            validation_results['final_obs'] = 0
            return validation_results

    # Step 3: Studies filter
    if args.studies:
        print(f"\n[STEP 3] Studies Filter")
        print(f"  Query: {args.studies}")

        study_list = [s.strip() for s in args.studies.split(',')]
        print(f"  Filter type: Exact match")
        print(f"  Search terms: {study_list}")

        filtered_adata = adata[base_filter]
        unique_studies = filtered_adata.obs.study.unique()
        print(f"\n  Available studies in current dataset: {len(unique_studies)}")

        # Find exact matches
        matching_studies = [s for s in study_list if s in unique_studies]

        if matching_studies:
            print(f"\n  ✓ Found {len(matching_studies)} matching study(ies):")
            for study in matching_studies:
                count = (filtered_adata.obs.study == study).sum()
                print(f"    - {study}: {count:,} obs")

            study_filter = adata.obs.study.isin(matching_studies)
            base_filter = base_filter & study_filter
            current_obs = base_filter.sum()

            print(f"\n  → Observations after studies filter: {current_obs:,}")
            validation_results['filters_applied'].append({
                'step': 3,
                'filter': 'studies',
                'query': args.studies,
                'matches': len(matching_studies),
                'obs_remaining': current_obs,
                'status': 'success'
            })

            if current_obs == 0:
                print(f"\n  ✗ Filter resulted in 0 observations!")
                validation_results['success'] = False
                validation_results['failed_filter'] = 'studies'
                validation_results['suggestions'].append(
                    f"Study filter resulted in 0 observations. Check if studies {study_list} exist in the filtered dataset."
                )
                validation_results['final_obs'] = 0
                return validation_results
        else:
            print(f"\n  ✗ No matching studies found!")
            print(f"\n  Available studies (showing first 20):")
            for i, study in enumerate(sorted(unique_studies)[:20]):
                print(f"    - {study}")

            validation_results['success'] = False
            validation_results['failed_filter'] = 'studies'
            validation_results['suggestions'].append(
                f"No studies match {study_list}. Please check available studies in schema viewer."
            )
            validation_results['final_obs'] = 0
            return validation_results

    # Step 4: Comparison category filter
    if args.comparison_category:
        print(f"\n[STEP 4] Comparison Category Filter")
        print(f"  Query: {args.comparison_category}")

        category_list = [c.strip() for c in args.comparison_category.split(',')]
        print(f"  Filter type: Exact match")
        print(f"  Search terms: {category_list}")

        filtered_adata = adata[base_filter]
        unique_categories = filtered_adata.obs.comparison_category.unique()
        print(f"\n  Available comparison categories in current dataset: {len(unique_categories)}")
        for cat in sorted(unique_categories):
            count = (filtered_adata.obs.comparison_category == cat).sum()
            print(f"    - {cat}: {count:,} obs")

        # Find exact matches
        matching_categories = [c for c in category_list if c in unique_categories]

        if matching_categories:
            print(f"\n  ✓ Found {len(matching_categories)} matching category(ies):")
            for cat in matching_categories:
                count = (filtered_adata.obs.comparison_category == cat).sum()
                print(f"    - {cat}: {count:,} obs")

            category_filter = adata.obs.comparison_category.isin(matching_categories)
            base_filter = base_filter & category_filter
            current_obs = base_filter.sum()

            print(f"\n  → Observations after comparison category filter: {current_obs:,}")
            validation_results['filters_applied'].append({
                'step': 4,
                'filter': 'comparison_category',
                'query': args.comparison_category,
                'matches': len(matching_categories),
                'obs_remaining': current_obs,
                'status': 'success'
            })

            if current_obs == 0:
                print(f"\n  ✗ Filter resulted in 0 observations!")
                validation_results['success'] = False
                validation_results['failed_filter'] = 'comparison_category'
                validation_results['suggestions'].append(
                    f"Comparison category filter resulted in 0 observations. "
                    f"Not all studies have all comparison categories."
                )
                validation_results['final_obs'] = 0
                return validation_results
        else:
            print(f"\n  ✗ No matching comparison categories found!")
            print(f"  Note: Exact match required. Available categories:")
            for cat in sorted(unique_categories):
                print(f"    - {cat}")

            validation_results['success'] = False
            validation_results['failed_filter'] = 'comparison_category'
            validation_results['suggestions'].append(
                f"No comparison categories match {category_list}. Use exact category names from schema."
            )
            validation_results['final_obs'] = 0
            return validation_results

    # Step 5: Case treatment filter
    if args.case_treatment:
        print(f"\n[STEP 5] Case Treatment Filter")
        print(f"  Query: {args.case_treatment}")

        treatment_list = [t.strip() for t in args.case_treatment.split(',')]
        print(f"  Filter type: Exact match")
        print(f"  Search terms: {treatment_list}")

        filtered_adata = adata[base_filter]
        unique_treatments = filtered_adata.obs.case_treatment.unique()
        print(f"\n  Available case treatments in current dataset: {len(unique_treatments)}")

        matching_treatments = [t for t in treatment_list if t in unique_treatments]

        if matching_treatments:
            print(f"\n  ✓ Found {len(matching_treatments)} matching treatment(s):")
            for treatment in matching_treatments:
                count = (filtered_adata.obs.case_treatment == treatment).sum()
                print(f"    - {treatment}: {count:,} obs")

            treatment_filter = adata.obs.case_treatment.isin(matching_treatments)
            base_filter = base_filter & treatment_filter
            current_obs = base_filter.sum()

            print(f"\n  → Observations after case treatment filter: {current_obs:,}")
            validation_results['filters_applied'].append({
                'step': 5,
                'filter': 'case_treatment',
                'query': args.case_treatment,
                'matches': len(matching_treatments),
                'obs_remaining': current_obs,
                'status': 'success'
            })

            if current_obs == 0:
                print(f"\n  ✗ Filter resulted in 0 observations!")
                validation_results['success'] = False
                validation_results['failed_filter'] = 'case_treatment'
                validation_results['final_obs'] = 0
                return validation_results
        else:
            print(f"\n  ✗ No matching case treatments found!")
            validation_results['success'] = False
            validation_results['failed_filter'] = 'case_treatment'
            validation_results['final_obs'] = 0
            return validation_results

    # Step 6: Comparison filter
    if args.comparison:
        print(f"\n[STEP 6] Comparison Filter")
        print(f"  Query: {args.comparison}")

        comparison_terms = [c.strip() for c in args.comparison.split(',')]
        print(f"  Filter type: Substring match (case-insensitive, fuzzy)")
        print(f"  Search terms: {comparison_terms}")

        filtered_adata = adata[base_filter]
        unique_comparisons = filtered_adata.obs.comparison.unique()
        print(f"\n  Available comparisons in current dataset: {len(unique_comparisons)}")

        matching_comparisons = []
        for comparison in unique_comparisons:
            for term in comparison_terms:
                if term.lower() in comparison.lower():
                    matching_comparisons.append(comparison)
                    break

        if matching_comparisons:
            print(f"\n  ✓ Found {len(matching_comparisons)} matching comparison(s):")
            for comp in sorted(set(matching_comparisons))[:10]:
                count = (filtered_adata.obs.comparison == comp).sum()
                print(f"    - {comp}: {count:,} obs")
            if len(set(matching_comparisons)) > 10:
                print(f"    ... and {len(set(matching_comparisons)) - 10} more")

            comparison_pattern = '|'.join(comparison_terms)
            comparison_filter = adata.obs.comparison.str.contains(comparison_pattern, case=False, na=False)
            base_filter = base_filter & comparison_filter
            current_obs = base_filter.sum()

            print(f"\n  → Observations after comparison filter: {current_obs:,}")
            validation_results['filters_applied'].append({
                'step': 6,
                'filter': 'comparison',
                'query': args.comparison,
                'matches': len(set(matching_comparisons)),
                'obs_remaining': current_obs,
                'status': 'success'
            })

            if current_obs == 0:
                print(f"\n  ✗ Filter resulted in 0 observations!")
                validation_results['success'] = False
                validation_results['failed_filter'] = 'comparison'
                validation_results['final_obs'] = 0
                return validation_results
        else:
            print(f"\n  ✗ No matching comparisons found!")
            validation_results['success'] = False
            validation_results['failed_filter'] = 'comparison'
            validation_results['final_obs'] = 0
            return validation_results

    # Final summary
    validation_results['final_obs'] = current_obs

    print(f"\n" + "=" * 80)
    print(f"FILTER VALIDATION SUMMARY")
    print(f"=" * 80)
    print(f"\n✓ All filters validated successfully!")
    print(f"\n  Initial observations: {adata.n_obs:,}")
    print(f"  Final observations: {current_obs:,}")
    print(f"  Reduction: {((adata.n_obs - current_obs) / adata.n_obs * 100):.1f}%")

    print(f"\n  Filters applied: {len(validation_results['filters_applied'])}")
    for filter_info in validation_results['filters_applied']:
        print(f"    {filter_info['step']}. {filter_info['filter']}: {filter_info['obs_remaining']:,} obs")

    print(f"\n✓ Ready to proceed with full analysis")
    print(f"  Suggested command:")
    cmd_parts = [
        "conda run -n <env_name> python scripts/deg_analysis.py",
        f"--file {args.file}",
        f"--target-name {args.target_name}",
        f"--signatures {args.signatures}"
    ]
    if args.diseases:
        cmd_parts.append(f'--diseases "{args.diseases}"')
    if args.exclude_diseases:
        cmd_parts.append(f'--exclude-diseases "{args.exclude_diseases}"')
    if args.tissues:
        cmd_parts.append(f'--tissues "{args.tissues}"')
    if args.studies:
        cmd_parts.append(f'--studies "{args.studies}"')
    if args.comparison_category:
        cmd_parts.append(f'--comparison-category "{args.comparison_category}"')
    if args.case_treatment:
        cmd_parts.append(f'--case-treatment "{args.case_treatment}"')
    if args.comparison:
        cmd_parts.append(f'--comparison "{args.comparison}"')
    if args.targets:
        cmd_parts.append(f'--targets "{args.targets}"')
    cmd_parts.append(f"--lfc-threshold {args.lfc_threshold}")
    cmd_parts.append(f"--padj-threshold {args.padj_threshold}")
    if not args.run_gsea:
        cmd_parts.append("--no-gsea")
    # GSEA is ON by default, so no need to add --run-gsea

    print(f"\n{' '.join(cmd_parts)}")

    return validation_results


def main():
    parser = argparse.ArgumentParser(description='Validate filters before running full DEG analysis')
    parser.add_argument('--file', required=True, help='Path to h5ad file')
    parser.add_argument('--target-name', required=True, help='Name of target for analysis')
    parser.add_argument('--signatures', required=True, help='Signatures in format "Name:Gene1,Gene2"')

    # Target genes for plotting and leading edge tracking
    parser.add_argument('--targets', default=None, help='Optional: Target genes for leading edge tracking')

    # Optional filtering parameters
    parser.add_argument('--diseases', default=None, help='Disease keywords (substring match)')
    parser.add_argument('--exclude-diseases', default=None, help='Disease keywords to exclude (substring match)')
    parser.add_argument('--studies', default=None, help='Study names (exact match)')
    parser.add_argument('--tissues', default=None, help='Tissue keywords (substring match)')
    parser.add_argument('--comparison-category', default=None, help='Comparison categories (exact match)')
    parser.add_argument('--case-treatment', default=None, help='Case treatment (exact match)')
    parser.add_argument('--comparison', default=None, help='Comparison keywords (fuzzy match)')

    # Thresholds
    parser.add_argument('--lfc-threshold', type=float, default=0.0, help='Log2FC threshold')
    parser.add_argument('--padj-threshold', type=float, default=0.05, help='Adjusted p-value threshold')

    # Analysis options
    parser.add_argument('--run-gsea', action='store_true', default=True, help='Run GSEA analysis (ON by default)')
    parser.add_argument('--no-gsea', action='store_true', help='Disable GSEA analysis')

    args = parser.parse_args()

    # Handle GSEA flag logic
    if args.no_gsea:
        args.run_gsea = False

    validation_results = validate_filters_stepwise(args)

    if not validation_results['success']:
        print(f"\n" + "=" * 80)
        print(f"VALIDATION FAILED")
        print(f"=" * 80)
        print(f"\n✗ Failed at filter: {validation_results['failed_filter']}")
        print(f"\nSuggestions:")
        for suggestion in validation_results['suggestions']:
            print(f"  - {suggestion}")
        print(f"\nPlease run schema exploration to find valid filter values:")
        print(f"  conda run -n <env> python scripts/explore_h5ad_schema.py --file {args.file} --format json --output schema.json")
        print(f"  conda run -n <env> python scripts/generate_schema_viewer.py --json schema.json --output schema_viewer.html")
        sys.exit(1)
    else:
        print(f"\n✓ Validation complete - ready to run full analysis")
        sys.exit(0)


if __name__ == "__main__":
    main()
