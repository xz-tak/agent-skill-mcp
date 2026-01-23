"""
NCBI Paper Query - Validation Module

Functions for validating and proofreading retrieval results across
multiple rounds, detecting and correcting common data issues.
"""

import re
import logging
from typing import Dict, TYPE_CHECKING

import pandas as pd

from ..config import logger

if TYPE_CHECKING:
    from ..metadata import ImpactFactorLookup


def validate_results(df: pd.DataFrame, round_num: int, if_lookup: 'ImpactFactorLookup') -> Dict:
    """
    Validate and proofread results from a retrieval round.

    Args:
        df: DataFrame with results
        round_num: Current round number
        if_lookup: ImpactFactorLookup instance for verification

    Returns:
        Dictionary with validation results
    """
    validation_results = {
        "round": round_num,
        "total_papers": len(df),
        "issues": [],
        "warnings": [],
        "corrections": [],
        "verified": 0,
        "valid": True
    }

    if len(df) == 0:
        validation_results["warnings"].append("No papers to validate")
        return validation_results

    # Check for duplicate PMIDs
    if 'PMID' in df.columns:
        duplicates = df['PMID'].duplicated().sum()
        if duplicates > 0:
            validation_results["issues"].append(f"Found {duplicates} duplicate PMIDs")
            validation_results["valid"] = False

    # Validate IF values
    if 'Impact_Factor' in df.columns and 'Journal' in df.columns:
        for idx, row in df.iterrows():
            journal = row.get('Journal', '')
            if_val = row.get('Impact_Factor')

            # Skip if no IF
            if if_val is None or (isinstance(if_val, float) and pd.isna(if_val)):
                continue

            # Check for year-like values and AUTO-CORRECT
            if 1900 <= float(if_val) <= 2100:
                # Clear cache for this journal to avoid using cached bad value
                if hasattr(if_lookup.get_impact_factor, 'cache_clear'):
                    if_lookup.get_impact_factor.cache_clear()

                # Try to get correct IF from lookup (skip cache by clearing first)
                corrected_if = if_lookup.get_impact_factor(journal)

                # Check if corrected IF is also year-like (cache might have bad value)
                if corrected_if is not None and 1900 <= corrected_if <= 2100:
                    # Corrected value is also year-like, set to None
                    df.at[idx, 'Impact_Factor'] = None
                    validation_results["corrections"].append(
                        f"PMID {row.get('PMID')}: IF {if_val} (year-like) set to None (lookup also returned year-like value)"
                    )
                elif corrected_if is not None:
                    df.at[idx, 'Impact_Factor'] = corrected_if
                    validation_results["corrections"].append(
                        f"PMID {row.get('PMID')}: IF {if_val} (year-like) corrected to {corrected_if}"
                    )
                else:
                    # Set to None if can't find correct IF
                    df.at[idx, 'Impact_Factor'] = None
                    validation_results["corrections"].append(
                        f"PMID {row.get('PMID')}: IF {if_val} (year-like) set to None (journal not in IF database)"
                    )
                continue  # Skip further validation for this row

            # Verify IF against lookup
            verified_if = if_lookup.get_impact_factor(journal)
            if verified_if is not None and abs(verified_if - float(if_val)) > 0.1:
                validation_results["warnings"].append(
                    f"PMID {row.get('PMID')}: IF mismatch ({if_val} vs {verified_if})"
                )

            validation_results["verified"] += 1

    # Check for required columns
    required_cols = ['PMID', 'Title', 'Journal', 'Publication_Year']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        validation_results["issues"].append(f"Missing columns: {missing_cols}")
        validation_results["valid"] = False

    # Validate omics data
    if 'Omics_Count' in df.columns and 'Omics_Accessions' in df.columns:
        omics_issues = []
        omics_warnings = []

        for idx, row in df.iterrows():
            pmid = row.get('PMID', 'unknown')
            omics_count = row.get('Omics_Count', 0)
            accessions = row.get('Omics_Accessions', '')
            condition_counts = row.get('Omics_Condition_Counts', '')
            sample_count = row.get('Omics_Sample_Count', '')

            # Check for inconsistent omics count
            if omics_count > 0 and not accessions:
                omics_issues.append(f"PMID {pmid}: Omics_Count={omics_count} but no accessions")

            # Validate accession formats
            if accessions:
                valid_patterns = [
                    r'GSE\d+', r'GSM\d+', r'GPL\d+',  # GEO
                    r'PRJNA\d+', r'SRP\d+', r'SRR\d+',  # SRA
                    r'E-[A-Z]+-\d+',  # ArrayExpress
                    r'PXD\d+',  # PRIDE
                    r'MTBLS\d+',  # MetaboLights
                    r'CRA\d+', r'HRA\d+', r'PRJCA\d+',  # NGDC/GSA
                ]
                for acc in accessions.split('; '):
                    acc = acc.strip()
                    if acc and not any(re.match(p, acc) for p in valid_patterns):
                        omics_warnings.append(f"PMID {pmid}: Unrecognized accession format: {acc}")

            # Check for missing condition counts when we have GEO data
            if accessions and 'GSE' in accessions and not condition_counts:
                omics_warnings.append(f"PMID {pmid}: GEO dataset found but no condition counts extracted")

            # Validate sample counts are reasonable
            if sample_count:
                for count_str in sample_count.split('; '):
                    try:
                        count = int(count_str.strip())
                        if count > 10000:
                            omics_warnings.append(f"PMID {pmid}: Unusually high sample count: {count}")
                        elif count < 1:
                            omics_issues.append(f"PMID {pmid}: Invalid sample count: {count}")
                    except ValueError:
                        pass

        validation_results["issues"].extend(omics_issues)
        validation_results["warnings"].extend(omics_warnings)
        if omics_issues:
            validation_results["valid"] = False

        # Summary stats
        papers_with_omics = (df['Omics_Count'] > 0).sum()
        validation_results["papers_with_omics"] = int(papers_with_omics)

    # Log validation results
    logger.info(f"Round {round_num} validation: {validation_results['verified']} papers verified, "
                f"{len(validation_results['issues'])} issues, {len(validation_results['warnings'])} warnings, "
                f"{len(validation_results['corrections'])} corrections")

    return validation_results


def compare_rounds(df1: pd.DataFrame, df2: pd.DataFrame, round1: int, round2: int) -> Dict:
    """
    Compare results between two retrieval rounds.

    Args:
        df1: DataFrame from first round
        df2: DataFrame from second round
        round1: First round number
        round2: Second round number

    Returns:
        Dictionary with comparison results
    """
    comparison = {
        "rounds": (round1, round2),
        "consistent": True,
        "differences": [],
        "common_pmids": 0,
        "unique_to_round1": 0,
        "unique_to_round2": 0
    }

    if 'PMID' not in df1.columns or 'PMID' not in df2.columns:
        comparison["differences"].append("Cannot compare: missing PMID column")
        comparison["consistent"] = False
        return comparison

    pmids1 = set(df1['PMID'].astype(str))
    pmids2 = set(df2['PMID'].astype(str))

    comparison["common_pmids"] = len(pmids1 & pmids2)
    comparison["unique_to_round1"] = len(pmids1 - pmids2)
    comparison["unique_to_round2"] = len(pmids2 - pmids1)

    # Check for significant differences
    if comparison["unique_to_round1"] > len(pmids1) * 0.1:
        comparison["differences"].append(
            f"Round {round1} has {comparison['unique_to_round1']} unique papers (>10%)"
        )
        comparison["consistent"] = False

    if comparison["unique_to_round2"] > len(pmids2) * 0.1:
        comparison["differences"].append(
            f"Round {round2} has {comparison['unique_to_round2']} unique papers (>10%)"
        )
        comparison["consistent"] = False

    logger.info(f"Round comparison {round1} vs {round2}: {comparison['common_pmids']} common, "
                f"{comparison['unique_to_round1']} unique to R{round1}, "
                f"{comparison['unique_to_round2']} unique to R{round2}")

    return comparison
