#!/usr/bin/env python3
"""
Convert manually downloaded Cortellis Excel files to JSON format matching API schema.
This is a fallback when API calls timeout or don't provide comprehensive records.

Usage:
    python convert_excel_to_json.py <excel_file> [gene_name] [--output output.json]

Example:
    python convert_excel_to_json.py Drugs___Biologics_Dec_19_2025_tyk2.xlsx TYK2
"""

import pandas as pd
import json
import sys
import os
import re
from pathlib import Path
from datetime import datetime

def parse_phase_from_excel(phase_str):
    """
    Parse phase from Excel 'Highest Phase' field to match API format.
    API format: 'Phase 1 Clinical', 'Phase 2 Clinical', 'Phase 3 Clinical', 'Launched', 'Preclinical', 'Discontinued'
    """
    if not phase_str or pd.isna(phase_str):
        return {'@id': 'Unknown', '$': 'Unknown'}

    phase_lower = str(phase_str).lower()

    # Map to API format
    if 'launched' in phase_lower:
        # Extract year if present: "Launched - 2022" → "Launched"
        return {'@id': 'LA', '$': 'Launched'}
    elif 'marketed' in phase_lower or 'approved' in phase_lower:
        return {'@id': 'LA', '$': 'Launched'}
    elif 'phase iii' in phase_lower or 'phase 3' in phase_lower:
        return {'@id': 'C3', '$': 'Phase 3 Clinical'}
    elif 'phase ii' in phase_lower or 'phase 2' in phase_lower:
        return {'@id': 'C2', '$': 'Phase 2 Clinical'}
    elif 'phase i' in phase_lower or 'phase 1' in phase_lower:
        return {'@id': 'C1', '$': 'Phase 1 Clinical'}
    elif 'preclinical' in phase_lower:
        return {'@id': 'PC', '$': 'Preclinical'}
    elif 'discontinued' in phase_lower or 'suspended' in phase_lower:
        return {'@id': 'DI', '$': 'Discontinued'}
    elif 'no development' in phase_lower:
        return {'@id': 'ND', '$': 'No Development Reported'}
    else:
        return {'@id': phase_str[:10], '$': phase_str}

def split_multiline_field(field_value, delimiter='\n'):
    """Split a field by delimiter and return list of non-empty strings."""
    if pd.isna(field_value) or not field_value:
        return []

    values = str(field_value).split(delimiter)
    return [v.strip() for v in values if v.strip()]

def extract_gene_name_from_filename(filename):
    """
    Extract gene name from filename like 'Drugs___Biologics_Dec_19_2025_tyk2.xlsx'
    Returns 'TYK2' (uppercase)
    """
    basename = os.path.basename(filename)
    name_without_ext = os.path.splitext(basename)[0]

    # Try to extract gene name from various patterns
    # Pattern 1: ends with _GENENAME
    parts = name_without_ext.split('_')
    if len(parts) > 0:
        # Get last non-empty part
        for part in reversed(parts):
            if part and not part.isdigit() and len(part) >= 2:
                return part.upper()

    return 'UNKNOWN'

def clean_drug_name(name):
    """Clean drug name by removing annotations in parentheses."""
    if pd.isna(name) or not name:
        return ''

    # Remove text in parentheses like "(Rec INN; USAN)"
    cleaned = re.sub(r'\s*\([^)]*\)', '', str(name)).strip()
    return cleaned

def convert_excel_to_json(excel_file, output_json=None, gene_name=None):
    """
    Convert Cortellis Excel file to JSON format matching API schema.

    Args:
        excel_file: Path to downloaded Cortellis Excel file
        output_json: Optional output JSON filename
        gene_name: Optional gene name (extracted from filename if not provided)

    Returns:
        Path to generated JSON file
    """

    if not os.path.exists(excel_file):
        print(f"ERROR: File not found: {excel_file}")
        return None

    print(f"Converting: {excel_file}")

    # Extract gene name from filename if not provided
    if not gene_name:
        gene_name = extract_gene_name_from_filename(excel_file)
        print(f"Detected gene: {gene_name}")

    # Set default output filename
    if not output_json:
        output_json = f"{gene_name}_cortellis_data.json"

    # Read all sheets
    try:
        excel_data = pd.ExcelFile(excel_file)
        sheets = excel_data.sheet_names
        print(f"Found sheets: {', '.join(sheets)}")
    except Exception as e:
        print(f"ERROR reading Excel file: {e}")
        return None

    # Initialize JSON structure matching API format
    json_data = {
        'annotation': {
            '@Id': '',
            '@namemain': gene_name,
            'Symbol': gene_name,
            'GeneId': '',
            'UniprotId': '',
            'TargetType': 'Protein',
            'Description': f'Manually downloaded data for {gene_name}',
            'Organism': {'$': 'Homo sapiens'},
            '_source': 'Manual Excel Download',
            '_conversion_date': datetime.now().isoformat()
        },
        'drug': {
            'Target': {
                '@namemain': gene_name,
                '@id': '',
            },
            'Drug': [],  # Basic drug list from Targets API
            'Trial': [],  # Clinical trials
            'DrugRecord': {}  # Comprehensive drug records (main data)
        }
    }

    # Read Product List sheet (main data source)
    if 'Product List' in sheets:
        df_products = pd.read_excel(excel_file, sheet_name='Product List')
        print(f"Processing {len(df_products)} drugs from Product List...")

        # Process each drug
        for idx, row in df_products.iterrows():
            # Get Entry Number as drug ID
            drug_id = str(row.get('Entry Number', idx))

            # Get drug name (prefer Generic Name, fallback to Code Name)
            generic_name = row.get('Generic Name', '')
            code_name = row.get('Code Name', '')
            brand_name = row.get('Brand Name', '')

            drug_name = clean_drug_name(generic_name)
            if not drug_name:
                drug_name = clean_drug_name(code_name)
            if not drug_name:
                drug_name = clean_drug_name(brand_name)
            if not drug_name:
                drug_name = f"Drug_{drug_id}"

            # Parse phase
            phase_highest = parse_phase_from_excel(row.get('Highest Phase', ''))

            # Get all drug names/synonyms
            all_names = split_multiline_field(row.get('Drug Name (All)', ''))

            # Build synonym structure
            drug_synonyms = {'Name': []}
            for name in all_names[:10]:  # Limit to 10
                drug_synonyms['Name'].append({'Value': name})

            # Parse conditions/indications
            conditions_str = row.get('Condition', '')
            indications_list = split_multiline_field(conditions_str)

            indications_primary = {'Indication': []}
            for ind_idx, indication in enumerate(indications_list[:10]):  # Limit to 10
                indications_primary['Indication'].append({
                    '@id': str(ind_idx),
                    '$': indication
                })

            # If only one indication, unwrap from list
            if len(indications_primary['Indication']) == 1:
                indications_primary['Indication'] = indications_primary['Indication'][0]

            # Parse organizations/companies
            org_str = row.get('Organization', '')
            companies = split_multiline_field(org_str)

            company_originator = {}
            companies_primary = {}
            if companies:
                company_originator = {
                    '@id': '0',
                    '$': companies[0]
                }

                # If multiple companies, add to primary
                if len(companies) > 1:
                    companies_primary = {
                        'Company': [{'@id': str(i), '$': comp} for i, comp in enumerate(companies[:5])]
                    }
                    if len(companies_primary['Company']) == 1:
                        companies_primary['Company'] = companies_primary['Company'][0]

            # Parse mechanism of action
            moa_str = row.get('Mechanism of Action', '')
            moa_list = split_multiline_field(moa_str)

            actions_primary = {'Action': []}
            for action_idx, action in enumerate(moa_list[:5]):  # Limit to 5
                actions_primary['Action'].append({
                    '@id': str(action_idx),
                    '$': action
                })

            if len(actions_primary['Action']) == 1:
                actions_primary['Action'] = actions_primary['Action'][0]

            # Parse therapeutic groups as therapy areas
            therapy_group_str = row.get('Therapeutic Group', '')
            therapy_areas_list = split_multiline_field(therapy_group_str)

            therapy_areas = {}
            if therapy_areas_list:
                if len(therapy_areas_list) == 1:
                    therapy_areas = {'TherapyArea': therapy_areas_list[0]}
                else:
                    therapy_areas = {'TherapyArea': therapy_areas_list[:5]}  # Limit to 5

            # Parse drug type as technology
            drug_type = row.get('Drug Type', '')
            technologies = {}
            if drug_type and not pd.isna(drug_type):
                technologies = {
                    'Technology': [{
                        '@id': '0',
                        '$': drug_type
                    }]
                }

            # Get chemical structure (SMILES)
            smiles = row.get('Smiles', '')
            if pd.isna(smiles) or smiles == 'nan':
                smiles = ''
            else:
                smiles = str(smiles)

            # Build comprehensive drug record matching API schema
            drug_record = {
                '@id': drug_id,
                'DrugName': drug_name,
                'DrugNamesKey': {
                    'Name': {
                        '@id': drug_id,
                        '$': drug_name
                    }
                },
                'DrugSynonyms': drug_synonyms if drug_synonyms['Name'] else {},
                'PhaseHighest': phase_highest,
                'CompanyOriginator': company_originator,
                'CompaniesPrimary': companies_primary,
                'IndicationsPrimary': indications_primary if indications_primary['Indication'] else {},
                'ActionsPrimary': actions_primary if actions_primary['Action'] else {},
                'TherapyAreas': therapy_areas,
                'Technologies': technologies,
                'StructureSmiles': smiles,
                '_source_drug_name': drug_name,
                '_source_display_name': drug_name,
                '_excel_row': idx
            }

            # Add optional fields
            if 'Molecular Formula' in row and not pd.isna(row['Molecular Formula']):
                drug_record['MolecularFormula'] = str(row['Molecular Formula'])

            if 'Molecular Weight' in row and not pd.isna(row['Molecular Weight']):
                drug_record['MolecularWeight'] = str(row['Molecular Weight'])

            if 'CAS Registry Number' in row and not pd.isna(row['CAS Registry Number']):
                drug_record['CASNumber'] = str(row['CAS Registry Number'])

            if 'Product Summary' in row and not pd.isna(row['Product Summary']):
                # Add to development profile
                drug_record['DevelopmentProfile'] = {
                    'Summary': {
                        'displayLabel': 'Summary',
                        'value': f'<Summary><para>{row["Product Summary"]}</para></Summary>'
                    }
                }

            # Add to DrugRecord dict
            json_data['drug']['DrugRecord'][drug_id] = drug_record

            # Also add basic drug entry to Drug list (simpler format)
            chem_desc = row.get('Chemical Name/Description', '')
            if pd.isna(chem_desc):
                chem_desc = ''

            basic_drug = {
                '@id': drug_id,
                '@namemain': drug_name,
                'NamesChemicalAndDescriptions': {
                    'Name': str(chem_desc) if chem_desc else ''
                },
                'NamesCode': {
                    'Name': code_name if code_name else drug_name
                },
                'MechanismsMolecular': {
                    'Mechanism': {
                        '@id': '0',
                        '$': moa_list[0] if moa_list else ''
                    }
                } if moa_list else {}
            }
            json_data['drug']['Drug'].append(basic_drug)

    # Read Development Status sheet (if available)
    if 'Development Status' in sheets:
        df_dev = pd.read_excel(excel_file, sheet_name='Development Status')
        print(f"Processing {len(df_dev)} development status entries...")

        # Group by Entry Number to add regional development info
        for entry_num, group in df_dev.groupby('Entry Number'):
            drug_id = str(entry_num)
            if drug_id in json_data['drug']['DrugRecord']:
                # Add regional development status
                dev_statuses = []
                for _, dev_row in group.iterrows():
                    dev_status = {
                        'Country': dev_row.get('Country/Region', ''),
                        'Phase': dev_row.get('Phase', ''),
                        'Organization': dev_row.get('Organization', ''),
                        'Indication': dev_row.get('Indication', ''),
                        'FormulationRoute': f"{dev_row.get('Formulation', '')} ({dev_row.get('Administration Route', '')})"
                    }
                    dev_statuses.append(dev_status)

                json_data['drug']['DrugRecord'][drug_id]['RegionalDevelopment'] = dev_statuses

    # Read Milestones sheet (if available)
    if 'Milestones' in sheets:
        df_milestones = pd.read_excel(excel_file, sheet_name='Milestones')
        print(f"Processing {len(df_milestones)} milestone entries...")

        # Group by Entry Number to add milestone info
        for entry_num, group in df_milestones.groupby('Entry Number'):
            drug_id = str(entry_num)
            if drug_id in json_data['drug']['DrugRecord']:
                milestones = []
                for _, ms_row in group.iterrows():
                    milestone = {
                        'Date': str(ms_row.get('Milestone Date', '')),
                        'Type': ms_row.get('Milestone', ''),
                        'Notes': ms_row.get('Notes', ''),
                        'Organization': ms_row.get('Organization', ''),
                        'Country': ms_row.get('Country/Region', '')
                    }
                    milestones.append(milestone)

                json_data['drug']['DrugRecord'][drug_id]['Milestones'] = milestones

    # Clean NaN values before saving (pandas NaN is not JSON-valid)
    def clean_nan(obj):
        """Recursively replace NaN with None (null in JSON)."""
        if isinstance(obj, dict):
            return {k: clean_nan(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [clean_nan(item) for item in obj]
        elif isinstance(obj, float):
            # Check for NaN
            if pd.isna(obj):
                return None
            return obj
        elif pd.isna(obj):
            return None
        else:
            return obj

    json_data_cleaned = clean_nan(json_data)

    # Save to JSON
    with open(output_json, 'w') as f:
        json.dump(json_data_cleaned, f, indent=2)

    print(f"\n✓ Conversion successful!")
    print(f"  Output: {output_json}")
    print(f"  Comprehensive Drug Records: {len(json_data['drug']['DrugRecord'])}")
    print(f"  Basic Drug Entries: {len(json_data['drug']['Drug'])}")
    print(f"  Trials: {len(json_data['drug']['Trial'])}")
    print()

    return output_json

def main():
    """Main function to handle command line arguments."""

    if len(sys.argv) < 2:
        print("Usage: python convert_excel_to_json.py <excel_file> [gene_name] [--output output.json]")
        print("\nExamples:")
        print("  python convert_excel_to_json.py Drugs___Biologics_Dec_19_2025_tyk2.xlsx")
        print("  python convert_excel_to_json.py cortellis_export.xlsx TYK2")
        print("  python convert_excel_to_json.py file.xlsx TYK2 --output custom_output.json")
        print("\nProcesses multiple files:")
        print("  python convert_excel_to_json.py file1.xlsx file2.xlsx file3.xlsx")
        sys.exit(1)

    # Parse arguments
    excel_files = []
    gene_name = None
    output_json = None

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]

        if arg == '--output' and i + 1 < len(sys.argv):
            output_json = sys.argv[i + 1]
            i += 2
        elif arg.endswith(('.xlsx', '.xls')):
            excel_files.append(arg)
            i += 1
        elif not arg.startswith('--'):
            # Assume it's gene name
            gene_name = arg
            i += 1
        else:
            i += 1

    if not excel_files:
        print("ERROR: No Excel files specified")
        sys.exit(1)

    print("=" * 80)
    print("Cortellis Excel to JSON Converter (API Schema)")
    print("=" * 80)
    print()

    converted_files = []

    for excel_file in excel_files:
        # Convert
        output_file = convert_excel_to_json(
            excel_file,
            output_json=output_json if len(excel_files) == 1 else None,
            gene_name=gene_name
        )
        if output_file:
            converted_files.append(output_file)

    print("=" * 80)
    print(f"Conversion Complete: {len(converted_files)} file(s) processed")
    print("=" * 80)
    print()

    if converted_files:
        print("Generated files:")
        for f in converted_files:
            print(f"  - {f}")
        print()
        print("These files now match the API JSON schema format.")

if __name__ == '__main__':
    main()
