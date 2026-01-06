# Cortellis Excel to JSON Conversion - Complete Guide

## Overview

The conversion script transforms manually downloaded Cortellis Excel files into JSON format that matches the API schema exactly. This allows you to use Excel downloads as a fallback when API calls timeout or for offline analysis.

## What Was Fixed

### ✅ Issues Resolved

1. **NaN Values in JSON** - Pandas NaN values were being written as `NaN` (not JSON-valid)
   - **Solution:** Added `clean_nan()` function to recursively replace NaN with `null`
   - **Result:** JSON is now 100% valid and can be opened in any JSON viewer

2. **Schema Compatibility** - Original script didn't match API JSON structure
   - **Solution:** Complete rewrite to match API schema fields exactly
   - **Result:** Generated JSON works with existing analysis scripts

3. **Missing Fields** - Many Excel columns weren't being extracted
   - **Solution:** Added extraction for synonyms, mechanisms, therapy areas, etc.
   - **Result:** Comprehensive data extraction with 15+ fields per drug

## Script Location

```
/home/sagemaker-user/.claude/skills/cortellis/cortellis_api/scripts/convert_excel_to_json.py
```

## Usage

### Basic Usage

```bash
# Auto-detect gene name from filename
python convert_excel_to_json.py Drugs___Biologics_Dec_19_2025_tyk2.xlsx

# Specify gene name explicitly
python convert_excel_to_json.py Drugs___Biologics_Dec_19_2025_tyk2.xlsx TYK2

# Custom output filename
python convert_excel_to_json.py file.xlsx TYK2 --output custom_output.json

# Process multiple files
python convert_excel_to_json.py file1.xlsx file2.xlsx file3.xlsx
```

### Example Output

```
================================================================================
Cortellis Excel to JSON Converter (API Schema)
================================================================================

Converting: Drugs___Biologics_Dec_19_2025_tyk2.xlsx
Found sheets: Product List, Development Status, Milestones, ...
Processing 1765 drugs from Product List...
Processing 170 development status entries...
Processing 303 milestone entries...

✓ Conversion successful!
  Output: TYK2_cortellis_data.json
  Comprehensive Drug Records: 1765
  Basic Drug Entries: 1765
  Trials: 0
```

## Generated JSON Structure

### Top-Level Schema

```json
{
  "annotation": {
    "@Id": "",
    "@namemain": "TYK2",
    "Symbol": "TYK2",
    "GeneId": "",
    "UniprotId": "",
    "TargetType": "Protein",
    "Description": "Manually downloaded data for TYK2",
    "Organism": {"$": "Homo sapiens"},
    "_source": "Manual Excel Download",
    "_conversion_date": "2025-12-20T..."
  },
  "drug": {
    "Target": {
      "@namemain": "TYK2",
      "@id": ""
    },
    "Drug": [...],
    "Trial": [],
    "DrugRecord": {
      "896086": {...},
      "10464": {...},
      ...
    }
  }
}
```

### DrugRecord Schema (Per Drug)

```json
{
  "@id": "896086",
  "DrugName": "deucravacitinib",
  "DrugNamesKey": {
    "Name": {
      "@id": "896086",
      "$": "deucravacitinib"
    }
  },
  "DrugSynonyms": {
    "Name": [
      {"Value": "deucravacitinib (Rec INN; USAN) (Generic)"},
      {"Value": "Sotyktu (Brand)"},
      {"Value": "BMS-986165 (Code)"},
      ...
    ]
  },
  "PhaseHighest": {
    "@id": "LA",
    "$": "Launched"
  },
  "CompanyOriginator": {
    "@id": "0",
    "$": "Bristol-Myers Squibb (Originator)"
  },
  "CompaniesPrimary": {...},
  "IndicationsPrimary": {
    "Indication": [
      {"@id": "0", "$": "Alopecia areata"},
      {"@id": "1", "$": "Arthritis, psoriatic"},
      ...
    ]
  },
  "ActionsPrimary": {
    "Action": [
      {"@id": "0", "$": "Non-Receptor Tyrosine-Protein Kinase TYK2 (JH2 Domain) Allosteric Inhibitors"},
      ...
    ]
  },
  "TherapyAreas": {
    "TherapyArea": [
      "Autoimmune Diseases, Treatment of",
      "Dermatologic Drugs",
      ...
    ]
  },
  "Technologies": {
    "Technology": [
      {"@id": "0", "$": "Small Molecules (>350 - 500 Da)"}
    ]
  },
  "StructureSmiles": "[2H]C([2H])([2H])NC(=O)c1c(cc(nn1)NC(=O)C2CC2)Nc3cccc(c3OC)c4ncn(n4)C",
  "MolecularFormula": "C20 H19 D3 N8 O3",
  "MolecularWeight": "425.459",
  "CASNumber": "1609392-27-9",
  "DevelopmentProfile": {
    "Summary": {
      "displayLabel": "Summary",
      "value": "<Summary><para>...</para></Summary>"
    }
  },
  "RegionalDevelopment": [
    {
      "Country": "European Union",
      "Phase": "Launched - 2023",
      "Organization": "Bristol-Myers Squibb",
      "Indication": "Treatment of adults with moderate-to-severe plaque psoriasis",
      "FormulationRoute": "Tablets, 6 mg (Oral)"
    },
    ...
  ],
  "Milestones": [
    {
      "Date": "Dec 10, 2025",
      "Type": "Licensed",
      "Notes": "Licensed to...",
      "Organization": "...",
      "Country": "Worldwide"
    },
    ...
  ],
  "_source_drug_name": "deucravacitinib",
  "_source_display_name": "deucravacitinib",
  "_excel_row": 0
}
```

## Field Mapping (Excel → JSON)

| Excel Column | JSON Field | Notes |
|--------------|------------|-------|
| Entry Number | @id | Drug identifier |
| Generic Name | DrugName | Preferred name (cleaned) |
| Drug Name (All) | DrugSynonyms | All alternative names |
| Highest Phase | PhaseHighest | Mapped to API format |
| Organization | CompanyOriginator, CompaniesPrimary | Split by newline |
| Condition | IndicationsPrimary | Split by newline |
| Mechanism of Action | ActionsPrimary | Split by newline |
| Therapeutic Group | TherapyAreas | Split by newline |
| Drug Type | Technologies | Drug modality |
| Smiles | StructureSmiles | Chemical structure |
| Molecular Formula | MolecularFormula | Chemical formula |
| Molecular Weight | MolecularWeight | MW in Da |
| CAS Registry Number | CASNumber | CAS identifier |
| Product Summary | DevelopmentProfile.Summary | Full narrative |
| Chemical Name/Description | NamesChemicalAndDescriptions | For Drug list |
| Code Name | NamesCode | Research code |

### Phase Mapping

| Excel Value | JSON @id | JSON $ |
|-------------|----------|--------|
| Launched - YYYY | LA | Launched |
| Phase III | C3 | Phase 3 Clinical |
| Phase II | C2 | Phase 2 Clinical |
| Phase I | C1 | Phase 1 Clinical |
| Preclinical | PC | Preclinical |
| Discontinued | DI | Discontinued |
| No Development | ND | No Development Reported |

## Key Features

### ✅ API Schema Compatible

- Matches all core fields from API JSON
- Compatible with existing analysis scripts
- Proper nested structure with '@id' and '$' keys

### ✅ Enhanced Data

Excel provides additional information not in API:
- **Regional Development:** Country-specific phase info from "Development Status" sheet
- **Milestones:** Detailed development timeline from "Milestones" sheet
- **Chemical Properties:** Formula, MW, CAS, SMILES
- **Full Summaries:** Complete product development narratives

### ✅ Robust Handling

- **NaN Cleanup:** All pandas NaN values converted to null
- **Multiline Fields:** Properly splits conditions, companies, mechanisms
- **Nested Structures:** Maintains API's nested dict/list format
- **List Unwrapping:** Single items unwrapped from lists (matches API behavior)

## Testing

### Verify JSON Validity

```bash
python3 -c "import json; data = json.load(open('TYK2_cortellis_data.json')); print(f'Valid JSON with {len(data[\"drug\"][\"DrugRecord\"])} drugs')"
```

Expected output:
```
Valid JSON with 1765 drugs
```

### Compare with API Schema

```bash
python3 -c "
import json
tyk2 = json.load(open('TYK2_cortellis_data.json'))
itga4 = json.load(open('ITGA4_cortellis_data.json'))
print('TYK2 keys:', list(tyk2.keys()))
print('ITGA4 keys:', list(itga4.keys()))
print('Match:', set(tyk2.keys()) == set(itga4.keys()))
"
```

Expected output:
```
TYK2 keys: ['annotation', 'drug']
ITGA4 keys: ['annotation', 'drug']
Match: True
```

## Use Cases

### 1. API Fallback

When API calls timeout or are unavailable:

```bash
# API query
python cortellis_gene_query.py TYK2 --excel
# If this times out...

# Use Excel fallback
python convert_excel_to_json.py Drugs___Biologics_Dec_19_2025_tyk2.xlsx
```

### 2. Bulk Processing

Process large datasets offline:

```bash
# Convert multiple Excel files
python convert_excel_to_json.py \
  Drugs_TYK2.xlsx \
  Drugs_JAK1.xlsx \
  Drugs_JAK2.xlsx \
  Drugs_JAK3.xlsx
```

### 3. Historical Analysis

Analyze historical snapshots:

```bash
# Convert historical exports
python convert_excel_to_json.py \
  Dec_2024_TYK2.xlsx \
  --output TYK2_2024_Q4.json
```

## Limitations

The following API fields are NOT available in Excel exports:

- `CrossReferences` - External database links
- `RegulatoryDesignations` - Special regulatory status
- `ActionsSecondary` - Secondary mechanisms
- `EphmraCodes` - Therapeutic classification codes
- `Targets` - Detailed target protein information
- `Deals` - Licensing/partnership details
- `PatentFamilies` - Patent information
- `Trial` - Detailed clinical trial records

## Troubleshooting

### JSON Won't Open

**Error:** "Unexpected token 'N', ... is not valid JSON"

**Solution:** This was the NaN issue - now fixed. Regenerate the JSON:

```bash
rm old_file.json
python convert_excel_to_json.py input.xlsx GENE
```

### Empty Fields

**Issue:** Many fields showing empty or null

**Solution:** Excel export may not include all data. Check:
1. Download complete Excel export from Cortellis
2. Verify "Product List" sheet exists
3. Check if columns like "Condition", "Organization" have data

### Gene Name Not Detected

**Issue:** Output is "UNKNOWN_cortellis_data.json"

**Solution:** Specify gene name explicitly:

```bash
python convert_excel_to_json.py file.xlsx GENENAME
```

## Summary

✅ **Script fixed and tested**
✅ **Valid JSON output (no NaN issues)**
✅ **API schema compatible**
✅ **Enhanced with regional & milestone data**
✅ **1,765 drugs converted successfully for TYK2**

The conversion script now provides a robust fallback for Cortellis data acquisition and enriches the data with regional development and milestone information not available through the API.
