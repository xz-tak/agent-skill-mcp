# Cortellis API Analysis Examples

This directory contains reference implementations demonstrating end-to-end target analysis workflows using Cortellis API v2.0 data.

## Overview

These scripts demonstrate how to:
1. Analyze gene targets and calculate clinical validation scores
2. Generate comprehensive clinical intelligence reports
3. Work with the updated Cortellis API v2.0 JSON schema

## Files

### Analysis Scripts

| File | Description | Input | Output |
|------|-------------|-------|--------|
| **`analyze_targets_example.py`** | Target scoring and ranking analysis | `{GENE}_cortellis_data.json` files | `ibd_analysis_results.json` |
| **`generate_report_example.py`** | Clinical intelligence report generation | `ibd_analysis_results.json` | `IBD_Target_Clinical_Intelligence_Report.md` |

### Supporting Files

| File | Description |
|------|-------------|
| **`EXCEL_TO_JSON_CONVERSION_GUIDE.md`** | Comprehensive guide for Excel to JSON conversion |
| **`IBD_Target_Clinical_Intelligence_Report.md`** | Sample output report |
| **`README.md`** | This file |

## Quick Start

### Step 1: Query Cortellis API

First, query the Cortellis API for your target genes:

```bash
# Query multiple genes
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_api/scripts/cortellis_gene_query.py \
  TYK2 GREM1 ITGA4 ITGB7 --fields drug --excel

# Files created: TYK2_cortellis_data.json, GREM1_cortellis_data.json, etc.
```

**Alternative:** If API timeouts occur, use Excel conversion:
```bash
# Download Excel manually from Cortellis web interface
# Convert to JSON
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_api/scripts/convert_excel_to_json.py \
  Drugs___Biologics_Dec_19_2025_jak1.xlsx JAK1
```

### Step 2: Analyze Targets

Run the analysis script to calculate scores and rankings:

```bash
cd /home/sagemaker-user/.claude/skills/cortellis/cortellis_api/examples/

python analyze_targets_example.py
```

**Output:**
```
================================================================================
IBD Target Drug Development Analysis
================================================================================

Analyzing CDKN2D...
  ✓ Found 0 drugs, Total Score: 0.0
Analyzing GREM1...
  ✓ Found 2 drugs, Total Score: 1.6
Analyzing ITGA4...
  ✓ Found 12 drugs, Total Score: 8.5
...

================================================================================
Gene Rankings by Score
================================================================================

Rank   Gene         Total Score  Total Drugs  Key Phases
--------------------------------------------------------------------------------
1      TYK2         45.2         67           Approved_IBD: 2, Phase 3: 5, ...
2      ITGA4        8.5          12           Approved: 1, Phase 2: 2, ...
...
```

**Generated File:** `ibd_analysis_results.json`

### Step 3: Generate Report

Generate a comprehensive markdown report:

```bash
python generate_report_example.py
```

**Output:**
```
Generating comprehensive IBD target analysis report...

================================================================================
Report generated successfully!
================================================================================

Output file: IBD_Target_Clinical_Intelligence_Report.md
Report sections: Executive Summary, Target Rankings, Individual Analysis,
                 Combination Analysis, Strategic Recommendations, Methodology
```

**Generated File:** `IBD_Target_Clinical_Intelligence_Report.md`

## Customization

### Modify Gene Lists

Edit `analyze_targets_example.py` to change the gene lists:

```python
GENE_LISTS = {
    'List 1': ['TYK2', 'JAK1'],
    'List 2': ['YOUR', 'GENES'],
    'List 3': ['MORE', 'GENES'],
}
```

### Change Disease Context

Update the indication keywords for your therapeutic area:

```python
# IBD example (current)
IBD_INDICATIONS = [
    'inflammatory bowel disease', 'ibd', "crohn's disease",
    'crohn', 'ulcerative colitis', 'uc', 'colitis'
]

# Example: Oncology
ONCOLOGY_INDICATIONS = [
    'cancer', 'carcinoma', 'tumor', 'melanoma',
    'leukemia', 'lymphoma', 'sarcoma'
]
```

### Adjust Scoring Weights

Modify the scoring matrix in `calculate_drug_score()`:

```python
def calculate_drug_score(phase, has_disease_indication):
    """Calculate score for a single drug based on phase and indication."""
    # Scoring matrix - adjust weights as needed
    if phase == 'Approved':
        return 7 if has_disease_indication else 4  # ← Adjust these
    elif phase == 'Phase 3':
        return 3 if has_disease_indication else 2  # ← Adjust these
    # ... etc
```

**See:** `../references/scoring_framework.md` for detailed guidance on weight selection.

## Schema Compatibility

These scripts are compatible with **Cortellis API v2.0** JSON structure:

### Key Features

✅ **DrugRecord is a dictionary** keyed by drug ID (not drug name)
```python
# Correct usage
for drug_id, drug_info in drug_records.items():
    drug_name = drug_info.get('DrugName')
```

✅ **Phase extraction** uses `$` field
```python
phase_data = drug.get('PhaseHighest', {})
phase = extract_value(phase_data)  # Gets "$" value
```

✅ **Single/Array pattern** handled with helper functions
```python
# Indications can be single object or array
ind_list = as_list(ind_data.get('Indication', []))
```

✅ **Backwards compatible** with legacy field names
```python
# Tries IndicationsPrimary (v2.0) then IndicationsSecondary (legacy)
ind_data = drug.get('IndicationsPrimary') or drug.get('IndicationsSecondary', {})
```

### Helper Functions

Both scripts include helper functions from `json_schema.md`:

```python
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
```

## Data Requirements

### Input Files

The scripts expect Cortellis JSON files in the **current working directory**:

```
working_directory/
├── TYK2_cortellis_data.json
├── JAK1_cortellis_data.json
├── GREM1_cortellis_data.json
├── ITGA4_cortellis_data.json
└── ... (other gene files)
```

### JSON Structure

Each file must follow the Cortellis API v2.0 schema:

```json
{
  "annotation": {
    "@Id": "...",
    "Symbol": "TYK2",
    ...
  },
  "drug": {
    "Target": {...},
    "Drug": [...],
    "DrugRecord": {
      "drug_id_1": {
        "DrugName": "...",
        "PhaseHighest": {"@id": "LA", "$": "Launched"},
        "IndicationsPrimary": {
          "Indication": [
            {"@id": "213", "$": "Rheumatoid arthritis"}
          ]
        },
        ...
      }
    }
  }
}
```

**See:** `../references/json_schema.md` for complete schema documentation.

## Output Files

### ibd_analysis_results.json

Structured analysis data:

```json
{
  "gene_results": {
    "TYK2": {
      "gene": "TYK2",
      "total_drugs": 67,
      "drugs": [
        {
          "name": "tofacitinib",
          "phase": "Approved",
          "indications": ["Rheumatoid arthritis", "Ulcerative colitis"],
          "ibd_indication": true,
          "score": 7.0
        }
      ],
      "score_breakdown": {...},
      "total_score": 45.2
    }
  },
  "gene_rankings": [
    ["TYK2", 45.2],
    ["ITGA4", 8.5],
    ...
  ],
  "list_scores": {...}
}
```

### IBD_Target_Clinical_Intelligence_Report.md

Comprehensive markdown report with:
- Executive Summary
- Target Rankings
- Individual Target Analysis (per gene)
- Gene Combination Analysis
- Strategic Recommendations
- Methodology

**Sample:** See `IBD_Target_Clinical_Intelligence_Report.md` in this directory.

## Troubleshooting

### Issue: "No data available" for genes

**Cause:** JSON file missing or incorrectly named

**Solution:**
```bash
# Check files exist
ls *_cortellis_data.json

# Ensure filenames match: {GENE}_cortellis_data.json
# Example: TYK2_cortellis_data.json (not tyk2_cortellis_data.json)
```

### Issue: "DrugRecord is not a dict"

**Cause:** Incorrect JSON schema (old format or corrupted)

**Solution:** Re-query API or verify Excel conversion:
```bash
# Re-query from API
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_api/scripts/cortellis_gene_query.py \
  GENE_NAME --fields drug

# Or reconvert from Excel
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_api/scripts/convert_excel_to_json.py \
  excel_file.xlsx GENE_NAME
```

### Issue: Empty indications or phase data

**Cause:** Normal - some drugs lack comprehensive data

**Behavior:** Script handles gracefully:
- Missing phase → 'Unknown' → score 0
- Missing indications → empty list → non-specific scoring

**No action needed** - this is expected for research compounds.

### Issue: ibd_analysis_results.json not found

**Cause:** `analyze_targets_example.py` not run yet

**Solution:**
```bash
# Run analysis first
python analyze_targets_example.py

# Then generate report
python generate_report_example.py
```

## Performance

**Typical Runtime:**
- Analysis (10 genes, ~500 drugs): ~2-5 seconds
- Report generation: <1 second
- **Total workflow**: ~3-6 seconds (after API queries complete)

**Memory Usage:**
- JSON files: ~100KB - 5MB per gene (depending on drug count)
- Analysis: ~50-100MB RAM
- Report generation: ~20MB RAM

## Related Documentation

- **JSON Schema:** [`../references/json_schema.md`](../references/json_schema.md) - Complete API v2.0 data structure
- **API Reference:** [`../references/api_reference.md`](../references/api_reference.md) - API endpoints and performance
- **Scoring Framework:** [`../references/scoring_framework.md`](../references/scoring_framework.md) - Scoring methodology
- **API Access Guide:** [`../API_ACCESS.md`](../API_ACCESS.md) - Complete API usage documentation
- **Excel Conversion:** [`EXCEL_TO_JSON_CONVERSION_GUIDE.md`](./EXCEL_TO_JSON_CONVERSION_GUIDE.md) - Excel to JSON conversion

## Support

For issues or questions:
1. Check `../references/json_schema.md` for data structure questions
2. Review `../references/scoring_framework.md` for scoring methodology
3. See `EXCEL_TO_JSON_CONVERSION_GUIDE.md` for conversion issues
4. Consult `../API_ACCESS.md` for API query problems

---

**Version:** 2.0 (Updated for Cortellis API v2.0)
**Last Updated:** 2025-12-21
