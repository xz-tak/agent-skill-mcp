# Schema Exploration Guide

This document provides detailed guidance on using the schema exploration tools to discover available metadata values in Omicsoft h5ad files before running DEG analysis.

## Try the Example First

**New to schema exploration?** Start by opening the included example:

1. **Open** `assets/schema_viewer_example.html` in your browser (just double-click)
2. **Explore** real data from IBD/MASH/Fibro/Derm/Rheum dataset:
   - 7,604 observations (comparisons)
   - 60,071 variables (genes)
   - 40 metadata columns fully documented
   - 178 diseases, 131 tissues, 70+ treatment statuses
3. **Search** for diseases, tissues, or treatments
4. **Copy** exact filter values for your own analysis

**Example files included in this skill**:
- `references/schema_report_example.json`: Complete schema data (3.6 MB JSON)
- `assets/schema_viewer_example.html`: Interactive viewer (2.4 MB) - works offline
- `references/schema_quickstart.md`: Quick start guide

## Overview

The schema exploration tools help users:
- **Discover available filter values** before running analysis
- **Explore metadata comprehensively** including demographics, diseases, tissues, treatments
- **Find exact strings** needed for filtering parameters
- **Understand data structure** and column relationships
- **Validate data contents** before analysis

## Tools

### 1. explore_h5ad_schema.py

Generates comprehensive reports of h5ad file metadata in JSON or text format.

**Purpose**: Extract and document all unique values from every metadata column with counts and percentages.

**Usage**:
```bash
# Generate JSON format (recommended for programmatic access)
conda run -n <env> python scripts/explore_h5ad_schema.py \
  --file <path_to_h5ad> \
  --output schema_report.json \
  --format json

# Generate text format (human-readable summary)
conda run -n <env> python scripts/explore_h5ad_schema.py \
  --file <path_to_h5ad> \
  --output schema_report.txt \
  --format text

# Search for specific terms (text format only)
conda run -n <env> python scripts/explore_h5ad_schema.py \
  --file <path_to_h5ad> \
  --search scleroderma ibd skin \
  --format text
```

**Output**:
- **JSON format**: Complete structured data with all values, no truncation (typically 2-5 MB)
- **Text format**: Human-readable summary with top values (may truncate long lists)

### 2. generate_schema_viewer.py

Creates standalone HTML viewers with embedded JSON data for interactive exploration.

**Purpose**: Generate self-contained HTML files that work without web servers or external dependencies.

**Usage**:
```bash
# Create standalone HTML viewer
conda run -n <env> python scripts/generate_schema_viewer.py \
  --json schema_report.json \
  --output schema_viewer_standalone.html
```

**Output**:
- **Standalone HTML file** (2-4 MB) with all data embedded
- Works by double-clicking (no web server needed)
- Fully functional search, filtering, and navigation

## JSON Schema Structure

The JSON output from `explore_h5ad_schema.py` has the following structure:

```json
{
  "file_info": {
    "file_path": "path/to/file.h5ad",
    "n_observations": 7604,
    "n_variables": 60071,
    "observation_description": "Each observation represents a differential expression comparison",
    "variable_description": "Each variable represents a gene"
  },

  "obs_columns": {
    "comparison_category": {
      "total_values": 7604,
      "unique_values": 8,
      "null_count": 0,
      "dtype": "category",
      "category": "key_filtering",
      "all_values": [
        {"value": "Treatment vs. Control", "count": 1620, "percentage": 21.3},
        {"value": "Disease vs. Normal", "count": 1063, "percentage": 13.98},
        ...
      ]
    },
    "case_ethnicity": {
      "total_values": 7604,
      "unique_values": 8,
      "null_count": 45,
      "dtype": "category",
      "category": "demographic",
      "all_values": [
        {"value": "NA", "count": 7516, "percentage": 98.8},
        {"value": "European American", "count": 19, "percentage": 0.2},
        ...
      ]
    },
    ...
  },

  "column_categories": {
    "key_filtering": ["comparison_category", "tissue", "disease", ...],
    "demographic": ["case_gender", "case_ethnicity", ...],
    "case_disease": ["case_disease_state", ...],
    "treatment": ["case_dosage", "case_treatment", ...],
    "response": ["case_response", "control_response"],
    "sample": ["case_tissue", "sample", ...],
    "comparison_details": ["comparison_id", "comparison_contrast"]
  },

  "filtering_guide": {
    "exact_match_filters": ["comparison_category", "study", "case_treatment", ...],
    "partial_match_filters": ["tissue", "disease", "comparison"],
    "note": "Exact match filters require exact string matching. Partial match filters are case-insensitive substring matching."
  },

  "var_info": {
    "n_genes": 60071,
    "gene_list_sample": ["A1BG", "A1CF", "A2M", ...],
    "all_genes_available": true
  },

  "layers": {
    "log2fc": {"shape": [7604, 60071], "type": "csr_matrix", ...},
    "padj": {"shape": [7604, 60071], "type": "csr_matrix", ...},
    "sig_score": {"shape": [7604, 60071], "type": "csr_matrix", ...}
  }
}
```

### Key Field Descriptions

**obs_columns**: Dictionary of all metadata columns
- Each column contains `all_values` array with every unique value
- `total_values`: Number of observations
- `unique_values`: Number of distinct values
- `null_count`: Missing values
- `category`: Column category (key_filtering, demographic, etc.)

**column_categories**: Organizes columns by purpose
- `key_filtering`: Columns used for deg_analysis.py filtering
- `demographic`: Patient/control demographics (gender, ethnicity, age)
- `case_disease`/`control_disease`: Disease information
- `treatment`: Treatment details
- `response`: Treatment response data
- `sample`: Sample metadata
- `comparison_details`: Comparison identifiers

**filtering_guide**: Instructions for using filter values
- `exact_match_filters`: Must match string exactly (e.g., comparison_category)
- `partial_match_filters`: Case-insensitive substring match (e.g., disease, tissue)

## Using the HTML Viewer

### Opening the Viewer

**Standalone HTML (Recommended)**:
1. Double-click `schema_viewer_standalone.html`
2. Opens directly in default browser
3. No setup required

**Alternative Methods**:
- Right-click → "Open with" → Select browser
- Drag file onto browser window
- Use `file://` URL in browser address bar

### Viewer Features

**Search Functionality**:
- Type in search box to filter across all columns and values
- Searches both column names and value strings
- Real-time filtering with highlighting

**Category Filtering**:
- Dropdown menu to filter by category
- Categories: Demographics, Key Filtering, Treatment, etc.
- Shows only relevant columns for selected category

**Interactive Navigation**:
- Click category headers to expand/collapse sections
- Click column headers to see all unique values
- "Expand All" / "Collapse All" buttons for bulk operations

**Value Display**:
- Each value shows: exact string, count, percentage
- Visual progress bars show relative frequencies
- Copy exact strings for use in filter parameters

**Color Coding**:
- Blue tags: Key filtering columns (used in deg_analysis.py)
- Green tags: Demographic columns

## Finding Filter Values for deg_analysis.py

### Workflow

1. **Generate schema report**:
   ```bash
   python scripts/explore_h5ad_schema.py --file data.h5ad --output schema.json --format json
   python scripts/generate_schema_viewer.py --json schema.json --output viewer.html
   ```

2. **Open viewer** (`viewer.html`) in browser

3. **Search for values**:
   - Search "scleroderma" → find all scleroderma-related diseases
   - Search "skin" → find skin-related tissues
   - Search "infliximab" → find treatment statuses

4. **Copy exact strings** from viewer

5. **Validate filters before full analysis** (Recommended):
   ```bash
   python scripts/validate_filters.py \
     --file data.h5ad \
     --target-name MyAnalysis \
     --diseases "systemic scleroderma,diffuse scleroderma" \
     --tissues "skin,epidermis,dermis" \
     --comparison-category "Disease vs. Normal" \
     --signatures "Sig1:Gene1,Gene2,Gene3" \
     ...
   ```

   This step-by-step validation:
   - Tests each filter incrementally
   - Shows matching values and observation counts
   - Detects zero-result queries before running full analysis
   - Provides suggestions for fixing problematic filters

6. **Run full analysis** (after validation succeeds):
   ```bash
   python scripts/deg_analysis.py \
     --diseases "systemic scleroderma,diffuse scleroderma" \
     --tissues "skin,epidermis,dermis" \
     --comparison-category "Disease vs. Normal" \
     ...
   ```

### Filter Types

**Exact Match** (comma-separated, exact strings):
- `--comparison-category`: e.g., "Disease vs. Normal,Treatment vs. Control"
- `--studies`: e.g., "GSE130955,GSE181549"
- `--case-treatment`: e.g., "none,NA"
- `--control-treatment`: e.g., "none,culture medium"
- `--case-treatment-status`: e.g., "infliximab,adalimumab"
- `--control-treatment-status`: e.g., "placebo,none"

**Partial Match** (comma-separated, case-insensitive substring):
- `--diseases`: e.g., "scleroderma,sclerosis" (matches "systemic scleroderma", "diffuse scleroderma", etc.)
- `--tissues`: e.g., "skin" (matches "skin", "skin biopsy", etc.)
- `--comparison`: e.g., "response vs no response" (substring match)

### Examples

**Example 1: Find all scleroderma diseases**
1. Open viewer, search "sclero"
2. Expand "disease" column
3. See: "systemic scleroderma" (80 obs), "diffuse scleroderma" (40 obs), "limited scleroderma" (18 obs)
4. Use: `--diseases "scleroderma,sclerosis"`

**Example 2: Find treatment response comparisons**
1. Search "response"
2. Look at "comparison" column
3. Find: "response vs no response" (20 obs), "responder vs non-responder" variations
4. Use: `--comparison "response vs no response"`

**Example 3: Filter by exact comparison category**
1. Click "comparison_category" column
2. See all 8 options with counts
3. Copy exact strings: "Disease vs. Normal", "Treatment vs. Control"
4. Use: `--comparison-category "Disease vs. Normal,Treatment vs. Control"`

## Demographic Information

The schema includes demographic metadata when available:

**Columns**:
- `case_gender`, `control_gender`: male, female, NA
- `case_ethnicity`, `control_ethnicity`: European American, African American, Asian Chinese, Hispanic, White, Mixed Ethnicity, etc.
- `case_age_category`, `control_age_category`: adult, child, fetus, NA

**Availability**:
Most demographic fields are "NA" (>95%) because:
- Not all studies report demographics
- Privacy concerns in public datasets
- Variability in reporting standards

**Usage**:
While demographics are sparse, they can be useful for:
- Understanding study populations
- Identifying studies with demographic data
- Subset analyses when available

**Note**: Country/nation information is not stored in h5ad metadata. Geographic information may be inferred from study accessions (e.g., GEO database records) but is not directly available in the file.

## Programmatic Access

### Python Example

```python
import json

# Load schema
with open('schema_report.json', 'r') as f:
    schema = json.load(f)

# Get all comparison categories
categories = schema['obs_columns']['comparison_category']['all_values']
for cat in categories:
    print(f"{cat['value']}: {cat['count']} ({cat['percentage']:.1f}%)")

# Find scleroderma-related diseases
diseases = schema['obs_columns']['disease']['all_values']
scleroderma = [d for d in diseases if 'sclero' in d['value'].lower()]
print(f"Found {len(scleroderma)} scleroderma-related diseases")
for d in scleroderma:
    print(f"  - {d['value']}: {d['count']} obs")

# Get all key filtering columns
key_filters = schema['column_categories']['key_filtering']
print(f"Key filtering columns: {', '.join(key_filters)}")

# Check ethnicity distribution
ethnicity = schema['obs_columns']['case_ethnicity']['all_values']
for eth in ethnicity:
    if eth['value'] != 'NA':
        print(f"{eth['value']}: {eth['count']} ({eth['percentage']:.1f}%)")
```

### JavaScript Example (for web apps)

```javascript
// Load and parse JSON
fetch('schema_report.json')
  .then(r => r.json())
  .then(schema => {
    // Get all tissue types
    const tissues = schema.obs_columns.tissue.all_values;
    console.log(`Total tissues: ${tissues.length}`);

    // Filter for skin-related
    const skinTissues = tissues.filter(t =>
      t.value.toLowerCase().includes('skin')
    );
    console.log('Skin tissues:', skinTissues);

    // Get filtering guide
    const guide = schema.filtering_guide;
    console.log('Exact match filters:', guide.exact_match_filters);
    console.log('Partial match filters:', guide.partial_match_filters);
  });
```

## Column Categories Explained

### Key Filtering (10 columns)
Primary columns for filtering in deg_analysis.py:
- `comparison_category`: Type of comparison (Disease vs. Normal, etc.)
- `tissue`: Tissue type (peripheral blood, skin, etc.)
- `disease`: Disease name (crohn's disease, RA, etc.)
- `case_treatment`, `control_treatment`: Specific treatments applied
- `case_treatment_status`, `control_treatment_status`: Treatment/drug names
- `comparison`: Specific comparison description
- `study`: Study identifier (GSE accessions)
- `database`: Source database

### Demographics (6 columns)
Patient and control demographic information:
- `case_age_category`, `control_age_category`: Age groups
- `case_gender`, `control_gender`: Sex/gender
- `case_ethnicity`, `control_ethnicity`: Race/ethnicity

### Disease Information (8 columns)
Detailed disease characteristics:
- `case_disease_state`, `control_disease_state`: Disease status
- `case_disease_subtype`, `control_disease_subtype`: Disease subtypes
- `case_disease_group`, `control_disease_group`: Disease groupings
- `case_disease_location`, `control_disease_location`: Anatomical location

### Treatment (6 columns)
Treatment details:
- `case_dosage`, `control_dosage`: Dosage information
- `case_treatment_group`, `control_treatment_group`: Treatment arms
- `case_treat_time`, `control_treat_time`: Treatment duration/timing

### Response (2 columns)
Treatment response data:
- `case_response`: Patient response (response, no response, partial, etc.)
- `control_response`: Control response

### Sample (6 columns)
Sample metadata:
- `case_tissue`, `control_tissue`: Tissue sampled
- `case_sample_material`, `control_sample_material`: Sample type
- `sample`: Sample identifier
- `project_id`: Project/study identifier

### Comparison Details (2 columns)
Comparison identifiers:
- `comparison_id`: Unique comparison identifier
- `comparison_contrast`: Comparison description

## Troubleshooting

### Schema Generation Issues

**Problem**: "File not found" error
**Solution**: Verify h5ad file path is correct and file exists

**Problem**: "Out of memory" error
**Solution**: H5ad file may be very large. The script automatically removes embedded expression data, but ensure sufficient RAM (16+ GB recommended for large files)

**Problem**: S3 access errors
**Solution**: Ensure `s3fs` and `boto3` are installed, AWS credentials configured, and IAM permissions allow S3 access

### HTML Viewer Issues

**Problem**: Viewer opens but shows no data
**Solution**: Make sure you're using `schema_viewer_standalone.html` (not `schema_viewer.html` which requires web server)

**Problem**: Search not working
**Solution**: Try expanding categories first or use "Expand All" button

**Problem**: Can't open HTML file
**Solution**: Right-click → "Open with" → Choose browser, or drag file onto browser window

### Filter Value Issues

**Problem**: No observations after filtering
**Solution**:
- Check filter values match available values in schema
- Use schema viewer to verify exact strings
- Try broader filter terms (partial matching for diseases/tissues)
- Verify organism-specific gene symbols

**Problem**: Too many results
**Solution**: Add more specific filters or increase thresholds (lfc, padj)

## Best Practices

1. **Always generate schema first** when working with new h5ad files
2. **Use standalone HTML viewer** for easiest exploration
3. **Copy exact strings** from viewer for exact-match filters
4. **Use partial matching wisely** for diseases and tissues (broader terms find more data)
5. **Check demographic availability** before planning demographic subgroup analyses
6. **Save schema reports** for reference and sharing with collaborators
7. **Update schema** if h5ad file is modified or regenerated

## Integration with deg_analysis.py

The schema exploration workflow integrates seamlessly with DEG analysis:

```bash
# Step 1: Explore schema
python scripts/explore_h5ad_schema.py --file data.h5ad --output schema.json --format json
python scripts/generate_schema_viewer.py --json schema.json --output viewer.html

# Step 2: Open viewer.html, find filter values

# Step 3: Run analysis with discovered values
python scripts/deg_analysis.py \
  --file data.h5ad \
  --target-name MyAnalysis \
  --diseases "systemic scleroderma,diffuse scleroderma" \
  --tissues "skin,epidermis" \
  --comparison-category "Disease vs. Normal" \
  --signatures "MySignature:GENE1,GENE2,GENE3" \
  --lfc-threshold 0.5 \
  --padj-threshold 0.05 \
  --run-gsea
```

This ensures you use correct filter values and understand your data before running analysis.
