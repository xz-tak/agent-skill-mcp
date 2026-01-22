# Cortellis API Access - Detailed Guide

Comprehensive guide for querying the Cortellis API to retrieve structured gene annotations and related data.

## Table of Contents

- [Prerequisites](#prerequisites)
- [When to Use API Access](#when-to-use-api-access)
- [Core Functionality](#core-functionality)
- [Query Options](#query-options)
- [Available Fields](#available-fields)
- [Output Files](#output-files)
- [Workflow Examples](#workflow-examples)
- [Data Structure](#data-structure)
- [Excel to JSON Conversion](#excel-to-json-conversion)
- [Data Analysis and Scoring](#data-analysis-and-scoring)
- [Example Analysis Workflows](#example-analysis-workflows)
- [Reference Documentation](#reference-documentation)
- [Error Handling](#error-handling)

## Prerequisites

**IMPORTANT: Cortellis API Credentials Required**

To use API access, valid Cortellis API credentials are required. The script reads credentials from:
1. A `.env` file in the working directory (preferred)
2. Environment variables `CORTELLIS_API_KEY` and `CORTELLIS_API_SECRET`

**Setup:**
1. Create a `.env` file in your working directory:
   ```bash
   cd /your/working/directory
   cat > .env << 'EOF'
   CORTELLIS_API_KEY=your_key_here
   CORTELLIS_API_SECRET=your_secret_here
   EOF
   ```

2. Verify credentials are configured:
   ```bash
   cat .env
   ```

If credentials are missing, the script will exit with an error message.

## When to Use API Access

Use API access when users request:
- Gene or protein target information from Cortellis
- Drug-target associations and development status
- Biomarker data for genes
- Protein-protein interactions
- Disease-gene associations
- Pharmaceutical research data for specific genes
- Batch queries for multiple genes
- Structured JSON data for programmatic analysis

## Core Functionality

### Main Script: `cortellis_gene_query.py`

The primary tool for querying Cortellis API. Execute this script to retrieve gene data.

**Location:** `/home/sagemaker-user/.claude/skills/cortellis/scripts/cortellis_gene_query.py`

**Basic Usage:**
```bash
python /home/sagemaker-user/.claude/skills/cortellis/scripts/cortellis_gene_query.py <GENE> [OPTIONS]
```

## Query Options

### Single Gene Query
```bash
python scripts/cortellis_gene_query.py BRCA1
```

### Multiple Genes
```bash
python scripts/cortellis_gene_query.py BRCA1 TP53 EGFR
```

### Specify Fields
```bash
python scripts/cortellis_gene_query.py BRCA1 --fields drug biomarker
```

### Query All Fields
```bash
python scripts/cortellis_gene_query.py BRCA1 --all
```

### Generate Excel Output
```bash
python scripts/cortellis_gene_query.py BRCA1 --excel
```

### Custom Output Directory
```bash
python scripts/cortellis_gene_query.py BRCA1 --output-dir ./results
```

### Skip Summary Output
```bash
python scripts/cortellis_gene_query.py BRCA1 --no-summary
```

## Available Fields

Reference [`references/api_fields.md`](../references/api_fields.md) for detailed field descriptions.

**Key fields:**
- **annotation** (always included): Target name, gene ID, UniProt ID, description
- **drug**: Drug associations, development phase, indications, mechanisms
- **biomarker**: Biomarker uses, roles, validity, studied drugs
- **interaction**: Protein-protein interactions with evidence
- **association**: Disease-gene associations and variants

**Default fields:** `drug` and `biomarker` (if not specified)

## Output Files

For each queried gene, the script generates:

### 1. JSON File: `{GENE}_cortellis_data.json`
- Complete API response data
- Nested JSON structure with all fields
- Contains both `Drug` (basic from Targets API) and `DrugRecord` (comprehensive from Investigational Drugs API)

### 2. Markdown Summary: `{GENE}_summary.md`
- Markdown-formatted summary
- Sections for each queried field
- Top results with key information
- Shows counts from both Targets API and Investigational Drugs API

### 3. Excel Workbook: `{GENE}_cortellis_data.xlsx` (with `--excel` flag)
- Excel workbook with multiple sheets
- Separate sheet per data field
- Flattened tabular data
- Auto-adjusted column widths
- Drug sheet uses comprehensive records from Investigational Drugs API

**Output Location:** Files save to the current working directory by default, or to `--output-dir` if specified.

## Workflow Examples

### Example 1: Comprehensive Gene Data
When a user asks: "Get drug and biomarker data for BRCA1 and EGFR, and create Excel files"

1. **Check for credentials**: Verify `.env` file exists
2. **Execute query**:
   ```bash
   python scripts/cortellis_gene_query.py BRCA1 EGFR --fields drug biomarker --excel
   ```
3. **Review outputs**:
   - Check generated JSON files for completeness
   - Review markdown summaries for key findings
   - Verify Excel files contain expected data

### Example 2: Research Analysis
```bash
python scripts/cortellis_gene_query.py MYC --all --excel --output-dir ./myc_analysis
```

### Example 3: Drug Target Investigation
```bash
python scripts/cortellis_gene_query.py EGFR KRAS --fields drug --excel
```

### Example 4: Batch Gene Query
```bash
python scripts/cortellis_gene_query.py BRCA1 BRCA2 TP53 --fields drug biomarker --excel
```

### Example 5: Quick Lookup (JSON only)
```bash
python scripts/cortellis_gene_query.py GREM1 --no-summary
```

## Data Structure

### Comprehensive Drug Records

The script queries two drug data sources:

#### 1. Targets API → `Drug` field
- Basic drug associations
- Molecular mechanisms
- Condition-drug associations

#### 2. Investigational Drugs API → `DrugRecord` field
- Development phase (Preclinical, Phase 1/2/3, Approved)
- Company information (Originator, Primary Companies)
- Therapeutic indications
- Mechanisms of action
- Therapy areas

Both fields are preserved in JSON output. The Excel Drug sheet uses comprehensive `DrugRecord` data when available.

**Search Process:**
1. Query Targets API for associated drugs
2. Extract drug `@namemain` identifier
3. Search Investigational Drugs API: `drugNamesAll:<@namemain>`
4. Fetch full drug record: `/drug/<drug_id>`
5. Store in separate `DrugRecord` field

**Note:** Not all drugs will have comprehensive records (e.g., research compounds may only have basic data).

### JSON Data Structure

All Cortellis data follows a consistent schema documented in [`references/json_schema.md`](../references/json_schema.md).

**Key Structure:**
```json
{
  "annotation": { ... },
  "drug": {
    "Target": { ... },
    "DrugRecord": {
      "drug_name": {
        "DrugName": "...",
        "PhaseHighest": {"$": "Approved"},
        "IndicationsSecondary": {
          "Indication": [{"$": "Disease name"}]
        }
      }
    }
  }
}
```

**Important Notes:**
- `DrugRecord` is nested: `data['drug']['DrugRecord']`
- `DrugRecord` is a dictionary keyed by drug name, not a list
- Important values often in `'$'` field
- Indications can be single object or array - always check type

## Excel to JSON Conversion

### Script: `convert_excel_to_json.py`

When API queries timeout (common for targets with extensive drug portfolios like JAK1, EGFR), convert manually downloaded Cortellis Excel files to compatible JSON format.

**📖 Comprehensive Guide:** See [`examples/EXCEL_TO_JSON_CONVERSION_GUIDE.md`](../examples/EXCEL_TO_JSON_CONVERSION_GUIDE.md) for detailed instructions, troubleshooting, and schema mapping.

**Location:** `/home/sagemaker-user/.claude/skills/cortellis/cortellis_api/scripts/convert_excel_to_json.py`

**Basic Usage:**
```bash
python scripts/convert_excel_to_json.py <excel_file> [gene_name]
```

**Examples:**
```bash
# Single file - gene name auto-detected from filename
python scripts/convert_excel_to_json.py Drugs___Biologics_Dec_19_2025_tyk2.xlsx

# Specify gene name explicitly
python scripts/convert_excel_to_json.py cortellis_export.xlsx TYK2

# Multiple files (batch processing)
python scripts/convert_excel_to_json.py file1.xlsx file2.xlsx file3.xlsx

# Custom output filename
python scripts/convert_excel_to_json.py file.xlsx TYK2 --output custom_output.json
```

**When to Use:**
- API timeouts for genes with >500 drug associations
- Manual exports from Cortellis web interface already available
- Batch processing of multiple manual downloads
- Offline analysis scenarios
- Need for additional Excel-specific data (regional development, milestones)

**Requirements:**
- Excel file must contain "Product List" sheet (required)
- Optional sheets: "Development Status", "Milestones"
- Key columns: Entry Number, Drug Name (All), Generic Name, Highest Phase, Condition, Organization, Mechanism of Action
- Output JSON format matches API structure exactly

**Output:**
- Generates `{GENE}_cortellis_data.json` in current directory
- Compatible with all analysis scripts using API JSON
- Includes phase normalization to API format (e.g., "Launched - 2022" → `{"@id": "LA", "$": "Launched"}`)
- Automatic indication and mechanism parsing (newline-delimited fields → structured arrays)
- Additional fields from Excel: `RegionalDevelopment`, `Milestones`, `_excel_row`
- Valid JSON (NaN values automatically converted to null)

**Key Features:**
- **Schema Compatibility:** Output matches Cortellis API JSON schema exactly
- **Phase Normalization:** Excel phase strings converted to API format with @id and $ keys
- **Multiline Fields:** Automatic parsing of newline-delimited conditions, companies, mechanisms
- **NaN Handling:** Pandas NaN values cleaned to JSON-valid null
- **Multiple Sheets:** Processes Product List, Development Status, and Milestones sheets
- **Batch Processing:** Convert multiple Excel files in one command

## Data Analysis and Scoring

### Flexible Scoring Framework

The skill supports flexible target scoring for clinical validation and competitive intelligence. Scoring weights should be adjusted based on therapeutic area, analysis purpose, and market dynamics.

See [`references/scoring_framework.md`](../references/scoring_framework.md) for detailed guidance.

#### Default Scoring Matrix (Clinical Validation)

| Phase | Disease-Specific | Non-Specific | Rationale |
|-------|------------------|--------------|-----------|
| FDA Approved / On Market | 7 | 4 | Highest validation; disease-specific proves efficacy |
| Phase 3 Clinical | 3 | 2 | Late-stage validation; high approval probability |
| Phase 2 Clinical | 2 | 1 | Mid-stage validation; proof-of-concept established |
| Phase 1 Clinical | 1 | 0.5 | Early validation; safety and PK/PD characterized |
| Preclinical | 0.1 | 0.1 | Minimal validation; target engagement shown |
| Discontinued | 0 | 0 | No active development |

**Disease-Specific:** Drug indication matches target disease (e.g., IBD drug for IBD target analysis)

#### Context-Specific Adjustments

Modify weights based on:

**High-Risk Areas (Oncology, Rare Diseases):**
- Increase early-phase weights (Phase 1/2)
- Any clinical data is valuable

**Crowded Markets (Competitive Intelligence):**
- Decrease early-phase weights
- Emphasize approved drugs only

**Repurposing Opportunities:**
- Increase "non-specific" weights for approved drugs
- Approved drugs in adjacent indications easier to repurpose

**Competitive Intensity Measurement:**
- Set all weights to 1 (equal counting)
- Measures total competitive activity

#### Score Interpretation

| Score Range | Interpretation | Strategic Guidance |
|-------------|----------------|-------------------|
| > 100 | Extremely validated, high competition | Need strong differentiation strategy |
| 50-100 | Well-validated, competitive | Proven biology; focus on novelty |
| 20-50 | Moderate validation | Balanced risk-reward opportunity |
| 5-20 | Emerging target | Early stage; first-mover advantage |
| 1-5 | Minimal activity | High risk; potential breakthrough |
| 0 | No activity | Blue ocean or undruggable target |

## Example Analysis Workflows

Complete reference implementations demonstrating end-to-end target analysis pipelines using Cortellis API v2.0 data.

**📖 Comprehensive Guide:** See [`examples/README.md`](../examples/README.md) for:
- Quick start guide (3-step workflow)
- Customization instructions (gene lists, disease context, scoring weights)
- Schema compatibility notes and helper functions
- Troubleshooting and performance metrics
- Complete usage documentation

**Example Scripts:**
- [`examples/analyze_targets_example.py`](../examples/analyze_targets_example.py) - Target scoring, ranking, and combination analysis
- [`examples/generate_report_example.py`](../examples/generate_report_example.py) - Clinical intelligence report generation
- [`examples/IBD_Target_Clinical_Intelligence_Report.md`](../examples/IBD_Target_Clinical_Intelligence_Report.md) - Sample output report

**Quick Workflow:**
```bash
# 1. Query API for genes
python scripts/cortellis_gene_query.py TYK2 GREM1 ITGA4 --fields drug

# 2. Analyze targets
cd examples/
python analyze_targets_example.py

# 3. Generate report
python generate_report_example.py
```

For detailed instructions, customization options, and troubleshooting, see [`examples/README.md`](../examples/README.md).

## Reference Documentation

Comprehensive documentation for data structures, APIs, and methodologies:

### API and Data Structures
- [`references/api_fields.md`](../references/api_fields.md) - API field descriptions, data types, output formats
- [`references/api_reference.md`](../references/api_reference.md) - API endpoints, authentication, request/response structures, error handling
- [`references/json_schema.md`](../references/json_schema.md) - JSON data structure, schema guide, accessing data in Python, common pitfalls

### Analysis Methodology
- [`references/scoring_framework.md`](../references/scoring_framework.md) - Flexible scoring methodology, context-specific weight adjustments, interpretation guidelines, configuration examples

### Example Implementations
- **[`examples/README.md`](../examples/README.md)** - Comprehensive guide to example workflows, customization, and usage
- [`examples/analyze_targets_example.py`](../examples/analyze_targets_example.py) - Target scoring and ranking (compatible with API v2.0)
- [`examples/generate_report_example.py`](../examples/generate_report_example.py) - Clinical intelligence report generation
- [`examples/EXCEL_TO_JSON_CONVERSION_GUIDE.md`](../examples/EXCEL_TO_JSON_CONVERSION_GUIDE.md) - Excel to JSON conversion guide
- [`examples/IBD_Target_Clinical_Intelligence_Report.md`](../examples/IBD_Target_Clinical_Intelligence_Report.md) - Sample report output

## Error Handling

The script handles:
- **Missing credentials** → Exits with error message
- **Gene not found** → Displays error and continues with remaining genes
- **API timeouts** → Reports warning and continues
- **Network errors** → Displays error message
- **Invalid field names** → Shows available options

If errors occur for one gene in a multi-gene query, the script continues processing remaining genes.

## Dependencies

Required Python packages (automatically installed in user environment):
- requests
- python-dotenv
- pandas
- openpyxl

## Common Issues and Solutions

### Issue: API Timeout for Large Datasets
**Solution:** Use Excel to JSON conversion
```bash
# Download Excel manually from Cortellis web interface
# Convert to JSON
python scripts/convert_excel_to_json.py downloaded_file.xlsx GENE_NAME
```

### Issue: Missing Credentials
**Solution:** Create `.env` file
```bash
cd /your/working/directory
echo "CORTELLIS_API_KEY=your_key" > .env
echo "CORTELLIS_API_SECRET=your_secret" >> .env
```

### Issue: Gene Not Found
**Solution:** Check gene symbol spelling and try alternatives (e.g., EGFR vs egfr)

---

**For additional help, see:**
- **Main Documentation:** [`../SKILL.md`](../SKILL.md)
- **API Reference:** [`../references/api_reference.md`](../references/api_reference.md)
- **Complete Navigation:** [`../DOCUMENTATION_MAP.md`](../DOCUMENTATION_MAP.md)
