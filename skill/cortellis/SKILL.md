---
name: cortellis
description: Unified toolkit for Cortellis Drug Discovery Intelligence + OFF-X: (1) build Excel→interactive CI dashboards with competition/opportunity scoring, (2) query targets/drugs via API for structured analysis, (3) automate web exports across Cortellis categories, and (4) export OFF-X safety evidence and adverse events.
---

# Cortellis Intelligence Skill

Access Cortellis Drug Discovery Intelligence and OFF-X through four complementary subskills for competitive intelligence, program tracking, and safety evaluation.

## What This Skill Does

- `cortellis_ci_web`: Read user Excel exports, normalize assets/targets/MoA/mechanisms/phases/routes, score competition vs blue-ocean opportunity, and generate a single-file interactive HTML dashboard.
- `cortellis_api`: Programmatic API queries for genes/targets/drugs/biomarkers; produces JSON/Excel/MD for analysis and reporting.
- `cortellis_targetdrug_web`: Playwright automation to download Cortellis category exports (Drugs & Biologics, Clinical Studies, Patents, etc.) to Excel/CSV.
- `offx_web`: Playwright automation to export OFF-X safety tables/adverse events for targets and drugs.

## 📁 Skill Structure

```
/home/sagemaker-user/.claude/skills/cortellis/
├── SKILL.md                          # This file - High-level overview
├── DOCUMENTATION_MAP.md              # Complete navigation guide
│
├── cortellis_api/                    # API Access Subskill
│   ├── API_ACCESS.md                 # Detailed API documentation
│   ├── scripts/                      # Python scripts
│   │   ├── cortellis_gene_query.py   # Main API query script
│   │   └── convert_excel_to_json.py  # Excel to JSON converter
│   ├── references/                   # Reference documentation
│   │   ├── api_reference.md          # API endpoints and authentication
│   │   ├── api_fields.md             # API field descriptions
│   │   ├── json_schema.md            # JSON data structure guide
│   │   └── scoring_framework.md      # Target scoring methodology
│   ├── examples/                     # Example analysis scripts
│   │   ├── analyze_targets_example.py
│   │   ├── generate_report_example.py
│   │   └── IBD_Target_Clinical_Intelligence_Report.md
│   └── assets/                       # Additional resources
│
├── cortellis_ci_web/                 # CI Web Dashboard Subskill (Excel → interactive HTML)
│   ├── CI_WEB.md                     # Dashboard generation guide
│   ├── README.md                     # Quick reference
│   ├── scripts/                      # Excel → dashboard generator scripts
│   ├── references/                   # Taxonomy + scoring notes
│   └── assets/                       # HTML template
│
├── cortellis_targetdrug_web/         # Cortellis Target-Drug Web Automation Subskill
│   ├── WEB_AUTOMATION.md             # Detailed web automation documentation
│   ├── SETUP_COMPLETE.md             # Setup guide with test results
│   ├── README.md                     # Quick reference
│   ├── okta_auth_setup.py            # One-time authentication setup
│   ├── cortellis-automation.js       # Single-category export script
│   ├── cortellis-download.js         # Multi-category download script
│   ├── run-cortellis.sh              # Wrapper for automation.js
│   ├── run-cortellis-download.sh     # Wrapper for download.js
│   └── example/                      # Example config files
│       ├── example_clinical_studies.json
│       ├── example_comprehensive.json
│       └── example_patents.json
│
└── offx_web/                         # OFF-X Safety Web Automation Subskill
    ├── OFFX_AUTOMATION.md            # Detailed OFF-X automation documentation
    ├── README.md                     # Quick reference
    ├── okta_auth_setup.py            # One-time authentication setup
    ├── offx-download.js              # Multi-entity safety download script
    ├── run-offx-download.sh          # Wrapper for offx-download.js
    └── example/                      # Example config files
        ├── example_targets.json
        ├── example_drugs_biologics.json
        ├── example_adverse_events.json
        └── example_comprehensive.json
```

## 📚 Documentation Map

**Complete Navigation:** See [`DOCUMENTATION_MAP.md`](DOCUMENTATION_MAP.md) for full documentation structure

**Quick Start:**
- **For CI dashboards from Excel exports:** See [`cortellis_ci_web/CI_WEB.md`](cortellis_ci_web/CI_WEB.md)
- **For API Access:** See [`cortellis_api/API_ACCESS.md`](cortellis_api/API_ACCESS.md)
- **For Cortellis Target-Drug Web Automation:** See [`cortellis_targetdrug_web/WEB_AUTOMATION.md`](cortellis_targetdrug_web/WEB_AUTOMATION.md)
- **For OFF-X Safety Assessment:** See [`offx_web/OFFX_AUTOMATION.md`](offx_web/OFFX_AUTOMATION.md)

## 🔀 Choose Your Access Method

### CI Dashboard from Excel (Python) - For Competitive Intelligence

Turn Cortellis exports or curated landscape spreadsheets into a single-file interactive dashboard with competitive intensity scoring and opportunity ranking.

**Best for:**
- Visualizing crowded vs blue-ocean spaces across targets/target types/MoA/routes
- Fast iteration with taxonomy rule updates (YAML)
- Sharing a portable HTML with filtering + export

**Documentation:** [`cortellis_ci_web/CI_WEB.md`](cortellis_ci_web/CI_WEB.md)

**Main Script:** `cortellis_ci_web/scripts/generate_ci_dashboard.py`

---

### API Access (Python) - For Programmatic Queries

Query the Cortellis API to retrieve structured gene annotations and related data with optimized async performance.

**Best for:**
- Structured JSON data for analysis
- Python data analysis workflows
- Gene/target research with biomarker data
- Protein-protein interactions
- Integration with analysis pipelines
- High-performance batch queries (1700+ drugs in ~18 seconds)

**Main Script:** `cortellis_api/scripts/cortellis_gene_query.py`

**Features:**
- ⚡ **Async optimization**: 11.8x faster than sequential with adaptive concurrency
- 🎯 **99.5%+ data completeness** with intelligent retry logic
- 🔗 **Comprehensive coverage**: Automatically queries related targets (isoforms, variants)
- 📋 **Excel conversion**: Fallback for API timeouts with schema-compatible output

**Quick Example:**
```bash
# Query gene data
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_api/scripts/cortellis_gene_query.py \
  BRCA1 --fields drug biomarker --excel
```

**Detailed Documentation:** [`cortellis_api/API_ACCESS.md`](cortellis_api/API_ACCESS.md)

---

### Cortellis Target-Drug Web Automation (Playwright) - For File Downloads

Browser-based automation to search, navigate, and download target-drug relationship data from specific Cortellis categories. Supports queries for targets with associated drugs and drugs with associated targets.

**Best for:**
- Target-to-Drug queries (e.g., "What drugs target EGFR?")
- Drug-to-Target queries (e.g., "What targets does imatinib affect?")
- Downloading Excel/CSV files with target-drug relationships
- Bulk data exports
- Category-specific downloads (Genes & Targets, Drugs & Biologics, Clinical Studies, Patents, Literature, etc.)
- Multiple target or drug queries with file downloads

**Main Scripts:**
- `cortellis_targetdrug_web/run-cortellis.sh` - Single-category export
- `cortellis_targetdrug_web/run-cortellis-download.sh` - Multi-category downloads

**Quick Example:**
```bash
# Download target data with associated drugs
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh \
  --queries "EGFR,TYK2" --categories "Genes & Targets,Drugs & Biologics"

# Download drug data with associated targets
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh \
  --queries "imatinib,gefitinib" --categories "Drugs & Biologics,Clinical Studies"
```

**Detailed Documentation:** [`cortellis_targetdrug_web/WEB_AUTOMATION.md`](cortellis_targetdrug_web/WEB_AUTOMATION.md)

---

### OFF-X Safety Web Automation (Playwright) - For Safety Assessment

Browser-based automation to evaluate target and drug safety, identify adverse events, and download Master view safety data from OFF-X (TargetSafety.info).

**Best for:**
- Target safety profile evaluation
- Drug safety assessment
- Adverse event analysis
- Master view data export
- Comprehensive safety tables in Excel format
- Bulk safety assessments for multiple entities

**Main Script:**
- `offx_web/run-offx-download.sh` - Multi-entity safety downloads

**Quick Example:**
```bash
# Assess target safety for ITGA4
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh itga4

# Drug safety assessment
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh \
  imatinib --field "Drugs and biologics"

# Multiple targets safety
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh \
  --queries "tyk2,jak1,egfr"
```

**Detailed Documentation:** [`offx_web/OFFX_AUTOMATION.md`](offx_web/OFFX_AUTOMATION.md)

## When to Use Which Method

Use this decision matrix to choose the appropriate method:

| User Request | Recommended Method | Script to Use | Reason |
|-------------|-------------------|---------------|---------|
| "What drugs target EGFR?" | API Access | `cortellis_gene_query.py` | Structured data, programmatic analysis |
| "Get biomarker data for BRCA1" | API Access | `cortellis_gene_query.py` | Gene-specific structured data |
| "Download clinical trials for imatinib" | Target-Drug Web | `run-cortellis.sh` | Need Excel file export |
| "Get patents for 5 kinase inhibitors" | Target-Drug Web | `run-cortellis-download.sh` | Bulk file downloads |
| "Analyze drug development phases" | API Access | `cortellis_gene_query.py` | JSON data for analysis |
| "Export all data types for gefitinib" | Target-Drug Web | `run-cortellis-download.sh` | Multiple category downloads |
| "Assess safety profile for ITGA4" | OFF-X Web | `run-offx-download.sh` | Target safety assessment |
| "Get adverse events for imatinib" | OFF-X Web | `run-offx-download.sh` | Drug safety data |
| "Evaluate safety for multiple targets" | OFF-X Web | `run-offx-download.sh` | Bulk safety assessments |
| "Master view data for TYK2" | OFF-X Web | `run-offx-download.sh` | Safety Master view export |

## Prerequisites

### For API Access
- **Required:** Cortellis API credentials
- **Setup:** Create `.env` file in working directory with:
  ```
  CORTELLIS_API_KEY=your_key
  CORTELLIS_API_SECRET=your_secret
  ```
- **Documentation:** See [`cortellis_api/API_ACCESS.md`](cortellis_api/API_ACCESS.md#prerequisites)

### For Cortellis Target-Drug Web Automation
- **Required:** Okta authentication session
- **Setup:** Run once from working directory:
  ```bash
  cd /your/working/directory
  python /home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/okta_auth_setup.py
  ```
- Creates `okta_auth_state.json` (valid for weeks/months)
- **Documentation:** See [`cortellis_targetdrug_web/WEB_AUTOMATION.md`](cortellis_targetdrug_web/WEB_AUTOMATION.md#prerequisites)

### For OFF-X Safety Web Automation
- **Required:** Okta authentication session (same as Cortellis)
- **Setup:** Run once from working directory:
  ```bash
  cd /your/working/directory
  python /home/sagemaker-user/.claude/skills/cortellis/offx_web/okta_auth_setup.py
  ```
- Creates `okta_auth_state.json` (valid for weeks/months)
- **Documentation:** See [`offx_web/OFFX_AUTOMATION.md`](offx_web/OFFX_AUTOMATION.md#prerequisites)

## Common Workflows

### Workflow 1: Gene Target Research (API)
```bash
# Use: cortellis_api/scripts/cortellis_gene_query.py
# Get comprehensive gene data
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_api/scripts/cortellis_gene_query.py \
  EGFR KRAS --fields drug biomarker --excel

# Output: JSON + Markdown summary + Excel files
# Location: Current working directory
```

### Workflow 2: Clinical Data Export (Web)
```bash
# Use: cortellis_targetdrug_web/run-cortellis-download.sh
# Download clinical studies for multiple drugs
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh \
  --queries "drug1,drug2,drug3" --categories "Clinical Studies"

# Output: Excel files
# Location: cortellis_playwright_result/ in working directory
```

### Workflow 3: Multi-Category Analysis (Web)
```bash
# Use: cortellis_targetdrug_web/run-cortellis-download.sh
# Get all data types for comprehensive analysis
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh \
  imatinib --categories "Drugs & Biologics,Clinical Studies,Patents,Literature"

# Output: 4 Excel files, one per category
# Location: cortellis_playwright_result/{query}_{categories}_{datetime}/
```

### Workflow 4: Single Category Export (Web)
```bash
# Use: cortellis_targetdrug_web/run-cortellis.sh
# Export one category with detailed UI interaction
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis.sh \
  erlotinib --category "Clinical Studies"

# Output: Excel file with screenshots
# Location: cortellis_playwright_result/{query}_{category}_{datetime}/
```

### Workflow 5: API Timeout Handling
```bash
# Use: cortellis_api/scripts/cortellis_gene_query.py + convert_excel_to_json.py
# When API times out for large datasets
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_api/scripts/cortellis_gene_query.py JAK1 --fields drug
# If timeout occurs, download Excel manually from Cortellis web interface

# Convert Excel to JSON format
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_api/scripts/convert_excel_to_json.py \
  Dec_18_2025_JAK1.xlsx

# Continue with analysis using JSON
```

## Output Formats

### API Access Output
For each queried gene:
- `{GENE}_cortellis_data.json` - Complete API response
- `{GENE}_summary.md` - Markdown summary
- `{GENE}_cortellis_data.xlsx` - Excel workbook (with --excel flag)

**Location:** Current working directory (or `--output-dir` if specified)

### Web Automation Output
```
working_directory/cortellis_playwright_result/
└── {query}_{categories}_{datetime}/
    ├── search_overview.png                   # Search results page
    ├── {Category}_page.png                   # Category results page
    ├── {Category}_after_button_click.png     # Export menu opened
    ├── {Category}_Dec_19_2025_{query}.xlsx   # Downloaded Excel file
    └── metadata.json                         # Download metadata
```

## Comparison: API vs Web Automation

| Feature | API Access | Web Automation |
|---------|-----------|----------------|
| **Data Format** | JSON (+ Excel/MD via --excel) | Excel/CSV |
| **Speed** | Very Fast (~18s for 1700+ drugs) | Slower (~1-3 min/query) |
| **Performance** | 11.8x faster with async | Sequential downloads |
| **Bulk Downloads** | Excellent (with Excel fallback) | Excellent |
| **Category Selection** | No | Yes (13 categories) |
| **File Downloads** | Excel output via --excel flag | Native Excel/CSV downloads |
| **Rate Limits** | API limits (adaptive handling) | Web limits |
| **Setup** | API credentials | Okta session |
| **Scripts** | `cortellis_gene_query.py` | `run-cortellis.sh` / `run-cortellis-download.sh` |
| **Best Use Case** | Data analysis, batch queries | Multi-category file exports |

## Available Categories (Web Automation Only)

The following categories are available for web automation downloads:

- **Drugs & Biologics** (default)
- **Clinical Studies** (most common for research)
- **Patents** (IP landscape analysis)
- **Literature** (publications)
- Genes & Targets
- Organic Synthesis
- Experimental Pharmacology
- Experimental Models
- Pharmacokinetics
- Drug Metabolism
- Drug-Drug Interactions
- Organizations
- Disease Briefings

## Quick Reference by Script

### API Access Scripts

**Main Query Script:** `cortellis_api/scripts/cortellis_gene_query.py`
```bash
# Single gene
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_api/scripts/cortellis_gene_query.py \
  BRCA1 --fields drug biomarker --excel

# Multiple genes
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_api/scripts/cortellis_gene_query.py \
  BRCA1 TP53 EGFR --fields drug --excel
```

**Excel Converter Script:** `cortellis_api/scripts/convert_excel_to_json.py`
```bash
# Convert downloaded Excel to JSON
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_api/scripts/convert_excel_to_json.py \
  downloaded_file.xlsx GENE_NAME
```

### Web Automation Scripts

**Single-Category Export:** `cortellis_targetdrug_web/run-cortellis.sh`
```bash
# Default category (first available)
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis.sh imatinib

# Specific category
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis.sh \
  imatinib --category "Clinical Studies"
```

**Multi-Category Download:** `cortellis_targetdrug_web/run-cortellis-download.sh`
```bash
# Multiple categories
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh \
  imatinib --categories "Clinical Studies,Patents"

# Multiple queries
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh \
  --queries "drug1,drug2" --categories "Clinical Studies"

# Config file
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh \
  --config my_config.json
```

## Detailed Documentation

### API Access
- **Main Documentation:** [`cortellis_api/API_ACCESS.md`](cortellis_api/API_ACCESS.md) - Complete API guide
- **API Reference:** [`cortellis_api/references/api_reference.md`](cortellis_api/references/api_reference.md) - Endpoints, authentication, and performance optimization
- **Field Descriptions:** [`cortellis_api/references/api_fields.md`](cortellis_api/references/api_fields.md) - Available data fields
- **JSON Schema:** [`cortellis_api/references/json_schema.md`](cortellis_api/references/json_schema.md) - Data structure guide (API v2.0)
- **Scoring Framework:** [`cortellis_api/references/scoring_framework.md`](cortellis_api/references/scoring_framework.md) - Target analysis methodology
- **Excel Conversion Guide:** [`cortellis_api/examples/EXCEL_TO_JSON_CONVERSION_GUIDE.md`](cortellis_api/examples/EXCEL_TO_JSON_CONVERSION_GUIDE.md) - Comprehensive Excel to JSON conversion guide
- **Examples:** [`cortellis_api/examples/`](cortellis_api/examples/) - Complete analysis workflows

### Web Automation
- **Main Documentation:** [`cortellis_targetdrug_web/WEB_AUTOMATION.md`](cortellis_targetdrug_web/WEB_AUTOMATION.md) - Complete web automation guide
- **Setup Guide:** [`cortellis_targetdrug_web/SETUP_COMPLETE.md`](cortellis_targetdrug_web/SETUP_COMPLETE.md) - Authentication setup with test results
- **Quick Reference:** [`cortellis_targetdrug_web/README.md`](cortellis_targetdrug_web/README.md) - Command examples and troubleshooting
- **Example Configs:** [`cortellis_targetdrug_web/example/`](cortellis_targetdrug_web/example/) - Config file templates

## Support

### Re-authentication
- **API:** Update `.env` file when credentials change
- **Web:** Rerun `okta_auth_setup.py` when session expires (after weeks/months)
  ```bash
  cd /your/working/directory
  python /home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/okta_auth_setup.py
  ```

### Troubleshooting
- **API:** Check error messages for credential issues or network errors
- **Web:** Check output screenshots in result directories (`{category}_page.png`, `{category}_after_button_click.png`)
- Review metadata JSON files for error details
- See detailed documentation for specific issues

## Validation

Both methods have been comprehensively tested and validated:

**API Access:**
- Tested with multiple genes
- All output formats verified (JSON, Markdown, Excel)
- Timeout handling with Excel conversion validated

**Web Automation:**
- **Test Results:** 100% success rate across 4 test scenarios (7/7 downloads successful)
- **Performance:** Averages 2-3 minutes per query, 1.5 minutes per category
- **Features Validated:**
  - Single and multiple query processing
  - Single and multiple category downloads
  - "..." menu button export workflow
  - Environment-agnostic execution
  - Screenshot capture and metadata generation

---

**For detailed instructions, examples, and reference material, see the subskill documentation:**
- **API Access Details:** [`cortellis_api/API_ACCESS.md`](cortellis_api/API_ACCESS.md)
- **Web Automation Details:** [`cortellis_targetdrug_web/WEB_AUTOMATION.md`](cortellis_targetdrug_web/WEB_AUTOMATION.md)
- **Complete Navigation:** [`DOCUMENTATION_MAP.md`](DOCUMENTATION_MAP.md)
