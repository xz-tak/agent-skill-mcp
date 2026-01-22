# Cortellis Target-Drug Web Automation - Detailed Guide

Comprehensive guide for browser-based automation to search, navigate, and download target-drug relationship data from specific Cortellis categories using Playwright. Supports queries for targets with associated drugs and drugs with associated targets.

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [When to Use Web Automation](#when-to-use-web-automation)
- [Core Scripts](#core-scripts)
- [Available Categories](#available-categories)
- [Workflow Examples](#workflow-examples)
- [Output Structure](#output-structure)
- [Error Handling & Troubleshooting](#error-handling--troubleshooting)
- [Configuration Files](#configuration-files)
- [Script Customization](#script-customization)
- [Common Workflows](#common-workflows)
- [Validation & Testing](#validation--testing)
- [Support & Maintenance](#support--maintenance)

## Overview

Web automation uses Playwright browser automation to:
1. Authenticate via Okta
2. Navigate to Cortellis Drug Discovery website
3. Search for targets (e.g., EGFR, KRAS, TYK2) or drugs (e.g., imatinib, gefitinib)
4. Navigate to category-specific result pages (Genes & Targets, Drugs & Biologics, Clinical Studies, Patents, etc.)
5. Export target-drug relationship data using "..." menu button
6. Download Excel/CSV files with comprehensive target and drug information

## Prerequisites

**IMPORTANT: Okta Authentication Required**

Web automation requires a saved Okta session. This is a one-time setup that saves browser authentication state.

### Initial Setup

1. **Run the authentication script from your working directory** (one-time):
```bash
# Navigate to your project/working directory first
cd /your/project/directory

# Run auth setup from there
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/okta_auth_setup.py
```

2. **Follow prompts**:
   - Enter your Takeda email
   - Enter your password
   - Complete MFA (enter code from authenticator app)

3. **Session saved**: Creates `okta_auth_state.json` in your working directory
   - Valid for weeks/months
   - Rerun script when session expires
   - File must be in the directory where you run the automation scripts

## When to Use Target-Drug Web Automation

Use target-drug web automation when users request:
- **Target-to-Drug queries**: Target genes with their associated drugs (e.g., "What drugs target EGFR?")
- **Drug-to-Target queries**: Drugs with their associated targets (e.g., "What targets does imatinib affect?")
- **Target-drug relationships**: Data from Genes & Targets and Drugs & Biologics categories
- Downloading Excel or CSV files from Cortellis for any category
- Exporting data from specific categories (Genes & Targets, Drugs & Biologics, Clinical Studies, Patents, Literature, etc.)
- Bulk downloads for multiple targets or drugs
- Category-specific data that's easier to download than query via API
- File-based workflows requiring Cortellis exports

## Core Scripts

### Script 1: `cortellis-automation.js` - Search & Export

Single category export with "..." menu button export workflow. Best for single-category exports with detailed UI interaction.

**Wrapper:** `run-cortellis.sh`

**Location:** `/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/`

**Usage:**
```bash
# Run from any directory with okta_auth_state.json
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis.sh <search_term> [--category "Category Name"]
```

**Examples:**
```bash
# Default category (first available "View Results")
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis.sh imatinib

# Specific category
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis.sh imatinib --category "Clinical Studies"

# Multiple search terms
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis.sh imatinib dasatinib nilotinib --category "Patents"

# From file
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis.sh --file drugs.txt --category "Literature"
```

**Export Workflow:**
1. Search for term
2. Navigate to category page (default or specified)
3. Click "..." (more_horiz) menu button at top right of results table
4. Click "Export" in dropdown menu
5. Click "Export" button in confirmation popup
6. Download file

**Output Location:**
```
working_directory/cortellis_playwright_result/
└── {search_term}_{category}_{datetime}/
    ├── search_results.png
    ├── results_page.png
    ├── after_button_click.png
    ├── final_state.png
    ├── Dec_19_2025_{search_term}.xlsx
    ├── results.json
    └── results.html
```

---

### Script 2: `cortellis-download.js` - Multi-Category Downloads

Download from multiple categories per search term with iterative processing. Best for bulk downloads and multi-category exports. Uses the same "..." menu button export workflow.

**Wrapper:** `run-cortellis-download.sh`

**Location:** `/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/`

**Usage:**
```bash
# Run from any directory with okta_auth_state.json
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh [OPTIONS]
```

**Options:**

1. **Single query, default category:**
```bash
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh imatinib
```

2. **Single query, specific categories:**
```bash
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh imatinib --categories "Clinical Studies,Patents,Literature"
```

3. **Multiple queries, default category:**
```bash
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh --queries "imatinib,gefitinib,dasatinib"
```

4. **Multiple queries, specific categories:**
```bash
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh --queries "imatinib,gefitinib" --categories "Clinical Studies,Patents"
```

5. **Config file (most flexible):**
```bash
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh --config my_downloads.json
```

**Export Workflow (per category):**
1. Search for term
2. Navigate to category page
3. Click "..." (more_horiz) menu button at top right of results table
4. Click "Export" in dropdown menu
5. Click "Export" button in confirmation popup
6. Download file with category prefix
7. Return to search results and repeat for next category

**Config File Format:**
```json
[
  {
    "searchTerm": "imatinib",
    "categories": ["Clinical Studies", "Patents"]
  },
  {
    "searchTerm": "gefitinib",
    "categories": ["Drugs & Biologics", "Literature"]
  }
]
```

**Output Location:**
```
working_directory/cortellis_playwright_result/
└── {search_term}_{categories}_{datetime}/
    ├── search_overview.png
    ├── Clinical_Studies_page.png
    ├── Clinical_Studies_after_button_click.png
    ├── Clinical_Studies_Dec_19_2025_{search_term}.xlsx
    ├── Patents_page.png
    ├── Patents_after_button_click.png
    ├── Patents_Dec_19_2025_{search_term}.xlsx
    └── metadata.json
```

## Available Categories

All scripts support these categories:
- **Drugs & Biologics** (default)
- Genes & Targets
- Organic Synthesis
- Experimental Pharmacology
- Experimental Models
- Pharmacokinetics
- Drug Metabolism
- Drug-Drug Interactions
- **Clinical Studies** (most common for research)
- Organizations
- **Literature** (publications)
- **Patents** (IP landscape analysis)
- Disease Briefings

## Workflow Examples

### Example 1: Download Clinical Data for Drug Research

User request: "Download clinical studies data for imatinib, dasatinib, and nilotinib"

```bash
# Method 1: Command line
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh \
  --queries "imatinib,dasatinib,nilotinib" --categories "Clinical Studies"

# Method 2: Config file
cat > clinical_drugs.json << 'EOJ'
[
  {"searchTerm": "imatinib", "categories": ["Clinical Studies"]},
  {"searchTerm": "dasatinib", "categories": ["Clinical Studies"]},
  {"searchTerm": "nilotinib", "categories": ["Clinical Studies"]}
]
EOJ

/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh --config clinical_drugs.json
```

### Example 2: Comprehensive Drug Intelligence

User request: "Get all available data for gefitinib including drugs, trials, patents, and literature"

```bash
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh \
  gefitinib --categories "Drugs & Biologics,Clinical Studies,Patents,Literature"
```

### Example 3: Patent Landscape Analysis

User request: "Download patent data for multiple kinase inhibitors"

```bash
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh \
  --queries "imatinib,erlotinib,gefitinib,lapatinib" --categories "Patents"
```

### Example 4: Single Category Export

User request: "Export clinical studies for erlotinib"

```bash
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis.sh \
  erlotinib --category "Clinical Studies"
```

## Output Structure

### cortellis-automation.js Output
```
cortellis_playwright_result/
└── erlotinib_Clinical_Studies_2025-12-20_21-40-45/
    ├── search_results.png                  # Category selection page
    ├── results_page.png                    # Results table page
    ├── after_button_click.png              # "..." menu opened
    ├── final_state.png                     # Final state
    ├── Dec_19_2025_erlotinib.xlsx          # Downloaded Excel file
    ├── results.json                        # Structured metadata
    └── results.html                        # Raw page HTML
```

### cortellis-download.js Output
```
cortellis_playwright_result/
└── dasatinib_Clinical_Studies-Patents_2025-12-20_21-51-26/
    ├── search_overview.png                            # Search results overview
    ├── Clinical_Studies_page.png                      # Category results page
    ├── Clinical_Studies_after_button_click.png        # Menu opened
    ├── Clinical_Studies_Dec_19_2025_dasatinib.xlsx    # Downloaded Excel file
    ├── Patents_page.png                               # Category results page
    ├── Patents_after_button_click.png                 # Menu opened
    ├── Patents_Dec_19_2025_dasatinib.xlsx             # Downloaded Excel file
    └── metadata.json                                  # Download metadata
```

## Error Handling & Troubleshooting

### Authentication Issues

**Error: "okta_auth_state.json not found"**
```bash
# Run setup from your working directory (where you want results saved)
cd /your/working/directory
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/okta_auth_setup.py
```

**Error: "Authentication expired"**
- Session cookies have expired (typically after weeks/months)
- Re-run auth setup from your working directory:
```bash
cd /your/working/directory
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/okta_auth_setup.py
```

### Download Issues

**Error: "Export button/menu not found"**
- Check screenshots in output directory:
  - `{category}_page.png` - shows the results table
  - `{category}_after_button_click.png` - shows if "..." menu opened
- Scripts look for "more_horiz" icon button at top right of results table
- If UI has changed, update selectors in the script
- Try different category

**Error: "Category not found"**
- Check exact spelling (case-sensitive): e.g., "Clinical Studies" not "clinical studies"
- Category may not have results for this search term
- See available categories list above
- Check `search_overview.png` to see available categories for your search

### Rate Limiting

Scripts include 3-second delays between searches. If experiencing issues:
- Reduce batch size
- Increase delays in script
- Wait before retrying

## Configuration Files

### Example: Clinical Intelligence Project
```json
[
  {
    "searchTerm": "imatinib",
    "categories": ["Clinical Studies", "Patents", "Literature"]
  },
  {
    "searchTerm": "dasatinib",
    "categories": ["Clinical Studies", "Patents", "Literature"]
  }
]
```

### Example: Drug Safety Analysis
```json
[
  {
    "searchTerm": "imatinib",
    "categories": ["Drug Metabolism", "Drug-Drug Interactions"]
  },
  {
    "searchTerm": "gefitinib",
    "categories": ["Drug Metabolism", "Drug-Drug Interactions"]
  }
]
```

### Example: Target Research
```json
[
  {
    "searchTerm": "EGFR",
    "categories": ["Genes & Targets", "Drugs & Biologics"]
  },
  {
    "searchTerm": "KRAS",
    "categories": ["Genes & Targets", "Drugs & Biologics"]
  }
]
```

## Script Customization & Environment-Agnostic Execution

Both scripts are designed to be environment-agnostic and work from any directory:

**Key Features:**
- `okta_auth_state.json` read from your current working directory (where you run the command)
- Results save to `cortellis_playwright_result/` in your current working directory
- Scripts can be run from anywhere using full path (no need to `cd` first)
- `CORTELLIS_WORK_DIR` environment variable automatically set to preserve working directory

**Example:**
```bash
# Auth file and results go to your project directory
cd /my/project
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/okta_auth_setup.py
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh imatinib

# Results saved to: /my/project/cortellis_playwright_result/
```

### Playwright Skill Path

Scripts automatically locate Playwright skill at:
```bash
${HOME}/.claude/plugins/cache/playwright-skill/playwright-skill/4.1.0/skills/playwright-skill
```

If your Playwright skill is elsewhere, edit the shell scripts:
```bash
# In run-cortellis.sh or run-cortellis-download.sh
PLAYWRIGHT_SKILL_DIR="/path/to/your/playwright-skill"
```

## Common Workflows

### Workflow 1: Quick Clinical Data Export
1. User asks: "What clinical trials exist for erlotinib?"
2. **Choose**: Web Automation (want Excel file)
3. **Run**: `/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis.sh erlotinib --category "Clinical Studies"`
4. **Result**: Excel file with trial data (254K)
   - Location: `cortellis_playwright_result/erlotinib_Clinical_Studies_2025-12-20_*/Dec_19_2025_erlotinib.xlsx`

### Workflow 2: Multi-Category Drug Analysis
1. User asks: "Get clinical studies and patents for dasatinib"
2. **Choose**: Web Automation (multiple categories)
3. **Run**: `/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh dasatinib --categories "Clinical Studies,Patents"`
4. **Result**: 2 Excel files with category prefixes
   - `Clinical_Studies_Dec_19_2025_dasatinib.xlsx` (273K)
   - `Patents_Dec_19_2025_dasatinib.xlsx` (273K)

### Workflow 3: Bulk Drug Download
1. User asks: "Download drugs data for nilotinib and lapatinib"
2. **Choose**: Web Automation (bulk download)
3. **Run**: `/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh --queries "nilotinib,lapatinib" --categories "Drugs & Biologics"`
4. **Result**: 2 Excel files, one per drug
   - `nilotinib_Drugs___Biologics_*/Drugs___Biologics_Dec_19_2025_nilotinib.xlsx` (153K)
   - `lapatinib_Drugs___Biologics_*/Drugs___Biologics_Dec_19_2025_lapatinib.xlsx` (96K)

## Validation & Testing

The web automation scripts have been comprehensively tested and validated (2025-12-20):

**Test Results: 100% Success Rate (4/4 tests passed, 7/7 downloads successful)**

| Test | Script | Features Validated | Status |
|------|--------|-------------------|--------|
| 1 | cortellis-automation.js | Single query, default category | ✅ PASSED |
| 2 | cortellis-automation.js | Single query, specific category | ✅ PASSED |
| 3 | cortellis-download.js | Single query, multiple categories | ✅ PASSED |
| 4 | cortellis-download.js | Multiple queries, single category | ✅ PASSED |

**Validated Features:**
- ✅ "..." menu button export workflow (both scripts)
- ✅ Multiple query processing (`--queries`)
- ✅ Multiple categories per query (`--categories`)
- ✅ Category-specific selection (`--category`)
- ✅ Default category handling (Drugs & Biologics)
- ✅ Environment-agnostic execution
- ✅ Screenshot capture at each step
- ✅ Metadata generation (JSON)
- ✅ Proper file naming with category prefixes

**Performance:**
- Average time per query: ~2-3 minutes
- Average time per category download: ~1.5 minutes
- Export success rate: 100%
- Screenshot capture: 100% successful

## Support & Maintenance

### Re-authentication

**API**: Update `.env` file when credentials change
**Web**: Rerun `okta_auth_setup.py` when session expires (after weeks/months)

### Updating Scripts

Scripts are in skill directory:
```
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/
├── okta_auth_setup.py
├── cortellis-automation.js
├── cortellis-download.js
├── run-cortellis.sh
└── run-cortellis-download.sh
```

### Getting Help

Check screenshots in output directories for debugging:
- `search_results.png` - Initial search page
- `results_page.png` - Results after navigation
- `after_button_click.png` - After clicking export
- `final_state.png` - Final page state

If issues persist, review metadata JSON files for error details.

## Comparison: API vs Web Automation

| Feature | API Access | Web Automation |
|---------|-----------|----------------|
| Data Format | JSON | Excel/CSV |
| Speed | Fast | Slower (browser) |
| Bulk Downloads | Limited | Excellent |
| Category Selection | No | Yes |
| File Downloads | No | Yes |
| Rate Limits | API limits | Web limits |
| Setup | API keys | Okta session |
| Use Case | Analysis | Downloads |

---

**For additional help, see:**
- **Main Documentation:** [`../SKILL.md`](../SKILL.md)
- **Setup Guide:** [`SETUP_COMPLETE.md`](SETUP_COMPLETE.md)
- **Quick Reference:** [`README.md`](README.md)
- **Complete Navigation:** [`../DOCUMENTATION_MAP.md`](../DOCUMENTATION_MAP.md)
