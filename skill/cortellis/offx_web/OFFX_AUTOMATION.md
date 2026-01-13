# OFF-X Web Automation - Detailed Guide

Comprehensive guide for browser-based automation to search, evaluate target and drug safety, identify adverse events, and download Master view data from OFF-X database using Playwright.

## Table of Contents

- [Overview](#overview)
- [Purpose](#purpose)
- [Prerequisites](#prerequisites)
- [When to Use OFF-X Automation](#when-to-use-off-x-automation)
- [Core Scripts](#core-scripts)
- [Available Fields](#available-fields)
- [Workflow Examples](#workflow-examples)
- [Output Structure](#output-structure)
- [Error Handling & Troubleshooting](#error-handling--troubleshooting)
- [Configuration Files](#configuration-files)
- [Common Workflows](#common-workflows)
- [Validation & Testing](#validation--testing)
- [Support & Maintenance](#support--maintenance)

## Overview

OFF-X web automation uses Playwright browser automation to:
1. Authenticate via Okta
2. Navigate to OFF-X (TargetSafety.info) website
3. Select field (Targets, Drugs and biologics, Drug combinations, Adverse events)
4. Search for specified entities
5. Process all dropdown matches (or exact match only)
6. Navigate to Target safety profile → Master view
7. Export and download Master view Excel files

## Purpose

**OFF-X Web Automation** is designed to evaluate target and drug safety and identify adverse events using the OFF-X database (TargetSafety.info) via browser automation. This tool enables:

- **Target Safety Assessment**: Evaluate safety profiles of drug targets
- **Drug Safety Analysis**: Assess safety data for drugs and biologics
- **Adverse Event Monitoring**: Identify and analyze adverse events
- **Master View Exports**: Download comprehensive Master view tables

## Prerequisites

**IMPORTANT: Okta Authentication Required**

OFF-X automation requires a saved Okta session at `~/.okta/auth_state.json`. Use the centralized `ai-sci:okta-sso` skill for authentication.

### Initial Setup

1. **Check if you have a valid session:**
```bash
~/ai-sci-claude-skills/ai-sci/skills/okta-sso/run-okta-login.sh --status
```

2. **If no valid session, authenticate using the okta-sso skill:**

   **Note:** Claude cannot handle secrets interactively. Run this command yourself in your terminal:
   ```bash
   OKTA_EMAIL="your.email@takeda.com" OKTA_PASSWORD="your-password" \
     ~/ai-sci-claude-skills/ai-sci/skills/okta-sso/run-okta-login.sh
   ```

3. **Complete MFA**: Approve the push notification on your phone (match the verification number displayed)

4. **Session saved**: Creates `~/.okta/auth_state.json`
   - Valid for weeks/months
   - Centralized location shared by all Cortellis and OFF-X scripts
   - Rerun okta-sso skill when session expires

## When to Use OFF-X Automation

Use OFF-X automation when users request:
- Target safety profile evaluation
- Drug safety assessment
- Adverse event analysis
- Master view data export
- Safety data for specific targets, drugs, or conditions
- Comprehensive safety tables in Excel format
- Bulk safety assessments for multiple entities

## Core Scripts

### Script: `offx-download.js` - Multi-Entity Safety Downloads

Download Master view safety data for multiple entities across different fields with iterative processing. Best for bulk downloads and comprehensive safety assessments.

**Wrapper:** `run-offx-download.sh`

**Location:** `/home/sagemaker-user/.claude/skills/cortellis/offx_web/`

**Usage:**
```bash
# Run from any directory with okta_auth_state.json
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh [OPTIONS]
```

**Options:**

1. **Single entity, default field (Targets), all matches:**
```bash
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh itga4
```

2. **Single entity, specific field, all matches:**
```bash
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh imatinib --field "Drugs and biologics"
```

3. **Single entity, exact match only:**
```bash
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh egfr --exact
```

4. **Multiple entities, default field:**
```bash
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh --queries "tyk2,jak1,itga4"
```

5. **Multiple entities, specific field:**
```bash
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh --queries "imatinib,gefitinib" --field "Drugs and biologics"
```

6. **Config file (most flexible):**
```bash
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh --config my_safety_assessment.json
```

**Export Workflow (per entity match):**
1. Select field (Targets, Drugs and biologics, etc.)
2. Enter entity name in search box
3. Capture all dropdown matches (ignore "select multiple targets" and ML/Pathway results)
4. For each match:
   a. Navigate to match page
   b. Click "Target safety profile" (left side)
   c. Select "Master view" from dropdown
   d. Click Export icon (far right of table)
   e. Download Master view Excel file
5. Return to home page and repeat for next match/entity

**Config File Format:**
```json
[
  {
    "entity": "tyk2",
    "field": "Targets",
    "exactMatch": false
  },
  {
    "entity": "imatinib",
    "field": "Drugs and biologics",
    "exactMatch": false
  },
  {
    "entity": "cardiac failure",
    "field": "Adverse events",
    "exactMatch": true
  }
]
```

**Output Location:**
```
working_directory/offx_playwright_result/
└── {entity}_{field}_{datetime}/
    ├── dropdown_matches.png
    ├── {match1}_page.png
    ├── {match1}_after_safety_profile.png
    ├── {match1}_master_view.png
    ├── OFFX targets {match1} master-view {datetime}.xlsx
    ├── {match2}_page.png
    ├── {match2}_after_safety_profile.png
    ├── {match2}_master_view.png
    ├── OFFX targets {match2} master-view {datetime}.xlsx
    └── metadata.json
```

## Available Fields

The following fields are available for OFF-X safety assessment:

- **Targets** (default) - Evaluate target safety profiles
- **Drugs and biologics** - Assess drug and biologic safety
- **Drug combinations** - Analyze combination therapy safety
- **Adverse events** - Investigate specific adverse events

## Workflow Examples

### Example 1: Target Safety Assessment (All Matches)

Evaluate safety for all ITGA4-related targets:

```bash
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh itga4
```

**Output:**
- 5 Excel files (one for each ITGA4 match: alpha4beta1, alpha4beta7, alpha 4 subunit, etc.)
- Screenshots at each step
- Metadata JSON

### Example 2: Drug Safety Assessment (Exact Match)

Get safety data for specific drug:

```bash
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh imatinib --field "Drugs and biologics" --exact
```

**Output:**
- 1 Excel file (exact match only)
- Master view table with safety data

### Example 3: Multiple Target Safety Assessment

Bulk safety assessment for multiple targets:

```bash
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh --queries "tyk2,jak1,egfr"
```

**Output:**
- Multiple Excel files for each target
- All dropdown matches processed

### Example 4: Adverse Event Analysis

Investigate specific adverse events:

```bash
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh "cardiac failure" --field "Adverse events"
```

**Output:**
- Master view Excel files for cardiac failure-related entries

### Example 5: Comprehensive Safety Config

Use config file for mixed field analysis:

```bash
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh --config comprehensive_safety.json
```

**Config content:**
```json
[
  {
    "entity": "egfr",
    "field": "Targets",
    "exactMatch": false
  },
  {
    "entity": "erlotinib",
    "field": "Drugs and biologics",
    "exactMatch": false
  },
  {
    "entity": "rash",
    "field": "Adverse events",
    "exactMatch": false
  }
]
```

## Safety Score Analysis Workflow

### Overview

After downloading Master view Excel files, analyze safety scores using the "OFF-X Target/Class Score label" column from the "Data" sheet. A comprehensive analysis script is provided for automated scoring and report generation.

### Scoring Methodology

**Negative Scoring System** (Lower/more negative = worse safety):

| Severity Label | Score | Interpretation |
|----------------|-------|----------------|
| **Very high** | **-10** | Most severe adverse events |
| **High** | **-8** | Serious adverse events |
| **Medium** | **-6** | Moderate adverse events |
| **Low** | **-2** | Minor adverse events |
| **Very low** | **-1** | Minimal adverse events |
| **Not associated** | **+2** | No causal association |
| **Empty/NA** | **0** | No scoring data |

**Key Metrics:**
- **Average Score**: Sum of all scores / Total adverse events
- **Medium-High %**: (Very high + High + Medium) / Total AEs × 100
- **NA %**: Empty/NA cases / Total AEs × 100

### Analysis Script

**Location:** `/home/sagemaker-user/.claude/skills/cortellis/offx_web/analyze_offx_safety.py`

**Features:**
- Reads metadata.json files to accurately extract gene names from search terms
- Individual gene and combination safety analysis
- Automated markdown report generation with rankings
- JSON export for programmatic use
- Severity breakdowns and clinical recommendations
- Handles multiple files per gene (isoforms, variants, etc.)

**Basic Usage:**

```bash
# Run from directory containing offx_playwright_result/ folder
cd /your/working/directory
python3 /home/sagemaker-user/.claude/skills/cortellis/offx_web/analyze_offx_safety.py
```

**Output Files:**
- `offx_safety_analysis_report.md` - Comprehensive markdown report with:
  - Gene safety rankings (best to worst)
  - Detailed individual gene analysis
  - Combination/list analysis for drug development
  - Clinical validation status and competitive intelligence
- `offx_safety_analysis_data.json` - Raw data for further analysis

**Script Features:**
- **Metadata-based gene extraction**: Reads entity names from metadata.json files for accurate gene identification
- **Score aggregation**: Combines multiple files per gene (e.g., IL17A across multiple isoforms/variants)
- **Combination analysis**: Automatically analyzes pre-defined gene lists
- **Clinical interpretation**: Provides actionable insights for target selection and development strategy

**Key Outputs:**
- Average safety score per gene
- % Medium-High severity adverse events
- % NA/uncharacterized events
- Total adverse events analyzed
- Interpretation and recommendations

### Interpreting Results

**Score Ranges:**
- **< -1.0**: High concern (significant safety liabilities)
- **-1.0 to -0.5**: Moderate concern (notable adverse events)
- **-0.5 to 0**: Low concern (manageable safety profile)
- **> 0**: Excellent (minimal safety concerns)

**Important:** Score of 0.00 is ambiguous:
- Could mean 100% NA (uncharacterized) ⚠️ OR perfect safety ✅
- **Always check NA % to distinguish**

**Examples:**
- JAK1: -1.17 | 4.0% Med-High | 9.0% NA → Significant concerns, well-characterized
- GREM1: -0.23 | 0.0% Med-High | 76.9% NA → Minimal concerns, limited data
- CDKN2D: 0.00 | 0.0% Med-High | 100% NA → Uncharacterized (unknown, not safe)

### Known Issues & Fixes

**Dropdown Detection Bug (Fixed Dec 2024):**

Original script incorrectly matched `[class*="autocomplete"]` (0 items) before `[role="listbox"]` (actual dropdown), causing "No matches found" for valid genes.

**Fix:** Modified `getDropdownMatches()` (lines 268-334 in offx-download.js) to:
- Wait for `[role="listbox"] [role="option"]` to exist
- Directly query `page.$$('[role="option"]')` bypassing parent containers
- Increased wait from 5s to 8s for dropdown loading

**Result:** All genes (GREM1, ITGA4, ITGB7, etc.) now correctly detected and processed

**Exact Matching Feature (Added Jan 2026):**

The script now implements exact matching for dropdown items to ensure only relevant targets are processed.

**Matching Logic:**
- Extracts text before "[" bracket in dropdown items
- Normalizes both search term and dropdown text by removing spaces, hyphens, and underscores
- Performs case-insensitive comparison
- Only processes exact matches

**Examples:**
- Searching "IL17F" matches: "IL17F [...]", "IL-17F [...]", "IL 17F [...]"
- Searching "IL17F" excludes: "IL17A [...]", "IL-17 [...]", "INTERLEUKIN 17F [...]"
- Searching "TNFSF13" matches: "TNFSF13 [BAFF]" but excludes "TNFSF13B [APRIL]"
- Searching "CD3G" matches: "CD3G [...]", "CD3g [...]" but excludes "CD3D [...]", "CD3E [...]"

**Implementation:** Modified `getDropdownMatches()` function with filtering logic after dropdown detection:
```javascript
const entityNormalized = entity.replace(/[\s\-_]+/g, '').toLowerCase();
const exactMatches = dropdownMatches.filter(match => {
  const textBeforeBracket = match.text.split('[')[0].trim();
  const textNormalized = textBeforeBracket.replace(/[\s\-_]+/g, '').toLowerCase();
  return textNormalized === entityNormalized;
});
```

**Cookie Overlay Fix (Added Jan 2026):**

Added cookie consent overlay handling to prevent timeout errors during dropdown clicks.

**Fix:** Uses `{ force: true }` option for dropdown clicks to bypass overlays, and adds `dismissCookieOverlay()` helper function for explicit overlay dismissal when needed.

**Result:** Eliminates "OneTrust privacy center overlay intercepts pointer events" errors

## Output Structure

### Main Result Directory

```
working_directory/offx_playwright_result/
├── {entity}_{field}_{datetime}/
│   ├── dropdown_matches.png              # Dropdown screenshot
│   ├── {match1}_page.png                 # Entity page
│   ├── {match1}_after_safety_profile.png # After clicking safety profile
│   ├── {match1}_master_view.png          # Master view table
│   ├── OFFX {field} {match1} master-view {datetime}.xlsx  # Excel export
│   ├── {match2}_page.png
│   ├── {match2}_after_safety_profile.png
│   ├── {match2}_master_view.png
│   ├── OFFX {field} {match2} master-view {datetime}.xlsx
│   └── metadata.json                     # Download metadata
├── summary_{datetime}.json               # Overall summary
└── ...
```

### Metadata JSON Structure

```json
{
  "entity": "itga4",
  "field": "Targets",
  "exactMatch": false,
  "downloads": [
    {
      "match": "ITGA4 [Integrin alpha4beta1]",
      "success": true,
      "filePath": "/path/to/OFFX targets Integrin alpha4beta1 master-view.xlsx",
      "fileName": "OFFX targets Integrin alpha4beta1 master-view.xlsx",
      "fileSize": 17155,
      "fileSizeMB": "0.02"
    },
    {
      "match": "ITGA4 [Integrin alpha4beta7]",
      "success": true,
      "filePath": "/path/to/OFFX targets Integrin alpha4beta7 master-view.xlsx",
      "fileName": "OFFX targets Integrin alpha4beta7 master-view.xlsx",
      "fileSize": 93839,
      "fileSizeMB": "0.09"
    }
  ],
  "timestamp": "2025-12-21T03:51:28.929Z"
}
```

## Error Handling & Troubleshooting

### Common Issues

1. **Authentication Expired**
   - **Symptom:** Script redirects to login page
   - **Solution:** Rerun `okta_auth_setup.py` to refresh session

2. **Entity Not Found**
   - **Symptom:** "No matches found in dropdown"
   - **Solution:** Verify entity name spelling, try alternative names

3. **Export Failed**
   - **Symptom:** "Export button not found"
   - **Solution:** Check screenshots in output directory, page structure may have changed

4. **Headless Mode Issues**
   - **Symptom:** "Missing X server" error
   - **Solution:** Script defaults to headless mode; if needed, edit `offx-download.js` line 723 to set `headless: false`

### Debugging

All runs generate detailed screenshots at each step:
- `dropdown_matches.png` - Search results
- `{match}_page.png` - Entity page
- `{match}_after_safety_profile.png` - Safety profile dropdown
- `{match}_master_view.png` - Master view table

Check these screenshots to identify where the automation failed.

## Configuration Files

### Example Config Files

Located in `/home/sagemaker-user/.claude/skills/cortellis/offx_web/example/`:

1. **example_targets.json** - Target safety assessment
2. **example_drugs_biologics.json** - Drug safety assessment
3. **example_adverse_events.json** - Adverse event analysis
4. **example_comprehensive.json** - Mixed field analysis

### Creating Custom Configs

```json
[
  {
    "entity": "your_entity_name",
    "field": "Targets|Drugs and biologics|Drug combinations|Adverse events",
    "exactMatch": false  // true for exact match only, false for all matches
  }
]
```

## Common Workflows

### Workflow 1: Single Target Safety Profile

```bash
cd /your/working/directory
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh itga4
```

**Use case:** Quick safety assessment for one target with all isoforms/variants

### Workflow 2: Bulk Drug Safety Assessment

```bash
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh \
  --queries "imatinib,gefitinib,erlotinib" \
  --field "Drugs and biologics"
```

**Use case:** Comparative safety analysis across multiple drugs

### Workflow 3: Adverse Event Investigation

```bash
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh \
  --config example/example_adverse_events.json
```

**Use case:** Systematic adverse event analysis for safety monitoring

### Workflow 4: Comprehensive Safety Report

```bash
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh \
  --config my_comprehensive_safety.json
```

**Config includes:**
- Multiple targets
- Related drugs
- Known adverse events

**Use case:** Complete safety dossier for regulatory submission or safety review

## Validation & Testing

### Test Results

**Test Entity:** itga4 (Targets field)

**Test Results:**
- ✅ 5/5 downloads successful (100% success rate)
- ✅ All dropdown matches processed
- ✅ Master view exports completed
- ✅ Metadata generated correctly

**Downloaded Files:**
1. ITGA4 [Integrin alpha4beta1] - 17 KB
2. ITGA4 [Integrin alpha4beta7] - 92 KB
3. ITGA4 [Integrin, alpha 4 subunit] - 71 KB
4. ITGA4L [Integrin, alpha 9 subunit] - 11 KB
5. ITGA4 [Integrin, alpha 4/Paxillin] - 11 KB

**Performance:**
- Average time: 1-2 minutes per entity
- Per match: 30-45 seconds

### Features Validated

- ✅ Okta authentication
- ✅ OFF-X application launch
- ✅ Field selection (Targets tested)
- ✅ Entity search with autocomplete
- ✅ Dropdown match detection (ignoring ML/Pathway results)
- ✅ Iterative match processing
- ✅ Target safety profile navigation
- ✅ Master view selection
- ✅ Export icon detection and click
- ✅ File download and naming
- ✅ Screenshot capture at each step
- ✅ Metadata generation
- ✅ Summary JSON creation

## Support & Maintenance

### Re-authentication

When Okta session expires (after weeks/months):

```bash
cd /your/working/directory
python /home/sagemaker-user/.claude/skills/cortellis/offx_web/okta_auth_setup.py
```

### Script Location

All scripts are located in:
```
/home/sagemaker-user/.claude/skills/cortellis/offx_web/
```

**Main Scripts:**
- `run-offx-download.sh` - Bash wrapper
- `offx-download.js` - Playwright automation
- `okta_auth_setup.py` - Authentication setup

**Example Configs:**
- `example/example_targets.json`
- `example/example_drugs_biologics.json`
- `example/example_adverse_events.json`
- `example/example_comprehensive.json`

### Environment Requirements

- Playwright installed (automatic via skill setup)
- Node.js (for Playwright execution)
- Python 3.6+ (for auth setup)
- Valid Okta credentials

---

**For parent skill documentation, see:**
- [`../SKILL.md`](../SKILL.md) - Cortellis parent skill overview
- [`../DOCUMENTATION_MAP.md`](../DOCUMENTATION_MAP.md) - Complete documentation map
