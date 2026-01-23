---
name: ncbi-paper-query
description: Query NCBI PubMed for publications matching disease/tissue/organism criteria and extract omics accessions (GEO, SRA, ArrayExpress) with metadata enrichment. This skill should be used when users request literature searches for omics studies, gene expression datasets, or publication lists for specific diseases and tissues.
---

# NCBI Paper Query

## Overview

Query PubMed for disease/tissue-specific publications and extract omics dataset accessions (GEO, SRA, ArrayExpress) with full metadata enrichment. The workflow interviews the user to collect search parameters, confirms settings, then executes the query.

## Quick Start

**Run from your project directory** (outputs saved to `./output/`):

```bash
# Abstract-only mode (fastest, no full-text access needed)
PYTHONPATH=~/.claude/skills/ncbi-paper-query/scripts python -m ncbi_paper_query \
    --disease "Crohn's disease" \
    --tissue intestine colon \
    --organism human \
    --if-cutoff 7.0 \
    --year-cutoff 2020 \
    --max-results 100 \
    --abstract-only \
    --output my_study

# Web-scrape mode (extracts accessions from article HTML)
PYTHONPATH=~/.claude/skills/ncbi-paper-query/scripts python -m ncbi_paper_query \
    --disease "IBD" \
    --tissue intestine \
    --organism human \
    --max-results 50 \
    --web-scrape \
    --output ibd_webscrape

# Subscription-download mode (web-scrape free papers, download subscription only)
PYTHONPATH=~/.claude/skills/ncbi-paper-query/scripts python -m ncbi_paper_query \
    --disease "fibrosis" \
    --tissue lung liver \
    --organism human \
    --max-results 100 \
    --subscription-download \
    --output fibrosis_hybrid

# PDF download mode (default - downloads all PDFs)
PYTHONPATH=~/.claude/skills/ncbi-paper-query/scripts python -m ncbi_paper_query \
    --disease "ulcerative colitis" \
    --tissue colon \
    --organism human \
    --max-results 20 \
    --output uc_full
```

**Output Structure:**
```
./output/{study_name}/
├── results.csv        # Main results with paper metadata
├── accessions.csv     # Omics accessions (one row per accession)
├── failed_studies.csv # Papers that couldn't be accessed
└── downloads/         # Downloaded PDFs (PDF mode only)
```

## Skill Directory

The skill directory is **auto-detected** - no manual setup required. The script automatically finds resources in `~/.claude/skills/ncbi-paper-query`.

To override, set the `SKILL_DIR` environment variable before running.

**Skill Structure:**
```
ncbi-paper-query/
├── SKILL.md                    # This file
├── references/
│   └── parameters.md           # Parameter documentation
└── scripts/
    ├── extraction_template.py  # Extraction template for Codex
    └── ncbi_paper_query/       # Main Python package
        ├── __init__.py         # Package exports
        ├── __main__.py         # CLI entry point
        ├── config.py           # Configuration & logging
        ├── models.py           # Data models (Publication, OmicsDataset)
        ├── retrieval.py        # PubMedSearcher class
        ├── matching.py         # EntityMatcher class
        ├── data/
        │   └── journal_if.json # Impact factor reference
        ├── metadata/
        │   └── impact_factor.py # ImpactFactorLookup class
        ├── download/
        │   ├── downloader.py   # PaperDownloader class
        │   └── omics_extractor.py # OmicsExtractor class
        ├── export/
        │   └── exporter.py     # ResultsExporter class
        └── core/
            ├── retrieval.py    # retrieve_publications() function
            └── validation.py   # validate_results(), compare_rounds()
```

**Working Directory Files:**
```
{current_working_dir}/
├── .env                        # Credentials (TAK_ACCOUNT, TAK_KEY)
└── output/                     # Output directory (auto-created)
    └── {study_name}/           # Per-study output
        ├── results.csv
        ├── accessions.csv
        └── downloads/
```

## Workflow

Execute these phases in sequence:

### Phase 1: User Interview

Use `AskUserQuestion` to collect ALL required parameters. Always ask all questions - do not skip any.

**Question 1: Disease/Indication**
```
header: "Disease"
question: "What disease(s) or indication(s) are you searching for? (Comma-separated for multiple, uses OR logic)"
options:
  - label: "IBD (Crohn's/UC)"
    description: "Inflammatory bowel disease - will expand to Crohn's disease, ulcerative colitis"
  - label: "Fibrosis"
    description: "Fibrotic conditions - will expand to pulmonary fibrosis, hepatic fibrosis, etc."
  - label: "Autoimmune"
    description: "Autoimmune conditions - will expand to RA, lupus, MS, etc."
```
This is a required field with no default - user must specify.

**Multi-Disease Interpretation:**
When user selects "Other" and enters custom text, Claude will interpret and expand the input:
- Parse comma-separated diseases (OR logic)
- Expand each disease to related terms, synonyms, subtypes
- Example: "IBD" → ["Crohn's disease", "ulcerative colitis", "inflammatory bowel disease", "colitis"]
- Example: "Crohns, UC" → ["Crohn's disease", "ulcerative colitis"]
- Fuzzy matching corrects typos when Claude CLI is unavailable

**Question 2: Tissue(s)**
```
header: "Tissue"
question: "Which tissue(s) or organ(s) should be included? (Multiple selections use OR logic)"
multiSelect: true
options:
  - label: "Intestine/Colon"
    description: "GI tract including ileum, colon, rectum"
  - label: "Lung"
    description: "Pulmonary tissue"
  - label: "Liver"
    description: "Hepatic tissue"
  - label: "Skin"
    description: "Dermal tissue"
```
This is a required field with no default - user must specify.

**Tissue Term Expansion (Claude-based):**
When user specifies a tissue, Claude automatically expands it to include anatomically related terms:
- Parse comma-separated tissues (OR logic)
- Expand each tissue to anatomical synonyms, subregions, and related terms
- Example: "intestine" → ["intestine", "colon", "ileum", "jejunum", "duodenum", "cecum", "rectum", "gut", "bowel", "gastrointestinal", "enteric", "colonic", "mucosa"]
- Example: "lung" → ["lung", "pulmonary", "bronchial", "alveolar", "airway", "respiratory"]
- Example: "colon" → ["colon", "colonic", "colorectal", "large intestine", "sigmoid", "cecum", "rectum"]
- Fallback to predefined anatomical mappings when Claude CLI is unavailable

This ensures papers using different anatomical terminology (e.g., "colon" vs "intestine") are included in search results.

**Question 2.5: Exact Match Mode (Optional)**
```
header: "Term Expansion"
question: "How should disease and tissue terms be handled?"
options:
  - label: "Auto-expand terms (Recommended)"
    description: "Claude/fuzzy expansion to related terms (e.g., 'intestine' → 'colon', 'ileum', etc.)"
  - label: "Exact match only"
    description: "Use terms exactly as provided - no expansion (use --exact flag)"
```
Default: Auto-expand

**Exact Match Mode:**
When `--exact` flag is used:
- Disease and tissue terms are used exactly as provided
- No Claude-based or fuzzy expansion is performed
- OR logic still applies for comma-separated or multiple terms
- Useful when you want precise control over search terms

Example:
```bash
# With expansion (default): "intestine" expands to colon, ileum, etc.
PYTHONPATH=~/.claude/skills/ncbi-paper-query/scripts python -m ncbi_paper_query \
    --disease "IBD" --tissue "intestine" ...

# Exact match: only "intestine" is searched
PYTHONPATH=~/.claude/skills/ncbi-paper-query/scripts python -m ncbi_paper_query \
    --disease "IBD" --tissue "intestine" --exact ...

# Exact with multiple terms (OR logic): searches "colon" OR "ileum"
PYTHONPATH=~/.claude/skills/ncbi-paper-query/scripts python -m ncbi_paper_query \
    --disease "IBD" --tissue "colon,ileum" --exact ...
```

**Question 3: Organism**
```
header: "Organism"
question: "Which organism(s) should be included?"
multiSelect: true
options:
  - label: "Human (Recommended)"
    description: "Homo sapiens studies"
  - label: "Mouse"
    description: "Mus musculus studies"
  - label: "Both"
    description: "Include both human and mouse"
```
Default: human

**Question 3.5: Keywords Filter (NEW)**
```
header: "Keywords"
question: "Which technology/method keywords should be used to filter results?"
multiSelect: true
options:
  - label: "Omics technologies (Recommended)"
    description: "single-cell, CITE-seq, spatial, omics, CyTOF, RNA-seq, etc."
  - label: "Custom keywords"
    description: "I'll specify my own keywords"
  - label: "No keyword filter"
    description: "Broader search without technology filtering (may return many more results)"
```
Default: Omics technologies

**Keywords Filter (OR logic):**
By default, queries include technology-focused keywords to prioritize omics-related papers:
- Default keywords: `single-cell`, `single cell`, `scRNA-seq`, `CITE-seq`, `spatial`, `omics`, `transcriptomics`, `proteomics`, `CyTOF`, `RNA-seq`, `gene expression`, `sequencing`
- These are combined with OR logic: papers matching ANY keyword are included
- Use `--no-keywords` to disable filtering (broader search, more results)
- Use `--keywords "term1" "term2"` for custom keywords

Example: "intestine" + "IBD" + default keywords narrows 6,000+ papers to ~500 focused omics papers

**Question 4: Impact Factor Cutoff**
```
header: "IF Cutoff"
question: "What minimum journal impact factor should be required?"
options:
  - label: "7.0 (Recommended)"
    description: "High-impact journals only"
  - label: "5.0"
    description: "Moderate-impact and above"
  - label: "3.0"
    description: "Broader journal coverage"
  - label: "No cutoff"
    description: "Include all journals (use --include-unknown-if)"
```
Default: 7.0

**Question 5: Year Cutoff**
```
header: "Year"
question: "What is the earliest publication year to include?"
options:
  - label: "2015 (Recommended)"
    description: "Last ~10 years of research"
  - label: "2018"
    description: "Last ~7 years"
  - label: "2020"
    description: "Last ~5 years"
  - label: "2010"
    description: "Extended historical coverage"
```
Default: 2015

**Question 6: Extraction Mode**
```
header: "Mode"
question: "How should omics accessions be extracted from papers?"
options:
  - label: "Abstract only (Recommended)"
    description: "Extract from abstract/title only - fastest, no authentication needed"
  - label: "Web scrape"
    description: "Scrape HTML using Playwright - automatically explores access buttons (Check Access, Institution, PDF) before falling back"
  - label: "PDF download"
    description: "Download and parse full PDFs - slowest, best coverage, may need credentials"
```
Default: Abstract only

**Question 7: Max Results**
```
header: "Max Results"
question: "Maximum number of papers to retrieve?"
options:
  - label: "1000 (Recommended)"
    description: "Comprehensive coverage"
  - label: "500"
    description: "Moderate coverage"
  - label: "200"
    description: "Quick overview"
  - label: "2000"
    description: "Extensive search"
```
Default: 1000

**Question 8: Study Name (Output Directory)**
```
header: "Study Name"
question: "What name should be used for the output directory?"
options:
  - label: "Custom name"
    description: "I'll specify a study name (e.g., 'ibd_study', 'crohns_2024')"
  - label: "Auto-generate"
    description: "Use timestamped filename in output/ directory"
```

When a study name is provided, the script creates a structured output directory:
```
output/{study_name}/
├── results.csv       # Main results with all paper metadata
├── accessions.csv    # Omics accessions in flat format (one row per accession)
└── downloads/        # Downloaded PDFs and scraped content
```

### Phase 2: Credential Check

**ALWAYS ask about credentials** when user selects **PDF download** or **Web scrape** mode:

1. First, search for .env file in working directory and parent directories:
   ```bash
   # Search up to 3 parent levels for .env
   ENV_FILE=$(find . .. ../.. ../../.. -maxdepth 1 -name ".env" 2>/dev/null | head -1)
   if [ -n "$ENV_FILE" ]; then
       echo "Found .env at: $ENV_FILE"
       cat "$ENV_FILE" | grep -E "^TAK_" | sed 's/=.*/=***/'
   else
       echo "No .env file found"
   fi
   ```

2. **If credentials found**, use `AskUserQuestion` to confirm:
   ```
   header: "Credentials"
   question: "How should institutional authentication be handled?"
   options:
     - label: "Use existing credentials (Recommended)"
       description: "Use TAK_ACCOUNT/TAK_KEY from .env file for subscription papers"
     - label: "Skip institutional auth"
       description: "Only access freely available papers (PMC open access)"
     - label: "Update credentials"
       description: "I want to provide different credentials"
   ```

3. **If NO .env found**, use `AskUserQuestion` to get credentials directly:
   ```
   header: "Credentials"
   question: "No .env file found. How would you like to proceed?"
   options:
     - label: "Enter credentials now"
       description: "I'll provide my institutional email and SSO key"
     - label: "Skip institutional auth"
       description: "Only access freely available papers (PMC open access)"
   ```

4. **If user selects "Enter credentials now"**, ask for the values:
   ```
   header: "Email"
   question: "Enter your institutional email (TAK_ACCOUNT):"
   options:
     - label: "I'll type my email"
       description: "e.g., your.name@institution.com"
   ```

   Then ask:
   ```
   header: "SSO Key"
   question: "Enter your SSO key (TAK_KEY):"
   options:
     - label: "I'll type my SSO key"
       description: "Your institutional single sign-on key"
   ```

5. **Ask user to confirm .env creation**:
   ```
   header: "Save Credentials"
   question: "Save credentials to .env file in current working directory ({pwd})?"
   options:
     - label: "Yes, save to .env (Recommended)"
       description: "Credentials will be saved for future queries"
     - label: "No, use for this session only"
       description: "Credentials will not be saved to disk"
   ```

6. **If user agrees**, create .env file:
   ```bash
   cat > .env << 'EOF'
   TAK_ACCOUNT={user_provided_email}
   TAK_KEY={user_provided_key}
   EOF
   echo "Created .env file in $(pwd)"
   ```

   **If user declines**, pass credentials via environment variables for this session only:
   ```bash
   export TAK_ACCOUNT="{user_provided_email}"
   export TAK_KEY="{user_provided_key}"
   ```

7. Map user choice to CLI flag:
   - "Use existing credentials" or credentials entered → (no flag, default behavior)
   - "Skip institutional auth" → `--no-institutional-auth`

**Note:** TAK_ACCOUNT is also used as the NCBI email for PubMed queries.

### Phase 3: Parameter Confirmation

Display all collected parameters in a table format and show the PubMed query that will be executed:

```markdown
## Query Parameters

| Parameter | Value |
|-----------|-------|
| Disease | [user input] |
| Tissue(s) | [user input] |
| Organism | [user input or "human"] |
| Impact Factor Cutoff | [user input or "7.0"] |
| Year Cutoff | [user input or "2015"] |
| Extraction Mode | [user input or "abstract-only"] |
| Max Results | [user input or "200"] |
| Study Name | [user input or auto-generated] |

## Command to Execute
```bash
PYTHONPATH=~/.claude/skills/ncbi-paper-query/scripts python -m ncbi_paper_query \
    --disease "[disease]" \
    --tissue [tissues] \
    --organism [organism] \
    --if-cutoff [if_cutoff] \
    --year-cutoff [year_cutoff] \
    [--exact] \
    [--abstract-only | --web-scrape | --subscription-download] \
    [--keywords "term1" "term2" | --no-keywords] \
    --max-results [max_results] \
    --output [study_name] \
    --format csv
```

**Note:** Run from your project directory. The PYTHONPATH ensures the package is found while outputs are saved to `./output/` in your current working directory.

**Note on Term Expansion:**
- By default, disease and tissue terms are auto-expanded (Claude/fuzzy)
- Use `--exact` for exact matching (no expansion, OR logic still applies)

**Note on Keywords:**
- Default keywords (omics-focused) are applied automatically
- Use `--no-keywords` for broader search without technology filtering
- Use `--keywords "custom1" "custom2"` to specify your own keywords

Proceed with this query?
```

Wait for user approval before executing.

### Phase 4: Execute Query

Run the query from the user's project directory (outputs will be saved there):

```bash
PYTHONPATH=~/.claude/skills/ncbi-paper-query/scripts python -m ncbi_paper_query \
    --disease "[disease]" \
    --tissue [tissues] \
    --organism [organism] \
    --if-cutoff [if_cutoff] \
    --year-cutoff [year_cutoff] \
    [extraction_mode_flag] \
    [keyword_flag] \
    --max-results [max_results] \
    --output [study_name] \
    --format csv
```

**Map extraction mode to flag:**
- "Abstract only" → `--abstract-only`
- "Web scrape" → `--web-scrape`
- "Subscription download" → `--subscription-download` (web-scrape free papers, download subscription only)
- "PDF download" → (no flag, default behavior - downloads all PDFs)

**Map keyword options to flag:**
- Default (omics keywords) → (no flag)
- Custom keywords → `--keywords "custom1" "custom2"`
- No filtering → `--no-keywords`

Add `--include-unknown-if` if user selected "No cutoff" for impact factor.

### Phase 5: Results & Export

After execution completes, present results:

1. **Summary Statistics**
   ```markdown
   ## Results Summary

   | Metric | Value |
   |--------|-------|
   | Total papers retrieved | [count] |
   | Papers with omics data | [count] |
   | GEO accessions found | [count] |
   | SRA accessions found | [count] |
   | ArrayExpress accessions found | [count] |
   ```

2. **Output Location**
   ```markdown
   ## Output Files

   Study directory: `./output/{study_name}/`
   ├── results.csv      # [N] papers with full metadata
   ├── accessions.csv   # [N] omics accessions (one row per accession)
   └── downloads/       # [N] PDFs downloaded
   ```

3. **Top Papers** (show first 5-10 with highest impact factor)
   ```markdown
   ## Top Papers by Impact Factor

   | PMID | Title | Journal | IF | Year | Omics |
   |------|-------|---------|----|----- |-------|
   ```

4. **Omics Datasets Found** (if any)
   ```markdown
   ## Omics Datasets Identified

   | Accession | Type | Paper PMID | Description |
   |-----------|------|------------|-------------|
   ```

### Phase 6: Extraction Interview (After Successful Access)

**TRIGGER:** When papers have been successfully accessed (Download_Status = "Downloaded" or "Full access without download" in the results).

**Step 1: Ask User for Extraction Requirements**

Use `AskUserQuestion`:
```
header: "Extract Info"
question: "What information would you like to extract from these papers?"
multiSelect: true
options:
  - label: "Study design & methods"
    description: "Sample size, inclusion/exclusion criteria, methodology"
  - label: "Key findings & results"
    description: "Main outcomes, statistical results, effect sizes"
  - label: "Biomarkers & targets"
    description: "Genes, proteins, pathways mentioned"
  - label: "Dataset details"
    description: "Omics data types, platforms, sample characteristics"
```

**Step 2: Generate Extraction Plan**

Based on user selections, generate an extraction plan:

```markdown
## Extraction Plan

### Papers to Process
- [N] papers with full access (Downloaded or web-scraped)

### Information to Extract
1. [Selected field 1] - extraction criteria
2. [Selected field 2] - extraction criteria
...

### Output Format
- Structured CSV/Excel with columns for each field
- One row per paper

### Processing Method
- Iterative extraction using Codex or Claude Code
- Each paper processed individually to ensure accuracy
```

**Step 3: User Alignment**

Show extraction plan and use `AskUserQuestion`:
```
header: "Confirm Plan"
question: "Does this extraction plan look correct?"
options:
  - label: "Yes, proceed"
    description: "Start extracting information from papers"
  - label: "Modify fields"
    description: "I want to change what's extracted"
  - label: "Add custom fields"
    description: "I have specific fields not listed"
  - label: "Skip extraction"
    description: "I don't need extraction now"
```

**Step 4: Execute Extraction**

Once user agrees, use the extraction template from the user's project directory:

```bash
python ~/.claude/skills/ncbi-paper-query/scripts/extraction_template.py \
  --study-name {study_name} \
  --fields "study_design,key_findings,biomarkers,genes_mentioned,sample_size,methodology"
```

**Available extraction fields:**
- `study_design` - Study type (RCT, cohort, scRNA-seq, etc.)
- `sample_size` - Number of patients/subjects
- `methodology` - Methods used
- `key_findings` - Main results and conclusions
- `statistical_results` - Statistical analyses
- `biomarkers` - Biomarkers identified
- `genes_mentioned` - Gene names/symbols
- `pathways` - Pathways analyzed
- `dataset_types` - Types of data generated
- `omics_platforms` - Technologies used
- `sample_characteristics` - Sample demographics
- `therapeutic_targets` - Potential drug targets
- `drug_candidates` - Drugs mentioned
- `cell_types_analyzed` - Cell populations studied
- `conclusions` - Author conclusions

For more sophisticated extraction, use Claude Code or Codex to process the PDFs in `./output/{study_name}/downloads/`.

**Step 5: Report Results**

After extraction completes:
```markdown
## Extraction Complete

### Summary
- Papers processed: [N]
- Successful extractions: [N]
- Failed/partial: [N]

### Output File
Results saved to: `./output/{study_name}/extraction_results.csv`

### Preview
| PMID | Title | [Field1] | [Field2] |
|------|-------|----------|----------|
| [First 5 rows] |

### Next Steps
- Review extraction results
- Codex session can be resumed with: `codex exec resume --last`
```

## Parameter Reference

For detailed parameter documentation, see [references/parameters.md](references/parameters.md).

## Supported Omics Databases

The script extracts accessions from these databases:
- **GEO** - GSE (series), GSM (sample), GPL (platform)
- **SRA** - SRP, SRR, SRX, SRS, PRJNA, PRJEB (BioProject)
- **ArrayExpress** - E-MTAB, E-GEOD
- **Single Cell Portal** - SCP (Broad Institute single-cell studies)
- **PRIDE** - PXD (proteomics)
- **MetaboLights** - MTBLS (metabolomics)
- **dbGaP** - phs (controlled access studies)
- **ENA** - ERP, ERR (European Nucleotide Archive)
- **NGDC/GSA** - CRA, HRA, PRJCA (China National Genomics Data Center)

## Impact Factor Lookup

Impact factors are retrieved in this order:
1. **`journal_if.json`** - User-maintainable reference file (~70 common journals with 2024 IF values)
   - Location: `~/.claude/skills/ncbi-paper-query/scripts/ncbi_paper_query/data/journal_if.json`
2. **`.if_cache.json`** - Cached web search results (created in working directory)
3. **Web search** - bioxbio.com lookup for journals not in either file

To update IF values, edit `journal_if.json` in the package data directory. No code changes required.

## Validation & Auto-Correction

The script runs multiple validation rounds (default: 3) that:
- Detect and **reject** year-like impact factor values (1900-2100)
- Validate IF values from web search before caching
- Check omics accession format validity
- Log all corrections, warnings, and issues

## Output Structure

When `--output {study_name}` is provided, the script creates the following structure **in the current working directory**:

```
./output/{study_name}/
├── results.csv        # Main results with all paper metadata
├── accessions.csv     # Omics accessions in flat format (one row per accession)
├── failed_studies.csv # Papers that couldn't be accessed (if any)
└── downloads/         # Downloaded PDFs (when using PDF download mode)
```

**Note:** Output is always created relative to where the command is executed, not in the skill directory.

### failed_studies.csv Columns

Created when papers cannot be downloaded (paywall, access issues):

- `PMID` - PubMed ID of inaccessible paper
- `Title` - Publication title
- `Journal` - Journal name
- `DOI` - Digital Object Identifier
- `Failure_Reason` - Why access failed (e.g., "Paywall", "No PDF available")

### results.csv Columns

- `PMID` - PubMed ID
- `Title` - Publication title
- `Authors` - Author list
- `Journal` - Journal name
- `Publication_Year` - Publication year
- `Impact_Factor` - Journal impact factor (auto-corrected if year-like value detected)
- `DOI` - Digital Object Identifier
- `Abstract` - Paper abstract (truncated to 500 chars)
- `Omics_Accessions` - Semicolon-separated accessions
- `Omics_Count` - Number of accessions found
- `Download_Status` - Free/Downloaded/Full access without download/Not Downloaded
- `Download_Path` - Path to downloaded PDF

### accessions.csv Columns

One row per accession for easy filtering and analysis:

- `PMID` - Source paper PubMed ID
- `Title` - Source paper title
- `Journal` - Source journal
- `Publication_Year` - Publication year
- `Accession` - Dataset accession (e.g., GSE134809)
- `Database` - Database type (GEO, SRA, ArrayExpress, etc.)
- `Technology` - Omics technology (RNA-seq, scRNA-seq, etc.)
- `Organism` - Species
- `Sample_Count` - Number of samples
- `Condition_Counts` - Condition breakdown
- `Data_URL` - Direct link to dataset

## Download Status Values

| Status | Description |
|--------|-------------|
| Free | Open access paper (PMC, preprint, etc.) |
| Downloaded | Full PDF downloaded to local disk |
| Full access without download | Content accessed via web-scraping (HTML) without PDF download |
| Not Downloaded | Abstract-only mode used, or download/web-scrape failed |

## Access Exploration (Playwright)

In Web scrape and PDF download modes, the system uses Playwright to:
- Navigate to journal/publisher pages
- Automatically click "Check Access", "Get PDF", or institution buttons
- Attempt institutional SSO login (Takeda) when configured
- Try multiple access strategies before falling back

## Supported Publishers (12 Total)

### Subscription Publishers (Require Institutional Auth)

| Publisher | Domains | SSO Method | Example Journals |
|-----------|---------|------------|------------------|
| **Elsevier** | sciencedirect.com, cell.com, thelancet.com, linkinghub.elsevier.com | OpenAthens | Gastroenterology, Cell, Immunity, Lancet |
| **Springer/Nature** | nature.com, springer.com, link.springer.com | WAYF | Nature, Cell & Mol Immunology, Mucosal Immunology |
| **Wiley** | wiley.com, onlinelibrary.wiley.com | Shibboleth | Immunology, Neurogastroenterology & Motility |
| **Oxford** | academic.oup.com, oup.com | Shibboleth | Inflammatory Bowel Diseases |
| **Taylor & Francis** | tandfonline.com | Shibboleth | Expert Review journals |
| **ACS** | pubs.acs.org | OpenAthens | Journal of Proteome Research |
| **AAAS/Science** | science.org, sciencemag.org | OpenAthens | Science, Science Signaling, Science Immunology |
| **Wolters Kluwer** | journals.lww.com, lww.com | Institution Login | Current Opinion in Gastroenterology, European J Gastro |
| **BMJ** | bmj.com, gut.bmj.com | BMJ SSO | Gut, Thorax |
| **Portland Press** | portlandpress.com | Institution Login | Clinical Science |

### Free Access Publishers (No Auth Required)

| Publisher | Domains | Notes |
|-----------|---------|-------|
| **bioRxiv/medRxiv** | biorxiv.org, medrxiv.org | Preprint servers - always free |
| **PMC Open Access** | ncbi.nlm.nih.gov/pmc | Open access subset |

### Download Strategy by Publisher

The system attempts downloads in this order:

1. **PMC Open Access** - Direct PDF download (no auth)
2. **bioRxiv/medRxiv** - Direct PDF via selectors (no auth)
3. **Publisher DOI** - Resolve DOI, detect publisher, attempt:
   - Direct PDF link (some open access)
   - Institutional access button → SSO login → PDF
   - Alternative PDF URL patterns

### Institutional Access Selectors

Each publisher has specific selectors for finding institutional access:

```python
# Elsevier
'a:has-text("Access through your institution")'
'a:has-text("Institutional access")'

# Springer/Nature
'a:has-text("Access through your institution")'
'a:has-text("Log in via an institution")'

# Wiley
'a:has-text("Institutional Login")'
'a:has-text("Shibboleth")'

# Wolters Kluwer
'a:has-text("Log In")'
'a:has-text("Institutional Access")'

# BMJ
'a:has-text("Log in")'
'a:has-text("Institutional access")'

# bioRxiv (free - PDF selectors)
'a:has-text("Download PDF")'
'a[href*=".full.pdf"]'
```

### Troubleshooting Publisher Access

| Issue | Possible Cause | Solution |
|-------|---------------|----------|
| "Not Downloaded" for supported publisher | SSO timeout or network issue | Retry, check VPN connection |
| Publisher not detected | New domain pattern | Check DOI resolution URL |
| PDF is HTML error page | Access denied | Verify institutional subscription |
| bioRxiv not downloading | Missing PDF selectors | Should be free - check URL |

## Example Usage

User: "Find papers about Crohn's disease with gene expression data"

Claude executes:
1. Interview for disease (Crohn's disease), tissues, organism, etc.
2. Check credentials if needed
3. Show parameter confirmation
4. Run query after approval
5. Display results with omics datasets found
