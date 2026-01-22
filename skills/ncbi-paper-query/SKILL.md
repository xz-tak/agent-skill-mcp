---
name: ncbi-paper-query
description: Query NCBI PubMed for publications matching disease/tissue/organism criteria and extract omics accessions (GEO, SRA, ArrayExpress) with metadata enrichment. This skill should be used when users request literature searches for omics studies, gene expression datasets, or publication lists for specific diseases and tissues.
---

# NCBI Paper Query

## Overview

Query PubMed for disease/tissue-specific publications and extract omics dataset accessions (GEO, SRA, ArrayExpress) with full metadata enrichment. The workflow interviews the user to collect search parameters, confirms settings, then executes the query.

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
  - label: "200 (Recommended)"
    description: "Comprehensive coverage"
  - label: "100"
    description: "Moderate coverage"
  - label: "50"
    description: "Quick overview"
  - label: "500"
    description: "Extensive search"
```
Default: 200

### Phase 2: Credential Check

**ALWAYS ask about credentials** when user selects **PDF download** or **Web scrape** mode:

1. First, check current credential status:
   ```bash
   cat /home/sagemaker-user/claude_code/ncbi_paper/.env 2>/dev/null | grep -E "^TAK_" | sed 's/=.*/=***/' || echo "No .env file found"
   ```

2. Display credential status to user and ALWAYS use `AskUserQuestion` to confirm:
   ```
   header: "Credentials"
   question: "How should institutional authentication be handled?"
   options:
     - label: "Use existing credentials (Recommended)"
       description: "Use TAK_ACCOUNT/TAK_KEY from .env file for subscription papers"
     - label: "Skip institutional auth"
       description: "Only access freely available papers (PMC open access)"
     - label: "Update credentials"
       description: "I'll provide instructions to update .env file"
   ```

   If no credentials found, adjust options:
   ```
   header: "Credentials"
   question: "No institutional credentials configured. How would you like to proceed?"
   options:
     - label: "Continue without (Recommended)"
       description: "Use only freely available papers (PMC open access)"
     - label: "Set credentials now"
       description: "I'll provide instructions to create .env file"
   ```

3. If user wants to set/update credentials, provide instructions:
   ```
   Create or update .env file in /home/sagemaker-user/claude_code/ncbi_paper/:

   TAK_ACCOUNT=your_email@institution.com
   TAK_KEY=your_sso_key

   Then re-run the query.
   ```

4. Map user choice to CLI flag:
   - "Use existing credentials" → (no flag, default behavior)
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

## Command to Execute
```bash
python /home/sagemaker-user/claude_code/ncbi_paper/ncbi_paper_retriever.py \
    --disease "[disease]" \
    --tissue [tissues] \
    --organism [organism] \
    --if-cutoff [if_cutoff] \
    --year-cutoff [year_cutoff] \
    [--abstract-only | --web-scrape] \
    --max-results [max_results] \
    --format csv
```

Proceed with this query?
```

Wait for user approval before executing.

### Phase 4: Execute Query

Run the query using Bash:

```bash
python /home/sagemaker-user/claude_code/ncbi_paper/ncbi_paper_retriever.py \
    --disease "[disease]" \
    --tissue [tissues] \
    --organism [organism] \
    --if-cutoff [if_cutoff] \
    --year-cutoff [year_cutoff] \
    [extraction_mode_flag] \
    --max-results [max_results] \
    --format csv
```

Map extraction mode to flag:
- "Abstract only" → `--abstract-only`
- "Web scrape" → `--web-scrape`
- "PDF download" → (no flag, default behavior)

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
   ## Output File

   Results saved to: `./output/results_YYYYMMDD_HHMMSS.csv`
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
- Iterative extraction using Codex (gpt-5.2 Xhigh full-auto)
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

**Step 4: Execute Extraction with Codex**

Once user agrees, invoke Codex using the skill-codex pattern:

```bash
codex exec -m gpt-5.2 \
  --config model_reasoning_effort="xhigh" \
  --sandbox workspace-write \
  --full-auto \
  --skip-git-repo-check \
  -C /home/sagemaker-user/claude_code/ncbi_paper \
  "Extract the following information from papers in ./output/downloads/:

   EXTRACTION REQUIREMENTS:
   [list extraction fields from user selection]

   INSTRUCTIONS:
   1. Process each PDF/HTML file iteratively
   2. For each paper, extract the requested fields
   3. Output results to ./output/extraction_results.csv
   4. Include columns: PMID, Title, [field1], [field2], ...
   5. Use 'N/A' for missing information
   6. Be thorough and accurate

   OUTPUT FORMAT:
   CSV with headers matching extraction fields" 2>/dev/null
```

**Step 5: Report Results**

After extraction completes:
```markdown
## Extraction Complete

### Summary
- Papers processed: [N]
- Successful extractions: [N]
- Failed/partial: [N]

### Output File
Results saved to: `./output/extraction_results.csv`

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
2. **`.if_cache.json`** - Cached web search results
3. **Web search** - bioxbio.com lookup for journals not in either file

To update IF values, edit `journal_if.json` in the script directory. No code changes required.

## Validation & Auto-Correction

The script runs multiple validation rounds (default: 3) that:
- Detect and **reject** year-like impact factor values (1900-2100)
- Validate IF values from web search before caching
- Check omics accession format validity
- Log all corrections, warnings, and issues

## Output Format

Results are saved as CSV with columns:
- `pmid` - PubMed ID
- `title` - Publication title
- `authors` - Author list
- `journal` - Journal name
- `year` - Publication year
- `impact_factor` - Journal impact factor (auto-corrected if year-like value detected)
- `doi` - Digital Object Identifier
- `abstract` - Paper abstract
- `omics_accessions` - Extracted accessions from all supported databases
- `accession_type` - Type of omics data
- `omics_metadata` - Enriched metadata from omics databases

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
