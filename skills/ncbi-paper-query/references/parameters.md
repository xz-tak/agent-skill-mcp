# NCBI Paper Query - Parameter Reference

## Required Parameters

| Parameter | CLI Flag | Description |
|-----------|----------|-------------|
| Disease | `--disease`, `-d` | Target disease or indication (e.g., "Crohn's disease", "ulcerative colitis", "fibrosis") |
| Tissue(s) | `--tissue`, `-t` | Target tissue(s) or organ(s). Multiple values use OR logic (e.g., "intestine colon ileum") |

## Optional Parameters

| Parameter | CLI Flag | Default | Description |
|-----------|----------|---------|-------------|
| Organism | `--organism`, `-o` | human | Organisms to include. Can specify multiple (e.g., "human mouse") |
| Cell Type | `--cell-type`, `-c` | None | Specific cell type to filter for (e.g., "fibroblast", "epithelial") |
| Impact Factor Cutoff | `--if-cutoff` | 7.0 | Minimum journal impact factor. Papers below this are excluded |
| Year Cutoff | `--year-cutoff` | 2015 | Earliest publication year. Papers before this year are excluded |
| Max Results | `--max-results` | 200 | Maximum number of papers to retrieve |
| Validation Rounds | `--rounds`, `-r` | 3 | Number of validation/proofreading rounds for omics extraction |
| Output Path | `--output`, `-O` | Auto-generated | Custom output file path |
| Format | `--format` | csv | Output format: `csv` or `excel` |

## Extraction Mode Flags

| Mode | CLI Flag | Description |
|------|----------|-------------|
| Abstract Only | `--abstract-only` | Extract omics accessions from abstract/title only (no full-text access required) |
| Web Scrape | `--web-scrape` | Extract omics from web HTML using Playwright. Auto-fallback to PDF if no accessions found |
| (Default) | (none) | Full PDF download mode - downloads and parses PDF full text |

## Filter Flags

| Flag | CLI Option | Description |
|------|------------|-------------|
| No Preprints | `--no-preprints` | Exclude preprint publications |
| Include Unknown IF | `--include-unknown-if` | Include papers with unknown impact factors (excluded by default) |

## Authentication Flags

| Flag | CLI Option | Description |
|------|------------|-------------|
| Subscription Only | `--subscription-only` | Only download subscription papers (skip free PMC papers). Requires VPN |
| No Institutional Auth | `--no-institutional-auth` | Disable institutional authentication |

## Credentials (.env file)

For institutional/subscription paper access, credentials are loaded from `.env` file in the script directory (`/home/sagemaker-user/claude_code/ncbi_paper/.env`), NOT from shell environment variables.

| Variable | Description |
|----------|-------------|
| `TAK_ACCOUNT` | Institutional email (also used as NCBI email for queries) |
| `TAK_KEY` | Institutional SSO authentication key |

Example `.env` file:
```
TAK_ACCOUNT=user@institution.com
TAK_KEY=your_sso_key
```

## Example Commands

### Basic Query (Abstract Only)
```bash
python ncbi_paper_retriever.py \
    --disease "Crohn's disease" \
    --tissue intestine colon ileum \
    --organism human \
    --if-cutoff 7.0 \
    --year-cutoff 2015 \
    --max-results 200 \
    --abstract-only
```

### Web Scraping Mode
```bash
python ncbi_paper_retriever.py \
    --disease "ulcerative colitis" \
    --tissue colon rectum \
    --organism human mouse \
    --if-cutoff 5.0 \
    --year-cutoff 2018 \
    --web-scrape
```

### Full PDF Mode (with institutional access)
```bash
python ncbi_paper_retriever.py \
    --disease "pulmonary fibrosis" \
    --tissue lung \
    --organism human \
    --if-cutoff 10.0 \
    --year-cutoff 2020 \
    --max-results 100
```

## Output

Results are saved to CSV format in the `./output/` directory with timestamp:
- `results_YYYYMMDD_HHMMSS.csv`

### Output Columns

| Column | Description |
|--------|-------------|
| pmid | PubMed ID |
| title | Publication title |
| authors | Author list |
| journal | Journal name |
| year | Publication year |
| impact_factor | Journal impact factor |
| doi | Digital Object Identifier |
| abstract | Paper abstract |
| omics_accessions | Extracted accessions from all supported databases |
| accession_type | Type of omics data |
| omics_metadata | Enriched metadata from omics databases |
| Download_Status | How content was accessed: Free, Downloaded, Full access without download, Not Downloaded |

### Download Status Values

| Status | Description |
|--------|-------------|
| Free | Open access paper (PMC, preprint, etc.) |
| Downloaded | Full PDF downloaded to local disk |
| Full access without download | Content accessed via web-scraping (HTML) without PDF download |
| Not Downloaded | Abstract-only mode used, or download/web-scrape failed |

## Supported Omics Databases

| Database | Accession Pattern | Description |
|----------|-------------------|-------------|
| GEO | GSE, GSM, GPL | Gene Expression Omnibus |
| SRA | SRP, SRR, SRX, SRS, PRJNA, PRJEB | Sequence Read Archive / BioProject |
| ArrayExpress | E-MTAB, E-GEOD | EBI ArrayExpress |
| Single Cell Portal | SCP | Broad Institute single-cell studies |
| PRIDE | PXD | Proteomics database |
| MetaboLights | MTBLS | Metabolomics database |
| dbGaP | phs | Database of Genotypes and Phenotypes |
| ENA | ERP, ERR | European Nucleotide Archive |
| NGDC/GSA | CRA, HRA, PRJCA, SAMC, CRR, CRX | China National Genomics Data Center |

## Supported Publishers (12 Total)

### Subscription Publishers

| Publisher | Domains | DOI Prefix | SSO URL |
|-----------|---------|------------|---------|
| Elsevier | sciencedirect.com, cell.com, thelancet.com, linkinghub.elsevier.com | 10.1016, 10.1053 | auth.elsevier.com |
| Springer/Nature | nature.com, springer.com, link.springer.com | 10.1038, 10.1007 | wayf.springernature.com |
| Wiley | wiley.com, onlinelibrary.wiley.com | 10.1111, 10.1002 | onlinelibrary.wiley.com/action/ssostart |
| Oxford | academic.oup.com, oup.com | 10.1093 | academic.oup.com/my-account/oauth/shibboleth |
| Taylor & Francis | tandfonline.com | 10.1080 | tandfonline.com/action/ssostart |
| ACS | pubs.acs.org | 10.1021 | pubs.acs.org/action/ssostart |
| AAAS/Science | science.org, sciencemag.org | 10.1126 | science.org/action/ssostart |
| Wolters Kluwer | journals.lww.com, lww.com | 10.1097 | journals.lww.com/pages/login.aspx |
| BMJ | bmj.com, gut.bmj.com | 10.1136 | sso.bmj.com |
| Portland Press | portlandpress.com | 10.1042 | (institution login) |

### Free Access Publishers

| Publisher | Domains | DOI Prefix | Notes |
|-----------|---------|------------|-------|
| bioRxiv/medRxiv | biorxiv.org, medrxiv.org | 10.1101 | Preprints - always free |
| PMC Open Access | ncbi.nlm.nih.gov/pmc | N/A | Open access subset |

### Publisher Detection

Publishers are detected by:
1. DOI prefix (e.g., 10.1038 → Springer/Nature)
2. Final URL domain after DOI resolution

### Not Supported

Papers from these publishers cannot be downloaded:
- Edizioni Minerva (minervamedica.it) - Italian publisher
- Papers without DOI

## Impact Factor Lookup

Impact factors are retrieved via web search (not hardcoded) in this order:

| Source | File | Description |
|--------|------|-------------|
| Reference File | `journal_if.json` | User-maintainable JSON with ~70 common journals (2024 IF values) |
| Cache | `.if_cache.json` | Cached results from web searches |
| Web Search | bioxbio.com | Live lookup for journals not in reference or cache |

To update IF values, edit `journal_if.json` in the script directory (`/home/sagemaker-user/claude_code/ncbi_paper/`).

## Validation & Auto-Correction

The script performs multi-round validation (default: 3 rounds):

| Feature | Description |
|---------|-------------|
| Year-like IF rejection | Rejects impact factors that look like years (1900-2100) |
| IF validation | Validates IF values from web search before caching |
| Cache cleanup | Removes invalid cached IF values automatically |
| Accession validation | Verifies omics accession format validity |
| Logging | All corrections, warnings, and issues are logged |
