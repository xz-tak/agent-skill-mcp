# Cortellis Target-Drug Web Automation - Quick Reference

Browser-based automation for downloading target-drug relationship data from Cortellis Drug Discovery Intelligence. Supports queries for targets with associated drugs and drugs with associated targets.

## 📚 Documentation

- **👉 First Time Setup:** See [`SETUP_COMPLETE.md`](SETUP_COMPLETE.md) for complete setup guide with test results
- **📖 Full Documentation:** See [`../SKILL.md`](../SKILL.md) for comprehensive guide including API access options
- **⚡ This File:** Quick command reference and examples

## Quick Start

### 1. Setup Authentication (One-time)

```bash
~/ai-sci-claude-skills/ai-sci/okta-sso/run-okta-login.sh
```

Follow prompts to save your Okta session. Creates `~/.okta/auth_state.json` (shared by all Cortellis/OFF-X automation).

### 2. Run Searches

**Single query, single category:**
```bash
./run-cortellis.sh imatinib --category "Clinical Studies"
```

**Multiple queries, multiple categories:**
```bash
./run-cortellis-download.sh --queries "imatinib,gefitinib" --categories "Clinical Studies,Patents"
```

**Using config file:**
```bash
./run-cortellis-download.sh --config example_clinical_studies.json
```

## Which Script to Use?

### Use `run-cortellis.sh` when:
- Searching single or multiple terms
- Exporting ONE category per term
- Want advanced export options (via "..." menu)

### Use `run-cortellis-download.sh` when:
- Need MULTIPLE categories per term
- Batch downloading for analysis
- Using config files for complex queries

## Available Scripts

### `run-cortellis.sh` → `cortellis-automation.js`
Single category export with full UI navigation.

**Examples:**
```bash
# Default category (Drugs & Biologics)
./run-cortellis.sh imatinib

# Specific category
./run-cortellis.sh imatinib --category "Clinical Studies"

# Multiple terms, same category
./run-cortellis.sh imatinib dasatinib --category "Patents"

# From file
echo -e "imatinib\ndasatinib\nnilotinib" > drugs.txt
./run-cortellis.sh --file drugs.txt --category "Clinical Studies"
```

### `run-cortellis-download.sh` → `cortellis-download.js`
Multiple categories per term, bulk downloads.

**Examples:**
```bash
# Single term, default category
./run-cortellis-download.sh imatinib

# Single term, multiple categories
./run-cortellis-download.sh imatinib --categories "Clinical Studies,Patents,Literature"

# Multiple terms, default category
./run-cortellis-download.sh --queries "imatinib,gefitinib,dasatinib"

# Multiple terms, multiple categories
./run-cortellis-download.sh --queries "imatinib,gefitinib" --categories "Clinical Studies,Patents"

# Config file (most flexible)
./run-cortellis-download.sh --config example_comprehensive.json
```

## Available Categories

All scripts support these categories (case-sensitive):
- **Drugs & Biologics** (default)
- **Clinical Studies** (most common for drug research)
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

## Example Config Files

Three example configs are provided:

### `example_clinical_studies.json`
Clinical trial downloads for multiple drugs:
```json
[
  {"searchTerm": "imatinib", "categories": ["Clinical Studies"]},
  {"searchTerm": "dasatinib", "categories": ["Clinical Studies"]},
  {"searchTerm": "nilotinib", "categories": ["Clinical Studies"]}
]
```

### `example_comprehensive.json`
Multiple categories per drug:
```json
[
  {
    "searchTerm": "imatinib",
    "categories": ["Drugs & Biologics", "Clinical Studies", "Patents", "Literature"]
  }
]
```

### `example_patents.json`
Patent landscape for kinase inhibitors:
```json
[
  {"searchTerm": "imatinib", "categories": ["Patents"]},
  {"searchTerm": "erlotinib", "categories": ["Patents"]},
  {"searchTerm": "gefitinib", "categories": ["Patents"]},
  {"searchTerm": "lapatinib", "categories": ["Patents"]}
]
```

## Output Location

All results save to working directory:
```
working_directory/cortellis_playwright_result/
└── {searchTerm}_{category}_{datetime}/
    ├── {category}_page.png
    ├── {category}_{filename}.xlsx
    └── metadata.json
```

## Common Use Cases

### Drug Development Research
```bash
# Get all development data for a drug
./run-cortellis-download.sh imatinib --categories "Drugs & Biologics,Clinical Studies,Patents"
```

### Clinical Intelligence
```bash
# Clinical trials for multiple drugs
./run-cortellis-download.sh --queries "drug1,drug2,drug3" --categories "Clinical Studies"
```

### Patent Analysis
```bash
# Patent landscape
./run-cortellis-download.sh --queries "drug1,drug2,drug3,drug4" --categories "Patents"
```

### Literature Review
```bash
# Publications
./run-cortellis.sh imatinib --category "Literature"
```

## Troubleshooting

### Authentication Errors

**Error: "Okta session not found"**
```bash
# Run the okta-sso skill to authenticate
~/ai-sci-claude-skills/ai-sci/okta-sso/run-okta-login.sh
```

**Error: "Authentication expired"**
```bash
# Session expired (happens after weeks/months)
~/ai-sci-claude-skills/ai-sci/okta-sso/run-okta-login.sh  # Re-authenticate
```

### Download Issues

**Error: "Category not found"**
- Check exact spelling (case-sensitive)
- Category may have no results for this search
- View `search_results.png` to see available categories

**Error: "Download button not found"**
- Check screenshots in output folder
- UI may have changed
- Try different category

### Script Issues

**Scripts not executable:**
```bash
chmod +x *.sh
```

**Playwright not found:**
```bash
# Edit shell scripts to update PLAYWRIGHT_SKILL_DIR path
# Default: ${HOME}/.claude/plugins/cache/playwright-skill/playwright-skill/4.1.0/skills/playwright-skill
```

## Tips

1. **Start Small**: Test with one drug and one category first
2. **Check Screenshots**: Always review output screenshots to verify results
3. **Use Config Files**: Easier to manage complex multi-drug/category queries
4. **Default Category**: If unsure, "Drugs & Biologics" is a safe default
5. **Rate Limiting**: Scripts include 3-second delays between searches

## Need More Help?

- **Full Documentation**: See `../SKILL.md`
- **API Access**: For programmatic queries, see `../scripts/`
- **Examples**: Check the `example_*.json` files
- **Debug**: Check screenshots and metadata.json in output folders

## File Structure

```
cortellis_targetdrug_web/
├── README.md                          # This file
├── cortellis-automation.js            # Single-category script
├── cortellis-download.js              # Multi-category script
├── run-cortellis.sh                   # Wrapper for automation
├── run-cortellis-download.sh          # Wrapper for download
├── example_clinical_studies.json      # Example config
├── example_comprehensive.json         # Example config
└── example_patents.json               # Example config

# Okta auth handled by centralized skill:
~/ai-sci-claude-skills/ai-sci/okta-sso/
```

## Quick Command Reference

```bash
# Setup (one-time) - uses centralized okta-sso skill
~/ai-sci-claude-skills/ai-sci/okta-sso/run-okta-login.sh

# Single category
./run-cortellis.sh DRUG --category "Category"

# Multiple categories
./run-cortellis-download.sh DRUG --categories "Cat1,Cat2"

# Multiple drugs, multiple categories
./run-cortellis-download.sh --queries "Drug1,Drug2" --categories "Cat1,Cat2"

# Config file
./run-cortellis-download.sh --config config.json
```
