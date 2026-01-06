# ✅ Cortellis Target-Drug Web Automation - Setup & Quick Start

Complete guide for setting up and using Cortellis target-drug web automation to download target-drug relationship data files.

## 📚 Documentation Navigation

- **⚡ This File:** Setup guide with tested examples and quick start
- **📖 Quick Reference:** [`README.md`](README.md) - Command examples and troubleshooting
- **📚 Full Documentation:** [`../SKILL.md`](../SKILL.md) - Comprehensive guide including API access

---

The Cortellis skill has been successfully updated with comprehensive web automation capabilities!

## 🎯 What's New

### 1. Complete Skill Organization
```
/home/sagemaker-user/.claude/skills/cortellis/
├── SKILL.md                           # Main skill documentation
├── cortellis_api/                     # API access (Python)
│   └── scripts/
│       ├── cortellis_gene_query.py
│       └── convert_excel_to_json.py
└── cortellis_targetdrug_web/          # Target-Drug Web automation (Playwright)
    ├── README.md                      # Quick start guide
    ├── okta_auth_setup.py             # Authentication setup
    ├── cortellis-automation.js        # Single-category export
    ├── cortellis-download.js          # Multi-category downloads
    ├── run-cortellis.sh               # Wrapper for automation
    ├── run-cortellis-download.sh      # Wrapper for downloads
    ├── example_clinical_studies.json  # Example config
    ├── example_comprehensive.json     # Example config
    └── example_patents.json           # Example config
```

### 2. Environment-Agnostic Execution

✅ Scripts now work from any directory
✅ Auth file (`okta_auth_state.json`) read from working directory
✅ Results save to working directory in `cortellis_playwright_result/`

### 3. Enhanced Features

**cortellis-automation.js** (via `run-cortellis.sh`):
- Single category per search
- Export via "..." menu button
- Advanced UI navigation

**cortellis-download.js** (via `run-cortellis-download.sh`):
- Multiple queries: `--queries "drug1,drug2,drug3"`
- Multiple categories per query
- Default category: "Drugs & Biologics"
- Config file support

### 4. Folder Structure

Results automatically organize:
```
working_directory/
├── okta_auth_state.json              # Your Okta session
└── cortellis_playwright_result/
    └── {query}_{category}_{datetime}/
        ├── search_overview.png
        ├── {category}_page.png
        ├── {category}_data.xlsx
        └── metadata.json
```

## 🚀 Quick Start

### Step 1: Setup Authentication (One-Time)

```bash
cd /your/project/directory
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/okta_auth_setup.py
```

This creates `okta_auth_state.json` in your project directory.

### Step 2: Run Searches

**From any directory with the auth file:**

```bash
# Single category export
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis.sh imatinib --category "Clinical Studies"

# Multiple categories
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh imatinib --categories "Clinical Studies,Patents"

# Multiple drugs
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh --queries "imatinib,gefitinib" --categories "Clinical Studies"
```

## ✅ Test Results

Successfully tested on 2025-12-20:
- ✅ Authentication from working directory
- ✅ Results saved to working directory with proper naming
- ✅ Category navigation working
- ✅ Screenshot capture working
- ✅ Metadata generation working

Test output:
```
Results: /home/sagemaker-user/claude_code/playwright/cortellis_playwright_result/
Folder: imatinib_Clinical_Studies_2025-12-20_20-32-13/
Files: search_overview.png, Clinical_Studies_page.png, metadata.json
```

## 📖 Documentation

- **Main Documentation**: `/home/sagemaker-user/.claude/skills/cortellis/SKILL.md`
- **Quick Start**: `/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/README.md`
- **API Reference**: `references/api_reference.md`
- **Examples**: `examples/*.py`

## 🔧 Available Categories

All scripts support:
- **Drugs & Biologics** (default)
- **Clinical Studies** (most common)
- **Patents** (IP analysis)
- **Literature** (publications)
- Genes & Targets
- Experimental Pharmacology
- Drug Metabolism
- Drug-Drug Interactions
- Pharmacokinetics
- Experimental Models
- Organizations
- Disease Briefings

## 💡 Use Cases

### Clinical Intelligence
```bash
cd /your/project
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh \
  --queries "drug1,drug2,drug3" \
  --categories "Clinical Studies"
```

### Patent Landscape
```bash
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh \
  --config example_patents.json
```

### Comprehensive Drug Analysis
```bash
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh \
  imatinib \
  --categories "Drugs & Biologics,Clinical Studies,Patents,Literature"
```

## 🎓 When to Use What

### Use API Access (Python) When:
- Need structured JSON data
- Doing programmatic analysis
- Querying gene/target information
- Integrating with Python workflows

### Use Web Automation (Playwright) When:
- Need Excel/CSV file downloads
- Downloading from specific categories
- Bulk file exports
- Web-only features

## 📞 Support

### Re-authentication
If session expires (after weeks/months):
```bash
cd /your/project
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/okta_auth_setup.py
```

### Troubleshooting
- Check screenshots in output folders for debugging
- Review `metadata.json` for error details
- Ensure `okta_auth_state.json` exists in working directory

## 🎉 Ready to Use!

The Cortellis skill is now fully integrated and ready for both API and web automation workflows. All scripts are environment-agnostic and will work from any directory where you have the auth file.

**Next Steps:**
1. Set up authentication in your project directory
2. Try example queries
3. Check output in `cortellis_playwright_result/`
4. Explore the example config files

For detailed usage, see `SKILL.md` or `web_automation/README.md`!
