# 📚 Cortellis Skill - Documentation Map

Quick guide to navigate the Cortellis skill documentation.

## 🎯 Start Here Based on Your Need

### Need to Build an Interactive CI Dashboard from Excel? (CI Web Dashboard)

1. **👉 START:** [`cortellis_ci_web/CI_WEB.md`](cortellis_ci_web/CI_WEB.md)
   - Excel → interactive single-file HTML dashboard
   - Competitive intensity + opportunity ranking
   - Target/target-type/MoA/route filters + CSV export

2. **Quick Reference:** [`cortellis_ci_web/README.md`](cortellis_ci_web/README.md)

### Need Safety Assessment from OFF-X? (OFF-X Web Automation)

1. **👉 START:** [`offx_web/OFFX_AUTOMATION.md`](offx_web/OFFX_AUTOMATION.md)
   - Complete OFF-X automation guide
   - Target and drug safety assessment
   - Adverse event analysis
   - Master view exports

2. **Quick Reference:** [`offx_web/README.md`](offx_web/README.md)
   - Command syntax
   - Field options
   - Examples

### Need to Download Target-Drug Data from Cortellis? (Target-Drug Web Automation)

1. **👉 START:** [`cortellis_targetdrug_web/WEB_AUTOMATION.md`](cortellis_targetdrug_web/WEB_AUTOMATION.md)
   - Complete web automation guide for target-drug relationships
   - Authentication instructions
   - All features and options
   - Workflow examples for targets and drugs

2. **Quick Reference:** [`cortellis_targetdrug_web/README.md`](cortellis_targetdrug_web/README.md)
   - Command syntax
   - Examples
   - Troubleshooting

3. **Setup Guide:** [`cortellis_targetdrug_web/SETUP_COMPLETE.md`](cortellis_targetdrug_web/SETUP_COMPLETE.md)
   - Setup validation
   - Test results
   - Quick start commands

### Need to Query Data Programmatically? (API Access)

1. **👉 START:** [`cortellis_api/API_ACCESS.md`](cortellis_api/API_ACCESS.md)
   - Complete API guide
   - Query examples
   - Output formats
   - Data structure

2. **Examples:** [`cortellis_api/examples/`](cortellis_api/examples/) directory
   - `analyze_targets_example.py`
   - `generate_report_example.py`
   - `IBD_Target_Clinical_Intelligence_Report.md`

3. **Reference:** [`cortellis_api/references/`](cortellis_api/references/) directory
   - `api_fields.md`
   - `api_reference.md`
   - `json_schema.md`
   - `scoring_framework.md`

## 📁 Documentation Structure

```
/home/sagemaker-user/.claude/skills/cortellis/
│
├── 📄 SKILL.md                              # High-level overview & decision guide
├── 📄 DOCUMENTATION_MAP.md                  # This file - Navigation guide
│
├── 📂 cortellis_ci_web/                     # CI Web Dashboard Subskill
│   ├── 📄 CI_WEB.md                         # ⭐ Complete dashboard guide
│   ├── 📄 README.md                         # Quick reference
│   ├── 📂 scripts/                          # Excel → dashboard scripts
│   ├── 📂 references/                       # Taxonomy defaults
│   └── 📂 assets/                           # HTML template
│
├── 📂 cortellis_api/                        # API Access Subskill
│   ├── 📄 API_ACCESS.md                     # ⭐ Complete API guide
│   ├── 📂 scripts/                          # Python scripts
│   │   ├── cortellis_gene_query.py          # Main API query script
│   │   └── convert_excel_to_json.py         # Excel to JSON converter
│   ├── 📂 references/                       # API reference docs
│   │   ├── api_fields.md                    # Field descriptions
│   │   ├── api_reference.md                 # Endpoints & auth
│   │   ├── json_schema.md                   # Data structure
│   │   └── scoring_framework.md             # Target scoring
│   ├── 📂 examples/                         # Usage examples
│   │   ├── analyze_targets_example.py
│   │   ├── generate_report_example.py
│   │   └── IBD_Target_Clinical_Intelligence_Report.md
│   └── 📂 assets/                           # Additional resources
│
└── 📂 cortellis_targetdrug_web/             # Target-Drug Web Automation Subskill
    ├── 📄 WEB_AUTOMATION.md                 # ⭐ Complete web automation guide
    ├── 📄 SETUP_COMPLETE.md                 # Setup validation & test results
    ├── 📄 README.md                         # Quick reference
    ├── 🔧 cortellis-automation.js           # Single-category script
    ├── 🔧 cortellis-download.js             # Multi-category script
    ├── 🔧 run-cortellis.sh                  # Wrapper for automation.js
    ├── 🔧 run-cortellis-download.sh         # Wrapper for download.js
    └── 📂 example/                          # Config examples
        ├── example_clinical_studies.json
        ├── example_comprehensive.json
        └── example_patents.json

# NOTE: Okta authentication handled by centralized skill:
# ~/ai-sci-claude-skills/ai-sci/okta-sso/
```

## 🚦 Navigation Flow

### For Target-Drug Web Automation Users

```
SKILL.md (High-level overview)
    ↓
cortellis_targetdrug_web/WEB_AUTOMATION.md (⭐ COMPLETE GUIDE)
    ↓
cortellis_targetdrug_web/README.md (Quick Reference)
    ↓
cortellis_targetdrug_web/SETUP_COMPLETE.md (Setup validation)
```

### For API Access Users

```
SKILL.md (High-level overview)
    ↓
cortellis_api/API_ACCESS.md (⭐ COMPLETE GUIDE)
    ↓
cortellis_api/examples/ (Working examples)
    ↓
cortellis_api/references/ (Detailed reference)
```

## 🔍 Finding What You Need

| I Want To... | Go To... |
|-------------|----------|
| Download files from Cortellis | [`cortellis_targetdrug_web/WEB_AUTOMATION.md`](cortellis_targetdrug_web/WEB_AUTOMATION.md) |
| Quick web automation command reference | [`cortellis_targetdrug_web/README.md`](cortellis_targetdrug_web/README.md) |
| Query gene/drug data via API | [`cortellis_api/API_ACCESS.md`](cortellis_api/API_ACCESS.md) |
| Understand web automation options | [`cortellis_targetdrug_web/WEB_AUTOMATION.md`](cortellis_targetdrug_web/WEB_AUTOMATION.md) |
| See API example workflows | [`cortellis_api/examples/`](cortellis_api/examples/) directory |
| API field reference | [`cortellis_api/references/api_fields.md`](cortellis_api/references/api_fields.md) |
| Troubleshoot web automation | [`cortellis_targetdrug_web/README.md`](cortellis_targetdrug_web/README.md) → Troubleshooting |
| Choose between API and Web | [`SKILL.md`](SKILL.md) → "Choose Your Access Method" |

## ⚡ Quick Commands

### Web Automation Setup

```bash
# 1. Setup auth (one-time) - uses centralized okta-sso skill
~/ai-sci-claude-skills/ai-sci/okta-sso/run-okta-login.sh

# 2. Run search
/home/sagemaker-user/.claude/skills/cortellis/cortellis_targetdrug_web/run-cortellis-download.sh \
  imatinib --categories "Clinical Studies"
```

### API Access Setup

```bash
# 1. Setup credentials (one-time)
cd /your/project/directory
echo "CORTELLIS_API_KEY=your_key" > .env
echo "CORTELLIS_API_SECRET=your_secret" >> .env

# 2. Run query
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_api/scripts/cortellis_gene_query.py \
  BRCA1 --fields drug biomarker --excel
```

## 💡 Tips

- **New to Cortellis Skill?** Start with [`SKILL.md`](SKILL.md) to choose your approach
- **Ready to download files?** Jump to [`cortellis_targetdrug_web/WEB_AUTOMATION.md`](cortellis_targetdrug_web/WEB_AUTOMATION.md)
- **Need quick commands?** Check [`cortellis_targetdrug_web/README.md`](cortellis_targetdrug_web/README.md)
- **API user?** [`cortellis_api/API_ACCESS.md`](cortellis_api/API_ACCESS.md) has everything you need
- **Stuck?** Check the Troubleshooting sections in respective documentation files

## 📞 Help & Support

Each documentation file includes:
- ✅ Prerequisites and setup instructions
- ✅ Usage examples
- ✅ Troubleshooting sections
- ✅ Cross-references to related docs

Follow the documentation structure above to find exactly what you need!
