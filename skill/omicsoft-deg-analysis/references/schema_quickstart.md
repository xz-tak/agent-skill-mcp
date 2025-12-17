# Quick Start Guide - Schema Viewer

## ✅ SOLUTION: "Error loading schema_report.json"

The error you encountered is due to browser CORS restrictions when opening HTML files locally. **I've created 3 solutions for you:**

---

## 🌟 BEST SOLUTION (Recommended)

### Open: `schema_viewer_standalone.html`

**Just double-click the file to open in your browser!**

✓ **Works immediately** - No setup required
✓ **No errors** - All data embedded in single file
✓ **Fully functional** - Search, filter, expand/collapse
✓ **Portable** - Share single file via email
✓ **File size:** 2.4 MB (includes all data)

---

## Alternative Solutions

### Option 2: File Picker Version

**File:** `schema_viewer_local.html`

1. Double-click to open in browser
2. Click "📁 Select schema_report.json" button
3. Browse and select the JSON file
4. Data loads and viewer works

**Use when:** You want smaller HTML file, don't mind clicking to load

---

### Option 3: Web Server (Advanced)

**File:** `schema_viewer.html`

```bash
# Start web server
python -m http.server 8000

# Open in browser
http://localhost:8000/schema_viewer.html
```

**Use when:** You're comfortable with command line and web servers

---

## What's Inside the Viewer?

### 📊 Complete Metadata Documentation

**40 Columns Fully Documented:**

#### 🔍 Key Filtering (10 columns)
- comparison_category - 8 values
- tissue - 131 values
- disease - 178 values
- case_treatment - 377 values
- case_treatment_status - 70 values
- control_treatment - 228 values
- control_treatment_status - 50 values
- comparison - 6,635 values
- study - 845 values
- database - 1 value

#### 👥 Demographics (6 columns)
- case_gender, case_age_category, case_ethnicity
- control_gender, control_age_category, control_ethnicity

#### 🏥 Disease Information (8 columns)
- Case and control disease states, subtypes, groups, locations

#### 💊 Treatment Information (6 columns)
- Dosages, treatment groups, treatment times

#### 📈 Response Information (2 columns)
- Case and control response data

#### 🔬 Sample Information (6 columns)
- Tissue types, sample materials, project IDs

#### 📋 Comparison Details (2 columns)
- Comparison IDs and contrasts

---

## Features You Get

### 🔍 Real-time Search
Type anything to search across ALL columns and values:
- Search "scleroderma" → finds all related diseases
- Search "skin" → finds all skin tissues
- Search "infliximab" → finds treatment statuses

### 📂 Category Filtering
Filter by category:
- Demographics only
- Key filtering columns only
- Treatment info only
- etc.

### 📊 Visual Progress Bars
See distribution of values at a glance

### ⚡ Expand/Collapse
- Click category headers to expand/collapse
- Click column headers to see all values
- "Expand All" / "Collapse All" buttons

### 🎯 Value Details
For each unique value, see:
- **Exact string** (for copying to filter commands)
- **Count** (number of observations)
- **Percentage** (% of total)
- **Visual bar** (relative frequency)

---

## How to Use with deg_analysis.py

### Recommended Workflow

**Step 1**: Explore schema to find filter values
**Step 2**: Validate filters before full analysis
**Step 3**: Run full analysis if validation succeeds

### Example 1: Find all scleroderma diseases

1. Open schema_viewer_standalone.html
2. Search "sclero"
3. Click on "disease" column
4. Copy exact disease names you want

**First, validate filters:**
```bash
python validate_filters.py \
  --file data.h5ad \
  --target-name Scleroderma_Analysis \
  --diseases "systemic scleroderma,diffuse scleroderma,limited scleroderma" \
  --signatures "MySignature:Gene1,Gene2,Gene3" \
  ...
```

**If validation succeeds, run full analysis:**
```bash
python deg_analysis.py \
  --diseases "systemic scleroderma,diffuse scleroderma,limited scleroderma" \
  ...
```

### Example 2: Find treatment response comparisons

1. Search "response"
2. Look at "comparison" column
3. Copy comparison strings

**Validate first:**
```bash
python validate_filters.py \
  --file data.h5ad \
  --target-name Response_Analysis \
  --comparison "response vs no response" \
  --signatures "ResponseGenes:Gene1,Gene2" \
  ...
```

**Then run analysis:**
```bash
python deg_analysis.py \
  --comparison "response vs no response" \
  ...
```

### Example 3: Filter by exact comparison category

1. Click "comparison_category" column
2. See all 8 options
3. Copy exact strings (must match exactly!)

**Validate filters step-by-step:**
```bash
python validate_filters.py \
  --file data.h5ad \
  --target-name Category_Analysis \
  --comparison-category "Disease vs. Normal,Treatment vs. Control" \
  --signatures "MyGenes:Gene1,Gene2" \
  ...
```

**Run full analysis after validation:**
```bash
python deg_analysis.py \
  --comparison-category "Disease vs. Normal,Treatment vs. Control" \
  ...
```

### Why Validate First?

The validation script:
- ✓ Tests each filter incrementally
- ✓ Shows matching values and observation counts
- ✓ Detects zero-result queries immediately
- ✓ Provides suggestions for fixing problematic filters
- ✓ Prevents wasting time on full analyses that return 0 observations

---

## Understanding the Data

### Filtering Types

**Exact Match** (must match string exactly):
- comparison_category
- study
- case_treatment
- control_treatment
- case_treatment_status
- control_treatment_status

**Partial Match** (case-insensitive substring):
- tissue
- disease
- comparison

### Example Values

**comparison_category (8 total):**
```
Treatment vs. Control        (1,620 obs, 21.3%)
Treatment1 vs. Treatment2    (1,460 obs, 19.2%)
Other Comparisons            (1,173 obs, 15.4%)
Disease vs. Normal           (1,063 obs, 14.0%)
Disease1 vs. Disease2        (  948 obs, 12.5%)
CellType1 vs. CellType2      (  839 obs, 11.0%)
Tissue1 vs. Tissue2          (  349 obs,  4.6%)
Responder vs. Non-Responder  (  152 obs,  2.0%)
```

**case_ethnicity (8 total):**
```
NA                       (7,516 obs, 98.8%)
European American        (   19 obs,  0.2%)
African American         (   15 obs,  0.2%)
Mixed Ethnicity          (    3 obs,  0.0%)
White                    (    3 obs,  0.0%)
Asian Chinese            (    1 obs,  0.0%)
Hispanic                 (    1 obs,  0.0%)
non-Hispanic Caucasian   (    1 obs,  0.0%)
```

---

## Files Reference

```
schema_viewer_standalone.html  ⭐ RECOMMENDED - Works immediately
schema_viewer_local.html       File picker version
schema_viewer.html             Web server version
schema_report.json             Raw JSON data (3.6 MB)
schema_report.txt              Text summary
README_schema.md               Full documentation
QUICKSTART.md                  This file
```

---

## Troubleshooting

### Browser won't open HTML file?
- Right-click → "Open with" → Choose browser
- Or drag file onto browser window

### Viewer loads but shows no data?
- You're using schema_viewer.html (needs web server)
- Solution: Use schema_viewer_standalone.html instead

### Search not working?
- Make sure you've expanded the category first
- Or use "Expand All" button to see everything

### Can't find a specific disease/tissue?
- Use the search box at the top
- It searches ALL values across ALL columns

---

## Questions?

- **Schema structure:** See `README_schema.md`
- **deg_analysis usage:** See skill documentation in `/home/sagemaker-user/.claude/skills/omicsoft-deg-analysis/SKILL.md`
- **AnnData format:** See `references/anndata_schema.md`

---

**Ready to explore? Just open `schema_viewer_standalone.html`!** 🚀
