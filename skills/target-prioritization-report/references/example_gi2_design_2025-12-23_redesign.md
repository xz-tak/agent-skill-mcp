# Target Prioritization System - Redesign Specification

**Date:** December 23, 2025
**Purpose:** Complete redesign of data loading and scoring system for IBD target prioritization
**Status:** Design Phase - Ready for Implementation

---

## Executive Summary

This document specifies the complete redesign of the target prioritization system with:
- **Dual-score architecture**: Raw scores preserved for traceability, normalized (0-100) for display/calculation
- **Report-first parsing**: Prefer structured sidecars (CSV/JSON) and fall back to Markdown/HTML reports
- **Enhanced competitive intelligence**: Blue Ocean categorization with Jaccard similarity
- **Risk-weighted safety scoring**: Amplified penalties for severe adverse events
- **Synergy-aware combination scoring**: Pathway and PPI synergy adjustments

---

## Core Principles

1. **Raw score preservation** - All original values stored in result dictionary
2. **Lazy normalization** - 0-100 mapping only when needed for calculations/display
3. **Graceful degradation** - Missing data sources default to 50.0 (neutral)
4. **Report-first parsing** - Prefer structured sidecars (CSV/JSON) when available; fall back to Markdown/HTML
5. **Single data loading pass** - Load all sources once via orchestrator
6. **Gene naming in reports** - Use "GENE1-GENE2" format, not "List 1, 2, 3..."

---

## Module Structure

```
target_prioritization/
├── data_loaders/
│   ├── __init__.py
│   ├── parsers.py          # Markdown/HTML/JSON parsers
│   ├── loaders.py          # Individual data source loaders
│   └── orchestrator.py     # Master load_all_data()
├── scoring/
│   ├── __init__.py
│   ├── normalizers.py      # Raw → 0-100 normalization
│   ├── subscores.py        # Individual subscore calculators
│   └── scoring.py          # Overall scoring logic
```

---

## Complete System Workflow

### End-to-End Execution Flow

```
1. USER INPUT
   ├─ Gene/Combination name (e.g., "JAK1" or "TYK2-JAK1")
   └─ Results directory path (default: "results/")

2. DATA LOADING PHASE (orchestrator.py)
   ├─ Load all data sources in parallel
   ├─ Parse markdown/HTML/JSON files
   ├─ Extract raw scores and metadata
   ├─ Store in unified data dictionary
   └─ Log any missing/failed sources

3. SCORING PHASE (scoring.py)
   ├─ Calculate 5 subscores (using normalizers + subscores modules)
   │  ├─ Clinical Validation
   │  ├─ Disease Association
   │  ├─ Safety
   │  ├─ Opportunity
   │  └─ Novelty
   ├─ Calculate synergy metrics (for combinations)
   ├─ Calculate overall weighted score
   └─ Perform Blue Ocean categorization

4. REPORT GENERATION PHASE
   ├─ Compile scores into structured format
   ├─ Embed visualizations (heatmaps, networks, plots)
   ├─ Generate traceability sections with raw scores
   ├─ Create summary tables and rankings
   └─ Output HTML report

5. OUTPUT
   ├─ HTML report file
   ├─ JSON score bundles (required for reproducibility)
   └─ Log file with processing details
```

---

## Detailed Data Transformation Workflows

### Workflow 1: Clinical Validation Score

```
SOURCE FILE: results/cortellis/IBD_Target_Analysis_Report.md

STEP 1 - PARSE MARKDOWN
  ├─ Read file line by line
  ├─ Find target gene section (e.g., "### 1. JAK1")
  ├─ Extract using regex: r"Total Score:\s*(\d+\.?\d*)"
  └─ Store: raw_score = 134.5

STEP 2 - STORE RAW SCORE
  └─ data_dict["clinical_validation_raw"] = 134.5

STEP 3 - NORMALIZE (when needed for calculation)
  ├─ max_clinical_score = 134.5  # JAK1 is max
  ├─ normalized = (134.5 / 134.5) × 100
  └─ data_dict["clinical_validation"] = 100.0

STEP 4 - DISPLAY IN REPORT
  ├─ Summary Table: Show "100.0"
  └─ Traceability Section: Show "134.5 (Raw), 100.0 (Normalized)"

FOR COMBINATIONS:
  ├─ Load individual gene scores: JAK1=100.0, TYK2=62.2
  ├─ Average: (100.0 + 62.2) / 2 = 81.1
  └─ Store: combo_data["clinical_validation"] = 81.1
```

### Workflow 2: Safety Score

```
SOURCE FILES (preferred → fallback):
  ├─ results/offx/gene_safety_scores_detailed.json
  └─ results/offx/OFF-X_Safety_Analysis_Report.md

STEP 1 - LOAD SEVERITY COUNTS
  ├─ If JSON exists:
  │   ├─ Read gene_results[GENE]
  │   └─ Extract counts: very_high/high/medium/low/very_low/not_associated/na (+total_aes)
  └─ Else parse markdown:
      ├─ Find gene section: "### 1. JAK1"
      ├─ Extract severity breakdown table
      └─ Parse counts:
          ├─ Very high: 39
          ├─ High: 118
          ├─ Medium: 404
          ├─ Low: 669
          ├─ Very low: 11,601
          ├─ Not associated: 43
          └─ NA: 1,279

STEP 2 - STORE RAW BREAKDOWN
  └─ data_dict["safety_breakdown_raw"] = {
        "very_high": 39,
        "high": 118,
        "medium": 404,
        "low": 669,
        "very_low": 11601,
        "not_associated": 43,
        "na": 1279
      }

STEP 3 - CALCULATE RISK-WEIGHTED SCORE
  ├─ weighted_sum = (39×0×10) + (118×10×5) + (404×20×2) +
  │                 (669×60×1) + (11601×80×1) + (43×100×1)
  │               = 0 + 5900 + 16160 + 40140 + 928080 + 4300
  │               = 994,580
  │
  ├─ weighted_total = (39×10) + (118×5) + (404×2) +
  │                   (669×1) + (11601×1) + (43×1)
  │                 = 390 + 590 + 808 + 669 + 11601 + 43
  │                 = 14,101
  │
  └─ safety_score = 994,580 / 14,101 = 70.5

STEP 4 - STORE CALCULATED SCORE
  └─ data_dict["safety"] = 70.5  # Already 0-100, no normalization

STEP 5 - DISPLAY IN REPORT
  ├─ Summary Table: Show "70.5"
  └─ Traceability Section: Show breakdown table + formula

FOR COMBINATIONS:
  ├─ Check if OFF-X combo data exists
  ├─ If yes: Use combo breakdown directly
  └─ If no: Average individual gene scores
```

### Workflow 3: Disease Association Score

```
SOURCE FILES:
  ├─ results/deg_results/target_scores.csv (preferred)
  ├─ results/deg_results/IBD_Targets_Summary_Report.md (fallback)
  ├─ results/biobridge/individual_report/{gene}-IBD.md
  ├─ results/ultra/{gene}-predictions.json (optional)
  ├─ results/pathwaydb/centrality_gene_{GENE}_centrality.csv (preferred for centrality)
  ├─ results/pathwaydb/gene_{GENE}_pathways.csv (context/fallback)
  └─ results/primekg/{gene}-IBD.md

STEP 1 - LOAD DEG SCORE (30% weight)
  ├─ If CSV exists: parse row where Gene == {GENE} from target_scores.csv
  ├─ Else: parse markdown table row from IBD_Targets_Summary_Report.md
  ├─ Extract: JAK1 score = 9.0 (raw)
  ├─ Store raw: data_dict["deg_score_raw"] = 9.0
  ├─ Normalize: max_deg = 15.0 (PCOLCE)
  │             normalized = (9.0 / 15.0) × 100 = 60.0
  └─ Store: deg_normalized = 60.0

STEP 2 - LOAD BIOBRIDGE SCORE (25% weight)
  ├─ Parse markdown: "Percentile Rank: 90.3%"
  ├─ Store raw: data_dict["biobridge_percentile_raw"] = 90.3
  ├─ Already 0-100, no normalization
  └─ Store: biobridge_normalized = 90.3

STEP 3 - LOAD ULTRA SCORE (20% weight)
  ├─ Try to load JSON file
  ├─ If exists:
  │   ├─ Extract: percentile_rank = 0.85 (raw 0-1)
  │   ├─ Store raw: data_dict["ultra_percentile_raw"] = 0.85
  │   ├─ Normalize: 0.85 × 100 = 85.0
  │   └─ Store: ultra_normalized = 85.0
  └─ If missing:
      └─ Use default: ultra_normalized = 50.0

STEP 4 - LOAD PATHWAY CENTRALITY (15% weight)
  ├─ Parse pathway centrality CSV for centrality metrics
  ├─ If available:
  │   ├─ Extract centrality score
  │   ├─ Normalize to 0-100
  │   └─ Store: pathway_centrality_normalized = X
  └─ If not available:
      └─ Use default: pathway_centrality_normalized = 50.0

STEP 5 - LOAD PRIMEKG SCORE (10% weight)
  ├─ Parse markdown: "X connection(s)"
  ├─ Store raw: data_dict["primekg_connections_raw"] = X
  ├─ Find max connections across all genes: max_conn = 50
  ├─ Normalize: (X / 50) × 100
  └─ Store: primekg_normalized = Y

STEP 6 - CALCULATE WEIGHTED DISEASE ASSOCIATION
  └─ disease_assoc = (60.0 × 0.30) + (90.3 × 0.25) +
                     (85.0 × 0.20) + (50.0 × 0.15) + (Y × 0.10)
                   = 18.0 + 22.575 + 17.0 + 7.5 + (Y × 0.10)
                   = 65.075 + (Y × 0.10)

STEP 7 - STORE AND DISPLAY
  ├─ Store: data_dict["disease_association"] = 73.2
  ├─ Summary Table: Show "73.2"
  └─ Traceability Section: Show component breakdown
      ├─ DEG: 9.0 (raw) → 60.0 (30%)
      ├─ BioBridge: 90.3% (25%)
      ├─ ULTRA: 0.85 → 85.0 (20%)
      ├─ Pathway: 50.0 (15%, default)
      └─ PrimeKG: X → Y (10%)

FOR COMBINATIONS:
  ├─ DEG: Average individual genes
  ├─ BioBridge: USE COMBO PERCENTILE (from hardcoded table)
  │             e.g., TYK2-JAK1 combo = 99.8%
  ├─ ULTRA: Average individual genes
  ├─ Pathway: Average individual genes
  └─ PrimeKG: Average individual genes
```

### Workflow 4: Opportunity Score

```
DEPENDENCIES:
  ├─ Disease Association score (calculated above)
  ├─ Clinical Validation score (calculated above)
  └─ Competitive Intelligence data (from CI HTML)

STEP 1 - LOAD COMPETITIVE INTELLIGENCE
  ├─ Parse results/ci/ibd_dashboard.html
  ├─ Extract JSON from <script id="data"> tag
  ├─ Filter entries where target = "JAK1" AND ibdTags exist
  └─ Count programs by phase:
      ├─ Marketed: 5 drugs
      ├─ Phase III: 3 drugs
      ├─ Phase II: 2 drugs
      ├─ Phase I: 4 drugs
      └─ Preclinical: 8 drugs

STEP 2 - CALCULATE WEIGHTED COMPETITION
  ├─ weighted_comp = (5 × 1.0) + (3 × 0.7) + (2 × 0.4) + (4 × 0.2)
  │                = 5.0 + 2.1 + 0.8 + 0.8 = 8.7
  ├─ total_programs = 5 + 3 + 2 + 4 + 8 = 22
  ├─ pct_marketed = (5 / 22) × 100 = 22.7%
  └─ pct_market_left = 100 - 22.7 = 77.3

STEP 3 - CALCULATE CI SCORE
  └─ CI_score = pct_market_left = 77.3

STEP 4 - CALCULATE CLINICAL NOVELTY
  └─ clinical_novelty = 100 - clinical_validation_normalized
                      = 100 - 100.0 = 0.0

STEP 5 - CALCULATE OPPORTUNITY SCORE
  └─ opportunity = (disease_assoc × 0.5) + (clinical_novelty × 0.3) +
                   (CI_score × 0.2)
                 = (73.2 × 0.5) + (0.0 × 0.3) + (77.3 × 0.2)
                 = 36.6 + 0.0 + 15.46
                 = 52.06

STEP 6 - STORE AND DISPLAY
  ├─ Store: data_dict["opportunity"] = 52.06
  └─ Traceability Section:
      ├─ Disease Association: 73.2 (50%)
      ├─ Clinical Novelty: 0.0 (30%)
      └─ CI Score: 77.3% market left (20%)

FOR COMBINATIONS (e.g., TYK2-JAK1):

  STEP 1 - Calculate Individual Opportunities
    ├─ TYK2 opportunity = 62.5
    └─ JAK1 opportunity = 52.06

  STEP 2 - Calculate Combo CI with Jaccard
    ├─ drugs_targeting_TYK2 = {drug1, drug2, drug3, ...}  (10 drugs)
    ├─ drugs_targeting_JAK1 = {drug1, drug4, drug5, ...}  (22 drugs)
    ├─ intersection = {drug1}  (1 drug)
    ├─ union = {drug1, drug2, drug3, drug4, drug5, ...}  (31 drugs)
    ├─ jaccard = 1 / 31 = 0.032
    │
    ├─ average_comp = (TYK2_comp + JAK1_comp) / 2
    ├─ combo_comp_raw = average_comp × (1 - 0.032) + union_comp × 0.032
    │
    ├─ Check pathway_synergy (from synergy calculation)
    ├─ If pathway_synergy > 0.5:
    │   └─ adjustment = 1 - (pathway_synergy × 0.3)
    │       combo_comp_adjusted = combo_comp_raw × adjustment
    │
    ├─ Calculate % market left for union
    └─ CI_score_combo = pct_market_left_union

  STEP 3 - Calculate Combo Opportunity
    ├─ mean_individual_opps = (62.5 + 52.06) / 2 = 57.28
    ├─ novel_mechanism_bonus = 100 - (combined_synergy × 100)
    │                        = 100 - (0.65 × 100) = 35.0
    └─ opportunity_combo = (57.28 × 0.6) + (35.0 × 0.4)
                         = 34.37 + 14.0 = 48.37
```

### Workflow 5: Novelty Score

```
DEPENDENCIES:
  ├─ Clinical Validation score
  └─ PrimeKG connection count

STEP 1 - CALCULATE CLINICAL NOVELTY (70% weight)
  └─ clinical_novelty = 100 - clinical_validation_normalized
                      = 100 - 100.0 = 0.0

STEP 2 - CALCULATE LITERATURE NOVELTY (30% weight)
  ├─ connection_count = 45 (from PrimeKG)
  ├─ literature_percentile = min((45 / 50) × 100, 100) = 90.0
  └─ literature_novelty = 100 - 90.0 = 10.0

STEP 3 - CALCULATE NOVELTY SCORE
  └─ novelty = (0.0 × 0.7) + (10.0 × 0.3) = 3.0

STEP 4 - STORE AND DISPLAY
  ├─ Store: data_dict["novelty"] = 3.0
  └─ Traceability Section:
      ├─ Clinical Novelty: 0.0 (70%)
      └─ Literature Novelty: 10.0 (30%)

FOR COMBINATIONS:
  ├─ Calculate individual novelties: TYK2=18.5, JAK1=3.0
  ├─ mean_novelty = (18.5 + 3.0) / 2 = 10.75
  └─ combo_novelty = 10.75 × 1.10 = 11.83  # 10% premium
```

### Workflow 6: Overall Score Calculation

```
INPUT: All 5 subscores calculated above
  ├─ Clinical Validation: 100.0
  ├─ Disease Association: 73.2
  ├─ Safety: 70.5
  ├─ Opportunity: 52.06
  └─ Novelty: 3.0

STEP 1 - APPLY WEIGHTS
  └─ overall = (100.0 × 0.30) + (73.2 × 0.30) + (70.5 × 0.10) +
               (52.06 × 0.20) + (3.0 × 0.10)
             = 30.0 + 21.96 + 7.05 + 10.41 + 0.30
             = 69.72

STEP 2 - STORE OVERALL SCORE
  └─ data_dict["overall_score"] = 69.72

STEP 3 - BLUE OCEAN CATEGORIZATION
  ├─ Check Opportunity: 52.06 > 50 → High Opportunity
  ├─ Check Competition: pct_market_left = 77.3 > 25 → Low Competition
  └─ Category: "Blue Ocean"

STEP 4 - STORE CATEGORIZATION
  └─ data_dict["blue_ocean_category"] = "Blue Ocean"
      data_dict["blue_ocean_rationale"] = "High opportunity (52.1) with low competition (77.3% market left)"
```

### Workflow 7: Synergy Metrics (Combinations Only)

```
SOURCE FILES:
  ├─ results/pathwaydb/REPORT_list1_TYK2_JAK1.md
  └─ results/interactdb/shortest_paths_analysis_complete/COMPLETE_ANALYSIS_REPORT_ALL_DATABASES.md

STEP 1 - PATHWAY SYNERGY
  ├─ Parse pathway report
  ├─ Extract: "110 shared pathways"
  ├─ Calculate: pathway_synergy = 110 / 110 = 1.0  # Max normalized
  └─ Store: combo_data["pathway_synergy"] = 1.0

STEP 2 - PPI SYNERGY
  ├─ Parse PPI report for TYK2-JAK1 pair
  ├─ Extract:
  │   ├─ min_hops = 2
  │   └─ databases_found = 3 (STRING, BioGRID, IntAct)
  ├─ Calculate:
  │   ├─ hop_score = 1.0 / 2 = 0.5
  │   ├─ database_consensus = 3 / 3 = 1.0
  │   └─ ppi_synergy = 0.5 × 1.0 = 0.5
  └─ Store: combo_data["ppi_synergy"] = 0.5

STEP 3 - COMBINED SYNERGY
  └─ combined_synergy = (1.0 × 0.5) + (0.5 × 0.5) = 0.75

STEP 4 - STORE
  └─ combo_data["combined_synergy"] = 0.75
```

---

## Data Flow Architecture (Detailed)

```
┌─────────────────────────────────────────────────────────────────┐
│                      PHASE 1: DATA LOADING                       │
└─────────────────────────────────────────────────────────────────┘
                                │
                    ┌───────────┴───────────┐
                    │   orchestrator.py     │
                    │  load_all_data()      │
                    └───────────┬───────────┘
                                │
            ┌───────────────────┼───────────────────┐
            │                   │                   │
    ┌───────▼───────┐   ┌──────▼──────┐   ┌───────▼───────┐
    │  parsers.py   │   │ loaders.py  │   │File I/O Layer │
    │               │   │             │   │               │
    │parse_markdown │   │load_cortel  │   │Read MD/HTML   │
    │parse_html_json│   │load_offx    │   │Parse JSON     │
    │parse_table    │   │load_deg     │   │Extract regex  │
    └───────┬───────┘   │load_biobridge   │               │
            │           │load_primekg │   └───────────────┘
            │           │load_pathway │
            │           │load_ppi     │
            │           │load_coexpr  │
            │           │load_ci      │
            │           └──────┬──────┘
            │                  │
            └──────────┬───────┘
                       │
                ┌──────▼──────┐
                │ Raw Data    │
                │ Dictionary  │
                │             │
                │ {gene: ..., │
                │  clinical_  │
                │  validation_│
                │  raw: 134.5,│
                │  ...}       │
                └──────┬──────┘
                       │
┌──────────────────────┴──────────────────────────────────────────┐
│                    PHASE 2: SCORING                              │
└──────────────────────────────────────────────────────────────────┘
                       │
         ┌─────────────┴─────────────┐
         │      scoring.py           │
         │  score_target()           │
         └─────────────┬─────────────┘
                       │
        ┌──────────────┼──────────────┐
        │              │              │
┌───────▼───────┐ ┌───▼────────┐ ┌──▼─────────┐
│normalizers.py │ │subscores.py│ │scoring.py  │
│               │ │            │ │            │
│normalize_     │ │calc_clinical validate    │
│clinical()     │ │calc_disease│calc_overall │
│normalize_deg()│ │calc_safety │categorize   │
│normalize_bb() │ │calc_opp    │blue_ocean   │
│...            │ │calc_novelty│            │
└───────┬───────┘ └───┬────────┘ └──┬─────────┘
        │             │             │
        └─────────────┼─────────────┘
                      │
              ┌───────▼────────┐
              │ Scored Data    │
              │ Dictionary     │
              │                │
              │ {gene: ...,    │
              │  clinical: 100,│
              │  disease: 73.2,│
              │  safety: 70.5, │
              │  opportunity:  │
              │   52.06,       │
              │  novelty: 3.0, │
              │  overall: 69.7,│
              │  category:     │
              │   "Blue Ocean"}│
              └───────┬────────┘
                      │
┌─────────────────────┴────────────────────────────────────────────┐
│                  PHASE 3: REPORT GENERATION                       │
└───────────────────────────────────────────────────────────────────┘
                      │
        ┌─────────────┴─────────────┐
        │  report_generator.py      │
        │  generate_html_report()   │
        └─────────────┬─────────────┘
                      │
     ┌────────────────┼────────────────┐
     │                │                │
┌────▼─────┐   ┌─────▼──────┐  ┌──────▼──────┐
│Summary   │   │Traceability│  │Visualization│
│Tables    │   │Sections    │  │Embedding    │
│          │   │            │  │             │
│Overall   │   │Raw scores  │  │Heatmaps     │
│rankings  │   │Data sources│  │Networks     │
│Subscores │   │Calculation │  │Boxplots     │
│Blue Ocean│   │steps       │  │             │
└────┬─────┘   └─────┬──────┘  └──────┬──────┘
     │               │                │
     └───────────────┼────────────────┘
                     │
              ┌──────▼──────┐
              │  HTML       │
              │  Report     │
              │  File       │
              └─────────────┘
```

---

## Report Compilation Strategy

### Report Structure

```
IBD Target Prioritization Report
├── 1. Executive Summary
│   ├─ Overall Rankings Table (Top targets by overall score)
│   ├─ Blue Ocean Categorization Chart
│   └─ Key Findings Bullets
│
├── 2. Individual Gene Analysis
│   ├─ For each gene (JAK1, TYK2, GREM1, ...):
│   │   ├─ 2.1 Gene Overview Card
│   │   │   ├─ Overall Score (normalized, 0-100)
│   │   │   ├─ Blue Ocean Category Badge
│   │   │   └─ Quick stats (# drugs, trials, safety)
│   │   │
│   │   ├─ 2.2 Subscore Breakdown (Spider/Radar Chart)
│   │   │   ├─ Clinical Validation: 100.0
│   │   │   ├─ Disease Association: 73.2
│   │   │   ├─ Safety: 70.5
│   │   │   ├─ Opportunity: 52.1
│   │   │   └─ Novelty: 3.0
│   │   │
│   │   ├─ 2.3 Detailed Score Analysis
│   │   │   ├─ Clinical Validation Section
│   │   │   │   ├─ Raw: 134.5
│   │   │   │   ├─ Normalized: 100.0
│   │   │   │   ├─ Source: Cortellis Report
│   │   │   │   └─ Interpretation: "Highly validated..."
│   │   │   │
│   │   │   ├─ Disease Association Section
│   │   │   │   ├─ Component Breakdown Table:
│   │   │   │   │   ├─ DEG: 9.0 (raw) → 60.0 (norm) × 30% = 18.0
│   │   │   │   │   ├─ BioBridge: 90.3% × 25% = 22.6
│   │   │   │   │   ├─ ULTRA: 85.0 × 20% = 17.0
│   │   │   │   │   ├─ Pathway: 50.0 × 15% = 7.5
│   │   │   │   │   └─ PrimeKG: 70.0 × 10% = 7.0
│   │   │   │   ├─ Final: 72.1
│   │   │   │   └─ Data Sources Links
│   │   │   │
│   │   │   ├─ Safety Section
│   │   │   │   ├─ Risk-Weighted Score: 70.5
│   │   │   │   ├─ Severity Breakdown Table
│   │   │   │   ├─ Calculation Formula Display
│   │   │   │   └─ Source: OFF-X Report
│   │   │   │
│   │   │   ├─ Opportunity Section
│   │   │   │   ├─ Component Breakdown:
│   │   │   │   │   ├─ Disease Assoc: 72.1 × 50% = 36.1
│   │   │   │   │   ├─ Clinical Novelty: 0.0 × 30% = 0.0
│   │   │   │   │   └─ CI Score: 77.3 × 20% = 15.5
│   │   │   │   ├─ Final: 51.6
│   │   │   │   ├─ Competitive Landscape Table
│   │   │   │   └─ Blue Ocean Position Chart
│   │   │   │
│   │   │   └─ Novelty Section
│   │   │       ├─ Component Breakdown:
│   │   │       │   ├─ Clinical: 0.0 × 70% = 0.0
│   │   │       │   └─ Literature: 10.0 × 30% = 3.0
│   │   │       ├─ Final: 3.0
│   │   │       └─ Interpretation: "Well-studied target"
│   │   │
│   │   └─ 2.4 Data Traceability
│   │       ├─ All Source Files Table
│   │       ├─ Parse Timestamps
│   │       └─ Data Completeness (% sources loaded)
│   │
│   └─ Repeat for all genes...
│
├── 3. Combination Analysis
│   ├─ For each combination (TYK2-JAK1, TNFRSF25-GREM1, ...):
│   │   ├─ 3.1 Combination Overview Card
│   │   │   ├─ Combo Name (e.g., "TYK2-JAK1")
│   │   │   ├─ Overall Score
│   │   │   ├─ Blue Ocean Category
│   │   │   └─ Synergy Metrics Badge
│   │   │
│   │   ├─ 3.2 Subscore Breakdown
│   │   │   ├─ Same 5 subscores as individual
│   │   │   └─ Radar chart comparing to individuals
│   │   │
│   │   ├─ 3.3 Synergy Analysis
│   │   │   ├─ Pathway Synergy: 1.0
│   │   │   │   ├─ Shared pathways: 110
│   │   │   │   ├─ Pathway network PNG embedded
│   │   │   │   ├─ Pathway UpSet plot PNG embedded
│   │   │   │   └─ Link to pathway report MD
│   │   │   │
│   │   │   ├─ PPI Synergy: 0.5
│   │   │   │   ├─ Min hops: 2
│   │   │   │   ├─ Databases: 3/3
│   │   │   │   └─ Path visualization
│   │   │   │
│   │   │   └─ Combined Synergy: 0.75
│   │   │
│   │   ├─ 3.4 Coexpression Analysis
│   │   │   ├─ Bulk RNA-seq Correlation
│   │   │   │   ├─ Mean correlation: 0.82
│   │   │   │   ├─ P-value: < 0.001
│   │   │   │   ├─ Correlation heatmap (HTML embed)
│   │   │   │   └─ Expression boxplot (HTML embed)
│   │   │   │
│   │   │   └─ Single-Cell Correlation
│   │   │       ├─ By cell type table
│   │   │       └─ 5 heatmaps embedded (B, T, Myeloid, etc.)
│   │   │
│   │   ├─ 3.5 Competitive Intelligence (Combo-specific)
│   │   │   ├─ Jaccard Similarity: 0.032
│   │   │   ├─ Union drug set: 31 drugs
│   │   │   ├─ Adjusted competition score
│   │   │   └─ Synergy-adjusted CI score
│   │   │
│   │   └─ 3.6 Data Traceability
│   │       └─ All source files for combination
│   │
│   └─ Repeat for all combinations...
│
├── 4. Comparative Analysis
│   ├─ 4.1 Overall Rankings
│   │   ├─ Combined table (individuals + combos)
│   │   └─ Sortable by any subscore
│   │
│   ├─ 4.2 Blue Ocean Landscape
│   │   ├─ 2x2 scatter plot
│   │   ├─ X-axis: Opportunity (0-100)
│   │   ├─ Y-axis: Competition (% market left)
│   │   └─ Quadrant labels + targets plotted
│   │
│   ├─ 4.3 Subscore Heatmap
│   │   ├─ Rows: All targets
│   │   ├─ Columns: 5 subscores
│   │   └─ Color gradient: 0-100
│   │
│   └─ 4.4 Recommendation Summary
│       ├─ Top 3 Blue Ocean targets
│       ├─ Highest opportunity targets
│       └─ Most novel targets
│
├── 5. Methodology
│   ├─ 5.1 Data Sources Description
│   ├─ 5.2 Scoring Formulas
│   ├─ 5.3 Normalization Methods
│   ├─ 5.4 Weight Justifications
│   └─ 5.5 Limitations & Caveats
│
└── 6. Appendices
    ├─ A. Raw Data Tables (all raw scores)
    ├─ B. Calculation Examples (step-by-step)
    ├─ C. Data Source Files Index
    └─ D. Glossary of Terms
```

### Report Generation Steps

```
STEP 1: INITIALIZE REPORT STRUCTURE
  ├─ Create HTML template with styling
  ├─ Set up section placeholders
  └─ Initialize data containers

STEP 2: GENERATE EXECUTIVE SUMMARY
  ├─ Aggregate all scored targets
  ├─ Sort by overall_score (descending)
  ├─ Create rankings table (top 10)
  ├─ Generate Blue Ocean categorization chart
  └─ Write key findings (auto-generated insights)

STEP 3: GENERATE INDIVIDUAL GENE SECTIONS
  For each gene in [JAK1, TYK2, ITGA4, ...]:
    ├─ Create gene overview card with overall score
    ├─ Generate radar chart (5 subscores)
    ├─ Build detailed score analysis sections
    │   ├─ For each subscore:
    │   │   ├─ Show normalized score (display)
    │   │   ├─ Show raw score (traceability)
    │   │   ├─ Show component breakdown (if applicable)
    │   │   ├─ Display calculation formula
    │   │   └─ Link to source files
    │   └─ Add interpretation text
    ├─ Generate competitive intel visualization
    ├─ Create data traceability table
    └─ Append to report HTML

STEP 4: GENERATE COMBINATION SECTIONS
  For each combo in [TYK2-JAK1, TNFRSF25-GREM1, ...]:
    ├─ Create combo overview card
    ├─ Generate comparison radar chart (combo vs individuals)
    ├─ Build synergy analysis section
    │   ├─ Embed pathway network PNG
    │   ├─ Embed pathway UpSet plot PNG
    │   ├─ Display pathway synergy metrics
    │   ├─ Visualize PPI connections
    │   └─ Show combined synergy score
    ├─ Build coexpression section
    │   ├─ Embed bulk correlation heatmap (HTML)
    │   ├─ Embed bulk expression boxplot (HTML)
    │   ├─ Embed SC heatmaps (5 cell types, HTML)
    │   └─ Display correlation statistics
    ├─ Generate combo-specific CI analysis
    │   ├─ Show Jaccard similarity
    │   ├─ Display union/intersection drug sets
    │   └─ Explain synergy adjustment
    └─ Append to report HTML

STEP 5: GENERATE COMPARATIVE ANALYSIS
  ├─ Create sortable rankings table (all targets)
  ├─ Generate Blue Ocean 2x2 scatter plot
  │   ├─ Plot each target as point
  │   ├─ Color by category
  │   └─ Add quadrant boundaries
  ├─ Create subscore heatmap (all targets × 5 subscores)
  └─ Generate recommendation summary

STEP 6: ADD METHODOLOGY SECTION
  ├─ Document all data sources
  ├─ Show all formulas with examples
  ├─ Explain normalization methods
  ├─ Justify weight choices
  └─ List limitations

STEP 7: ADD APPENDICES
  ├─ Export raw data as tables
  ├─ Include calculation examples
  ├─ Index all source files
  └─ Add glossary

STEP 8: FINALIZE REPORT
  ├─ Add navigation menu
  ├─ Add table of contents with links
  ├─ Embed JavaScript for interactivity
  │   ├─ Sortable tables
  │   ├─ Expandable sections
  │   └─ Chart interactions
  ├─ Optimize embedded visualizations
  └─ Write to HTML file

STEP 9: VALIDATION
  ├─ Verify all scores present
  ├─ Check all visualizations loaded
  ├─ Validate internal links
  └─ Log report generation summary

STEP 10: OUTPUT
  ├─ Save HTML report:
  │   └─ results/reports/IBD_Target_Prioritization_Report_{timestamp}.html
  ├─ Save machine-readable score bundles (required for reproducibility):
  │   ├─ results/reports/IBD_Target_Prioritization_Scores_Raw_{timestamp}.json
  │   └─ results/reports/IBD_Target_Prioritization_Scores_Normalized_{timestamp}.json
  └─ Save input manifest (what files were read):
      └─ results/reports/IBD_Target_Prioritization_Inputs_{timestamp}.json
```

---

## Gene Combinations

**Fixed combinations (use genes in format like "TYK2-JAK1" in reports):**
- **TYK2-JAK1** (formerly List 1)
- **TNFRSF25-GREM1** (formerly List 2)
- **TNFRSF25-PCOLCE** (formerly List 3)
- **CDKN2D-ITGA4-ITGB7** (formerly List 4)
- **CDKN2D-PCOLCE** (formerly List 5)

---

## Data Sources & File Locations

| Data Source | File Location | Parser Type | Raw Score Format | Extracted Fields |
|-------------|---------------|-------------|------------------|------------------|
| **Cortellis** | `results/cortellis/IBD_Target_Analysis_Report.md` | Markdown regex | Numeric (0-134.5) | Total score, total drugs, total trials, IBD drugs |
| **OFF-X Safety (preferred)** | `results/offx/gene_safety_scores_detailed.json` | JSON | Severity breakdown | Very high, High, Medium, Low, Very low, Not assoc, NA counts (+ total AEs) |
| **OFF-X Safety (fallback)** | `results/offx/OFF-X_Safety_Analysis_Report.md` | Markdown table | Severity breakdown | Very high, High, Medium, Low, Very low, Not assoc, NA counts |
| **DEG (preferred)** | `results/deg_results/IBD_Targets_Summary_Report.md` | Markdown table | Numeric (0-15.0) | Score, Disease↑, Responder↓, Studies |
| **DEG (optional structured)** | `results/deg_results/target_scores.csv` | CSV | Numeric (0-15.0) | Score, Disease↑, Responder↓, Studies |
| **BioBridge Individual** | `results/biobridge/individual_report/{gene}-IBD.md` | Markdown regex | Percentile % (0-100) | Cosine similarity, percentile rank |
| **BioBridge Combo** | Hardcoded in combo table | Python dict | Percentile % (0-100) | Combo-specific percentile ranks |
| **ULTRA** | `results/ultra/{gene}-predictions.json` | JSON (optional) | Percentile (0-1) | Percentile rank, top predictions |
| **PrimeKG** | `results/primekg/{gene}-IBD.md` | Markdown regex | Connection count | Number of disease connections |
| **Pathway Overlap** | `results/pathwaydb/REPORT_list{N}_{genes}.md` | Markdown parsing | Pathway count, score | Shared pathways, synergy interpretation |
| **Pathway Centrality (preferred)** | `results/pathwaydb/centrality_gene_{GENE}_centrality.csv` | CSV | Centrality metrics | Degree/PageRank/Eigenvector centrality (used for “Pathway Centrality” component) |
| **Pathway Individual** | `results/pathwaydb/gene_{GENE}_pathways.csv` | CSV | Pathway list | Individual gene pathways for context |
| **PPI Connections** | `results/interactdb/shortest_paths_analysis_complete/COMPLETE_ANALYSIS_REPORT_ALL_DATABASES.md` | Markdown parsing | Hop count, databases | Min hops, database consensus, paths |
| **Bulk Coexpression** | `results/bulk_coexpression/ANALYSIS_REPORT.md` | Markdown parsing | Correlation values | Mean correlation, p-values, expression data |
| **SC Coexpression** | `results/sc_coexp/ibd_analysis_results/COMPREHENSIVE_REPORT.md` | Markdown parsing | Correlation by cell type | Cell-type specific correlation, p-values |
| **Competitive Intel** | `results/ci/ibd_dashboard.html` | JSON from `<script id="data" type="application/json">` | Program entries | Competing programs, IBD phase, MOA data, target families |

### Visualization Assets (for report embedding)

- **Pathway network**: `results/pathwaydb/list{N}_{genes}.png`
- **Pathway UpSet inputs**: `results/pathwaydb/list{N}_{genes}_upset_set_sizes.csv` (can be rendered into a plot during report compilation)
- **Pathway MD reports**: `results/pathwaydb/REPORT_list{N}_{genes}.md`
- **PPI paths**: From markdown report text
- **Bulk correlation heatmap**: `results/bulk_coexpression/results/list{N}_{genes}_correlation_heatmap.html`
- **Bulk expression boxplot**: `results/bulk_coexpression/results/list{N}_{genes}_expression_boxplot.html`
- **SC correlation heatmaps**: `results/sc_coexp/ibd_analysis_results/{cell_type}/...heatmap.html` (5 cell types: B lineage, T cell, Myeloid, Mesenchymal, Endothelial)

---

## Score Data Source & Transformation Map (Quick Reference)

This section answers: “I need to compute X — what file do I open, what field do I extract, and how do I transform it?”

| Output (normalized unless noted) | Entity | Retrieve from | Extract (raw) | Transform scenario |
|---|---|---|---|---|
| Clinical Validation | Gene | `results/cortellis/IBD_Target_Analysis_Report.md` | `Total Score` (0–134.5) | `score/134.5*100` |
| Clinical Validation | Combo | `results/cortellis/IBD_Target_Analysis_Report.md` | (optional evidence) list-level `Combined Score` | Scoring uses the mean of component gene clinical validation scores; store list-level score as contextual raw evidence only |
| Safety | Gene | `results/offx/gene_safety_scores_detailed.json` | severity counts (+ `na`, `total_aes`) | risk-weighted formula (0–100); if all NA/empty ⇒ 50 |
| Safety | Combo | `results/offx/OFF-X_Safety_Analysis_Report.md` | (only if present) combo-level severity breakdown | Prefer combo-level breakdown if present; otherwise average component gene safety scores |
| Disease Association: DEG component | Gene | `results/deg_results/IBD_Targets_Summary_Report.md` | `Score` (0–15) | `score/15*100` |
| Disease Association: BioBridge component | Gene | `results/biobridge/individual_report/{GENE}-IBD.md` | `Percentile Rank` | already 0–100 |
| Disease Association: ULTRA component | Gene | `results/ultra/{GENE}-predictions.json` | `percentile` (0–1) | `percentile*100`; missing ⇒ 50 |
| Disease Association: Pathway Centrality component | Gene | `results/pathwaydb/centrality_gene_{GENE}_centrality.csv` | chosen metric (recommend `PageRank`) | normalize to 0–100 within current run; missing ⇒ 50 |
| Disease Association: PrimeKG component | Gene | `results/primekg/{GENE}-IBD.md` | connection count (heuristic) | normalize to 0–100 within current run |
| Opportunity: Competitive Intensity | Gene/Combo | `results/ci/ibd_dashboard.html` | per-asset IBD phases for targets | compute `% market left` per rules below |
| Synergy: Pathway Synergy | Combo | `results/pathwaydb/REPORT_list{N}_{genes}.md` | shared pathway count | divide by chosen max (list1=110) and cap at 1.0 |
| Synergy: PPI Synergy | Combo | `results/interactdb/shortest_paths_analysis_complete/COMPLETE_ANALYSIS_REPORT_ALL_DATABASES.md` | `min_hops`, databases found | `ppi_synergy=(1/min_hops)*(db_found/3)` |

---

## Transformation Scenarios (Explicit Precedence and Fallback Rules)

To keep report generation deterministic and debuggable, apply these precedence rules consistently:

1. **Structured sidecar precedence**
   - If a structured sidecar exists for numeric values (CSV/JSON), treat it as canonical.
   - Use Markdown/HTML parsing only as fallback.

2. **Combination precedence**
   - Prefer list/combo-level values if they exist (Cortellis list table, OFF-X combo rankings, pathway list reports, PPI list sections).
   - If combo-level values are absent, fall back to aggregating component genes (usually mean).

3. **Missing/NA handling**
   - Missing optional source (ULTRA): set component to 50.0 and record `missing=true` in raw payload.
   - All-NA safety: set Safety to 50.0 and label “Uncharacterized”.

4. **Normalization set definition**
   - Any “max_in_dataset” normalization must use the **current run’s entity set** (the genes you are scoring).
   - Persist the `min/max` used so numbers can be reproduced from raw payloads.

5. **Gene symbol resolution (aliases and file stems)**
   - Canonical symbols for scoring/output should match the majority of numeric sources (Cortellis/DEG/OFF-X), e.g.:
     - `TNFRSF25` (canonical) ↔ `DR3` (BioBridge/PrimeKG file stems in current results set)
   - Resolution rule during loading:
     - Try canonical stem first (e.g., `TNFRSF25-IBD.md`)
     - If missing, try known aliases (e.g., `DR3-IBD.md`)
     - If neither exists, treat the component as missing and apply the standard fallback (`50.0`) with provenance noting the missing files.

6. **Competitive intel target matching (symbol extraction)**
   - CI entries store targets as strings like `"Janus kinase 1 (JAK1)"`.
   - Matching rule:
     - Extract gene symbols from parentheses when present (`(...)`), otherwise fall back to word-boundary search for the canonical symbol.
   - Persist in raw payload which CI target strings matched which gene(s) to keep the computation auditable.

## Scoring System Design

### Overall Score Weights (Both Individual & Combination)

- **Clinical Validation:** 30%
- **Disease Association:** 30%
- **Opportunity:** 20%
- **Novelty:** 10%
- **Safety:** 10%

**Total:** 100%

---

## Subscore Specifications

### 1. Clinical Validation Score

**Purpose:** Measure clinical evidence and market validation

**Raw Data Source:** Cortellis Report
**Raw Score Range:** 0-134.5 (JAK1 = max)

**Normalization to 0-100:**
```python
normalized_score = (raw_score / 134.5) × 100
```

**For Combinations:**
- Average of individual gene scores

**Display:**
- Raw score: Show original value (e.g., "134.5")
- Normalized: Show 0-100 for tables/plots (e.g., "100.0")

---

### 2. Disease Association Score

**Purpose:** Measure strength of gene-disease relationship

**Component Sources & Weights:**
- **DEG:** 30%
- **BioBridge:** 25%
- **ULTRA:** 20%
- **Pathway Centrality:** 15%
- **PrimeKG:** 10%

**Component Normalization (each to 0-100):**

1. **DEG:**
   - Raw range: 0-15.0 (PCOLCE = max in dataset)
   - Normalized: `(raw_score / 15.0) × 100`

2. **BioBridge:**
   - Raw: Percentile rank 0-100
   - Normalized: `percentile_rank × 1.0` (already 0-100)

3. **ULTRA:**
   - Raw: Percentile rank 0-1
   - Normalized: `percentile_rank × 100`
   - Fallback: If missing, use 50.0 in weighted calculation

4. **Pathway Centrality:**
   - Raw: Centrality metrics from `results/pathwaydb/centrality_gene_{GENE}_centrality.csv`
   - Transformation scenario (deterministic):
     - Choose a canonical metric (recommended: `PageRank`)
     - Compute `pathway_centrality_raw = mean(PageRank)` over rows where `Is_Seed == True`
       - If `Is_Seed` is missing/unusable, compute mean over the top `N=25` rows by `Degree`
   - Normalization:
     - Normalize within the current run’s gene set:
       - `norm = (raw - min_raw) / (max_raw - min_raw) * 100`
     - If `max_raw == min_raw`, set `norm = 50.0` for all genes
   - Default: `50.0` (neutral) if file missing or parse fails
   - Higher centrality = more central to disease pathways

5. **PrimeKG:**
   - Raw: Connection count (varies by gene)
   - Normalized: `(connection_count / max_in_dataset) × 100`
   - Max determined dynamically from current gene set

**Final Calculation:**
```python
disease_assoc = (DEG_norm × 0.30) + (BioBridge_norm × 0.25) +
                (ULTRA_norm × 0.20) + (PathwayCentrality_norm × 0.15) +
                (PrimeKG_norm × 0.10)
```

**For Combinations:**
- **DEG:** Average of individual gene scores
- **BioBridge:** Use combo percentile_rank directly (NOT averaged from individuals)
- **ULTRA:** Average of individual gene scores
- **Pathway Centrality:** Average of individual gene scores
- **PrimeKG:** Average of individual gene scores

---

### 3. Safety Score

**Purpose:** Assess safety profile based on adverse event severity

**Raw Data Source:** OFF-X Report severity breakdown
**Input:** Count of events by severity category

**Risk-Weighted Calculation with Amplified Penalties:**

```python
# Severity → Score × Amplification Factor
weighted_sum = (count_VeryHigh × 0 × 10) +    # Most dangerous
               (count_High × 10 × 5) +         # Serious
               (count_Medium × 20 × 2) +       # Moderate
               (count_Low × 60 × 1) +          # Minor
               (count_VeryLow × 80 × 1) +      # Minimal
               (count_NotAssoc × 100 × 1)      # Safe

weighted_total = (count_VeryHigh × 10) +
                 (count_High × 5) +
                 (count_Medium × 2) +
                 (count_Low × 1) +
                 (count_VeryLow × 1) +
                 (count_NotAssoc × 1)

safety_score = weighted_sum / weighted_total
```

**Severity Mapping:**
- Very high: 0 points (10× amplification)
- High: 10 points (5× amplification)
- Medium: 20 points (2× amplification)
- Low: 60 points (1× amplification)
- Very low: 80 points (1× amplification)
- Not associated: 100 points (1× amplification)

**Special Cases:**
- If NA/empty data: Use 50.0 (neutral)
- Already on 0-100 scale, no additional normalization needed

**For Combinations:**
- Prefer OFF-X combo data if available
- Fallback: Average of individual gene safety scores

---

### 4. Opportunity Score

**Purpose:** Identify targets with high disease relevance but low competitive intensity

**Components:**
- **Disease Association (normalized):** 50%
- **Clinical Novelty:** 30%
- **Competitive Intelligence:** 20%

**Calculation Steps:**

1. **Clinical Novelty:**
   ```python
   clinical_novelty = 100 - clinical_validation_normalized
   ```
   (Less validated = more opportunity)

2. **Competitive Intelligence Score:**

   **For Individual Genes:**
   ```python
   # Parse CI data from HTML to get drugs targeting this gene for IBD
   weighted_competition = (count_Marketed × 1.0) +
                         (count_PhaseIII × 0.7) +
                         (count_PhaseII × 0.4) +
                         (count_PhaseI × 0.2)

   total_programs = count_Marketed + count_PhaseIII + count_PhaseII +
                   count_PhaseI + count_Preclinical

   pct_marketed = (count_Marketed / total_programs) × 100
   pct_market_left = 100 - pct_marketed

   CI_score = pct_market_left  # Higher = more opportunity
   ```

   **For Combinations:**
   ```python
   # Step 1: Calculate individual gene scores
   gene_A_comp = weighted_competition(gene_A)
   gene_B_comp = weighted_competition(gene_B)

   # Step 2: Calculate union score
   union_drugs = drugs_targeting_A OR drugs_targeting_B (for IBD)
   union_comp = weighted_competition(union_drugs)

   # Step 3: Calculate Jaccard similarity
   intersection = drugs_targeting_A AND drugs_targeting_B
   union = drugs_targeting_A OR drugs_targeting_B
   jaccard = len(intersection) / len(union)

   # Step 4: Combine with Jaccard weighting
   average_comp = (gene_A_comp + gene_B_comp) / 2
   combo_comp_raw = average_comp × (1 - jaccard) + union_comp × jaccard

   # Step 5: Apply synergy adjustment (reduces competition for synergistic combos)
   if pathway_synergy > 0.5:
       adjustment_factor = 1 - (pathway_synergy × 0.3)  # Up to 30% reduction
       combo_comp_adjusted = combo_comp_raw × adjustment_factor

   # Step 6: Calculate % market left
   pct_market_left_combo = 100 - (marketed_in_union / total_union × 100)
   CI_score_combo = pct_market_left_combo
   ```

3. **Final Opportunity Score:**

   **For Individual Genes:**
   ```python
   opportunity = (disease_assoc_norm × 0.5) +
                (clinical_novelty × 0.3) +
                (CI_score × 0.2)
   ```

   **For Combinations:**
   ```python
   # Calculate individual opportunities
   mean_individual_opps = (opp_gene_A + opp_gene_B) / 2

   # Novel mechanism bonus
   novel_mechanism_bonus = (100 - combined_synergy × 100)

   # Combine
   opportunity_combo = (mean_individual_opps × 0.6) +
                      (novel_mechanism_bonus × 0.4)
   ```

---

### 5. Novelty Score

**Purpose:** Identify understudied targets with less literature coverage

**Components:**
- **Clinical novelty:** 70%
- **Literature novelty:** 30%

**Calculation:**

1. **Clinical Novelty:**
   ```python
   clinical_novelty = 100 - clinical_validation_normalized
   ```
   (Less clinically validated = more novel)

2. **Literature Novelty:**
   ```python
   # From PrimeKG connection count
   literature_percentile = min((connection_count / 50) × 100, 100)
   literature_novelty = 100 - literature_percentile
   ```
   (Fewer connections = less studied = more novel)

**Final Novelty Score:**

**For Individual Genes:**
```python
novelty = (clinical_novelty × 0.7) + (literature_novelty × 0.3)
```

**For Combinations:**
```python
mean_individual_novelty = (novelty_gene_A + novelty_gene_B) / 2
novelty_combo = mean_individual_novelty × 1.10  # 10% combination premium
```

---

## Synergy Metrics (Combinations Only)

### Pathway Synergy

**Purpose:** Measure shared regulatory mechanisms

**Data Source:** Pathway overlap report
**Calculation:**
```python
pathway_synergy = shared_pathways / 110  # Normalized to List 1 max
pathway_synergy = min(pathway_synergy, 1.0)  # Cap at 1.0
```

### PPI Synergy

**Purpose:** Measure protein-protein interaction connectivity

**Data Source:** PPI complete analysis report
**Calculation:**
```python
hop_score = 1.0 / min_hops  # Closer = higher
database_consensus = found_databases / 3  # Max 3: STRING, BioGRID, IntAct
ppi_synergy = hop_score × database_consensus
ppi_synergy = min(ppi_synergy, 1.0)  # Cap at 1.0
```

### Combined Synergy

```python
combined_synergy = (pathway_synergy × 0.5) + (ppi_synergy × 0.5)
```

---

## Competitive Intelligence - Blue Ocean Categorization

### 2x2 Matrix Framework

**Purpose:** Project targets onto opportunity/competitiveness landscape

**Axes:**
- **X-axis: Opportunity** (from Opportunity subscore)
- **Y-axis: Competitiveness** (from CI analysis)

**Category Definitions:**

1. **Blue Ocean** (Recommended)
   - High Opportunity (Opportunity subscore > 50)
   - Low Competition (% market left > 25%)
   - Interpretation: High disease relevance, low competitive intensity

2. **Crowded Leader** (Validated but competitive)
   - High Opportunity (Opportunity subscore > 50)
   - High Competition (% market left ≤ 25%)
   - Interpretation: Validated targets with many competitors

3. **Red Ocean** (Avoid)
   - Low Opportunity (Opportunity subscore ≤ 50)
   - High Competition (% market left ≤ 25%)
   - Interpretation: Crowded space with limited upside

4. **Avoid** (Low priority)
   - Low Opportunity (Opportunity subscore ≤ 50)
   - Low Competition (% market left > 25%)
   - Interpretation: Limited disease relevance despite open space

**Visualization:**
```
    High Opp
        ↑
        │  Blue Ocean   │ Crowded Leader
        │               │
  ──────┼───────────────┼────────→ High Comp
        │               │
        │    Avoid      │   Red Ocean
        │               │
```

---

## Data Output Format

### Raw Scores (for traceability)

Store in result dictionary with `_raw` suffix:
```python
{
    "gene": "JAK1",
    "clinical_validation_raw": 134.5,
    "deg_score_raw": 9.0,
    "biobridge_percentile_raw": 90.3,
    "safety_breakdown_raw": {
        "very_high": 39,
        "high": 118,
        "medium": 404,
        ...
    },
    "safety_score_raw": 76.8,
    ...
}
```

### Normalized Scores (for display)

Store with `_normalized` suffix or direct name for subscores:
```python
{
    "gene": "JAK1",
    "clinical_validation": 100.0,      # Displayed score
    "disease_association": 82.5,       # Displayed score
    "safety": 76.8,                    # Displayed score
    "opportunity": 65.2,               # Displayed score
    "novelty": 15.4,                   # Displayed score
    "overall_score": 78.3,             # Weighted combination
    ...
}
```

### Display Format

**Tables/Reports:** Always show normalized 0-100 scores
**Traceability sections:** Show raw scores with units
**Visualizations:** Use normalized scores for comparability

---

## Implementation Notes

### Error Handling

1. **Missing files:** Log warning, use default value 50.0
2. **Parse errors:** Log error with file path, skip that data source
3. **Missing genes:** Return empty result with all scores = 50.0
4. **Malformed data:** Validate and sanitize, use defaults if invalid

### Validation

1. **Score ranges:** Assert all normalized scores 0-100
2. **Weight sums:** Assert all weight sets sum to 1.0
3. **Data completeness:** Log % of data sources successfully loaded
4. **Consistency checks:** Verify raw scores can be denormalized

### Performance

1. **Single load pass:** Load all data once in orchestrator
2. **Lazy normalization:** Only normalize when needed for calculation
3. **Cache parsed results:** Avoid re-parsing same files
4. **Parallel loading:** Consider concurrent file reads for large datasets

---

## Next Steps

1. ✅ Design complete
2. ⏳ Implement parsers.py (markdown/HTML/JSON parsing utilities)
3. ⏳ Implement loaders.py (individual data source loaders)
4. ⏳ Implement orchestrator.py (load_all_data master function)
5. ⏳ Implement normalizers.py (raw → 0-100 conversion)
6. ⏳ Implement subscores.py (5 subscore calculators)
7. ⏳ Implement scoring.py (overall scoring + CI categorization)
8. ⏳ Add unit tests for all components
9. ⏳ Integration testing with actual data files
10. ⏳ Update HTML report generator to use new system

---

## Revision History

| Date | Version | Changes |
|------|---------|---------|
| 2025-12-23 | 1.0 | Initial design specification |
