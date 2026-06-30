# Commercial Risk Assessment Scoring Contract

## Overview

Commercial score assesses the market viability and strategic opportunity for a target combination in a given indication. The score is **0–100** (higher = better opportunity / lower risk).

## Scoring Formula

```
Commercial = 0.40 × Market_Opportunity + 0.40 × Competitive_Profile + 0.20 × Strategic_Fit
```

### Weights
| Checkpoint | Weight | Rationale |
|-----------|--------|-----------|
| Market Opportunity | 0.40 | Quantifiable market size drives commercial potential |
| Competitive Profile | 0.40 | Differentiation determines probability of commercial success |
| Strategic Fit | 0.20 | Alignment with therapeutic area focus |

## Checkpoint Scoring (Continuous 0–100, NOT binned)

Scores are continuous within each checkpoint. The guideline tiers below provide anchor points; the deep-research agent should interpolate based on evidence quality and specifics.

### 1. Market Opportunity (0–100)

Assesses the total addressable market (TAM) for the combination product in the target indication.

**Dynamic per-indication calibration:** Thresholds are NOT fixed — the research agent determines what constitutes "large" vs "small" market relative to the disease area. Example calibration for IBD (high-prevalence autoimmune): a $25B+ global market means even a niche combo can be commercially attractive.

**Scoring anchors (indication-adjusted):**
| Score Range | Evidence Pattern |
|------------|-----------------|
| 80–100 | Large unmet need in sizable patient population; limited effective options; growing incidence; clear reimbursement pathway |
| 60–79 | Moderate-to-large market with identifiable unmet need segments; some established SOC but room for improvement |
| 40–59 | Moderate market with partial unmet need; standard treatments exist but combo may offer incremental benefit |
| 20–39 | Small market or largely addressed unmet need; combo offers marginal improvement over existing treatments |
| 0–19 | Very small/saturated market; no clear unmet need for the combo mechanism |

**Evidence to gather:**
- Total market size (global, by geography)
- Patient population estimate (prevalent + incident)
- Current SOC and identified gaps
- Projected growth rate
- Reimbursement landscape

### 2. Competitive Profile (0–100)

Assesses differentiation potential and competitive intensity for the target combination.

**Scoring anchors:**
| Score Range | Evidence Pattern |
|------------|-----------------|
| 80–100 | First-in-class combo mechanism; no direct competitors; clear differentiation on efficacy/safety/convenience |
| 60–79 | Novel combo with few competitors; meaningful differentiation vs monotherapy SOC; limited pipeline overlap |
| 40–59 | Some competing combos in development; moderate differentiation; need to compete on specific endpoints |
| 20–39 | Crowded space with multiple similar combos; limited differentiation; biosimilar/generic risk for mono components |
| 0–19 | Highly saturated; direct combo competitors in late-stage; no credible differentiation path |

**Evidence to gather:**
- Competing programs (same targets, same combo concept)
- Phase advancement of competitors
- Patent landscape and exclusivity windows
- Differentiation opportunities (efficacy, safety, route, dosing)
- Biosimilar/generic timeline for component monotherapies

### 3. Strategic Fit (0–100)

Assesses alignment with therapeutic area (TA) focus. **Modality-agnostic** — scored on TA relevance only.

**Scoring anchors:**
| Score Range | Evidence Pattern |
|------------|-----------------|
| 80–100 | Core TA (GI, immunology, inflammation); active programs in the space; established KOL networks |
| 60–79 | Adjacent TA with clear biological connection to core areas; some internal expertise |
| 40–59 | Partially related TA; would require moderate capability build-up |
| 20–39 | Distant TA; limited internal expertise; significant capability gap |
| 0–19 | Completely outside current focus; no pathway to leverage existing capabilities |

**For IBD/GI-focused organizations:**
- Core: IBD, celiac, GI fibrosis, GI inflammation → 80–100
- Adjacent: other autoimmune (RA, psoriasis, dermatitis), liver disease → 60–79
- Partial: respiratory inflammation, renal → 40–59
- Distant: oncology, CNS, rare metabolic → 20–39

## Research Protocol

### Input Context (Hybrid: local + web)

Before deep-research, feed the agent:
1. **Cortellis data** (from `results/cortellis_<indication>/`): clinical phase, competing drugs, indication pipeline
2. **CI dashboard** (from `results/ci_<indication>/`): competitive intelligence, market landscape
3. **Gene approach** annotations (agonist/antagonist context)
4. **Combo mechanism** summary (what the combo targets mechanistically)

### Deep-Research Query Template

For each combo, the research agent should investigate:
```
Research the commercial viability of a [{COMBO_GENES joined by ' + '}] combination therapy 
for [{INDICATION}]. Specifically assess:

1. MARKET OPPORTUNITY: Total addressable market size, patient population, unmet need, 
   growth trajectory, reimbursement considerations for {INDICATION}.
   
2. COMPETITIVE PROFILE: Existing and pipeline combination therapies targeting similar 
   mechanisms ({GENE1}, {GENE2}, etc.), differentiation potential, patent landscape, 
   time-to-market considerations.
   
3. STRATEGIC FIT: Alignment with GI/immunology therapeutic focus, leverageability of 
   existing capabilities and programs.

Context from internal pipeline data:
- Clinical stages: {CORTELLIS_SUMMARY}
- Known competitors: {CI_SUMMARY}
```

### Output Schema (per combo)

```json
{
  "combo": "GENE1 + GENE2",
  "indication": "IBD",
  "commercial_score": 72.5,
  "confidence": "Medium",
  "components": {
    "market_opportunity": {
      "score": 85.0,
      "rationale": "IBD market >$25B globally with significant unmet need in moderate-severe patients...",
      "key_evidence": ["Market projected to reach $X by 20XX", "Y% patients inadequately treated"]
    },
    "competitive_profile": {
      "score": 62.0,
      "rationale": "Novel combo mechanism with limited direct competition, but established monotherapy competitors...",
      "key_evidence": ["No direct combo competitor in Phase 3", "Monotherapy biosimilar expected 20XX"]
    },
    "strategic_fit": {
      "score": 90.0,
      "rationale": "Core GI/immunology focus with active programs...",
      "key_evidence": ["Active IBD pipeline", "Established GI KOL network"]
    }
  },
  "citations": [
    {"source": "...", "relevance": "market size"},
    {"source": "...", "relevance": "competitor pipeline"}
  ],
  "narrative_summary": "One paragraph executive summary of commercial risk assessment..."
}
```

### Confidence Levels

| Level | Criteria |
|-------|----------|
| High | Multiple analyst reports, peer-reviewed market data, confirmed competitor pipeline from ClinicalTrials.gov |
| Medium | Some published market data, pipeline intelligence partially confirmed, analyst consensus available |
| Low | Limited public data, estimates extrapolated from related indications, early-stage competitive intelligence |

## Skill Integration Workflow

Commercial scoring is integrated into the `/target-prioritization-report` skill as follows:

### Step 1: Extract CI Context
- Parse CI HTML dashboard (`results/ci_<indication>/*.html`) for the embedded `<script id="data">` JSON
- Map each gene to its active drugs, phases, organizations, and modalities
- Build per-combo CI summary (genes → active pipeline programs)

### Step 2: Build Research Prompts
- For each combo, construct a research prompt using the template above
- Inject local context (Cortellis summary + CI pipeline data)
- Include indication-specific calibration hints

### Step 3: Execute Deep Research
- Use `/deep-research` skill (or workflow fan-out) per combo
- Each research agent scores Market/Competition/Strategic per the anchors
- Agent outputs structured JSON per the output schema

### Step 4: Aggregate & Save
- Collect all combo scores into `results/commercial_<indication>/commercial_scores.json`
- Compute gene-level commercial = mean of combos containing gene
- Report generator reads this JSON during report generation

### CI Dashboard Input Specification

The CI HTML dashboard must contain:
```html
<script id="data" type="application/json">
{
  "meta": {"sourceFile": "...", "count": N},
  "entries": [
    {
      "displayName": "drug name",
      "targets": ["target name (GENE)"],
      "ibdPhase": "Phase III",
      "overallPhase": "Marketed",
      "organization": "...",
      "underActiveDevelopment": "Yes",
      "modality": "..."
    }
  ]
}
</script>
```

Target name matching: case-insensitive substring match against gene symbols in the combo list.

## Integration with Report

- Commercial score flows directly into the Overall formula: `0.10 × Commercial`
- Score table shows Commercial column
- Each combo card includes a "Commercial" section with:
  - **Score bars:** 3 horizontal bars for Market/Competition/Strategic (color-coded 0-100)
  - **Bullet points:** 3-5 key findings from the research
  - **Narrative:** Executive summary paragraph
  - **Sub-tabs:** Overview | Market | Competition | Strategic Fit (detailed breakdown)
  - **Confidence badge:** High/Medium/Low with color coding
  - **Citations:** Source links used in scoring

## Score Differentiation (CRITICAL)

Scores MUST meaningfully differentiate between combos within the same indication:

- **Market Opportunity** will be similar across combos (same indication = same market). Minor differences arise from mechanism-specific patient segments.
- **Competitive Profile** is the PRIMARY differentiator:
  - Combos targeting crowded mechanisms (many Phase III+ drugs) → LOWER score (30-55)
  - Combos with ONE validated + ONE novel target → MODERATE (50-75)
  - Combos targeting novel mechanisms with no direct competitors → HIGH (75-95)
  - Use CI pipeline drug counts as quantitative grounding
- **Strategic Fit** varies if combos span different TAs or involve non-core modalities.

Example spread for IBD combos:
- TYK2 + JAK1 (both highly validated, crowded): Commercial ~45
- TNFRSF25 + GREM1 (both novel, little competition): Commercial ~78
- GLP2R + IL23A (one novel agonist + one validated): Commercial ~62

The agent must avoid "score compression" where all combos cluster around 50-60.

## Combo-Level vs Gene-Level

- **Combo:** Deep-research evaluates the combo AS A PRODUCT (not individual genes). Results stored in `test/commercial/<combo_slug>/commercial_result.json`.
- **Gene:** Individual target commercial assessment via dedicated deep research agent (Step 3b). Results stored in `test/commercial/genes/<GENE>/commercial_result.json`.
  - **Loading priority:** Per-gene result file → fallback: mean of all combo scores containing that gene.
  - **Same scoring formula:** 0.40×Market + 0.40×Competitive + 0.20×Strategic, but evaluates the single gene as a standalone target.
  - **Same output schema:** identical JSON structure with `"combo"` field containing the single gene name.

## Multi-Run Consensus (3 runs per entity)

Each combo AND each gene is assessed **3 independent times** by separate agents. Final merged result:

- **Scores:** Mean across 3 runs (commercial_score + each component score)
- **Evidence:** Union of all non-overlapping key_evidence items across runs (deduped by exact string)
- **Citations:** Union of all citations across runs (deduped by `source` field)
- **Narrative:** From the run with the highest commercial_score (most confident)
- **Confidence:** Majority vote across 3 runs

Per-run files: `commercial_result_run1.json`, `commercial_result_run2.json`, `commercial_result_run3.json`
Merged file: `commercial_result.json` (auto-merged by loader if runs exist but merged file doesn't)

Benefits:
- Reduces single-agent scoring bias
- Captures broader evidence from independent web searches
- `score_std` field in merged result indicates scoring agreement (low = high consensus)
