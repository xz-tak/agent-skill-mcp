# Commercial Scoring Deep-Research Prompt Template

This is the standard prompt template for scoring ANY combo in ANY indication. Variables in `{braces}` are filled dynamically from `pipeline_config.json` + CI dashboard data.

## Usage

When running commercial scoring for a new indication/combo set:
1. Load `pipeline_config.json` for combos, indication, diseases
2. Parse CI HTML dashboard for pipeline context per gene
3. Fill this template per combo
4. Run `/deep-research` with the filled prompt
5. Parse structured JSON output into `commercial_scores.json`

---

## Template

```
Research the commercial viability of a {combo_name} combination therapy for {indication_full}.

## Assessment Criteria

Score each dimension on a continuous 0-100 scale (higher = better opportunity / lower risk).

### 1. MARKET OPPORTUNITY (weight: 40%)
- Total addressable market size for {indication_full}
- Patient population (prevalent + incident), growth trajectory
- Unmet need: what current SOC fails to address that this combo could
- Reimbursement/payer landscape considerations

### 2. COMPETITIVE PROFILE (weight: 40%)
- Existing marketed therapies targeting {gene_list} (monotherapy or combo)
- Pipeline programs with similar mechanisms (phase, timeline)
- Patent/exclusivity landscape
- Differentiation: what makes this combo unique vs existing options
- Biosimilar/generic risk for component monotherapies

### 3. STRATEGIC FIT (weight: 20%)
- Alignment with GI/immunology therapeutic area focus
- Leverageability of existing programs and capabilities
- (Score modality-agnostic, only on TA relevance)

## Local Pipeline Context (from CI Dashboard)
{ci_context}

## Cortellis Summary
{cortellis_context}

## Required Output Format
Return a JSON object with this exact schema:
{
  "combo": "{combo_name}",
  "indication": "{indication}",
  "commercial_score": <float 0-100>,
  "confidence": "<High|Medium|Low>",
  "components": {
    "market_opportunity": {
      "score": <float 0-100>,
      "rationale": "<1-2 sentence justification>",
      "key_evidence": ["<evidence point 1>", "<evidence point 2>"]
    },
    "competitive_profile": {
      "score": <float 0-100>,
      "rationale": "<1-2 sentence justification>",
      "key_evidence": ["<evidence point 1>", "<evidence point 2>"]
    },
    "strategic_fit": {
      "score": <float 0-100>,
      "rationale": "<1-2 sentence justification>",
      "key_evidence": ["<evidence point 1>", "<evidence point 2>"]
    }
  },
  "citations": [
    {"source": "<url or reference>", "relevance": "<what it supports>"}
  ],
  "narrative_summary": "<one paragraph executive summary>"
}

## Scoring Guidance

- Market: {indication_full} is a {market_hint} market. Calibrate thresholds dynamically for the disease area.
- Competition: score LOW (0-30) if highly crowded with little differentiation, HIGH (70-100) if first-in-class with clear advantages.
- Strategic: GI/immunology core = 80-100, adjacent autoimmune = 60-79, distant TA = 20-59.
- commercial_score = 0.40*market + 0.40*competitive + 0.20*strategic

CRITICAL — Scores must DIFFERENTIATE between combos within the same indication:
- All combos share the same market/indication, so Market scores may be similar.
- Competitive Profile is the PRIMARY differentiator:
  * Novel mechanism with no direct competitors → HIGH (75-95)
  * ONE validated + ONE novel target → MODERATE (50-75)
  * Both targets highly validated and crowded → LOWER (30-55)
- Use CI pipeline drug counts as quantitative grounding for competition scoring.
- Combos with targets that have many Phase III+ programs = MORE competition = LOWER score.
- Combos with targets that have NO pipeline = HIGHER novelty = HIGHER score.
```

---

## Template Variables

| Variable | Source | Example |
|----------|--------|---------|
| `{combo_name}` | `" + ".join(combo)` from config | `TYK2 + JAK1` |
| `{indication}` | `config.diseases.primary` | `IBD` |
| `{indication_full}` | primary + subtypes | `IBD (Crohn's disease, Ulcerative colitis)` |
| `{gene_list}` | `", ".join(combo)` | `TYK2, JAK1` |
| `{ci_context}` | Extracted from CI HTML dashboard per gene | `TYK2: 3 active drugs\n  - zasocitinib (Phase II)...` |
| `{cortellis_context}` | From `results/cortellis_<ind>/<GENE>_summary.md` | First 10 lines of each gene's Cortellis summary |
| `{market_hint}` | Dynamic from indication (or "to be determined") | `large ($25B+ global)` |

## Indication Market Hints (dynamic, agent-determined)

The template includes a `{market_hint}` placeholder. The prep script provides a starting estimate, but the research agent should dynamically calibrate based on actual findings. Known anchors:

- IBD: large ($25B+ global)
- UC/CD: large (IBD subset, $15B+)
- SSc: moderate-to-small ($2-5B)
- IPF: moderate ($5-8B)
- NASH/MASH: large (projected $20B+)
- RA: large ($20B+)
- Any other: agent determines from research

For unknown indications, the agent should research market size as part of the Market Opportunity assessment and calibrate accordingly.
