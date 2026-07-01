# Scoring Contract (Abstract)

This is the documented contract your loaders + scorers should implement so the report remains traceable and consistent across indications.

## Key principles

- Preserve raw values in `data_dict` (`*_raw`) with a `source_file`.
- Normalize to 0–100 only at scoring time (or for display).
- Missing sources degrade gracefully to `50.0` (neutral), and should be visible in traceability output.

## Heatmap Color Specification (REQUIRED)

All scoring heatmaps MUST use a consistent color gradient:

| Score | Color | Description |
|-------|-------|-------------|
| 100 | Green (#22c55e) | High/Good |
| 50 | Yellow (#eab308) | Neutral |
| 0 | Red (#ef4444) | Low/Bad |

Use a continuous gradient interpolation: `red → yellow → green` (0 → 50 → 100).

CSS implementation:
```css
.score-cell {
  background: linear-gradient(90deg,
    rgb(239, 68, 68) 0%,    /* red at 0 */
    rgb(234, 179, 8) 50%,   /* yellow at 50 */
    rgb(34, 197, 94) 100%   /* green at 100 */
  );
}
```

JavaScript color interpolation (HSL hue-based, matching GI2 mockup):
```javascript
function scoreHue(score){
  let t = Math.max(0, Math.min(100, score)) / 100;
  if (t <= 0.5) { const tt = t / 0.5; return 0 + (60 * tt); }  // red (0°) → yellow (60°)
  const tt = (t - 0.5) / 0.5;
  return 60 + (60 * tt);  // yellow (60°) → green (120°)
}
// Apply: background = linear-gradient(180deg, hsla(hue, 88%, 86%, .92), hsla(hue, 88%, 74%, .78))
```

**Required implementation:** Linear HSL hue interpolation: 0→60→120 (red→yellow→green). No smoothstep easing. This matches the GI2 portfolio mockup exactly.

## Overall score

Default weighted sum (6 components):

`Overall = 0.25·Clinical + 0.25·Biology + 0.125·Safety + 0.125·Druggability + 0.125·Translation + 0.125·Commercial`

## Clinical Validation (0–100)

- Input: raw clinical score from a clinical report table/section.
- Normalize: `Clinical = clamp((raw / MAX_RAW)·100, 0, 100)`
- Combo: average member-gene raw, then normalize.

## Biology (0–100)

Components (example set; keep consistent in code + report):

- DEG evidence (**0.20**)
- Signature evidence (**0.20**)
- BioBridge percentile (**0.12**)
- ULTRA percentile (**0.12**)
- PrimeKG connectivity (path-based score) (**0.06**)
- GSP (Genetics Support Profiler) (**0.15**)
- UKBPPP (UKB Pharma Proteomics Project) (**0.15**)

Recommended biology score:

`Biology = 0.20·DEG + 0.20·Signature + 0.12·BioBridge + 0.12·ULTRA + 0.06·PrimeKG + 0.15·GSP + 0.15·UKBPPP`

Typical combo rule:
- **DEG (combo):** Average member-gene normalized DEG scores.
- **BioBridge (combo):** Use combo-level `pct_rank` from Part2_Combo_Analysis.md when available; fallback = mean of member-gene BioBridge scores.
- **ULTRA (combo):** Use `geo_pct_rank` from Part2_Combo_Analysis.md × 100 when available; fallback = mean of member-gene ULTRA scores. The `geo_pct_rank` is the geometric mean of individual gene ranks converted to percentile space — it reflects expected combo performance without relying on intersection queries that often show dilution.
- **PrimeKG (combo):** Use combo-level score from Part2/JSON when available; fallback = mean of member-gene PrimeKG scores.

**Important distinction (BioBridge vs ULTRA combo scoring):**
- BioBridge combo uses the **direct combo `pct_rank`** (from the intersection/mean-embedding query result), because BioBridge's mean-embedding approach typically gives meaningful combo-level scores.
- ULTRA combo uses **`geo_pct_rank`** (geometric mean of individual ranks), because ULTRA's intersection queries frequently show dilution artifacts (combo weaker than expected) due to query specificity constraints.

### DEG normalization (UPDATED)

For combo scoring, use `target_score_table.csv` total_score per target gene (not per-comparison weighted_score).

**Raw inputs:**
- `deg_score_raw`: differential expression score (can be positive or negative)
- `gene_approach`: therapeutic approach annotation (`"antagonist"` or `"agonist"`, default: `"antagonist"`)

**Approach-based transformation (Step 1):**
```python
if gene_approach == "agonist":
    adjusted_raw = -deg_score_raw  # negate for agonist targets
else:
    adjusted_raw = deg_score_raw   # antagonist (default): no change
```

**Dynamic normalization (Step 2):**
```python
# Compute max from all adjusted scores in the dataset (dynamic, not hardcoded)
max_deg_score = max(abs(s) for s in all_adjusted_scores)

# Normalize without clamping
DEG_normalized = (adjusted_raw / max_deg_score) * 100
```

**Key changes from previous implementation:**
- Use **dynamic max** (computed from actual data) instead of hardcoded `15.0`
- **No clamping** to [0, 100] - scores reflect actual distribution
- **Agonist approach** negates raw score before normalization
- **Default approach** is `"antagonist"` (no modification)

**Per-gene approach annotation format:**
```json
{
  "TYK2": {"approach": "antagonist"},
  "JAK1": {"approach": "antagonist"},
  "GLP2R": {"approach": "agonist"},
  ...
}
```

If no approach is specified for a gene, assume `"antagonist"`.

### Signature normalization (NEW)

**Raw inputs:**
- Source: `{sig_name}_score_table.csv` in each combo directory
- Per combo: for each target gene's signature, compute mean `total_score` across all signature genes
- `gene_approach`: therapeutic approach annotation (same as DEG)

**Approach-based transformation (Step 1):**
```python
if gene_approach == "agonist":
    adjusted_sig = -mean_sig_score  # negate for agonist targets
else:
    adjusted_sig = mean_sig_score
```

**Cross-combo normalization (Step 2):**
```python
# Compute max from all adjusted signature scores across ALL combos
max_sig_score = max(abs(s) for s in all_combo_adjusted_sig_scores)

# Normalize
Signature_normalized = (adjusted_sig / max_sig_score) * 100
```

**File format:** Same columns as `target_score_table.csv` (Gene, total_score, n_comparisons, n_significant, etc.) plus `associated_target` column linking each signature gene to its target.

### PrimeKG connectivity normalization

Preferred raw inputs from PrimeKG (from `results/kgpred_<indication>/primekg/...`):
- `primekg_connectivity_score_raw`: expected in `[0, 1]` where **higher = more connected** (shorter/stronger paths)
- optionally `primekg_path_length_raw`: integer path length (1 = direct)

Transform to 0–100 connectivity for scoring:

- If `primekg_connectivity_score_raw` is present:
  - `Connectivity = clamp(primekg_connectivity_score_raw * 100, 0, 100)`
- Else if only `primekg_path_length_raw` is present:
  - `score = 0.9^(path_length - 1)` (with `path_length >= 1`)
  - `Connectivity = clamp(score * 100, 0, 100)`
- Else:
  - `Connectivity = 50.0` (neutral)

This replaces older “connection-count / max-in-run” normalizations.

### GSP (Genetics Support Profiler)
- Human genetics evidence from GWAS, Mendelian, gene burden
- Scoring: High=100, Medium=50, Low=25, None/Unclear=0
- Combo: mean of member-gene GSP scores (averaged across all queried phenotypes)

### UKBPPP (UKB Pharma Proteomics Project)
- Source: `test/genetics_ukbppp/ukbppp_summary.md` (pre-computed Total Score)
- Scoring: FDR-significant model count. Net direction: +1→30, +2→60, +3→100. Average across K50/K51.
- Range: 0-100. Genes not on panel = 0.
- Combo: mean of member-gene scores
- Also used as Translation score (identical value)

## Safety (0–100)

Compute a severity-weighted average from OFF‑X breakdown counts:

- Base scores: `very_high=0`, `high=10`, `medium=20`, `low=60`, `very_low=80`, `not_assoc=100`
- NA rows are excluded from the denominator (unknown).

`Safety = Σ(count[class]·score[class]) / Σ(count[class])`

Combo rule:
- Prefer combo-level OFF‑X breakdown if you have it; otherwise **mean of member-gene safety scores** (fallback).
- Data source: `offx_safety_data.json` (`gene_stats` for genes, `combo_results` for combos)

## Druggability (0–100)

- Input: DrugnomeAI `composite_score` (0–1 scale)
- Normalize: `Druggability = composite_score × 100`
- Combo: prefer combo-level score if available; otherwise **mean of member-gene druggability scores** (fallback)
- Source: DrugnomeAI PU learning model (`drugnome_individual_score.tsv`, `drugnome_combo_score.tsv`)

## Translation (0–100)

**Composite formula:** `Translation = 0.5 × UKBPPP + 0.5 × Combo_Sig_Normalized`

- **UKBPPP (50%):** UKB Pharma Proteomics Project score (same source as Biology UKBPPP component, 0-100). Genes not on panel = 0.
- **Combo Signature (50%):** Agonist-corrected combo signature score, normalized by dynamic max across all combos (0-100).
  - Agonist genes (e.g., GLP2R): signature score negated (`*-1`) before combo mean calculation
  - Dynamic max: `max(abs(all_combo_sig_scores))` computed across all combos
  - Normalized: `(combo_sig / max_sig) * 100`
- **Combo:** `0.5 × mean-member-UKBPPP + 0.5 × combo_sig_normalized`
- **Gene:** `0.5 × gene_UKBPPP + 0.5 × best_normalized_sig_across_combos_containing_gene`

## Commercial (0–100)

**Composite formula:** `Commercial = 0.40 × Market_Opportunity + 0.40 × Competitive_Profile + 0.20 × Strategic_Fit`

Scored via deep-research per combo (web research + local CI data context). See `references/commercial-scoring-contract.md` for full schema.

- **Input:** CI HTML dashboard (`results/ci_<indication>/*.html`) + web research
- **Market Opportunity (40%):** TAM, unmet need, patient population, growth. Dynamic per-indication calibration.
- **Competitive Profile (40%):** Differentiation potential, pipeline competitors, patent landscape.
- **Strategic Fit (20%):** TA alignment (modality-agnostic). GI/immunology = core.
- **Output:** `results/commercial_<indication>/commercial_scores.json` (per-combo scores + confidence + citations)
- **Combo:** Deep-research evaluates the combo AS A PRODUCT
- **Gene:** mean of all combo scores containing that gene
- **Display:** Score bars + bullet points + narrative + confidence badge + citations

## Combo Fallback Rule (General)

When combo-level data is missing (loader returns `raw_value: None`), ALL scoring dimensions fall back to **mean of member-gene scores**. This applies to: Safety, Druggability, BioBridge, ULTRA, PrimeKG, GSP.
