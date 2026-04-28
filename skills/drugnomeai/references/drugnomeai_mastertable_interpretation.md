# DrugnomeAI Master Table Columns (44 total)

## Column 1: Identifier

| # | Column | Interpretation |
|---|--------|----------------|
| 1 | Gene_Name | HGNC gene symbol |

## Columns 2-34: Per-Run Triplets (11 runs x 3 metrics each)

Each run produces 3 columns. `_proba` = how druggable (0-1), `_perc` = rank among all 20K genes (0-100), `_known` = was this gene in the training set (1) or a novel prediction (0).

| Columns | Run | What it measures |
|---------|-----|------------------|
| 2-4 | tclin | Resemblance to approved drug targets (Tclin = has approved drug with known mechanism) |
| 5-7 | tchem | Resemblance to chemical tool targets (Tchem = has potent compound ≤30nM but no approved drug) |
| 8-10 | tclin_tchem | Combined Tclin + Tchem — resemblance to any pharmacologically active target |
| 11-13 | tier1 | Resemblance to Tier 1 targets — has clinical precedence (approved drug or clinical candidate) |
| 14-16 | tier12 | Resemblance to Tier 1+2 — clinical + discovery precedence |
| 17-19 | tier123A | Resemblance to Tier 1+2+3A — adds predicted tractable by structure (binding pocket) |
| 20-22 | tier123AB | Resemblance to Tier 1+2+3A+3B — broadest tractability (includes expression/pathway evidence) |
| 23-25 | tclin_tier1 | Intersection of Tclin and Tier 1 — highest-confidence validated targets |
| 26-28 | small_mol | Likelihood of being druggable by small molecule (intracellular: kinase, GPCR, enzyme) |
| 29-31 | antibody | Likelihood of being druggable by antibody (extracellular: receptor, cytokine, surface antigen) |
| 32-34 | protac | Likelihood of being druggable by PROTAC/degrader (undruggable by conventional means) |

## Columns 35-44: Derived/Composite

| # | Column | Interpretation |
|---|--------|----------------|
| 35 | pharos_mean_proba | Average of tclin, tchem, tclin_tchem probas — overall drug-target maturity |
| 36 | tier_mean_proba | Average of tier1, tier12, tier123A, tier123AB probas — overall structural tractability |
| 37 | best_modality | Which modality scored highest: Small Molecule, Antibody, or PROTAC |
| 38 | best_modality_proba | Probability of the best modality (0-1) |
| 39 | modality_specificity | Gap between best and 2nd-best modality. >0.15 = strong single-modality call, <0.05 = multi-modal |
| 40 | composite_score | Primary druggability score: 0.3 x pharos + 0.3 x tier + 0.4 x best_modality. >=0.8 strong, 0.6-0.8 moderate, <0.3 weak |
| 41 | n_runs_above_75perc | How many of the 11 runs ranked this gene in the top 25% (0-11) |
| 42 | consensus_tier | High (>=8 runs in top 10%), Moderate (>=5 runs in top 25%), Low |
| 43 | n_runs_top_decile | How many of the 11 runs ranked this gene in the top 10% (0-11) |
| 44 | is_novel_everywhere | 1 = gene was not in any training set (pure prediction), 0 = known in at least one run |

## Quick example — A1BG (first data row)

- `tclin_proba=0.715` — moderately resembles approved drug targets
- `tier123AB_known=1.0` — was in the Tier 1-3B training set
- `antibody_proba=0.760` — best modality is antibody
- `composite_score=0.702` — moderate druggability
- `consensus_tier=High` — 9/11 runs in top 10%, 10/11 in top 25%
- `is_novel_everywhere=0` — known target in at least one dimension
