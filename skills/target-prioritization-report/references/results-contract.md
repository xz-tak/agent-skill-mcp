# Results Tree Contract (Abstract)

This reference describes the *expected* artifact layout under an abstract working directory `<WORKDIR>`.

Principle: keep the **top-level `<WORKDIR>` configurable**, but treat **subdirectory names and key filenames as stable contracts**.

## Artifact root (preferred)

- `<WORKDIR>/results/`: clinical / DEG / safety / CI and optional analysis tabs **including** `kgpred_<indication>/` for KG model exports.

`<indication>` is typically a short, lowercased token used in folder suffixes (e.g., `hs`, `ibd`).

Legacy fallback (supported in some pipelines):
- `<WORKDIR>/kgpred_<indication>/` at the workdir root (instead of nested under `results/`).

## Minimum Directory Structure

- `<WORKDIR>/results/`
- `<WORKDIR>/reports/`

Recommended full tree (some sources are optional):

```
<WORKDIR>/
  results/
    cortellis_<indication>/     # markdown (clinical validation)
    offx_<indication>/          # markdown and/or json (optional; safety)
    deg_results_<indication>/   # markdown (or CSV sidecar)
    ci_<indication>/            # HTML dashboard with <script id="data"> JSON
    pathwaydb_<indication>/     # markdown + png (optional)
    interactdb_<indication>/    # markdown (optional)
    bulk_coexpression_<indication>/  # Plotly HTML exports + md narrative (optional)
    sc_coexp_<indication>/      # Plotly HTML exports + md narrative (optional)
    kgpred_<indication>/
      biobridge/                # markdown (Part1 individual, Part2 combo)
      ultra/                    # markdown (Part1 individual, Part2 combo)
      primekg/                  # markdown (Part1 individual, Part2 combo)
      data/                     # optional structured sidecars (json/parquet)
  reports/
```

## KG prediction reports: Step1 vs Step2

For scoring, treat these as the authoritative analysis pages:

- **Step1 (individual gene scoring)**
  - `results/kgpred_<indication>/biobridge/Part1_Individual_Analysis.md`
  - `results/kgpred_<indication>/ultra/Part1_Individual_Analysis.md`
  - `results/kgpred_<indication>/primekg/Part1_Individual_Analysis.md`

- **Step2 (combination scoring)**
  - `results/kgpred_<indication>/biobridge/Part2_Combo_Analysis.md`
  - `results/kgpred_<indication>/ultra/Part2_Combo_Analysis.md`
  - `results/kgpred_<indication>/primekg/Part2_Combo_Analysis.md`

Loaders should parse these and extract per-gene / per-combo raw values.

## Required vs Optional Inputs

**Required to compute the core score table (Clinical/Disease/Safety/Opportunity/Novelty):**
- Cortellis clinical report (`results/cortellis_<indication>/*.md`)
- DEG report or sidecar (`results/deg_results_<indication>/*.md` or CSV sidecar)
- CI dashboard export (`results/ci_<indication>/*.html` with `<script id="data">` JSON)
- KG prediction reports (Step1 and Step2 as above)

**Optional (report should still build, but sections may be empty or neutral):**
- OFF‑X safety (`results/offx_<indication>/*.md` and/or JSON sidecar). If missing, Safety should degrade to neutral.
- Pathways, PPI, bulk/single-cell coexpression tabs.

## Data sources (typical)

| Source | Typical path(s) under `<WORKDIR>/` | Used for |
|---|---|---|
| Cortellis clinical | `results/cortellis_<indication>/*Target_Analysis_Report.md` | Clinical Validation raw score |
| Cortellis gene XLSX (summary visuals) | `results/cortellis_<indication>/gene_cortellis_data.xlsx` or `results/cortellis_<indication>/<GENE>_cortellis_data.xlsx` | Gene card Summary stacked indication×phase bars |
| OFF‑X safety | `results/offx_<indication>/*.md` (and/or JSON sidecar) | Safety severity breakdown |
| DEG | `results/deg_results_<indication>/*Targets_Summary_Report.md` (and/or CSV sidecar) | DEG evidence score |
| BioBridge (kgpred) | `results/kgpred_<indication>/biobridge/Part1_Individual_Analysis.md` + `Part2_Combo_Analysis.md` | Disease evidence percentile |
| ULTRA (kgpred) | `results/kgpred_<indication>/ultra/Part1_Individual_Analysis.md` + `Part2_Combo_Analysis.md` | Model percentile |
| PrimeKG (kgpred) | `results/kgpred_<indication>/primekg/Part1_Individual_Analysis.md` + `Part2_Combo_Analysis.md` | Connectivity score / path length + narrative excerpt |
| Pathways | `results/pathwaydb_<indication>/REPORT_*.md` + `results/pathwaydb_<indication>/*.png` | Pathway overlap + interpretation |
| PPI | `results/interactdb_<indication>/**/*.md` | Cross-db shortest paths narrative |
| Bulk coexpression | `results/bulk_coexpression_<indication>/results/*.html` + `results/bulk_coexpression_<indication>/*.md` | Plotly exports + interpretation |
| Single-cell coexp | `results/sc_coexp_<indication>/**.html` + `results/sc_coexp_<indication>/**.md` | Plotly exports + interpretation |
| CI dashboard | `results/ci_<indication>/*.html` with `<script id="data">` JSON | Competitive Intelligence |

## Data Name → Path → Type (Suggested Contract)

Use this when implementing loaders; treat “Data name” as the canonical `data_dict` key.

| Data name | Path pattern (under `<WORKDIR>/`) | Type | Minimum extracted fields |
|---|---|---|---|
| `cortellis` | `results/cortellis_<indication>/*.md` | Markdown | `total_score_raw`, `total_drugs`, `total_trials` |
| `cortellis_xlsx` | `results/cortellis_<indication>/*cortellis_data.xlsx` or `gene_cortellis_data.xlsx` | XLSX | `highest_phase_counts_raw`, `top_indications_phase_counts_raw` |
| `deg` | `results/deg_results_<indication>/*.md` (or CSV sidecar) | Markdown/CSV | `deg_score_raw` (+ optional study counts) |
| `biobridge` | `results/kgpred_<indication>/biobridge/Part1_Individual_Analysis.md` and `Part2_Combo_Analysis.md` | Markdown | `biobridge_percentile_raw` (per gene or per combo) |
| `ultra` | `results/kgpred_<indication>/ultra/Part1_Individual_Analysis.md` and `Part2_Combo_Analysis.md` | Markdown | `ultra_percentile_raw` (per gene or per combo) |
| `primekg` | `results/kgpred_<indication>/primekg/Part1_Individual_Analysis.md` and `Part2_Combo_Analysis.md` | Markdown | `primekg_connectivity_score_raw` and/or `primekg_path_length_raw` (+ excerpt) |
| `ci` | `results/ci_<indication>/*.html` (`<script id="data">`) | HTML/JSON | `programs_by_phase`, `total_programs` |
| `offx` | `results/offx_<indication>/*.md` (or JSON sidecar) | Markdown/JSON | `safety_breakdown_raw` |
| `pathway` (combo) | `results/pathwaydb_<indication>/REPORT_*.md` + `results/pathwaydb_<indication>/*.png` | MD/PNG | `shared_pathways`, `distinguished_pathways` |
| `ppi` (combo) | `results/interactdb_<indication>/**.md` | Markdown | `min_hops`, `databases_found` |
| `coexpression` (combo) | `results/bulk_coexpression_<indication>/**/*.html` + `results/sc_coexp_<indication>/**/*.html` | Plotly HTML | keys to embed + summary metrics |

Legacy compatibility note:
- Some older pipelines store BioBridge/PrimeKG under `results/biobridge/` and `results/primekg/` with keys like `primekg_connections_raw`. If you must support those, keep the raw values under legacy keys **and** map them into the new scoring fields when computing Novelty.

## If your workspace does not match the contract

If your working directory does not contain a `results/` tree with the expected subfolders, you should:
- search for likely artifact roots and filenames (e.g., `cortellis_*`, `deg_*`, `ci_*`, `kgpred_*`)
- confirm with the user which artifacts are authoritative
- record the resolved mapping (source → path(s) → fields/figures) in `<WORKDIR>/PLAN.md`

## Naming patterns (recommendation)

For combos, keep a stable **combo slug** used in filenames:

- Display name: `GENE1-GENE2[-GENE3]`
- File slug: `GENE1_GENE2[_GENE3]` (underscores)
- If you use list numbering: prefix with `listN_` consistently across related artifacts.

The report generator can be written to either:
- rely on naming conventions (fast) or
- accept an explicit manifest/config mapping combo → artifact paths (robust).
