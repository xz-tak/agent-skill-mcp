---
name: target-prioritization-report
description: Create and maintain a traceable target-prioritization pipeline (loaders + scoring) and generate a single-file, fully offline interactive HTML report with the same style/layout as the GI2 combo prioritization report. Use when asked to (1) map redesign requirements into a concrete plan and implementation, (2) wire data loaders for `results/*` artifacts, (3) implement/adjust scoring (Clinical/Disease/Safety/Opportunity/Novelty; CI is folded into Opportunity), or (4) generate a standalone report that embeds Plotly HTML figures and markdown interpretation for any indication/target set.
---

# Target Prioritization Report

## Overview
Produce a **fully offline, single-HTML** prioritization dossier (tables + tabs + embedded Plotly + narrative excerpts) and the supporting **loader/scoring** pipeline. Keep raw values for traceability; normalize only for scoring/display; degrade gracefully when artifacts are missing.

This skill assumes a working directory `<WORKDIR>` with a `results/` tree containing precomputed artifacts (Markdown/HTML/JSON/PNG) for Cortellis/DEG/OFF-X/CI and optional analysis tabs, including `results/kgpred_<indication>/` for KG model outputs (BioBridge / ULTRA / PrimeKG) used for Disease + Literature Novelty.

**Default scope (unless the user overrides):**
- Build a report that includes **both**:
  - **Combinations (gene lists):** executive summary + rankings + one card per combo
  - **Single genes:** a single-gene score table **and** per-gene cards
- **Section order (default):** combos first, then single genes.
- **Combo display names (default):** use `GENE1-GENE2(-GENE3)`; do **not** prefix with labels like `list1:`, `list2:`, etc.
- Follow the **same order/sections** as the bundled example HTML (only the data/figures/narratives change).

**Default scoring (unless the user overrides):**
- Overall weights (CI is not a standalone overall component; it is part of Opportunity):
  - `Overall = 0.30·Clinical + 0.30·Disease + 0.10·Safety + 0.20·Opportunity + 0.10·Novelty`
- Disease subweights:
  - `Disease = 0.40·DEG + 0.25·BioBridge + 0.25·ULTRA + 0.10·PrimeKG`

## Required Workflow (Must Follow)
1) **Plan first (no questions yet):** propose a concrete plan that aligns on:
   - **layout/theme/style contract** (must match the bundled example report)
   - report summary structure (exec summary + tabs + what narrative excerpts to show)
   - figures to embed (which Plotly exports / PNGs; what gets omitted if missing)
   - scoring strategy (weights, normalizers, missing-data defaults, combo aggregation rules)
2) **Clarify after planning:** ask any questions needed to finalize the plan (paths, target lists, weights, which artifacts are authoritative).
3) **Write `<WORKDIR>/PLAN.md`:** once aligned, write/update `PLAN.md` in the working directory as the task tracker and proceed step-by-step based on it.

## Required Directory Structure (To Start Report Generation)
Minimum structure (top-level paths are abstract; subdirs are stable):

```
<WORKDIR>/
  results/                      # precomputed artifacts (inputs)
    cortellis_<indication>/     # markdown
    offx_<indication>/          # markdown (or JSON sidecar if your pipeline produces one)
    deg_results_<indication>/   # markdown (or CSV sidecar)
    ci_<indication>/            # HTML with <script id="data"> JSON
    bulk_coexpression_<indication>/  # Plotly-export HTML (optional but recommended)
    sc_coexp_<indication>/      # Plotly-export HTML (optional but recommended)
    pathwaydb_<indication>/     # PNG + markdown (optional but recommended)
      COMPREHENSIVE_PATHWAY_ANALYSIS_REPORT.md  # includes individual + list shared/distinguished summaries
    interactdb_<indication>/    # markdown (optional)
    kgpred_<indication>/        # KG prediction exports (inputs)
      biobridge/                # markdown (Part1 individual, Part2 combo)
      ultra/                    # markdown (Part1 individual, Part2 combo)
      primekg/                  # markdown (Part1 individual, Part2 combo)
  target_prioritization/        # loaders + scoring (or equivalent package)
  tools/                        # report generator script(s)
  reports/                      # generated output HTML
  PLAN.md                       # task plan + status (created/maintained by this skill)
```

Required inputs to generate a score table (minimum viable report):
- Clinical: `results/cortellis_<indication>/*.md` (or a JSON sidecar if that’s what your pipeline emits)
- Disease: `results/deg_results_<indication>/*.md` (or CSV sidecar) + `results/kgpred_<indication>/{biobridge,ultra,primekg}/Part1_Individual_Analysis.md`
- CI: `results/ci_<indication>/*.html` with `<script id="data">` JSON
- Safety: `results/offx_<indication>/*.md` (if absent, Safety will degrade to neutral)

**Important:** BioBridge / ULTRA / PrimeKG KG outputs are expected under `results/kgpred_<indication>/...` (Part1 for single genes, Part2 for combos). Treat `<WORKDIR>/kgpred_<indication>/*` or `results/biobridge/*` as **legacy** layouts only.

For combination scoring, Disease KG components should come from:
- `results/kgpred_<indication>/{biobridge,ultra,primekg}/Part2_Combo_Analysis.md`

Everything else (bulk/single-cell plots, pathways, PPI) is optional; the report should still build with those tabs empty or replaced by a “missing artifact” note.

### Flexible artifact discovery (default)

If an expected `results/<source>_<indication>/...` path is missing, **do not hard-fail**. Instead:
1) **Search recursively under `<WORKDIR>/<corresponding_section>_<indication>`** for the relevant `.md` / `.html` / `.png` artifacts.
2) **Confirm** the authoritative file(s) if multiple candidates exist.
3) **Record** the resolved mapping in `<WORKDIR>/PLAN.md` (source → resolved file path(s) → fields/figures used).

Optional (recommended) narrative inputs:
- PPI (InteractDB): `results/interactdb_<indication>/interactdb_results/SUMMARY_REPORT.md`
- Pathways (PathwayDB): `results/pathwaydb_<indication>/COMPREHENSIVE_PATHWAY_ANALYSIS_REPORT.md` (individual + list shared/distinguished/unique summary)

### Single-gene Summary (Cortellis XLSX, default)

For each gene card’s **Summary** section (`gene.overview_html` in the template), prefer a **Cortellis XLSX sidecar** when available:

- File: `results/cortellis_<indication>/gene_cortellis_data.xlsx`
- Alternate (also valid): `results/cortellis_<indication>/<GENE>_cortellis_data.xlsx`
- Sheet: `Drugs_Comprehensive`
- Columns:
  - `Highest Phase`
  - `Primary Indications` (may be multi-valued per row; split on `;` or `|` by default)

**REQUIRED presentation in the report (Summary section):**
- A **stacked bar plot** of **top indications** (x/y orientation optional), where:
  - each bar = an indication (only high-count indications; avoid long tail)
  - stacked segments = **Highest Phase** categories (consistent colors across genes)
  - **phase inclusion rule (REQUIRED):**
    - include **clinical and later** phases only (e.g., `Phase 1 Clinical`, `Phase 2 Clinical`, `Phase 3 Clinical`, `Phase I/II`, `Pre-registration`, `Registered`, `Launched/Marketed`)
    - **fallback:** if a gene has **no clinical-or-later** programs, then include `Preclinical`, `Discovery` (if present), and `Discontinued` as the next-best signal
  - **discard NaN/blank** indication labels and phase labels (do not render any `"nan"` bars/legend entries)
  - **phase label requirement:** include Cortellis phase labels like `Phase 1 Clinical`, `Phase 2 Clinical`, `Phase 3 Clinical` (and `Phase I/II` when present) so they reliably appear in the stacked bars/legend with consistent colors
- (Optional) a separate bar plot for overall Highest Phase distribution.

Helper script (bundled):
- `python scripts/summarize_cortellis_gene_xlsx.py --xlsx <path> --gene <GENE>` (markdown tables, including an indication×phase matrix)
- `python scripts/summarize_cortellis_gene_xlsx.py --xlsx <path> --gene <GENE> --format json` (structured counts for stacked-bar plotting)

## Core Outputs
- `PLAN.md`: the working plan + progress tracker for the run.
- `docs/plans/DESIGN_<date>_<indication>.md` (or equivalent): data sources + workflows + scoring formulas + traceability.
- `target_prioritization/data_loaders/*`: loaders that map artifacts → `data_dict` with raw fields + `source_file`.
- `target_prioritization/scoring/*`: normalizers + subscores + overall scoring.
- `reports/<IND>_Combo_Prioritization_Report_offline.html`: standalone report, offline Plotly embedded.

## Workflow (Decision Tree)
1) If the user already has precomputed artifacts under `<WORKDIR>/results/` (including `results/kgpred_<indication>/`):
   - Validate artifact existence (see `references/results-contract.md`).
   - Ensure loaders return raw values + `source_file` for traceability.
   - Ensure scoring matches the documented methodology (see `references/scoring-contract.md`).
   - Generate report using the template + offline embedding rules (see `references/report-generation.md`).

2) If artifacts are missing or the pipeline is not implemented:
   - Convert requirements into a design spec + extraction workflows (see `references/plan-and-design.md`).
   - Implement parsers → loaders → orchestrator → normalizers → subscores → scoring (same module boundaries).
   - Add a thin report generator that uses the shared HTML template and offline Plotly embedding.

## Minimal Inputs (Abstract)
You should ask the user for (after planning):
- Indication label (e.g., `IBD`, `UC`, `SSc`) and the resolved artifact suffix used in folder names (default: `<indication_lower>`):
  - `results/cortellis_<indication>/`, `results/deg_results_<indication>/`, `results/ci_<indication>/`, etc
  - `results/kgpred_<indication>/` for BioBridge/ULTRA/PrimeKG
- Target set:
  - single genes (e.g., `["TYK2", "JAK1"]`)
  - combinations (e.g., `["TYK2-JAK1", "TNFRSF25-GREM1"]`)
- Artifact roots:
  - `<WORKDIR>/results/` (default)
  - Report output: `<WORKDIR>/reports/`

## Style/Format Contract
To keep the same look & feel across indications:
- **REQUIRED:** follow the bundled example report’s **layout, theme, and styling**:
  - `assets/example_report_offline.html` (visual contract)
  - `assets/combo_prioritization_report_template.html` (template + JS behavior contract)
- Do **not** redesign the UI/section order; only swap in data/narratives/plots sourced from the configured `results/` directories.
- Use the bundled template in `assets/combo_prioritization_report_template.html` (copy into the project as needed).
- Keep CSS variables, typography scale, nav layout, tab behavior, and score color logic consistent.
- Embed Plotly fully offline by extracting a single Plotly bundle and reconstructing figures from exported Plotly HTML artifacts.

### Pathways / PPI summary requirements (default)

- **Combinations (Pathways tab) — REQUIRED:** summarize the “Shared vs Distinguished Pathways” portion for each combo as **biological functions/themes**:
  - clearly describe what is shared (convergent biology) vs what is distinguished/unique (differentiating biology)
  - include short interpretation of what those themes imply mechanistically
  - **DO NOT** list pathway items directly from the markdown report
- **Single genes (Pathways summary) — REQUIRED:** source from the PathwayDB comprehensive markdown report and summarize the gene’s dominant biological functions/pathway themes **with interpretation**:
  - preferred input: `results/pathwaydb_<indication>/COMPREHENSIVE_PATHWAY_ANALYSIS_REPORT.md`
  - extract the gene’s “top pathways” subsection (e.g., “Top 5 Central Pathways”) from the gene’s Individual Gene Analysis section
  - summarize into a few **biological themes/signatures** (e.g., cytokine signaling, immune activation, barrier remodeling), and add a short mechanistic interpretation
  - **DO NOT** list ranked pathway names/items directly from the markdown report
- **PPI tab:** keep PPI narrative focused on connectors/hops and mechanistic interpretation; do not repeat the pathway list.

## Input Data Contract (Names, Paths, Types)
Open `references/results-contract.md` for the exact *data name → path → type → fields extracted* mapping and which sources are required vs optional.

## Concrete Example Set (GI2 Reference Implementation)
Use these as “known-good” references for reproducing the same layout/style while swapping indication/targets/artifacts.

- **Design spec (example)**:
  - Project path: `<WORKDIR>/plans_xz/DESIGN_2025-12-23_redesign.md`
  - Bundled copy: `references/example_gi2_design_2025-12-23_redesign.md`
  - Use it to see how workflows are mapped end-to-end (sources → transforms → scoring → report compilation).
- **Report generator script (example)**:
  - Project path: `<WORKDIR>/tools/generate_combo_prioritization_report.py`
  - Bundled copy: `references/example_gi2_generate_combo_prioritization_report.py`
  - This is the concrete implementation of: loaders → scoring → narrative extraction → offline Plotly embedding → HTML rendering.
- **Key scoring source files (examples)**:
  - Project paths:
    - `<WORKDIR>/target_prioritization/data_loaders/orchestrator.py`
    - `<WORKDIR>/target_prioritization/scoring/subscores.py`
  - Bundled copies:
    - `references/example_gi2_orchestrator.py`
    - `references/example_gi2_subscores.py`
  - Use these to copy/adapt the loader→score data contract and the CI (gene+family blended) + Opportunity formula.
- **Final offline HTML output (example)**:
  - Project path: `<WORKDIR>/reports/IBD_Combo_Prioritization_Report_offline.html`
  - Bundled copy: `assets/example_report_offline.html`
  - Use this for visual QA to confirm color scale, typography, spacing, tabs, executive-summary formatting (top-3), and plot embedding behavior.

When adapting to a new indication:
- Start by cloning/scaffolding a new report tool (`scripts/clone_report_tool.py`) and update its combo manifest (targets, filenames, markdown section keys).
- Keep the HTML template stable (use `assets/combo_prioritization_report_template.html`).
- Ensure the scoring contract in `references/scoring-contract.md` matches the pipeline you implement (treat the GI2 design doc as an example snapshot).

## Bundled Resources
- `scripts/clone_report_tool.py`: generate a new report-tool scaffold for a new indication/target set (and initialize `PLAN.md`).
- `scripts/validate_results_tree.py`: validate required artifacts exist for a configured run (fast existence check).
- `references/plan-and-design.md`: how to turn requirements into an implementation-ready design + workflow mapping.
- `references/results-contract.md`: artifact naming/paths contract (abstract, configurable).
- `references/scoring-contract.md`: scoring formulas + combo aggregation rules + CI weighting.
- `references/report-generation.md`: offline Plotly embedding + HTML template integration.
- `assets/combo_prioritization_report_template.html`: the reference HTML template (layout/colors/JS).
- `assets/example_report_offline.html`: example output artifact for visual QA (do not treat as a data source).

## Quick Start (Existing Pipeline)
1) Confirm `<WORKDIR>/target_prioritization/` and `<WORKDIR>/tools/` exist.
2) Validate artifacts:
   - `python scripts/validate_results_tree.py --workdir <WORKDIR> --indication <IND>`
   - If validation fails, use the autodetect suggestions printed by the script to **map your actual artifact paths in `<WORKDIR>/PLAN.md`** (source → resolved path(s) → fields/figures), then proceed.
3) Initialize a working plan (created if missing):
   - `python scripts/clone_report_tool.py --workdir <WORKDIR> --indication <IND>`
   - Review/adjust `<WORKDIR>/PLAN.md` (scope + figures + scoring) before editing code.
4) Run the project’s report tool (or the scaffolded tool) to generate `<WORKDIR>/reports/<IND>_Combo_Prioritization_Report_offline.html`.
