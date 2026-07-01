---
name: target-prioritization-report
description: Create and maintain a traceable target-prioritization pipeline (loaders + scoring) and generate single-file, fully offline interactive HTML reports with the GI2 combo prioritization style/layout. Use when asked to (1) map redesign requirements into a concrete plan and implementation, (2) wire data loaders for `results/*` artifacts, (3) implement/adjust scoring (Clinical/Biology/Safety/Druggability/Translation/Commercial), (4) generate a standalone per‑indication report that embeds Plotly HTML figures and markdown interpretation, or (5) stitch multiple per‑indication offline reports into one portfolio summary HTML (final compilation step).
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
- **Combo display names (default):** use `GENE1 + GENE2 (+ GENE3)` (spaces around `+`); do **not** prefix with labels like `list1:`, `list2:`, etc.
- **Author attribution (default):** "Built by Xinghao Zhang" — applied to report headers and portfolio chrome.
- Follow the **same order/sections** as the bundled example HTML (only the data/figures/narratives change).

**Default scoring (unless the user overrides):**
- Overall weights (6 components):
  - `Overall = 0.25·Clinical + 0.25·Biology + 0.125·Safety + 0.125·Druggability + 0.125·Translation + 0.125·Commercial`
- Biology subweights:
  - `Biology = 0.20·DEG + 0.20·Signature + 0.12·BioBridge + 0.12·ULTRA + 0.06·PrimeKG + 0.15·GSP + 0.15·UKBPPP`
- DEG (0.20): Mean target gene total_score from target_score_table.csv, agonist-negated, dynamic-max normalized across all combos.
- Signature (0.20): Mean per-target signature total_score from {sig_name}_score_table.csv, agonist-negated, dynamic-max normalized across all combos.
- Druggability: DrugnomeAI composite_score × 100 (combo = mean of member genes)
- Translation: `0.5 × UKBPPP + 0.5 × Combo_Sig_Normalized` (agonist-corrected, dynamic-max normalized across all combos). Missing UKBPPP gene = 0.
- Commercial: `0.40 × Market_Opportunity + 0.40 × Competitive_Profile + 0.20 × Strategic_Fit` (deep-research scored per combo, 0-100). See `references/commercial-scoring-contract.md`.
- Combo biology scoring:
  - DEG: mean of member-gene normalized DEG scores
  - BioBridge: combo `pct_rank` from Part2 × 100 (fallback: mean of member-gene scores)
  - ULTRA: `geo_pct_rank` from Part2 × 100 (fallback: mean of member-gene scores)
  - PrimeKG: combo score from Part2/JSON × 100 (fallback: mean of member-gene scores)

**DEG normalization (dynamic max, approach support):**
- DEG uses **dynamic max** (computed from actual data, not hardcoded)
- **No clamping** - scores reflect actual distribution
- **Per-gene approach annotation** supported:
  - `"antagonist"` (default): use raw DEG score as-is
  - `"agonist"`: negate raw DEG score (-raw) before normalization
- Format: `{"GENE1": {"approach": "antagonist"}, "GLP2R": {"approach": "agonist"}, ...}`
- See `references/scoring-contract.md` for full DEG normalization specification

**Heatmap color specification (REQUIRED):**
- All scoring heatmaps use: **green (100, high) → yellow (50, neutral) → red (0, low)**
- See `references/scoring-contract.md` and `references/report-generation.md` for color spec

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
- CI: `results/ci_<indication>/*.html` with `<script id="data">` JSON (also used for Commercial scoring context)
- Safety: `results/offx_<indication>/*.md` (if absent, Safety will degrade to neutral)
- Commercial: `results/commercial_<indication>/commercial_scores.json` (if absent, Commercial = 0; run `scripts/run_commercial_scoring.py` + deep-research to generate)
- Translation: UKBPPP data (`test/genetics_ukbppp/`) + signature scores from DEG pipeline (auto-computed)

**Important:** BioBridge / ULTRA / PrimeKG KG outputs are expected under `results/kgpred_<indication>/...` (Part1 for single genes, Part2 for combos). Treat `<WORKDIR>/kgpred_<indication>/*` or `results/biobridge/*` as **legacy** layouts only.

For combination scoring, Disease KG components should come from:
- `results/kgpred_<indication>/{biobridge,ultra,primekg}/Part2_Combo_Analysis.md`

Everything else (bulk plots, pathways, PPI) is optional; the report should still build with those tabs empty or replaced by a “missing artifact” note. Single-cell results are pipeline-output-only (not embedded in the report HTML due to size).

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

3) For Commercial scoring (runs on-the-fly during this skill invocation):

   Follow the execution protocol in `references/commercial-deep-research-execution.md`.

   **Key behavior:**
   - Reads prompts from `<WORKDIR>/test/commercial/commercial_prompts.json`
   - Checks existing per-combo results in `test/commercial/<combo_slug>/commercial_result.json`
   - Checks existing per-gene results in `test/commercial/genes/<GENE>/commercial_result.json`
   - Fans out deep-research Agent() calls for BOTH combos AND individual genes without existing scores
   - Step 3: 3 agents per combo (parallel, independent runs) — evaluates combination as a product
   - Step 3.1: Merge 3 runs → mean scores + union of non-overlapping evidence/citations
   - Step 3b: 3 agents per unique gene (parallel, independent runs) — evaluates individual target
   - Same merge logic: mean scores, union evidence, narrative from highest-scoring run
   - Saves results to `test/commercial/<combo_slug>/` and `test/commercial/genes/<GENE>/`
   - Updates `test/commercial/commercial_scores.json` and `test/commercial/commercial_summary.md`
   - Validates scores via `scripts/run_commercial_scoring.py --mode validate`
   - Then proceeds to report generation with all scores ingested

   **Gene-level commercial loading priority:**
   1. Per-gene result file (`test/commercial/genes/<GENE>/commercial_result.json`) — preferred
   2. Fallback: mean of all combo scores containing that gene

   Results are saved to `test/commercial/` — alongside all other module outputs (`test/omicsoft_full/`, `test/cortellis/`, `test/offx/`, etc.). The orchestrator reads from `test/commercial/` like all other sources.

   **Supporting references:**
   - `references/commercial-deep-research-execution.md` — full step-by-step execution protocol (Steps 3 + 3b)
   - `references/commercial-research-prompt-template.md` — prompt template with variables
   - `references/commercial-scoring-contract.md` — scoring schema, weights, anchors

   When `commercial_scores.json` is absent and no research is run, loaders return `_missing()` and commercial stays 0.0 (graceful degradation).

4) For CI Dashboard embedding:

   If `results/ci_ibd/ibd_dashboard.html` exists, the report generator automatically embeds it as a lazy-loaded iframe in the "Competitive Intelligence" section. The dashboard HTML is stored inline in the report (keeping it fully offline) and loads on user click.

## Minimal Inputs (Abstract)
You should ask the user for (after planning):
- Indication label (e.g., `IBD`, `UC`, `SSc`) and the resolved artifact suffix used in folder names (default: `<indication_lower>`):
  - `results/cortellis_<indication>/`, `results/deg_results_<indication>/`, `results/ci_<indication>/`, etc
  - `results/kgpred_<indication>/` for BioBridge/ULTRA/PrimeKG
- Target set:
  - single genes (e.g., `["TYK2", "JAK1"]`)
  - combinations (e.g., `["TYK2 + JAK1", "TNFRSF25 + GREM1"]`)
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

### Translation section in combo cards (REQUIRED)

Each combo card must include a **Translation** section showing:
- **UKBPPP component:** Per-gene score table + bar chart showing UKBPPP scores for each member gene
  - If detailed UKBPPP TSV data is available (`test/genetics_ukbppp/<GENE>_<ICD10>_ukbppp.tsv`): show model counts, directions, FDR significance
  - Always show the final per-gene UKBPPP score (0-100)
- **Signature component:** Combo-level normalized signature score with agonist correction note
- **Combined Translation score:** `0.5 × UKBPPP + 0.5 × Sig` with visual comparison bar

### Commercial section in combo cards (REQUIRED)

Each combo card must include a **Commercial** section showing:
- **Score bars:** 3 horizontal progress bars for Market Opportunity / Competitive Profile / Strategic Fit (color-coded 0-100, green→yellow→red)
- **Bullet points:** 3-5 key findings from the deep-research assessment
- **Narrative:** Executive summary paragraph from the research
- **Sub-tabs:** Overview | Market | Competition | Strategic Fit (detailed breakdown per dimension)
- **Confidence badge:** High/Medium/Low with color coding (green/yellow/orange)
- **Citations:** Source links/references used in scoring
- **Overall Commercial score** prominently displayed

Data source: `results/commercial_<indication>/commercial_scores.json` (per combo entry matching output schema in `references/commercial-scoring-contract.md`).

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

## Portfolio Summary Report (Final Stitch Step)
Use this when the “basic reports” already exist and you want a **single standalone portfolio report** that knits them together:

- Inputs: multiple per‑indication `*_Combo_Prioritization_Report_offline.html` files (already fully offline)
- Output: `GI2_Combo_Portfolio_Report_offline.html` (single file, fully offline, no external reads)
- Requirements (high level):
  - Embed each per‑indication report inside the portfolio HTML (no runtime file reads)
  - Indication tabs use **page mode only**: auto-size the embedded report iframe height to content so the **portfolio page scrolls naturally** (no “constrained viewer” mode)
  - Rewrite “Built by …” to the configured name (default: **Xinghao Zhang**)
  - Enforce consistent titles: `<Indication> Combination Prioritization Report`
  - Rewrite each embedded report’s Executive Summary to include **only**:
    1) Scoring weights
    2) Top combinations (by Overall score)
  - Summary tab uses the reference layout: Executive Summary + overlap matrices + interactive pie hover tooltip
  - Heatmap scale is **fixed 0–100** with **0=red, 50=yellow, 100=green** and smooth gradients between
  - Heatmap QA: if any score-table column is constant (`max == min`), the heatmap fill for that column must be neutral (no background) to avoid misleading color

Primary docs + artifacts:
- `references/portfolio-summary-report.md` (what this is, required behaviors, required layout)
- `references/portfolio-summary-prompt.md` (copy/paste prompt for future projects)
- `scripts/build_combo_portfolio_report.py` (stitcher script to generate the portfolio report)
- `<WORKDIR>/GI2_Combo_Portfolio_Report_final_v4.html` (visual/layout contract for the portfolio report)

## Bundled Resources
- `scripts/clone_report_tool.py`: generate a new report-tool scaffold for a new indication/target set (and initialize `PLAN.md`).
- `scripts/validate_results_tree.py`: validate required artifacts exist for a configured run (fast existence check).
- `scripts/build_combo_portfolio_report.py`: build the **portfolio** report from multiple per‑indication offline HTML reports (final stitch step).
- `scripts/run_commercial_scoring.py`: extract CI context, generate research prompts, validate commercial scores JSON.
- `references/plan-and-design.md`: how to turn requirements into an implementation-ready design + workflow mapping.
- `references/results-contract.md`: artifact naming/paths contract (abstract, configurable).
- `references/scoring-contract.md`: scoring formulas + combo aggregation rules + CI weighting + Translation + Commercial.
- `references/commercial-scoring-contract.md`: full Commercial scoring schema (weights, research protocol, output JSON schema, differentiation guidance).
- `references/report-generation.md`: offline Plotly embedding + HTML template integration.
- `references/portfolio-summary-report.md`: portfolio stitcher requirements + layout contract.
- `references/portfolio-summary-prompt.md`: future-project prompt (explicit layout + constraints).
- `assets/combo_prioritization_report_template.html`: the reference HTML template (layout/colors/JS).
- `assets/example_report_offline.html`: example output artifact for visual QA (do not treat as a data source).
- `<WORKDIR>/GI2_Combo_Portfolio_Report_final_v4.html`: example portfolio output artifact for visual QA (layout + interactions contract).

## Quick Start (Existing Pipeline)
1) Confirm `<WORKDIR>/target_prioritization/` and `<WORKDIR>/tools/` exist.
2) Validate artifacts:
   - `python scripts/validate_results_tree.py --workdir <WORKDIR> --indication <IND>`
   - If validation fails, use the autodetect suggestions printed by the script to **map your actual artifact paths in `<WORKDIR>/PLAN.md`** (source → resolved path(s) → fields/figures), then proceed.
3) Initialize a working plan (created if missing):
   - `python scripts/clone_report_tool.py --workdir <WORKDIR> --indication <IND>`
   - Review/adjust `<WORKDIR>/PLAN.md` (scope + figures + scoring) before editing code.
4) Run the project’s report tool (or the scaffolded tool) to generate `<WORKDIR>/reports/<IND>_Combo_Prioritization_Report_offline.html`.
