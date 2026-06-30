# Plan + Design (Abstract Workflow)

Use this reference when asked to “map the pipeline end-to-end”, “make it reproducible”, “update loaders/scoring”, or “generate a traceable report”.

## 0) Always track the work in `<WORKDIR>/PLAN.md`

Before coding/parsing/scoring, create or update `<WORKDIR>/PLAN.md` with:
- agreed report scope (summary + tabs + figures)
- agreed scoring strategy (weights + normalizers + fallbacks)
- concrete steps with checkboxes and current status
- **resolved artifact paths** (exact files in this workspace for each score component/plot/tab)

Proceed by the plan and keep it updated as scope changes.

### Path resolution rule (important)

The contracts in `references/results-contract.md` describe **expected** locations. In practice, working directories often contain artifacts with different names (e.g., `results/cortellis_hs/` instead of `results/cortellis_<indication>/`).

If an expected file/path is missing:
1) **Search the working directory** for likely candidates (folders and filenames).
2) **Stop and clarify** which artifact is authoritative.
3) **Record the final mapping in `PLAN.md`** (source → resolved path(s) → fields extracted → figures embedded).

## 1) Requirements → Design Doc (1 file, concrete)

Create a design document in `<WORKDIR>/docs/plans/` (or `<WORKDIR>/plans_xz/`) that includes:

- **Data source mapping table**: source → file path(s) → fields extracted → storage key in `data_dict`
- **Per-subscore workflow**: Input → Parse → Store raw → Normalize → Weight → Output
- **Combination aggregation rules**: when to average member genes vs when to use a combo artifact directly
- **Fallback rules**: missing source → neutral default (typically `50.0`) + log warning
- **Report compilation steps**: how the HTML is assembled and where each plot/markdown comes from

Keep the doc **traceable**: every displayed score should be explainable with (a) source file and (b) transform.

## 2) Module boundaries (repeatable implementation order)

Implement bottom-up:

1. `target_prioritization/data_loaders/parsers.py`
   - Pure parsing utilities (markdown tables/sections, HTML `<script id="data">` JSON blobs, simple regex).
2. `target_prioritization/data_loaders/loaders.py`
   - One loader per source. Returns `{source_file, <field>_raw, ...}` only.
3. `target_prioritization/data_loaders/orchestrator.py`
   - `load_all_data(genes, results_dir, kgpred_dir, is_combo)` merges loaders into one `data_dict`.
4. `target_prioritization/scoring/normalizers.py`
   - Raw → 0–100 mappings only. No weights.
5. `target_prioritization/scoring/subscores.py`
   - Compute the 5 subscores; return `{score, components, raw}` as needed.
6. `target_prioritization/scoring/scoring.py`
   - Overall weighting, synergy metrics for combos, any categorization labels.
7. Report generator (e.g., `<WORKDIR>/tools/generate_*_report.py`)
   - Pulls scores + artifacts + narrative into one offline HTML file.

## 3) Data traceability contract

For each loader:
- always store `source_file` (workspace-relative path)
- store raw numeric values under `*_raw`
- do *not* overwrite raw values with normalized ones

For each subscore:
- include a `components` section that records normalized components used in scoring (so the report can render bar plots / tables).

## 4) Report build contract (format is stable; content varies)

The report should:
- Render the **Priority Rankings** table with column-wise green→red heat coloring.
- Provide **tabs per combination** (Overview / Disease / Bulk / Single-cell / Pathways / PPI, as available).
- Embed Plotly fully offline by extracting the Plotly bundle from one exported Plotly HTML and reconstructing each plot in-page.
- Pull interpretation text from markdown reports under `<WORKDIR>/results/` and render it as HTML.
- For pathway narratives, summarize “shared vs distinguished” *as biological functions/themes* (avoid dumping raw pathway name lists).
