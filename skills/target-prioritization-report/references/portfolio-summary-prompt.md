# Portfolio Summary Report (Final Stitch Step) — Future Project Prompt

Copy/paste the prompt below into a new Codex/ChatGPT session when you want to build a **single standalone portfolio report** from multiple already-generated per‑indication offline reports.

This is explicitly the **final compilation step**. The “basic reports” (per indication `*_Combo_Prioritization_Report_offline.html`) are assumed to already exist.

---

Use `target-prioritization-report`.

## Inputs (fill in)

- Workdir:
  - `<WORKDIR>`
  - Example: `/path/to/project`
- Input reports glob:
  - `<GLOB>`
  - Example: `results_*/reports/*_Combo_Prioritization_Report_offline.html`
- Output path:
  - `<OUT_HTML>`
  - Example: `GI2_Combo_Portfolio_Report_offline.html`
- Built-by (default):
  - `Xinghao Zhang`

## Reference examples (must follow)

- Treat `<WORKDIR>/GI2_Combo_Portfolio_Report_final_v4.html` as the **portfolio layout + interaction contract**:
  - top tabs
  - Summary layout (Shared/Unique left; Distribution+Strategy right)
  - hover tooltip behavior on pie
  - dropdown-driven red bar chart (overall scores)
  - lazy-loading embedded reports (must not be blank)
- Treat `assets/example_report_offline.html` as the **per‑indication styling contract**.

## Non-negotiable requirements

- Output is **ONE** fully offline HTML file.
- Each indication report tab must render an embedded copy of its source report:
  - do **not** read any external HTML at runtime
  - do **not** rely on CDN resources
  - must support very large reports (avoid `data:` URL size issues; stream into iframe via `document.write`)
  - **page mode only**: auto-size the iframe height to the embedded report’s content height so the portfolio page scrolls naturally
  - do **not** include a scroll/viewer-mode toggle or an “Open full report” action in indication tabs

## Required global rewrites

1) Built-by rewrite:
- Replace every “Built by …” in:
  - portfolio header + tabs
  - each embedded report HTML
to:
- `Built by Xinghao Zhang` (unless overridden)

2) Title normalization:
- Enforce `<Indication> Combination Prioritization Report` for all stitched panel titles and embedded report titles.

2) Embedded Executive Summary rewrite (each indication report):
- In each embedded report, the Executive Summary must contain **only**:
  - `1. Scoring weights`
  - `2. Top combinations (by Overall score)`
- Remove any other content inside the Executive Summary.

3) Embedded heatmap scale + QA:
- Heatmap scale is fixed **0–100** with **0=red, 50=yellow, 100=green** (smooth gradient).
- If a score-table column is constant (`max == min`), remove its heatmap fill (neutral / no background).

## Required Summary layout

Summary page structure:
- `Executive Summary` (full width)
- Lower row:
  - Left: `Shared / Unique`
  - Right: top `Indication Distribution`, bottom `Overall Score Strategy`

Inside `Executive Summary`:
- KPI tiles: Indications, Total Genes, Total Lists, Shared ≥2
- Indication overview table (left)
- Red bar chart “Overall scores (ranked lists)” (right) with indication dropdown
- Two matrices:
  - Gene × Indication (✓/✕)
  - Combo gene coverage × Indication (✓ if all genes present; else ✕)

## Implementation expectation

Prefer using the bundled stitcher script instead of writing a new one:

```bash
cd <WORKDIR>
python /home/sagemaker-user/.codex/skills/target-prioritization-report/scripts/build_combo_portfolio_report.py \
  --glob "<GLOB>" \
  --out "<OUT_HTML>" \
  --built-by "Xinghao Zhang"
```

If anything in the generated HTML differs materially from `<WORKDIR>/GI2_Combo_Portfolio_Report_final_v4.html` (layout/theme/behavior), fix the generator script rather than hand-editing the output.
