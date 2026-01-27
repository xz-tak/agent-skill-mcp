# Portfolio Summary Report (Multi‑Indication Stitcher)

This document describes the **final compilation step**: taking multiple already-generated, per‑indication offline reports (the “basic reports”) and **knitting them into one single, fully offline portfolio report**.

This is intentionally **post‑pipeline**: the underlying GI2 pipeline (loaders → scoring → per‑indication `*_Combo_Prioritization_Report_offline.html`) is assumed to already exist and produce valid standalone HTML reports.

## What this produces

- Output: one file, e.g. `GI2_Combo_Portfolio_Report_offline.html`
- Characteristics:
  - **single HTML**
  - **no external reads** (each indication report is embedded inside the portfolio HTML)
  - **same theme/style/layout contract** as the GI2 reports
  - Summary tab emphasizes **portfolio overview** (charts + overlap matrices), not deep per‑indication details

## Inputs (basic reports)

The stitcher expects inputs like:

- `results_*/reports/*_Combo_Prioritization_Report_offline.html`

Each file should be a GI2-style offline report generated from the template in `assets/combo_prioritization_report_template.html` (or a compatible derivative).

## Reference examples (visual + behavior contracts)

- **Per‑indication contract**: `assets/example_report_offline.html`
  - Use this to confirm typography, score coloring, nav behavior, etc.
- **Portfolio contract (this stitcher’s target)**: `<WORKDIR>/GI2_Combo_Portfolio_Report_final_v4.html`
  - Use this to confirm Summary layout + tabbed stitching behavior.

## Required behaviors

### 1) Embed reports (no external reads)

- The portfolio HTML must contain the full HTML content for each input report (commonly base64 in a `<script>` block).
- Each report is loaded lazily when its tab is selected.
- **Must work for very large HTML payloads** (Plotly-heavy):
  - Avoid giant `data:` URLs (can silently fail).
  - Recommended approach: decode bytes and stream into an iframe via:
    - `iframe.contentDocument.open(); iframe.contentDocument.write(...); iframe.contentDocument.close();`
- **Page mode (required):** the embedded report iframe should auto-size its height to the embedded report’s content height so the **outer portfolio page** scrolls naturally (no constrained iframe viewport).
- **UI constraint:** do **not** add a “Scroll mode” / “viewer mode” toggle and do **not** add an “Open full report” button/link in indication tabs (portfolio should remain single-page, page-mode only).

### 1b) Theme consistency (portfolio chrome)

- Keep the portfolio’s header, tabs, cards, and typography consistent with the GI2 per‑indication report theme (same CSS variables / palette: **red header**, red accent, paper cards).
- Embedded reports retain their own internal styling (they are rendered inside an iframe), so the portfolio should focus on consistent *outer* chrome rather than trying to restyle iframe contents.

### 2) Global “Built by” rewrite

Rewrite every occurrence of `Built by ...` to:

- `Built by Xinghao Zhang` (default; configurable)

Apply to:
- portfolio header
- indication tab headers
- each embedded report’s HTML

Also enforce consistent titles in the stitched report:
- `<Indication> Combination Prioritization Report`

### 3) Executive Summary rewrite (embedded reports)

In each embedded report, rewrite the `#executive` section so it includes **ONLY**:

1. Scoring weights
2. Top combinations (by Overall score)

If the report doesn’t contain those blocks, synthesize them from the rankings table + canonical defaults.

### 4) Heatmap QA (embedded reports)

GI2-style per‑indication reports apply per‑column heatmap fills to score tables.

- **Heatmap scale must be fixed to 0–100**:
  - `0 = red`, `50 = yellow`, `100 = green`
  - smooth gradient from red → yellow (0–50) and yellow → green (50–100)
- If a numeric score-table column is constant (`max == min`), remove the heatmap fill for that column (neutral / no background).
  - Rationale: constant columns contain no rank information and should not display a misleading “high/low” gradient.

## Required Summary tab layout (portfolio)

The Summary tab must follow the reference portfolio layout:

### Executive Summary (full width)

- KPI tiles
- Indication overview table
- Overall-score ranked-list bar chart with indication dropdown (single-color red bars)
- Two matrices:
  - Gene × Indication (✓/✕)
  - Combo gene coverage × Indication (✓ if all genes present; else ✕)

### Lower row (2 columns)

- Left: Shared / Unique
- Right: stacked:
  - Indication Distribution (interactive pie hover)
  - Overall Score Strategy

## Script usage (bundled)

Use the bundled stitcher (either run it directly from the skill path, or copy it into your project as `scripts/build_combo_portfolio_report.py`):

```bash
python /home/sagemaker-user/.codex/skills/target-prioritization-report/scripts/build_combo_portfolio_report.py \
  --glob "results_*/reports/*_Combo_Prioritization_Report_offline.html" \
  --out "GI2_Combo_Portfolio_Report_offline.html" \
  --built-by "Xinghao Zhang"
```

Optional controls:
- tab order + label overrides (see `--help` in the script)
