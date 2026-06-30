# Report Generation (Offline Standalone HTML)

Goal: produce **one HTML file** with no external dependencies (no CDNs, no filesystem links required at view time), while preserving interactive Plotly figures.

## 1) Use a stable HTML template

- Start from `assets/combo_prioritization_report_template.html` to preserve the format (colors, spacing, typography, tab logic, score heat coloring).
- Copy the template into the project at `<WORKDIR>/target_prioritization/templates/` (or point the report tool directly to the asset).
  - Example implementation: `references/example_gi2_generate_combo_prioritization_report.py`
- Treat `assets/example_report_offline.html` as the **layout/theme/style contract** for the report.
  - The report generator must keep the same section order and visual design; only the data/narratives/plots change.

## 2) Plotly offline embedding pattern

Preferred pattern when you already have Plotly-exported HTML artifacts under `<WORKDIR>/results/**`:

1. Pick a single Plotly export HTML that includes the Plotly.js bundle.
2. Extract that `plotly.js` bundle text and inline it once in the final report.
3. For each Plotly export to embed:
   - parse the `Plotly.newPlot(divId, data, layout, config)` call
   - extract `data/layout/config` payloads
   - recreate the plot in your report using a deterministic container ID

Outcome: all plots render offline, and the final report remains a single file.

## 3) Narrative / interpretation embedding

- Store interpretation in markdown under `<WORKDIR>/results/**.md`.
- Extract the relevant section (by heading or key phrase) and render markdown → HTML in the report.
- Avoid copying entire reports; include only the section needed for that tab.

### Pathway / PPI narrative expectations (recommended default)

- **Pathways (combos) — REQUIRED:** summarize the “Shared vs Distinguished Pathways” portion *as biology*, not as a long list of pathway names.
  - Include both: (a) shared (convergent) functions and (b) distinguished/unique (differentiating) functions, with a short interpretation of what that implies mechanistically.
- **Pathways (single genes) — REQUIRED:** summarize the *top biological functions/pathway themes* for the gene (e.g., immune activation, cytokine signaling, ECM remodeling), with interpretation.
  - Do **not** list pathway items directly from the markdown report.
- **PPI:** keep the PPI section focused on connectors/hops/mechanistic interpretation; do not duplicate the pathway list here.

## 4) Data traceability in report

For every score table shown:
- include subscore breakdown bars/tables where possible
- ensure each component can be traced to a loader output with `source_file` + raw value

## 4a) Cortellis gene summary visuals (required default)

For each gene card’s Summary section, render Cortellis XLSX-derived visuals (do not paste long lists from markdown):
- A **stacked bar plot** of **top primary indications**, stacked by **Highest Phase**.
- Phase inclusion rule (required default):
  - include **clinical and later** phases only (`Phase 1 Clinical`, `Phase 2 Clinical`, `Phase 3 Clinical`, `Phase I/II`, `Pre-registration`, `Registered`, `Launched/Marketed`)
  - fallback: if there are **no clinical-or-later** programs, include `Preclinical`, `Discovery` (if present), and `Discontinued`
- Discard NaN/blank indication labels and phase labels (do not render `"nan"`).
- Canonicalize phase labels so Cortellis-style labels (`Phase 1 Clinical`, `Phase 2 Clinical`, `Phase 3 Clinical`, `Phase I/II`) appear consistently in the stacked-bar legend.
- Source: `results/cortellis_<indication>/gene_cortellis_data.xlsx` or `results/cortellis_<indication>/<GENE>_cortellis_data.xlsx` (sheet `Drugs_Comprehensive`).
- Helper: `python scripts/summarize_cortellis_gene_xlsx.py --xlsx <path> --gene <GENE> --format json`

## 5) Heatmap Color Specification (REQUIRED)

All score heatmaps in the report MUST use a consistent color gradient:

| Score | Color | Hue | Description |
|-------|-------|-----|-------------|
| 100 | Green | 120° | High/Good |
| 50 | Yellow | 60° | Neutral |
| 0 | Red | 0° | Low/Bad |

**Implementation:** Linear HSL hue interpolation 0→60→120 (matching GI2 mockup).

```javascript
// REQUIRED heatmap color function (linear, no smoothstep)
function scoreHue(score){
  let t = Math.max(0, Math.min(100, score)) / 100;
  if (t <= 0.5) { const tt = t / 0.5; return 0 + (60 * tt); }  // red → yellow
  const tt = (t - 0.5) / 0.5;
  return 60 + (60 * tt);  // yellow → green
}

// Apply to cells:
cell.style.background = `linear-gradient(180deg, hsla(${hue}, 88%, 86%, .92), hsla(${hue}, 88%, 74%, .78))`;
cell.style.boxShadow = 'inset 0 0 0 1px rgba(12,19,35,.14), 0 10px 18px rgba(10,16,30,.08)';
cell.style.borderRadius = '12px';
```

**Important:** Ensure all scoring components use the same color scale:
- Overall score heatmap
- Subscore breakdown bars
- Ranking table cell backgrounds
- Individual gene score cells

## 6) Example QA artifact

`assets/example_report_offline.html` is included as a visual reference for layout/styling; use it for QA when porting to a new indication.
