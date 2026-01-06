# Cortellis CI Web Dashboard (Excel → Interactive HTML)

Convert Cortellis (or similar) Excel exports into a **single-file interactive HTML dashboard** for competitive intelligence.

## What this subskill is for

- Read and interpret user-provided Excel exports (often multi-sheet).
- Build a unified **asset table** with targets, MoA, mechanisms, phases, routes, and indications.
- Quantitatively score and rank **competition intensity vs blue-ocean opportunity** across:
  - Target
  - Target family (fine-grained)
  - Target type (broad: cytokine vs integrin vs kinase vs GPCR vs microbiome, etc.)
  - Modality / drug type
  - MoA / mechanism buckets
  - Route
- Generate an **interactive HTML** with the same layout/theme and core functions as the reference dashboard:
  - Global filters (clickable charts + checkboxes)
  - Click-to-open detail drawer
  - CSV export
  - Opportunity Lens selection plot (click bubbles to filter)

## Quick start

1) Inspect the workbook

```bash
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_ci_web/scripts/inspect_excel.py path/to/input.xlsx
```

2) Generate the dashboard

```bash
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_ci_web/scripts/generate_ci_dashboard.py   --input path/to/input.xlsx   --output ci_dashboard.html
```

3) (Optional) Provide a taxonomy config

```bash
python /home/sagemaker-user/.claude/skills/cortellis/cortellis_ci_web/scripts/generate_ci_dashboard.py   --input input.xlsx   --output ci_dashboard.html   --taxonomy /home/sagemaker-user/.claude/skills/cortellis/cortellis_ci_web/references/taxonomy_default.yaml
```

## Notes on “success rates”

The dashboard reports **observed maturity rates in the dataset** (e.g., fraction reaching clinical / late-stage / marketed within the filtered slice). These are **not** clinical PoS estimates.

## Key files

- Template: `assets/dashboard_template.html` (single-file, JS-driven)
- Generator: `scripts/generate_ci_dashboard.py`
- Taxonomy defaults: `references/taxonomy_default.yaml`
- Scoring notes: `references/opportunity_scoring.md`

## Opportunity Lens UX (current template)

- Opportunity subsections render as **scrollable tables** (not preformatted bullet blocks).
- Clicking a row opens a **detail drawer** with deeper stats (stage mix, dominant routes/MoA/modality/targets).
- The drawer includes **Toggle Filter**, which applies the selected group value globally (all charts + asset list).
