# Example Notes (GI2)

This skill bundles GI2 artifacts as *examples* to preserve formatting + workflow clarity.

Note: the bundled GI2 design/spec artifacts may reflect **legacy path conventions** (e.g., `results/biobridge/` or `results/primekg/`) and older PrimeKG “connection count” proxies. Treat the **contracts** as the source of truth:
- `references/results-contract.md` (includes `kgpred_<indication>/...` Step1/Step2)
- `references/scoring-contract.md` (PrimeKG connectivity-based novelty)

## Included example files

- Design spec (example snapshot):
  - `references/example_gi2_design_2025-12-23_redesign.md`
- Report generator (example implementation):
  - `references/example_gi2_generate_combo_prioritization_report.py`
- Offline HTML report (example output; visual QA only):
  - `assets/example_report_offline.html`
- HTML template (format contract):
  - `assets/combo_prioritization_report_template.html`

## How to use the examples

- Use the **template + report generator** as the canonical source for layout, styles, tab behavior, and Plotly embedding.
- Use the **design spec** as a blueprint for writing new design docs (source mapping + transformation workflows + scoring formulas + report compilation steps).
- Use the **contracts** to resolve path/scoring conflicts when examples differ from current requirements.

## Adapting the report generator

The two primary adaptation surfaces are:

1) **Target manifest** (your new genes/combos + any “list” numbering + any markdown keys).
2) **Artifact discovery rules** (how to match target → plotly html / png / md blocks under `results/` and `kgpred_<indication>/`).

Everything else should remain stable to preserve the report format.

## Updated helper scripts (recommended)

These scripts are bundled with the skill and match the current `kgpred_<indication>/...` + `PLAN.md` workflow:

- Validate required inputs early:
  - `python scripts/validate_results_tree.py --workdir <WORKDIR> --indication <IND>`
- Scaffold a new indication tool and initialize `<WORKDIR>/PLAN.md` (created only if missing):
  - `python scripts/clone_report_tool.py --workdir <WORKDIR> --indication <IND>`

