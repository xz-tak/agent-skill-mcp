#!/usr/bin/env python3
"""Clone a baseline offline report generator into a new, indication-specific tool file.

This is a lightweight scaffold generator meant for **human/agent editing** afterward.
It avoids hardcoding absolute paths: everything should be relative to <WORKDIR>.

In addition, it can initialize `<WORKDIR>/PLAN.md` (recommended) so work is tracked
explicitly and the report scope/scoring strategy are aligned before implementation.

Usage:
  python scripts/clone_report_tool.py --workdir <WORKDIR> --indication IBD
  python scripts/clone_report_tool.py --workdir <WORKDIR> --indication UC --out tools/generate_uc_combo_report.py

Optional:
  python scripts/clone_report_tool.py --workdir <WORKDIR> --indication IBD --no-init-plan
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re


def _replace_title_tokens(src: str, indication: str) -> str:
    # Conservative substitutions: adjust the visible title/subtitle strings only.
    src = re.sub(
        r'title="[^"]+"',
        f'title="{indication} Target Combination Prioritization Report"',
        src,
        count=1,
    )
    src = re.sub(
        r'nav_title="[^"]+"',
        f'nav_title="{indication} Target Combination Prioritization"',
        src,
    )
    src = re.sub(
        r'subtitle="[^"]+"',
        'subtitle="Interactive dossier with traceable narrative excerpts and embedded analysis."',
        src,
        count=1,
    )
    return src


def _replace_combos_block(src: str) -> str:
    """Replace the COMBOS constant with a placeholder list the user should edit.

    This assumes the baseline contains a `COMBOS: List[ComboSpec] = [` block.
    """

    start = src.find("COMBOS: List[ComboSpec] = [")
    if start < 0:
        return src

    end = src.find("]\n\n", start)
    if end < 0:
        return src

    block = "\n".join(
        [
            "COMBOS: List[ComboSpec] = [",
            "    # TODO: update combos for your indication/target set.",
            "    # Tips:",
            "    # - keep `genes` canonical symbols used in results artifacts",
            "    # - keep `sc_key` aligned with your sc_coexp markdown headings (if used)",
            "    # - list_num is optional, but useful when your artifacts are list-indexed",
            "    # Example:",
            "    # ComboSpec(1, (\"GENE1\", \"GENE2\"), \"list1_MECHANISM_KEY\"),",
            "]",
            "",
            "",
        ]
    )
    return src[:start] + block + src[end + 3 :]


def _init_plan_md(workdir: Path, indication: str) -> Path:
    plan_path = workdir / "PLAN.md"
    if plan_path.exists():
        return plan_path

    ind_lower = indication.lower()
    content = "\n".join(
        [
            "# PLAN",
            "",
            f"## Indication",
            f"- {indication}",
            "",
            "## Artifact Roots",
            "- `results/` (Cortellis/DEG/OFF-X/CI + optional analysis tabs + kgpred)",
            f"- `results/kgpred_{ind_lower}/` (BioBridge/ULTRA/PrimeKG)",
            f"  - `results/kgpred_{ind_lower}/biobridge/Part1_Individual_Analysis.md` (single genes)",
            f"  - `results/kgpred_{ind_lower}/biobridge/Part2_Combo_Analysis.md` (combos)",
            f"  - `results/kgpred_{ind_lower}/ultra/Part1_Individual_Analysis.md` (single genes; optional)",
            f"  - `results/kgpred_{ind_lower}/ultra/Part2_Combo_Analysis.md` (combos; optional)",
            f"  - `results/kgpred_{ind_lower}/primekg/Part1_Individual_Analysis.md` (single genes)",
            f"  - `results/kgpred_{ind_lower}/primekg/Part2_Combo_Analysis.md` (combos)",
            "",
            "## Resolved paths (edit if non-standard layout)",
            "- If an expected file is missing, search recursively under the corresponding section directory and record the resolved file(s) here.",
            f"- Clinical (Cortellis): `results/cortellis_{ind_lower}/...`",
            f"- Gene Summary (Cortellis XLSX): `results/cortellis_{ind_lower}/gene_cortellis_data.xlsx` or `results/cortellis_{ind_lower}/<GENE>_cortellis_data.xlsx`",
            f"  - REQUIRED visual: stacked bar for top indications, stacked by Highest Phase",
            f"  - helper: `python scripts/summarize_cortellis_gene_xlsx.py --xlsx <xlsx> --gene <GENE> --format json` (provides top_indications_by_phase)",
            f"- Disease (DEG): `results/deg_results_{ind_lower}/...` (or CSV sidecar)",
            f"- Disease (BioBridge): `results/kgpred_{ind_lower}/biobridge/Part1_Individual_Analysis.md` + `Part2_Combo_Analysis.md`",
            f"- Disease (ULTRA): `results/kgpred_{ind_lower}/ultra/Part1_Individual_Analysis.md` + `Part2_Combo_Analysis.md`",
            f"- Disease (PrimeKG): `results/kgpred_{ind_lower}/primekg/Part1_Individual_Analysis.md` + `Part2_Combo_Analysis.md`",
            f"- CI: `results/ci_{ind_lower}/*.html` with `<script id=\"data\">` JSON",
            f"- Safety: `results/offx_{ind_lower}/...` (optional)",
            f"- Pathways: `results/pathwaydb_{ind_lower}/COMPREHENSIVE_PATHWAY_ANALYSIS_REPORT.md` (optional; use for single-gene pathway theme summaries + list shared/distinguished summaries)",
            f"- PPI: `results/interactdb_{ind_lower}/interactdb_results/SUMMARY_REPORT.md` (optional)",
            "",
            "## Scoring Strategy (to confirm)",
            "- Default overall weights: Clinical=0.30, Disease=0.30, Safety=0.10, Opportunity=0.20, Novelty=0.10 (CI feeds Opportunity; no standalone CI weight)",
            "- Disease subweights (default): DEG=0.40, BioBridge=0.25, ULTRA=0.25, PrimeKG=0.10",
            "- Novelty: higher PrimeKG connectivity => lower literature novelty",
            "- Combo aggregation rules: mean across genes vs combo-level artifacts",
            "- Missing-data defaults (typically neutral 50)",
            "",
            "## Report Scope (to confirm)",
            "- Default scope: combos (ranked) + single-gene table + gene cards (same section order/layout as example HTML)",
            "- Executive summary: top 3 combos + rationale",
            "- Figures: which Plotly exports/PNGs to embed; omit gracefully if missing",
            "- Tabs: Overview / Disease / Bulk / Single-cell / Pathways / PPI (as available)",
            "",
            "## Steps",
            "- [ ] Validate required artifact directories/files",
            "- [ ] Implement/verify loaders (results + kgpred)",
            "- [ ] Implement/verify scoring (incl. novelty/connectivity)",
            "- [ ] Generate offline HTML report",
            "- [ ] QA: links/plots render offline; missing sections degrade cleanly",
            "",
        ]
    )

    plan_path.write_text(content, encoding="utf-8")
    return plan_path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True, help="Abstract project root (<WORKDIR>)")
    ap.add_argument("--indication", required=True, help="Indication label (e.g., IBD, UC, SSc)")
    ap.add_argument(
        "--baseline",
        default="tools/generate_combo_prioritization_report.py",
        help="Baseline report tool relative to <WORKDIR>",
    )
    ap.add_argument(
        "--out",
        default=None,
        help="Output tool path relative to <WORKDIR> (default: tools/generate_<indication>_combo_report.py)",
    )
    ap.add_argument(
        "--init-plan",
        action="store_true",
        help="Create <WORKDIR>/PLAN.md if missing (default behavior unless --no-init-plan)",
    )
    ap.add_argument(
        "--no-init-plan",
        action="store_true",
        help="Do not create <WORKDIR>/PLAN.md",
    )
    args = ap.parse_args()

    workdir = Path(args.workdir).resolve()

    init_plan = args.init_plan or (not args.no_init_plan)
    if init_plan:
        plan_path = _init_plan_md(workdir, args.indication)
        print(f"PLAN.md: {plan_path}")

    baseline = (workdir / args.baseline).resolve()
    if not baseline.exists():
        raise SystemExit(f"Baseline not found: {baseline}")

    out_rel = args.out or f"tools/generate_{args.indication.lower()}_combo_report.py"
    out_path = (workdir / out_rel).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    src = baseline.read_text(encoding="utf-8", errors="replace")
    src = _replace_title_tokens(src, args.indication)
    src = _replace_combos_block(src)

    out_path.write_text(src, encoding="utf-8")
    print(f"Wrote scaffold tool: {out_path}")
    print("Next: edit COMBOS and any indication-specific result file paths in the tool.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
