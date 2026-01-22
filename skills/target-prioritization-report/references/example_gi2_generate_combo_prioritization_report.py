#!/usr/bin/env python
"""
Generate an interactive IBD combination prioritization HTML report.

Key goals:
- Keep the report fully offline and standalone by extracting Plotly figures from
  the existing HTML artifacts under `results/` and inlining a single Plotly.js bundle.
- Pull interpretation text from the corresponding Markdown reports under `results/`.
- Provide a Disease Association tab with a horizontal bar visualization of the
  disease association subscore component values.
"""

from __future__ import annotations

import json
import os
import re
import sys
import base64
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import markdown as md
from jinja2 import Template
from bs4 import BeautifulSoup


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
RESULTS_DIR = REPO_ROOT / "results"
REPORTS_DIR = REPO_ROOT / "reports"
TEMPLATE_PATH = REPO_ROOT / "target_prioritization" / "templates" / "combo_prioritization_report_template.html"
INDICATION = "IBD"


def _resolve_kgpred_root(indication: str) -> Optional[Path]:
    ind = (indication or "").strip()
    if not ind:
        return None
    candidates = [
        RESULTS_DIR / f"kgpred_{ind.lower()}",
        RESULTS_DIR / f"kgpred_{ind}",
        RESULTS_DIR / f"kgpred_{ind.upper()}",
        REPO_ROOT / f"kgpred_{ind.lower()}",
        REPO_ROOT / f"kgpred_{ind}",
        REPO_ROOT / f"kgpred_{ind.upper()}",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


@dataclass(frozen=True)
class ComboSpec:
    list_num: int
    genes: Tuple[str, ...]
    sc_key: str  # key used in sc_coexp markdown sections

    @property
    def name(self) -> str:
        return "-".join(self.genes)

    @property
    def anchor(self) -> str:
        return f"list{self.list_num}-{self.name}".lower()

    @property
    def bulk_base(self) -> str:
        return f"list{self.list_num}_" + "_".join(self.genes)

    @property
    def pathway_base(self) -> str:
        return f"list{self.list_num}_" + "_".join(self.genes)


GENE_TITLES = {
    "TYK2": "Tyrosine kinase 2",
    "JAK1": "Janus kinase 1",
    "TNFRSF25": "TNF receptor superfamily member 25 (DR3)",
    "GREM1": "Gremlin 1",
    "PCOLCE": "Procollagen C-endopeptidase enhancer",
    "CDKN2D": "Cyclin dependent kinase inhibitor 2D",
    "ITGA4": "Integrin subunit alpha 4",
    "ITGB7": "Integrin subunit beta 7",
}

COMBOS: List[ComboSpec] = [
    ComboSpec(1, ("TYK2", "JAK1"), "list1_JAK_STAT"),
    ComboSpec(2, ("TNFRSF25", "GREM1"), "list2_TNFR_BMP"),
    ComboSpec(3, ("TNFRSF25", "PCOLCE"), "list3_TNFR_PCOLCE"),
    ComboSpec(4, ("CDKN2D", "ITGA4", "ITGB7"), "list4_cell_cycle_integrin"),
    ComboSpec(5, ("CDKN2D", "PCOLCE"), "list5_CDKN2D_PCOLCE"),
]


def _rel_from_reports(path: Optional[Path]) -> Optional[str]:
    if path is None:
        return None
    rel = os.path.relpath(path, REPORTS_DIR)
    return rel.replace(os.sep, "/")


def _read_text(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    return path.read_text(encoding="utf-8", errors="replace")


def _extract_heading_block(text: str, heading_prefix: str, startswith: str) -> Optional[str]:
    """
    Extract a markdown block starting with a specific heading.

    - heading_prefix: e.g. "##", "###", "####"
    - startswith: heading text (matched as line.startswith(f"{prefix} {startswith}"))
    """
    lines = text.splitlines()
    start_idx = None
    marker = f"{heading_prefix} {startswith}"
    for i, line in enumerate(lines):
        if line.startswith(marker):
            start_idx = i
            break
    if start_idx is None:
        return None

    next_idx = len(lines)
    for j in range(start_idx + 1, len(lines)):
        if lines[j].startswith(heading_prefix + " "):
            next_idx = j
            break
    return "\n".join(lines[start_idx:next_idx]).strip()


def _md_to_html(fragment: str) -> str:
    return md.markdown(
        fragment,
        extensions=["extra", "sane_lists", "tables"],
        output_format="html5",
    )


def _find_first_existing(paths: List[Path]) -> Optional[Path]:
    for p in paths:
        if p.exists():
            return p
    return None


def _bulk_variants(list_num: int) -> List[Dict[str, Optional[Path]]]:
    """
    Discover available bulk coexpression Plotly exports for a given list.

    Expects filenames like:
      - list1_<...>_<tissue>_correlation_heatmap.html
      - list1_<...>_<tissue>_expression_boxplot.html

    and supports both `results/bulk_coexpression/results/` and `results/bulk_coexpression/`.
    """

    roots = [RESULTS_DIR / "bulk_coexpression" / "results", RESULTS_DIR / "bulk_coexpression"]
    heatmap_suffix = "_correlation_heatmap"
    box_suffix = "_expression_boxplot"
    list_prefix = f"list{list_num}_"

    def _body_from_path(p: Path, suffix: str) -> Optional[str]:
        stem = p.stem
        if not stem.startswith(list_prefix) or not stem.endswith(suffix):
            return None
        return stem[len(list_prefix) : -len(suffix)]

    heatmaps: List[Path] = []
    boxplots: List[Path] = []
    for root in roots:
        if not root.exists():
            continue
        heatmaps.extend(list(root.glob(f"{list_prefix}*{heatmap_suffix}.html")))
        boxplots.extend(list(root.glob(f"{list_prefix}*{box_suffix}.html")))

    hm_by_body: Dict[str, Path] = {}
    bx_by_body: Dict[str, Path] = {}
    for p in heatmaps:
        body = _body_from_path(p, heatmap_suffix)
        if body:
            hm_by_body.setdefault(body, p)
    for p in boxplots:
        body = _body_from_path(p, box_suffix)
        if body:
            bx_by_body.setdefault(body, p)

    def _ctx_from_body(body: str) -> str:
        parts = [x for x in body.split("_") if x]
        return parts[-1] if parts else "default"

    by_ctx: Dict[str, Dict[str, Optional[Path]]] = {}
    for body, hm in hm_by_body.items():
        ctx = _ctx_from_body(body)
        by_ctx.setdefault(ctx, {"key": ctx, "heatmap_path": None, "boxplot_path": None})
        by_ctx[ctx]["heatmap_path"] = hm
        if body in bx_by_body:
            by_ctx[ctx]["boxplot_path"] = bx_by_body[body]
    for body, bx in bx_by_body.items():
        ctx = _ctx_from_body(body)
        by_ctx.setdefault(ctx, {"key": ctx, "heatmap_path": None, "boxplot_path": None})
        by_ctx[ctx]["boxplot_path"] = bx
        if body in hm_by_body:
            by_ctx[ctx]["heatmap_path"] = hm_by_body[body]

    prefs = ["lung", "colon", "ileum", "skin", "blood"]
    keys = list(by_ctx.keys())
    keys.sort(key=lambda k: (prefs.index(k) if k in prefs else 99, k))
    return [by_ctx[k] for k in keys]


def _sc_variant_key(path: Path, ct_key: str, sc_key: str, list_num: int) -> str:
    stem = path.stem
    if stem.endswith("_heatmap"):
        stem = stem[: -len("_heatmap")]
    prefix = f"{ct_key}_"
    if stem.startswith(prefix):
        stem = stem[len(prefix) :]
    if stem == sc_key:
        return "default"
    if stem.startswith(sc_key + "_"):
        return stem[len(sc_key) + 1 :] or "default"
    marker = f"list{list_num}_"
    if marker in stem:
        return stem.split(marker, 1)[1] or "default"
    return stem or "default"


def _png_to_data_uri(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    b = path.read_bytes()
    return "data:image/png;base64," + base64.b64encode(b).decode("ascii")


def _rel_from_repo_root(path: Optional[Path]) -> Optional[str]:
    if path is None:
        return None
    try:
        return str(path.relative_to(REPO_ROOT)).replace(os.sep, "/")
    except Exception:
        return str(path).replace(os.sep, "/")


def _file_to_b64(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    return base64.b64encode(path.read_bytes()).decode("ascii")


def _ci_dashboard_b64(path: Path, author_line: str) -> Optional[str]:
    """
    Load the offline CI dashboard HTML and inject an author line below the generated timestamp.
    This keeps the main report standalone (srcdoc iframe) without modifying the source dashboard file.
    """
    html = _read_text(path)
    if not html:
        return None

    marker = '<span class="pill">Generated <span id="genAt"></span></span>'
    if marker in html and author_line not in html:
        html = html.replace(
            marker,
            marker
            + '\n          <span style="flex-basis: 100%; height: 0;"></span>\n'
            + f'          <span class="pill">{author_line}</span>',
            1,
        )

    return base64.b64encode(html.encode("utf-8")).decode("ascii")


def _extract_plotly_bundle_js(plot_html: str) -> str:
    soup = BeautifulSoup(plot_html, "html.parser")
    for script in soup.find_all("script"):
        if script.string and "plotly.js v" in script.string:
            return script.string
    raise ValueError("Could not find Plotly bundle script in plot HTML.")


def _extract_plotly_newplot_script(plot_html: str) -> str:
    soup = BeautifulSoup(plot_html, "html.parser")
    for script in soup.find_all("script"):
        if script.string and "Plotly.newPlot" in script.string and "plotly.js v" not in script.string:
            return script.string
    raise ValueError("Could not find Plotly.newPlot script in plot HTML.")


def _skip_ws(s: str, i: int) -> int:
    while i < len(s) and s[i].isspace():
        i += 1
    return i


def _extract_balanced(s: str, i: int, open_ch: str, close_ch: str) -> Tuple[str, int]:
    """
    Return (substring_including_brackets, next_index_after_substring).
    Handles nested brackets and JS string literals.
    """
    if i >= len(s) or s[i] != open_ch:
        raise ValueError(f"Expected {open_ch!r} at {i}")
    depth = 0
    start = i
    in_str: Optional[str] = None
    escaped = False
    while i < len(s):
        ch = s[i]
        if in_str:
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == in_str:
                in_str = None
        else:
            if ch in ("'", '"', "`"):
                in_str = ch
            elif ch == open_ch:
                depth += 1
            elif ch == close_ch:
                depth -= 1
                if depth == 0:
                    i += 1
                    return s[start:i], i
        i += 1
    raise ValueError("Unterminated balanced section")


def _extract_plotly_figure_args(plot_html: str) -> Tuple[str, str, str]:
    """
    Extract (data_js, layout_js, config_js) from the Plotly.newPlot call.
    """
    script = _extract_plotly_newplot_script(plot_html)
    idx = script.find("Plotly.newPlot")
    if idx < 0:
        raise ValueError("Plotly.newPlot not found")
    idx = script.find("(", idx)
    if idx < 0:
        raise ValueError("Plotly.newPlot '(' not found")
    i = idx + 1
    i = _skip_ws(script, i)

    # arg1: div id string or variable; skip to first comma at top-level (no parens expected here).
    # Handle quoted strings.
    if script[i] in ("'", '"'):
        q = script[i]
        i += 1
        while i < len(script):
            if script[i] == "\\":
                i += 2
                continue
            if script[i] == q:
                i += 1
                break
            i += 1
    else:
        while i < len(script) and script[i] != ",":
            i += 1

    # comma after arg1
    while i < len(script) and script[i] != ",":
        i += 1
    i += 1
    i = _skip_ws(script, i)

    data_js, i = _extract_balanced(script, i, "[", "]")
    i = _skip_ws(script, i)
    if script[i] == ",":
        i += 1
    i = _skip_ws(script, i)
    layout_js, i = _extract_balanced(script, i, "{", "}")
    i = _skip_ws(script, i)

    config_js = "{responsive: true}"
    if i < len(script) and script[i] == ",":
        i += 1
        i = _skip_ws(script, i)
        if i < len(script) and script[i] == "{":
            config_js, i = _extract_balanced(script, i, "{", "}")

    return data_js.strip(), layout_js.strip(), config_js.strip()


def _js_obj_literal(obj: Dict[str, Dict[str, str]]) -> str:
    """
    Render a JS object literal with embedded raw JS for data/layout/config.
    """
    parts = ["{"]
    first = True
    for key, fig in obj.items():
        if not first:
            parts.append(",")
        first = False
        k = json.dumps(key)
        parts.append(
            f"\n  {k}: {{ data: {fig['data']}, layout: {fig['layout']}, config: {fig['config']} }}"
        )
    parts.append("\n}")
    return "".join(parts)


def _parse_deg_table_interpretations(deg_report_md: str) -> Dict[str, str]:
    """
    Extract the 'Interpretation' column from the 'Target Gene Rankings' table.
    """
    lines = deg_report_md.splitlines()
    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith("|") and "Interpretation" in line and "Gene" in line and "Score" in line:
            header_idx = i
            break

    if header_idx is None:
        return {}

    headers = [c.strip() for c in lines[header_idx].strip().strip("|").split("|")]
    interpretation_idx = None
    gene_idx = None
    for i, h in enumerate(headers):
        if h.lower().strip() == "interpretation":
            interpretation_idx = i
        if h.lower().strip() == "gene":
            gene_idx = i

    if interpretation_idx is None or gene_idx is None:
        return {}

    start = header_idx + 2
    out: Dict[str, str] = {}
    for line in lines[start:]:
        if not line.strip().startswith("|"):
            break
        row = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(row) <= max(interpretation_idx, gene_idx):
            continue
        gene = re.sub(r"[*_`]", "", row[gene_idx]).strip()
        interp = re.sub(r"[*_`]", "", row[interpretation_idx]).strip()
        if gene:
            out[gene] = interp
    return out


def _parse_biobridge_interpretation(gene: str) -> Optional[str]:
    # Preferred (current contract): results/kgpred_<indication>/biobridge/Part1_Individual_Analysis.md
    kgpred_root = _resolve_kgpred_root(INDICATION)
    if kgpred_root:
        p1 = kgpred_root / "biobridge" / "Part1_Individual_Analysis.md"
        content = _read_text(p1)
        if content:
            # Extract from the "Individual Gene Results" table if present.
            lines = content.splitlines()
            header_idx = None
            for i, line in enumerate(lines):
                if line.strip().startswith("|") and "Gene" in line and "pct_rank" in line:
                    header_idx = i
                    break
            if header_idx is not None:
                headers = [c.strip() for c in lines[header_idx].strip().strip("|").split("|")]
                try:
                    gi = [h.lower() for h in headers].index("gene")
                    pi = [h.lower() for h in headers].index("pct_rank")
                    ii = [h.lower() for h in headers].index("interpretation")
                except Exception:
                    gi = pi = ii = None  # type: ignore[assignment]
                if gi is not None and pi is not None:
                    for line in lines[header_idx + 2 :]:
                        if not line.strip().startswith("|"):
                            break
                        row = [c.strip() for c in line.strip().strip("|").split("|")]
                        if len(row) <= max(gi, pi):
                            continue
                        g = re.sub(r"[*_`]", "", row[gi]).strip()
                        if g != gene:
                            continue
                        pct = re.sub(r"[*_`]", "", row[pi]).strip()
                        interp = ""
                        if ii is not None and len(row) > ii:
                            interp = re.sub(r"[*_`]", "", row[ii]).strip()
                        return f"pct_rank {pct} — {interp}".strip(" —")

    # Legacy fallback: results/biobridge/individual_report/<gene>-<indication>.md
    candidates = [
        RESULTS_DIR / "biobridge" / "individual_report" / f"{gene}-{INDICATION}.md",
        RESULTS_DIR / "biobridge" / "individual_report" / f"DR3-{INDICATION}.md" if gene == "TNFRSF25" else None,
    ]
    candidates = [p for p in candidates if p is not None]
    for path in candidates:
        content = _read_text(path)
        if not content:
            continue
        m = re.search(r"\*\*Interpretation:\*\*\s*(.+)", content)
        if m:
            return m.group(1).strip()
    return None


def _parse_primekg_interpretation(gene: str) -> Optional[str]:
    # Preferred (current contract): results/kgpred_<indication>/primekg/Part1_Individual_Analysis.md
    kgpred_root = _resolve_kgpred_root(INDICATION)
    if kgpred_root:
        p1 = kgpred_root / "primekg" / "Part1_Individual_Analysis.md"
        content = _read_text(p1)
        if content:
            # Extract from the "Individual Gene Scores" table if present.
            lines = content.splitlines()
            header_idx = None
            for i, line in enumerate(lines):
                if line.strip().startswith("|") and "Gene" in line and "Path Length" in line and "Score" in line:
                    header_idx = i
                    break
            if header_idx is not None:
                headers = [c.strip() for c in lines[header_idx].strip().strip("|").split("|")]
                hl = [h.lower() for h in headers]
                try:
                    gi = hl.index("gene")
                    pli = hl.index("path length")
                    si = hl.index("score")
                except Exception:
                    gi = pli = si = None  # type: ignore[assignment]
                if gi is not None and pli is not None and si is not None:
                    for line in lines[header_idx + 2 :]:
                        if not line.strip().startswith("|"):
                            break
                        row = [c.strip() for c in line.strip().strip("|").split("|")]
                        if len(row) <= max(gi, pli, si):
                            continue
                        g = re.sub(r"[*_`]", "", row[gi]).strip()
                        if g != gene:
                            continue
                        pl = re.sub(r"[*_`]", "", row[pli]).strip()
                        score = re.sub(r"[*_`]", "", row[si]).strip()
                        return f"Shortest path length {pl}, score {score}."

    # Legacy fallback: results/primekg/<gene>-<indication>.md
    candidates = [
        RESULTS_DIR / "primekg" / f"{gene}-{INDICATION}.md",
        RESULTS_DIR / "primekg" / f"DR3-{INDICATION}.md" if gene == "TNFRSF25" else None,
    ]
    candidates = [p for p in candidates if p is not None]
    for path in candidates:
        content = _read_text(path)
        if not content:
            continue
        lines = [ln.strip() for ln in content.splitlines() if ln.strip()]
        excerpt = []
        for ln in lines[1:]:
            if ln.startswith("#"):
                continue
            if re.match(r"^:?\s*\|", ln):
                break
            if re.search(r"\|\s*drugbank\s*id\s*\|", ln, re.IGNORECASE):
                break
            excerpt.append(ln)
            if len(excerpt) >= 4:
                break
        if excerpt:
            return " ".join(excerpt)
    return None


def _parse_bulk_list_interpretation(list_num: int) -> Dict[str, str]:
    report_path = RESULTS_DIR / "bulk_coexpression" / "ANALYSIS_REPORT.md"
    content = _read_text(report_path)
    if not content:
        return {"html": "<p>Bulk report not found.</p>"}

    block = _extract_heading_block(content, "##", f"List {list_num}:")
    if not block:
        # fallback: some reports might use "## List 1" without colon
        block = _extract_heading_block(content, "##", f"List {list_num}")
    if not block:
        return {"html": "<p>No bulk section found.</p>"}

    # Keep only the Interpretation subsection if present; otherwise keep the whole block (trimmed).
    lines = block.splitlines()
    interp_idx = None
    for i, ln in enumerate(lines):
        if ln.strip().startswith("**Interpretation:**"):
            interp_idx = i
            break
    if interp_idx is not None:
        snippet = "\n".join(lines[interp_idx:]).strip()
    else:
        snippet = block.strip()

    return {"html": _md_to_html(snippet)}


def _parse_sc_celltype_interpretation(
    sc_md: str,
    cell_heading: str,
    list_key: str,
) -> Tuple[Optional[str], Optional[str]]:
    """
    Returns (interpretation_md, heatmap_rel_path_from_ibd_analysis_results).
    """
    cell_block = _extract_heading_block(sc_md, "###", cell_heading)
    if not cell_block:
        return None, None

    list_block = _extract_heading_block(cell_block, "####", list_key)
    if not list_block:
        return None, None

    lines = list_block.splitlines()

    # Heatmap path
    heatmap_rel = None
    for ln in lines:
        m = re.search(r"\*\*Interactive Heatmap:\*\*\s*`([^`]+)`", ln)
        if m:
            heatmap_rel = m.group(1).strip()
            break

    # Interpretation block
    interp_start = None
    for i, ln in enumerate(lines):
        if ln.strip().startswith("**Interpretation:**"):
            interp_start = i + 1
            break

    interpretation_lines: List[str] = []
    if interp_start is not None:
        for ln in lines[interp_start:]:
            if ln.strip().startswith("**Interactive Heatmap:**"):
                break
            interpretation_lines.append(ln)

    interpretation_md = "\n".join(interpretation_lines).strip() if interpretation_lines else None
    return interpretation_md, heatmap_rel


def _pathway_overview_html(combo: ComboSpec) -> str:
    path = RESULTS_DIR / "pathwaydb" / f"REPORT_{combo.pathway_base}.md"
    content = _read_text(path)
    if not content:
        return "<p>No pathway report found.</p>"

    overview = _extract_heading_block(content, "##", "Overview")
    if not overview:
        overview = content.split("\n\n", 1)[0]
    return _md_to_html(overview)


def _pathway_combo_summary_html(combo: ComboSpec) -> str:
    """
    Summarize shared vs distinguished pathways from the list-level pathway report.
    Do not dump ranked pathway lists; synthesize into biological themes + interpretation.
    """
    path = RESULTS_DIR / "pathwaydb" / f"REPORT_{combo.pathway_base}.md"
    content = _read_text(path)
    if not content:
        return "<p>No pathway report found.</p>"

    overview = _extract_heading_block(content, "##", "Overview") or ""
    total = shared = distinguished = None
    for ln in overview.splitlines():
        m = re.search(r"Total unique pathways:\s*(\d+)", ln)
        if m:
            total = int(m.group(1))
        m = re.search(r"Shared pathways.*:\s*(\d+)", ln)
        if m:
            shared = int(m.group(1))
        m = re.search(r"Distinguished pathways.*:\s*(\d+)", ln)
        if m:
            distinguished = int(m.group(1))

    # Extract top pathway names from each gene section and map to coarse themes.
    theme_map = {
        "immune": ["IMMUNE", "LYMPHOCYTE", "LEUKOCYTE", "T_CELL", "B_CELL", "INTERFERON", "INTERLEUKIN", "CYTOKINE"],
        "tnf_death": ["TNF", "APOPT", "DEATH_RECEPTOR", "NECRO", "NF_KB"],
        "ecm_fibrosis": ["COLLAGEN", "EXTRACELLULAR_MATRIX", "ECM", "MESENCHYM", "FIBRO"],
        "cell_cycle": ["CELL_CYCLE", "MITOT", "G1", "SENESC", "DNA_REPAIR", "HEMOPOIESIS"],
        "integrin_traffic": ["INTEGRIN", "ADHESION", "HOMING", "MIGRATION"],
        "bmp_wnt": ["BMP", "WNT", "TGF", "SMAD"],
    }

    def themes_for(names: List[str]) -> List[str]:
        out = set()
        for nm in names:
            up = nm.upper()
            for theme, keys in theme_map.items():
                if any(k in up for k in keys):
                    out.add(theme)
        return sorted(out)

    per_gene_themes = {}
    per_gene_names = {}
    for g in combo.genes:
        section = _extract_heading_block(content, "###", g)
        if not section:
            continue
        names = []
        for ln in section.splitlines():
            # lines like: "1. **GOBP_XYZ ...**"
            m = re.match(r"\d+\.\s+\*\*([^*]+)\*\*", ln.strip())
            if m:
                names.append(m.group(1).strip())
            if len(names) >= 10:
                break
        per_gene_names[g] = names
        per_gene_themes[g] = set(themes_for(names))

    if per_gene_themes:
        shared_themes = set.intersection(*per_gene_themes.values()) if len(per_gene_themes) > 1 else set()
        union_themes = set.union(*per_gene_themes.values())
        distinct_themes = union_themes - shared_themes
    else:
        shared_themes, distinct_themes = set(), set()

    def pretty(theme: str) -> str:
        return {
            "immune": "immune activation / cytokine signaling",
            "tnf_death": "TNF / cell-death programs",
            "ecm_fibrosis": "ECM remodeling / fibrosis",
            "cell_cycle": "cell-cycle / stress programs",
            "integrin_traffic": "integrin adhesion / trafficking",
            "bmp_wnt": "BMP/WNT/TGF signaling",
        }.get(theme, theme)

    lines = []
    counts = []
    if total is not None:
        counts.append(f"**{total}** total unique pathways")
    if shared is not None:
        counts.append(f"**{shared}** shared (multi-gene)")
    if distinguished is not None:
        counts.append(f"**{distinguished}** distinguished (single-gene)")
    if counts:
        lines.append("### Network overlap snapshot")
        lines.append("- " + " · ".join(counts))
        lines.append("")

    shared_pretty = [pretty(t) for t in sorted(shared_themes)]
    distinct_pretty = [pretty(t) for t in sorted(distinct_themes)]

    def interp_for(theme: str) -> str:
        return {
            "immune activation / cytokine signaling": "convergent inflammatory signaling consistent with immune-driven disease biology",
            "TNF / cell-death programs": "shared inflammatory stress and tissue-damage programs",
            "ECM remodeling / fibrosis": "tissue remodeling and chronic lesion architecture",
            "cell-cycle / stress programs": "proliferation/stress programs that can reflect tissue turnover and repair",
            "integrin adhesion / trafficking": "cell trafficking/adhesion programs relevant to immune infiltration",
            "BMP/WNT/TGF signaling": "developmental and differentiation cues linked to barrier remodeling",
        }.get(theme, "coherent pathway theme")

    if shared_pretty:
        lines.append("### Shared functional biology (convergent)")
        lines.append(
            "Shared themes indicate the combination converges on: "
            + "; ".join([f"**{t}** ({interp_for(t)})" for t in shared_pretty])
            + "."
        )
        lines.append("")

    if distinct_pretty:
        lines.append("### Distinguished/unique biology (differentiating)")
        lines.append(
            "Distinguished themes suggest additional leverage points contributed by individual genes: "
            + "; ".join([f"**{t}** ({interp_for(t)})" for t in distinct_pretty])
            + "."
        )
        lines.append("")

    lines.append("### Interpretation")
    lines.append(
        f"This list’s pathway profile supports a **{combo.sc_key.replace('_', ' ')}** mechanistic framing: "
        "shared biology reflects convergent disease-relevant processes, while distinguished biology highlights non-overlapping mechanism coverage."
    )

    return _md_to_html("\n".join(lines))


def _ppi_excerpt_html(combo: ComboSpec) -> str:
    path = RESULTS_DIR / "interactdb" / "shortest_paths_analysis_complete" / "COMPLETE_ANALYSIS_REPORT_ALL_DATABASES.md"
    # Common alternative in other workdirs:
    if not path.exists():
        path = RESULTS_DIR / "interactdb" / "interactdb_results" / "SUMMARY_REPORT.md"
    content = _read_text(path)
    if not content:
        return "<p>PPI report not found.</p>"

    # Identify list section by list number.
    section = _extract_heading_block(content, "###", f"List {combo.list_num}:")
    if not section:
        # fallback: match by gene arrows
        if len(combo.genes) >= 2:
            section = _extract_heading_block(content, "###", f"List {combo.list_num}: {combo.genes[0]} ↔ {combo.genes[1]}")

    if not section:
        return "<p>No matching PPI section found.</p>"

    # Keep the executive part of the section (first ~140 lines) so hops / connectors are retained.
    lines = section.splitlines()
    trimmed = "\n".join(lines[:140]).strip()
    return _md_to_html(trimmed)


def _gene_pathway_summary_html(gene: str) -> str:
    """
    Summarize a gene's top pathways from the list-level pathway centrality markdown
    (avoids dumping raw pathway lists).
    """
    # Find a list report that contains this gene (stable in this workspace).
    report = None
    for spec in COMBOS:
        if gene in spec.genes:
            candidate = RESULTS_DIR / "pathwaydb" / f"REPORT_{spec.pathway_base}.md"
            if candidate.exists():
                report = candidate
                break
    if report is None:
        return "<p>No pathway report found for this gene.</p>"

    content = _read_text(report)
    if not content:
        return "<p>No pathway report found for this gene.</p>"

    section = _extract_heading_block(content, "###", gene)
    if not section:
        return "<p>No gene pathway section found.</p>"

    stats = {}
    for ln in section.splitlines():
        m = re.search(r"Pathways in subnetwork:\s*(\d+)", ln)
        if m:
            stats["subnetwork"] = int(m.group(1))
        m = re.search(r"Hub pathways:\s*(\d+)", ln)
        if m:
            stats["hub"] = int(m.group(1))
        m = re.search(r"Bridge pathways:\s*(\d+)", ln)
        if m:
            stats["bridge"] = int(m.group(1))
        m = re.search(r"Leaf pathways:\s*(\d+)", ln)
        if m:
            stats["leaf"] = int(m.group(1))

    # Extract top central pathway names, then summarize into themes and a few exemplars.
    top_names = []
    for ln in section.splitlines():
        m = re.match(r"\d+\.\s+\*\*([^*]+)\*\*", ln.strip())
        if m:
            top_names.append(m.group(1).strip())
        if len(top_names) >= 10:
            break

    keywords = {
        "immune activation / cytokine signaling": ["IMMUNE", "CYTOKINE", "INTERFERON", "INTERLEUKIN", "LYMPHOCYTE", "T_CELL", "B_CELL"],
        "ECM remodeling / fibrosis": ["COLLAGEN", "EXTRACELLULAR_MATRIX", "ECM", "MESENCHYM", "FIBRO"],
        "cell-cycle / stress programs": ["CELL_CYCLE", "MITOT", "SENESC", "DNA_REPAIR", "HEMOPOIESIS"],
        "integrin adhesion / trafficking": ["INTEGRIN", "ADHESION", "MIGRATION", "HOMING"],
        "TNF / apoptosis programs": ["TNF", "APOPT", "DEATH_RECEPTOR", "NECRO"],
        "BMP/WNT/TGF signaling": ["BMP", "WNT", "TGF", "SMAD"],
    }

    theme_scores = {k: 0 for k in keywords}
    for nm in top_names:
        up = nm.upper()
        for theme, keys in keywords.items():
            if any(k in up for k in keys):
                theme_scores[theme] += 1

    top_themes = [t for t, s in sorted(theme_scores.items(), key=lambda kv: kv[1], reverse=True) if s > 0][:2]
    md_lines = []
    if stats:
        md_lines.append("**Centrality snapshot:** " + " · ".join([f"{k} {v}" for k, v in stats.items()]))
        md_lines.append("")
    if top_themes:
        md_lines.append("**Dominant themes:** " + "; ".join(top_themes))
        md_lines.append("")
        md_lines.append(
            "**Interpretation:** these themes summarize the gene’s most central biological functions in the pathway network and provide a mechanistic lens for prioritization (without listing ranked pathway names)."
        )

    return _md_to_html("\n".join(md_lines) if md_lines else "No pathway summary available.")


def _compute_scores() -> Tuple[List[Dict[str, Any]], int, Dict[str, Dict[str, Any]]]:
    """
    Returns:
      - combo payloads with scores
      - max_primekg used for normalization
      - scored per-gene payloads (for per-combo gene breakdown tables)
    """
    from target_prioritization.data_loaders.orchestrator import load_all_data
    from target_prioritization.scoring.scoring import score_target

    # Determine max_primekg across unique genes in this run.
    unique_genes = sorted({g for c in COMBOS for g in c.genes})
    max_primekg = 0
    gene_data: Dict[str, Dict[str, Any]] = {}
    for gene in unique_genes:
        data = load_all_data([gene], results_dir=str(RESULTS_DIR), is_combo=False)
        gene_data[gene] = data
        if "primekg" in data and isinstance(data["primekg"].get("primekg_connections_raw"), int):
            max_primekg = max(max_primekg, int(data["primekg"]["primekg_connections_raw"]))

    gene_scored: Dict[str, Dict[str, Any]] = {}
    for gene in unique_genes:
        gene_scored[gene] = score_target(gene_data[gene], max_primekg=max_primekg)

    combo_payloads: List[Dict[str, Any]] = []
    for combo in COMBOS:
        data = load_all_data(list(combo.genes), results_dir=str(RESULTS_DIR), is_combo=True)
        scored = score_target(data, max_primekg=max_primekg)
        combo_payloads.append({"spec": combo, "data": data, "scored": scored})

    return combo_payloads, max_primekg, gene_scored


def main() -> int:
    if not TEMPLATE_PATH.exists():
        raise FileNotFoundError(f"Template not found: {TEMPLATE_PATH}")

    REPORTS_DIR.mkdir(parents=True, exist_ok=True)

    sc_md_path = RESULTS_DIR / "sc_coexp" / "ibd_analysis_results" / "COMPREHENSIVE_REPORT.md"
    bulk_md_path = RESULTS_DIR / "bulk_coexpression" / "ANALYSIS_REPORT.md"
    deg_md_path = RESULTS_DIR / "deg_results" / "IBD_Targets_Summary_Report.md"

    sc_md = _read_text(sc_md_path) or ""
    bulk_md = _read_text(bulk_md_path) or ""
    deg_md = _read_text(deg_md_path) or ""

    deg_interps = _parse_deg_table_interpretations(deg_md)

    payloads, max_primekg, gene_scored = _compute_scores()

    # Sort by overall score descending
    payloads.sort(key=lambda p: float(p["scored"].get("overall_score", 0.0)), reverse=True)

    def fmt(x: Any) -> str:
        try:
            return f"{float(x):.1f}"
        except Exception:
            return "NA"

    def fmt_novelty_display(novelty_payload: Any) -> str:
        """
        Map novelty to 0-100 for display.

        Novelty scoring may apply a 10% combo premium, taking it up to 110.
        For display we map back to the base 0-100 range by removing the premium.
        """
        if novelty_payload is None:
            return "NA"

        if isinstance(novelty_payload, dict):
            raw = novelty_payload.get("score")
            has_premium = bool(novelty_payload.get("combo_premium"))
        else:
            raw = novelty_payload
            has_premium = False

        try:
            v = float(raw)
        except Exception:
            return "NA"

        if has_premium:
            v = v / 1.10

        return fmt(v)

    def novelty_display_value(novelty_payload: Any) -> float:
        s = fmt_novelty_display(novelty_payload)
        try:
            return float(s)
        except Exception:
            return 50.0

    def mechanism_phrase(spec: ComboSpec) -> str:
        mapping = {
            "list1_JAK_STAT": "JAK–STAT pathway",
            "list2_TNFR_BMP": "inflammation→fibrosis axis (TNFR–BMP)",
            "list3_TNFR_PCOLCE": "inflammation + ECM remodeling (TNFR–PCOLCE)",
            "list4_cell_cycle_integrin": "cell cycle + integrin trafficking",
            "list5_CDKN2D_PCOLCE": "cell cycle + ECM remodeling",
        }
        return mapping.get(spec.sc_key, "multi-pathway mechanism")

    ranking_rows = []
    combos_out = []
    sc_data: Dict[str, Dict[str, str]] = {}
    plot_figs: Dict[str, Dict[str, str]] = {}
    plot_targets: Dict[str, Dict[str, Dict[str, str]]] = {}

    def _float(x: Any, default: float = 50.0) -> float:
        try:
            return float(x)
        except Exception:
            return default

    # Radar plot (subscores for each combo).
    radar_theta = ["Clinical", "Disease", "Safety", "Opportunity", "Novelty"]
    radar_traces = []
    for payload in payloads:
        spec: ComboSpec = payload["spec"]
        scored: Dict[str, Any] = payload["scored"]
        subscores = scored.get("subscores", {})
        disease = subscores.get("disease_association", {})
        vals = [
            _float(subscores.get("clinical_validation", {}).get("normalized_score")),
            _float(disease.get("score")),
            _float(subscores.get("safety", {}).get("score")),
            _float(subscores.get("opportunity", {}).get("score")),
            _float(fmt_novelty_display(subscores.get("novelty"))),
        ]
        # close the shape
        radar_traces.append(
            {
                "type": "scatterpolar",
                "r": vals + [vals[0]],
                "theta": radar_theta + [radar_theta[0]],
                "fill": "toself",
                "name": spec.name,
                "opacity": 0.70,
            }
        )

    plot_figs["radar::rankings"] = {
        "data": json.dumps(radar_traces),
        "layout": json.dumps(
            {
                "polar": {"radialaxis": {"visible": True, "range": [0, 100]}},
                "showlegend": True,
                "legend": {"orientation": "h", "x": 0, "y": -0.12},
                "margin": {"l": 40, "r": 40, "t": 40, "b": 60},
                "paper_bgcolor": "rgba(0,0,0,0)",
                "plot_bgcolor": "rgba(0,0,0,0)",
            }
        ),
        "config": json.dumps({"displayModeBar": True, "displaylogo": False, "responsive": True}),
    }

    # cell type mapping in SC report
    cell_types = [
        ("mesenchymal", "Mesenchymal"),
        ("t_cell", "T Cell"),
        ("b_lineage", "B Lineage"),
        ("myeloid", "Myeloid"),
        ("endothelial", "Endothelial"),
    ]

    for rank, payload in enumerate(payloads, start=1):
        spec: ComboSpec = payload["spec"]
        scored: Dict[str, Any] = payload["scored"]

        subscores = scored.get("subscores", {})
        disease = subscores.get("disease_association", {})
        disease_components = disease.get("components", {})
        disease_components_norm = disease.get("components_normalized", {}) if isinstance(disease.get("components_normalized", {}), dict) else {}

        # Prefer component normalized values if present; fall back to overall.
        comp_bars = []
        weights = {
            "deg": 0.40,
            "biobridge": 0.25,
            "ultra": 0.25,
            "primekg": 0.10,
        }
        weight_sum = sum(weights.values()) or 1.0
        labels = {
            "deg": "DEG (expression evidence)",
            "biobridge": "BioBridge (combo percentile)",
            "ultra": "ULTRA (model percentile)",
            "primekg": "PrimeKG (connectivity)",
        }

        # Prefer the normalized component bundle produced by the scorer (best-effort and stable),
        # fall back to reconstructing from the structured component breakdown when possible.
        normalized_values: Dict[str, float] = {}
        for key in ["deg", "biobridge", "ultra", "primekg"]:
            if key in disease_components_norm:
                try:
                    normalized_values[key] = float(disease_components_norm[key])
                except Exception:
                    pass

        if "deg" not in normalized_values:
            # DEG average for combo
            if isinstance(disease_components.get("deg"), dict):
                vals = []
                for v in disease_components["deg"].values():
                    if isinstance(v, dict) and v.get("normalized") is not None:
                        vals.append(float(v["normalized"]))
                if vals:
                    normalized_values["deg"] = sum(vals) / len(vals)

        if "biobridge" not in normalized_values:
            if isinstance(disease_components.get("biobridge"), dict) and disease_components["biobridge"].get("normalized") is not None:
                normalized_values["biobridge"] = float(disease_components["biobridge"]["normalized"])

        if "ultra" not in normalized_values:
            if isinstance(disease_components.get("ultra"), dict):
                vals = []
                for v in disease_components["ultra"].values():
                    if isinstance(v, dict) and v.get("normalized") is not None:
                        vals.append(float(v["normalized"]))
                if vals:
                    normalized_values["ultra"] = sum(vals) / len(vals)

        if "primekg" not in normalized_values:
            if isinstance(disease_components.get("primekg"), dict):
                vals = []
                for v in disease_components["primekg"].values():
                    if isinstance(v, dict) and v.get("normalized") is not None:
                        vals.append(float(v["normalized"]))
                if vals:
                    normalized_values["primekg"] = sum(vals) / len(vals)

        for k in ["deg", "biobridge", "ultra", "primekg"]:
            v = normalized_values.get(k, 50.0)
            # Disease Association scoring normalizes weights to sum to 1.0; reflect that in the UI.
            w = weights[k] / weight_sum
            comp_bars.append({"key": k, "label": labels[k], "weight": f"{w:.2f}", "value": f"{v:.1f}"})

        # Plotly sources (extract fig payloads from exported Plotly HTML)
        bulk_variants = _bulk_variants(spec.list_num)
        if not bulk_variants:
            # Back-compat: explicit `bulk_base` paths (single view).
            bulk_heatmap = _find_first_existing(
                [
                    RESULTS_DIR / "bulk_coexpression" / "results" / f"{spec.bulk_base}_correlation_heatmap.html",
                    RESULTS_DIR / "bulk_coexpression" / f"{spec.bulk_base}_correlation_heatmap.html",
                ]
            )
            bulk_box = _find_first_existing(
                [
                    RESULTS_DIR / "bulk_coexpression" / "results" / f"{spec.bulk_base}_expression_boxplot.html",
                    RESULTS_DIR / "bulk_coexpression" / f"{spec.bulk_base}_expression_boxplot.html",
                ]
            )
            bulk_variants = [{"key": "default", "heatmap_path": bulk_heatmap, "boxplot_path": bulk_box}]

        def _label(ctx: str) -> str:
            return "Default" if ctx == "default" else ctx.replace("-", " ").replace("_", " ").title()

        def _default_key(keys: List[str]) -> str:
            prefs = ["lung", "colon", "ileum", "skin", "blood", "default"]
            for k in prefs:
                if k in keys:
                    return k
            return keys[0] if keys else "default"

        bulk_ctx_keys = [str(v.get("key") or "") for v in bulk_variants if v.get("key")]
        bulk_default_ctx = _default_key(bulk_ctx_keys)
        bulk_variant_targets = []
        bulk_heatmap_key = ""
        bulk_box_key = ""

        for v in bulk_variants:
            ctx = str(v.get("key") or "").strip() or "unknown"
            heatmap_path = v.get("heatmap_path")
            box_path = v.get("boxplot_path")
            heatmap_html = _read_text(heatmap_path) if heatmap_path else None
            box_html = _read_text(box_path) if box_path else None

            heatmap_fig_key = ""
            box_fig_key = ""

            if heatmap_html:
                fig_key = f"bulkHeatmap::{spec.anchor}::{ctx}"
                d, l, c = _extract_plotly_figure_args(heatmap_html)
                plot_figs[fig_key] = {"data": d, "layout": l, "config": c}
                heatmap_fig_key = fig_key

            if box_html:
                fig_key = f"bulkBox::{spec.anchor}::{ctx}"
                d, l, c = _extract_plotly_figure_args(box_html)
                plot_figs[fig_key] = {"data": d, "layout": l, "config": c}
                box_fig_key = fig_key

            bulk_variant_targets.append(
                {
                    "key": ctx,
                    "label": _label(ctx),
                    "heatmapFigKey": heatmap_fig_key,
                    "boxFigKey": box_fig_key,
                }
            )

            if ctx == bulk_default_ctx:
                bulk_heatmap_key = heatmap_fig_key
                bulk_box_key = box_fig_key

        # Ensure we have a default selection even if the preferred context is missing.
        if not bulk_heatmap_key:
            for v in bulk_variant_targets:
                if v.get("heatmapFigKey"):
                    bulk_heatmap_key = str(v["heatmapFigKey"])
                    break
        if not bulk_box_key:
            for v in bulk_variant_targets:
                if v.get("boxFigKey"):
                    bulk_box_key = str(v["boxFigKey"])
                    break

        # SC defaults
        sc_entries = []
        sc_default = None
        for ct_key, ct_heading in cell_types:
            interp_md, heatmap_rel = _parse_sc_celltype_interpretation(sc_md, ct_heading, spec.sc_key)
            expected_rel = f"{ct_key}/{ct_key}_{spec.sc_key}_heatmap.html"
            expected_abs = RESULTS_DIR / "sc_coexp" / "ibd_analysis_results" / expected_rel

            candidates: List[Path] = [expected_abs]

            # Prefer explicit markdown reference when it matches the list number (guard against stale links).
            if heatmap_rel and f"list{spec.list_num}_" in heatmap_rel:
                candidates.append(RESULTS_DIR / "sc_coexp" / "ibd_analysis_results" / heatmap_rel)

            # Fallback: any heatmap for this list number inside the cell-type folder.
            globbed = sorted(
                (RESULTS_DIR / "sc_coexp" / "ibd_analysis_results" / ct_key).glob(f"*list{spec.list_num}_*heatmap.html")
            )
            candidates.extend(globbed)

            interp_html = _md_to_html(interp_md) if interp_md else "<p>No interpretation found.</p>"

            # If multiple SC heatmap HTMLs exist for the same list + cell type, keep all and
            # expose them through a dropdown (template reads scData[*].variants).
            seen = set()
            existing: List[Path] = []
            for p in candidates:
                try:
                    rp = p.resolve()
                except Exception:
                    rp = p
                if not p.exists():
                    continue
                k = str(rp)
                if k in seen:
                    continue
                seen.add(k)
                existing.append(p)

            if not existing:
                continue

            default_path = expected_abs if expected_abs.exists() else existing[0]
            default_variant = _sc_variant_key(default_path, ct_key=ct_key, sc_key=spec.sc_key, list_num=spec.list_num)
            variants_out = []
            default_fig_key = ""

            for p in existing:
                vkey = _sc_variant_key(p, ct_key=ct_key, sc_key=spec.sc_key, list_num=spec.list_num)
                plot_html = p.read_text(encoding="utf-8", errors="replace")
                fig_key = f"sc::{spec.anchor}::{ct_key}::{vkey}"
                d, l, c = _extract_plotly_figure_args(plot_html)
                plot_figs[fig_key] = {"data": d, "layout": l, "config": c}

                if p == default_path:
                    default_fig_key = fig_key

                variants_out.append(
                    {
                        "key": vkey,
                        "label": _label(vkey),
                        "fig_key": fig_key,
                        "caption": f"{ct_heading} — {_label(vkey)} — List {spec.list_num}",
                    }
                )

            if not default_fig_key:
                default_fig_key = variants_out[0]["fig_key"]

            sc_entries.append({"key": ct_key, "label": ct_heading})
            if not sc_default:
                sc_default = {
                    "fig_key": default_fig_key,
                    "caption": f"{ct_heading} — {_label(default_variant)} — List {spec.list_num}",
                    "interpretation_html": interp_html,
                }

            sc_data[f"{spec.anchor}::{ct_key}"] = {
                "fig_key": default_fig_key,
                "caption": f"{ct_heading} — List {spec.list_num}",
                "interpretation_html": interp_html,
                "variants": variants_out,
                "default_variant_key": default_variant,
            }

        if sc_default is None:
            sc_default = {"fig_key": "", "caption": "Single-cell heatmap", "interpretation_html": "<p>No SC data.</p>"}

        # Interpretation blocks
        bulk_section = _extract_heading_block(bulk_md, "##", f"List {spec.list_num}:") or ""
        bulk_interp = _parse_bulk_list_interpretation(spec.list_num)["html"]

        # Disease interpretation: synthesize from per-gene markdown excerpts
        disease_bullets = []
        for gene in spec.genes:
            deg_interp = deg_interps.get(gene)
            if deg_interp:
                disease_bullets.append(f"- **DEG ({gene})**: {deg_interp}")
            bb_interp = _parse_biobridge_interpretation(gene)
            if bb_interp:
                disease_bullets.append(f"- **BioBridge ({gene})**: {bb_interp}")
            pk_excerpt = _parse_primekg_interpretation(gene)
            if pk_excerpt:
                disease_bullets.append(f"- **PrimeKG ({gene})**: {pk_excerpt}")
        disease_interpretation_html = _md_to_html("\n".join(disease_bullets) if disease_bullets else "No interpretation found.")

        # List-level pathway + PPI summaries.
        pathway_png = RESULTS_DIR / "pathwaydb" / f"{spec.pathway_base}.png"
        pathway_summary_html = _pathway_combo_summary_html(spec)
        ppi_html = _ppi_excerpt_html(spec)

        # Per-gene disease association subscores for the combo disease tab.
        disease_gene_rows = []
        for g in spec.genes:
            gs = gene_scored.get(g, {})
            gsub = gs.get("subscores", {}) if isinstance(gs, dict) else {}
            gd = gsub.get("disease_association", {}) if isinstance(gsub, dict) else {}
            comp_norm = gd.get("components_normalized", {}) if isinstance(gd.get("components_normalized", {}), dict) else {}
            disease_gene_rows.append(
                {
                    "gene": g,
                    "disease": fmt(gd.get("score")),
                    "deg": fmt(comp_norm.get("deg")),
                    "biobridge": fmt(comp_norm.get("biobridge")),
                    "ultra": fmt(comp_norm.get("ultra")),
                    "primekg": fmt(comp_norm.get("primekg")),
                }
            )

        ranking_rows.append(
            {
                "rank": rank,
                "anchor": spec.anchor,
                "name": spec.name,
                "genes": ", ".join(spec.genes),
                "overall": fmt(scored.get("overall_score")),
                "clinical": fmt(subscores.get("clinical_validation", {}).get("normalized_score")),
                "disease": fmt(disease.get("score")),
                "opportunity": fmt(subscores.get("opportunity", {}).get("score")),
                "novelty": fmt_novelty_display(subscores.get("novelty")),
                "safety": (
                    f"{fmt(subscores.get('safety', {}).get('score'))}*"
                    if subscores.get("safety", {}).get("breakdown") is None
                    and fmt(subscores.get("safety", {}).get("score")) == "50.0"
                    else fmt(subscores.get("safety", {}).get("score"))
                ),
            }
        )

        combos_out.append(
            {
                "name": spec.name,
                "anchor": spec.anchor,
                "genes": list(spec.genes),
                "long_genes": ", ".join(f"{g} — {GENE_TITLES.get(g, g)}" for g in spec.genes),
                "scores": {
                    "overall": fmt(scored.get("overall_score")),
                    "clinical": fmt(subscores.get("clinical_validation", {}).get("normalized_score")),
                    "disease": fmt(disease.get("score")),
                    "opportunity": fmt(subscores.get("opportunity", {}).get("score")),
                    "novelty": fmt_novelty_display(subscores.get("novelty")),
                    "safety": fmt(subscores.get("safety", {}).get("score")),
                },
                "overview_html": _md_to_html(
                    "\n".join(
                        [
                            f"**Blue Ocean category:** {scored.get('blue_ocean_category', 'Unknown')}",
                            "",
                            f"**Rationale:** {scored.get('blue_ocean_rationale', '')}",
                        ]
                    )
                ),
                "gene_rows": [
                    {
                        "gene": gene,
                        "overall": fmt(gene_scored.get(gene, {}).get("overall_score")),
                        "clinical": fmt(gene_scored.get(gene, {}).get("subscores", {}).get("clinical_validation", {}).get("normalized_score")),
                        "disease": fmt(gene_scored.get(gene, {}).get("subscores", {}).get("disease_association", {}).get("score")),
                        "opportunity": fmt(gene_scored.get(gene, {}).get("subscores", {}).get("opportunity", {}).get("score")),
                        "novelty": fmt_novelty_display(gene_scored.get(gene, {}).get("subscores", {}).get("novelty")),
                        "safety": fmt(gene_scored.get(gene, {}).get("subscores", {}).get("safety", {}).get("score")),
                    }
                    for gene in spec.genes
                ],
                "disease_components": comp_bars,
                "disease_interpretation_html": disease_interpretation_html,
                "disease_gene_rows": disease_gene_rows,
                "bulk_heatmap_fig_key": bulk_heatmap_key,
                "bulk_box_fig_key": bulk_box_key,
                "bulk_interpretation_html": bulk_interp,
                "sc_cell_types": sc_entries,
                "sc_default": sc_default,
                "pathway_png_path": _png_to_data_uri(pathway_png) if pathway_png.exists() else None,
                "pathway_summary_html": pathway_summary_html,
                "ppi_html": ppi_html,
            }
        )

        plot_targets[spec.anchor] = {
            "bulkHeatmap": {"containerId": f"bulk-heatmap-{spec.anchor}", "figKey": bulk_heatmap_key},
            "bulkBox": {"containerId": f"bulk-box-{spec.anchor}", "figKey": bulk_box_key},
            "bulkVariants": bulk_variant_targets,
            "bulkDefaultKey": bulk_default_ctx,
            "scDefault": {"containerId": f"sc-plot-{spec.anchor}", "figKey": sc_default.get("fig_key", "")},
        }

    # Single-gene score table (same scoring pipeline, shown after the last combo).
    gene_rows_all = []
    for gene, scored in gene_scored.items():
        subscores = scored.get("subscores", {}) if isinstance(scored, dict) else {}
        disease = subscores.get("disease_association", {}) if isinstance(subscores, dict) else {}
        gene_rows_all.append(
            {
                "gene": gene,
                "overall": fmt(scored.get("overall_score") if isinstance(scored, dict) else None),
                "clinical": fmt(subscores.get("clinical_validation", {}).get("normalized_score")),
                "disease": fmt(disease.get("score")),
                "opportunity": fmt(subscores.get("opportunity", {}).get("score")),
                "novelty": fmt_novelty_display(subscores.get("novelty")),
                "safety": fmt(subscores.get("safety", {}).get("score")),
            }
        )
    gene_rows_all.sort(key=lambda r: float(r["overall"]) if r["overall"] != "NA" else -1.0, reverse=True)

    # Gene cards (no coexpression tabs; just score + disease evidence + pathways).
    gene_cards = []
    for gene, scored in sorted(gene_scored.items(), key=lambda kv: float(kv[1].get("overall_score", 0.0)), reverse=True):
        subs = scored.get("subscores", {}) if isinstance(scored, dict) else {}
        disease_sub = subs.get("disease_association", {}) if isinstance(subs, dict) else {}
        comp_norm = disease_sub.get("components_normalized", {}) if isinstance(disease_sub.get("components_normalized", {}), dict) else {}

        # Build disease component bars for gene from normalized components and shared weights.
        weights = {"deg": 0.40, "biobridge": 0.25, "ultra": 0.25, "primekg": 0.10}
        weight_sum = sum(weights.values()) or 1.0
        labels = {
            "deg": "DEG (expression evidence)",
            "biobridge": "BioBridge (percentile)",
            "ultra": "ULTRA (model percentile)",
            "primekg": "PrimeKG (connectivity)",
        }
        gene_comp_bars = []
        for k in ["deg", "biobridge", "ultra", "primekg"]:
            v = comp_norm.get(k, 50.0)
            w = weights[k] / weight_sum
            gene_comp_bars.append({"key": k, "label": labels[k], "weight": f"{w:.2f}", "value": f"{float(v):.1f}"})

        # Reuse the same interpretation snippets used for combos.
        disease_bullets = []
        deg_interp = deg_interps.get(gene)
        if deg_interp:
            disease_bullets.append(f"- **DEG**: {deg_interp}")
        bb_interp = _parse_biobridge_interpretation(gene)
        if bb_interp:
            disease_bullets.append(f"- **BioBridge**: {bb_interp}")
        pk_excerpt = _parse_primekg_interpretation(gene)
        if pk_excerpt:
            disease_bullets.append(f"- **PrimeKG**: {pk_excerpt}")

        gene_cards.append(
            {
                "name": gene,
                "anchor": f"gene-{gene.lower()}",
                "long_name": GENE_TITLES.get(gene, gene),
                "scores": {
                    "overall": fmt(scored.get("overall_score")),
                    "clinical": fmt(subs.get("clinical_validation", {}).get("normalized_score")),
                    "disease": fmt(disease_sub.get("score")),
                    "opportunity": fmt(subs.get("opportunity", {}).get("score")),
                    "novelty": fmt_novelty_display(subs.get("novelty")),
                    "safety": fmt(subs.get("safety", {}).get("score")),
                },
                "overview_html": _md_to_html(
                    "\n".join(
                        [
                            f"**Blue Ocean category:** {scored.get('blue_ocean_category', 'Unknown')}",
                            "",
                            f"**Rationale:** {scored.get('blue_ocean_rationale', '')}",
                        ]
                    )
                ),
                "disease_components": gene_comp_bars,
                "disease_interpretation_html": _md_to_html("\n".join(disease_bullets) if disease_bullets else "No interpretation found."),
                "pathway_html": _gene_pathway_summary_html(gene),
            }
        )

    # Scoring methodology (end-of-report, tabbed in template)
    try:
        from target_prioritization.scoring.normalizers import MAX_CLINICAL_SCORE, MAX_DEG_SCORE
    except Exception:
        MAX_CLINICAL_SCORE, MAX_DEG_SCORE = 134.5, 15.0

    disease_base_weights = {"deg": 0.40, "biobridge": 0.25, "ultra": 0.25, "primekg": 0.10}
    disease_weight_sum = sum(disease_base_weights.values()) or 1.0
    disease_norm_weights = {k: (v / disease_weight_sum) for k, v in disease_base_weights.items()}

    methodology = {
        "overall_html": _md_to_html(
            "\n".join(
                [
                    "All subscores are intended to be on a **0–100** scale. Missing sources default to **50.0** (neutral).",
                    "",
                    "### Overall score (weighted sum)",
                    "- `Overall = 0.30·Clinical + 0.30·Disease + 0.10·Safety + 0.20·Opportunity + 0.10·Novelty`",
                    "",
                    "### Combo aggregation (how genes contribute)",
                    "- **Clinical:** average raw Cortellis scores across genes, then normalize to 0–100.",
                    "- **Disease:** DEG/ULTRA/PrimeKG are averaged across genes (on the normalized scale); BioBridge uses a combo percentile when available.",
                    "- **Safety:** average gene safety scores when no combo OFF‑X breakdown is available.",
                    "- **Opportunity (CI):** CI per gene is blended (0.5 gene-level, 0.5 family-level); combos use mean CI across member genes.",
                    "- **Novelty:** uses combo clinical score and average PrimeKG connectivity across genes; internal score includes a combo premium, but the report displays the base 0–100 scale.",
                ]
            )
        ),
        "clinical_html": _md_to_html(
            "\n".join(
                [
                    "### Data source",
                    "- `results/cortellis/IBD_Target_Analysis_Report.md`",
                    "",
                    "### Raw → 0–100 transform",
                    f"- `Clinical = clamp( (raw_total_score / {MAX_CLINICAL_SCORE}) · 100, 0, 100 )`",
                    "",
                    "### Combo rule",
                    "- For combinations, the raw Cortellis scores for each gene are averaged, then normalized using the same transform.",
                ]
            )
        ),
        "disease_html": _md_to_html(
            "\n".join(
                [
                    "### Components used (no pathway centrality)",
                    "- DEG (expression evidence)",
                    "- BioBridge percentile",
                    "- ULTRA model percentile",
                    "- PrimeKG connectivity (path score)",
                    "",
                    "### Component normalization (0–100)",
                    f"- **DEG:** `deg_norm = (deg_raw / {MAX_DEG_SCORE}) · 100` (clamped 0–100)",
                    "- **BioBridge:** percentile already 0–100",
                    "- **ULTRA:** percentile 0–1 converted to 0–100 when needed",
                    "- **PrimeKG:** `connectivity = clamp(primekg_score · 100, 0, 100)` (where `primekg_score` is the path-based connectivity score in [0,1])",
                    "",
                    "### Weighted combination (weights normalized to sum to 1.0)",
                    f"- DEG: `{disease_norm_weights['deg']:.2f}`",
                    f"- BioBridge: `{disease_norm_weights['biobridge']:.2f}`",
                    f"- ULTRA: `{disease_norm_weights['ultra']:.2f}`",
                    f"- PrimeKG: `{disease_norm_weights['primekg']:.2f}`",
                    "",
                    "### Combo rule",
                    "- DEG/ULTRA/PrimeKG are averaged across genes (after normalization).",
                    "- BioBridge uses a combo percentile when available; otherwise falls back to the first gene.",
                ]
            )
        ),
        "safety_html": _md_to_html(
            "\n".join(
                [
                    "### Data source",
                    "- `results/offx/OFF-X_Safety_Analysis_Report.md`",
                    "",
                    "### Severity scoring (risk-weighted average)",
                    "- Category base scores: `very_high=0`, `high=10`, `medium=20`, `low=60`, `very_low=80`, `not_assoc=100`",
                    "- **Amplification factors are all 1×** (equal weighting across classes).",
                    "- NA rows are excluded from the denominator (treated as unknown).",
                    "",
                    "Formula:",
                    "- `Safety = Σ(count[class]·score[class]) / Σ(count[class])`, clamped 0–100",
                    "",
                    "### Combo rule",
                    "- If no combo-level OFF‑X breakdown exists, the combo safety score is the average of the member gene safety scores.",
                    "",
                    "### Missing data",
                    "- If OFF‑X breakdown is missing, Safety defaults to **50.0*** in the Priority Rankings table.",
                ]
            )
        ),
        "opportunity_html": _md_to_html(
            "\n".join(
                [
                    "### Inputs",
                    "- Disease Association (0–100)",
                    "- Clinical Novelty = `100 - Clinical`",
                    "- Competitive Intelligence (CI) score (0–100, higher = more whitespace)",
                    "  - **Gene-level CI (phase-weighted):**",
                    "    - `Weighted_gene = 1.0·Marketed + 0.7·PhaseIII + 0.4·PhaseII + 0.2·PhaseI + 0.1·Preclinical`",
                    "    - `CI_gene = 100 - 100·(Weighted_gene / TotalPrograms_gene)`",
                    "  - **Family-level CI (crowding):**",
                    "    - Each CI entry has `targetFamilyPrimary` (e.g., JAK, WNT, IL-23).",
                    "    - `FamilyWeighted = Σ(weight(entry.ibdPhase))` over all IBD-tagged entries in that family.",
                    "    - `AllWeighted = Σ(weight(entry.ibdPhase))` over all IBD-tagged entries across families.",
                    "    - `CI_family = 100 - 100·(FamilyWeighted / AllWeighted)`",
                    "  - **Blended CI:** `CI = 0.5·CI_gene + 0.5·CI_family`",
                    "  - **Combinations:** `CI_combo = mean(CI_gene)` across member genes (no extra scaling).",
                    "  - If CI data is missing/empty, CI defaults to `50.0`.",
                    "- For combos: a synergy term influences the opportunity blend (via NovelMechanismBonus)",
                    "",
                    "### Individual gene formula",
                    "- `Opportunity = 0.40·Disease + 0.30·ClinicalNovelty + 0.30·CI`",
                    "",
                    "### Combo formula (as implemented)",
                    "- `MeanOpp = 0.40·Disease + 0.30·ClinicalNovelty + 0.30·CI`",
                    "- `NovelMechanismBonus = 100 - 100·Synergy`",
                    "- `Opportunity_combo = 0.60·MeanOpp + 0.40·NovelMechanismBonus`",
                    "",
                    "Notes:",
                    "- Synergy is derived from pathway overlap metrics (shared pathways normalized by a reference max).",
                ]
            )
        ),
        "novelty_html": _md_to_html(
            "\n".join(
                [
                    "### Inputs",
                    "- Clinical Novelty = `100 - Clinical`",
                    "- Literature Novelty (PrimeKG): higher connectivity → lower novelty",
                    "",
                    "### Individual gene formula",
                    "- `Connectivity = clamp(primekg_score · 100, 0, 100)`",
                    "- `LiteratureNovelty = 100 - Connectivity`",
                    "- `Novelty = 0.70·ClinicalNovelty + 0.30·LiteratureNovelty`",
                    "",
                    "### Combo rule",
                    "- PrimeKG connectivity scores are averaged across genes before computing LiteratureNovelty.",
                    "- Internally, combos apply a **+10% premium** to Novelty; the report displays the base 0–100 equivalent by removing this premium.",
                ]
            )
        ),
    }

    # Executive summary: top picks and what changed.
    top3 = combos_out[:3] if combos_out else []
    exec_md_lines = [
        "This report prioritizes **IBD target combinations** using the redesigned loader+scoring pipeline.",
        "",
        "## How combinations are scored",
        "- **Clinical Validation (30%)**: Cortellis clinical/asset evidence (normalized to 0–100).",
        "- **Disease Association (30%)**: weighted evidence from DEG, BioBridge, ULTRA, and PrimeKG (0–100).",
        "- **Opportunity (20%)**: whitespace signal blending disease strength, clinical novelty, and competitive intensity (0–100).",
        "- **Safety (10%)**: OFF‑X adverse event severity mix (0–100; missing defaults to 50.0*).",
        "- **Novelty (10%)**: inverse clinical + literature novelty proxy (0–100 display).",
        "",
    ]
    if top3:
        exec_md_lines += ["## Current 3 top-ranked combinations"]
        for row in top3:
            exec_md_lines.append(
                f"- **{row['name']}** (Overall **{row['scores']['overall']}**) with Disease **{row['scores']['disease']}** and Safety **{row['scores']['safety']}**."
            )
        exec_md_lines.append("")
    # Strategic portfolio recommendations (data-driven from this run).
    exec_md_lines += ["## Strategic Portfolio Recommendation"]
    if payloads:
        scored_by_anchor = {p["spec"].anchor: p for p in payloads}
        # Leader: highest overall (already sorted)
        leader = payloads[0]

        def overall(p):  # type: ignore[no-redef]
            try:
                return float(p["scored"].get("overall_score", 50.0))
            except Exception:
                return 50.0

        def sub(p, key, default=50.0):
            try:
                v = p["scored"]["subscores"].get(key, {}).get("score" if key != "clinical_validation" else "normalized_score")
            except Exception:
                v = None
            try:
                return float(v)
            except Exception:
                return default

        def disease_score(p):
            try:
                return float(p["scored"]["subscores"].get("disease_association", {}).get("score", 50.0))
            except Exception:
                return 50.0

        def safety_score(p):
            try:
                return float(p["scored"]["subscores"].get("safety", {}).get("score", 50.0))
            except Exception:
                return 50.0

        def opportunity_score(p):
            try:
                return float(p["scored"]["subscores"].get("opportunity", {}).get("score", 50.0))
            except Exception:
                return 50.0

        def novelty_score(p):
            try:
                return novelty_display_value(p["scored"]["subscores"].get("novelty"))
            except Exception:
                return 50.0

        def line(label: str, p, extra: str):
            spec = p["spec"]
            exec_md_lines.append(
                f"- **{label}: {spec.name}** (Overall: **{overall(p):.1f}**) — {extra}"
            )

        # High-confidence leader (top overall)
        leader_spec = leader["spec"]
        leader_reason = f"{mechanism_phrase(leader_spec)} with strong clinical precedent (Clinical: {sub(leader, 'clinical_validation'):.1f}, Disease: {disease_score(leader):.1f})."
        line("High-Confidence Leader", leader, leader_reason)

        remaining = payloads[1:]

        # Balanced opportunity: best combined overall + opportunity among remaining.
        if remaining:
            balanced = max(remaining, key=lambda p: (overall(p) + opportunity_score(p)) / 2.0)
            remaining = [p for p in remaining if p is not balanced]
            b_spec = balanced["spec"]
            b_reason = f"{mechanism_phrase(b_spec)} with favorable balance of opportunity (Opportunity: {opportunity_score(balanced):.1f}) and execution confidence (Safety: {safety_score(balanced):.1f})."
            line("Balanced Opportunity", balanced, b_reason)

        # Novel high-risk: highest novelty among remaining.
        if remaining:
            novel = max(remaining, key=lambda p: novelty_score(p))
            remaining = [p for p in remaining if p is not novel]
            n_spec = novel["spec"]
            n_reason = f"maximal novelty (Novelty: {novelty_score(novel):.1f}) via {mechanism_phrase(n_spec)}; de-risk with targeted validation (Safety: {safety_score(novel):.1f})."
            line("Novel High-Risk", novel, n_reason)

        # Exploratory: list the rest, highlight opportunity.
        if remaining:
            remaining_sorted = sorted(remaining, key=lambda p: opportunity_score(p), reverse=True)
            top_remaining = remaining_sorted[:2]
            names = " & ".join([f"**{p['spec'].name}** (Overall: **{overall(p):.1f}**)" for p in top_remaining])
            opps = ", ".join([f"{opportunity_score(p):.1f}" for p in top_remaining])
            exec_md_lines.append(
                f"- **Exploratory:** {names} — highest opportunity scores among remaining ({opps})."
            )
    else:
        exec_md_lines.append("- (No combinations scored in this run.)")
    executive_html = _md_to_html("\n".join(exec_md_lines))

    template = Template(TEMPLATE_PATH.read_text(encoding="utf-8"))

    # Extract Plotly bundle once from any available Plotly-export HTML.
    plotly_html = ""
    for spec in COMBOS:
        for v in _bulk_variants(spec.list_num):
            p = v.get("heatmap_path") or v.get("boxplot_path")
            plotly_html = _read_text(p) or ""
            if plotly_html:
                break
        if plotly_html:
            break

    if not plotly_html:
        # Back-compat: the baseline path used by earlier versions of this script.
        plotly_src = RESULTS_DIR / "bulk_coexpression" / "results" / f"{COMBOS[0].bulk_base}_correlation_heatmap.html"
        plotly_html = _read_text(plotly_src) or ""

    if not plotly_html:
        for p in (RESULTS_DIR / "sc_coexp" / "ibd_analysis_results").rglob("*.html"):
            plotly_html = _read_text(p) or ""
            if plotly_html:
                break

    plotly_bundle = _extract_plotly_bundle_js(plotly_html)
    # NOTE: this is raw JS embedded into the report's <script> tag (not a JS string),
    # so we must emit a real newline (not a literal "\\n" token).
    plotly_bundle_js = "window.PlotlyConfig = {MathJaxConfig: 'local'};\n" + plotly_bundle

    ci_dashboard_path = RESULTS_DIR / "ci" / "ibd_dashboard.html"
    ci_dashboard_b64 = _ci_dashboard_b64(ci_dashboard_path, "Xinghao Zhang")

    html_out = template.render(
        title="IBD Target Combination Prioritization Report",
        nav_title="IBD Target Combination Prioritization",
        subtitle="Interactive dossier with traceable narrative excerpts and embedded analysis.",
        build_id=f"{datetime.now().strftime('%Y-%m-%d')} / combo-report-v2",
        generation_date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        builders="Xinghao Zhang & Roger Tu",
        executive_html=executive_html,
        combos=combos_out,
        ranking_rows=ranking_rows,
        gene_rows_all=gene_rows_all,
        gene_cards=gene_cards,
        methodology=methodology,
        ci_dashboard_available=bool(ci_dashboard_b64),
        ci_dashboard_b64_json=json.dumps(ci_dashboard_b64 or ""),
        plotly_bundle_js=plotly_bundle_js,
        plot_figs_js=_js_obj_literal(plot_figs),
        plot_targets_json=json.dumps(plot_targets),
        sc_data_json=json.dumps(sc_data),
    )

    output_path = REPORTS_DIR / "IBD_Combo_Prioritization_Report_offline.html"
    output_path.write_text(html_out, encoding="utf-8")
    print(f"Wrote: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
