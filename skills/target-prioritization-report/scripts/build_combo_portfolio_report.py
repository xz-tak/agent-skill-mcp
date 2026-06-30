#!/usr/bin/env python3
from __future__ import annotations

import argparse
import glob
import json
import math
import re
from base64 import b64encode
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

try:
    from bs4 import BeautifulSoup  # type: ignore
except Exception:  # pragma: no cover
    BeautifulSoup = None


DEFAULT_GLOB = "results_*/reports/*_Combo_Prioritization_Report_offline.html"
DEFAULT_ORDER = ["IBD", "SSc", "ATD", "HS"]
DEFAULT_LABELS = {
    "SSc": "Systemic Sclerosis",
    "ATD": "Atopic Dermatitis",
    "HS": "Hidradenitis Suppurativa",
}
DEFAULT_BUILT_BY = "Xinghao Zhang"


@dataclass(frozen=True)
class ReportData:
    key: str
    title: str
    path: Path
    build_meta: str | None
    score_strategy_lines: list[str]
    combos: list[dict[str, Any]]
    genes: list[str]


def _safe_text(node: Any) -> str:
    if node is None:
        return ""
    return node.get_text(" ", strip=True)


def _slug(s: str) -> str:
    s = s.strip().lower()
    s = re.sub(r"[^a-z0-9]+", "-", s)
    return s.strip("-") or "x"


def _extract_key(path: Path, soup: Any | None) -> str:
    m = re.match(r"^([A-Za-z0-9]+)_Combo_Prioritization_Report_offline\.html$", path.name)
    if m:
        return m.group(1)
    if soup is not None and getattr(soup, "title", None) and getattr(soup.title, "string", None):
        title = (soup.title.string or "").strip()
        if title:
            return re.sub(r"[^A-Za-z0-9]+", "", title.split()[0]) or path.stem
    return path.stem


def _extract_build_meta(soup: Any) -> str | None:
    meta = _safe_text(soup.select_one(".brand .meta"))
    return meta or None


def _extract_genes(soup: Any) -> list[str]:
    genes: list[str] = []

    for a in soup.select("aside.nav a.navlink"):
        href = (a.get("href") or "").strip()
        if href.startswith("#gene-"):
            label = _safe_text(a)
            if label and label.upper() == label and label not in genes:
                genes.append(label)

    table = soup.find("table", class_="rank-table")
    if table and table.tbody:
        for tr in table.tbody.find_all("tr"):
            tds = tr.find_all("td")
            if len(tds) >= 3:
                cell = _safe_text(tds[2])
                for g in [x.strip() for x in cell.split(",") if x.strip()]:
                    if g.upper() == g and g not in genes:
                        genes.append(g)

    return sorted(genes)


def _extract_combos(soup: Any) -> list[dict[str, Any]]:
    table = soup.find("table", class_="rank-table")
    if not table or not table.tbody:
        return []

    combos: list[dict[str, Any]] = []
    for tr in table.tbody.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) < 9:
            continue

        def num(i: int) -> float | None:
            t = _safe_text(tds[i]).replace("*", "")
            try:
                return float(t)
            except Exception:
                return None

        combo_name = _safe_text(tds[1])
        genes = [g.strip() for g in _safe_text(tds[2]).split(",") if g.strip()]
        combos.append(
            {
                "rank": int(num(0) or (len(combos) + 1)),
                "name": combo_name,
                "genes": genes,
                "overall": num(3),
                "clinical": num(4),
                "disease": num(5),
                "opportunity": num(6),
                "novelty": num(7),
                "safety": num(8),
            }
        )
    return combos


def _extract_score_strategy_lines(soup: Any) -> list[str]:
    li = soup.find("li", string=re.compile(r"Clinical Validation\s*\(30%\)", re.I))
    if li and li.parent and li.parent.name == "ul":
        out = [_safe_text(x) for x in li.parent.find_all("li")]
        return [x for x in out if x]

    for ul in soup.find_all("ul"):
        lines = [_safe_text(x) for x in ul.find_all("li")]
        joined = " | ".join(lines)
        if re.search(r"\b30%\b", joined) and re.search(r"\b20%\b", joined) and re.search(r"\b10%\b", joined):
            return [x for x in lines if x]

    return []


def _parse_report(path: Path) -> ReportData:
    if BeautifulSoup is None:
        raise SystemExit(
            "Missing dependency: beautifulsoup4. Install it (e.g. `pip install beautifulsoup4`) "
            "or run this stitch step in an environment that already has it."
        )

    soup = BeautifulSoup(path.read_text(encoding="utf-8"), "html.parser")
    key = _extract_key(path, soup)
    title = (soup.title.string if soup.title else key).strip() or key
    return ReportData(
        key=key,
        title=title,
        path=path,
        build_meta=_extract_build_meta(soup),
        score_strategy_lines=_extract_score_strategy_lines(soup),
        combos=_extract_combos(soup),
        genes=_extract_genes(soup),
    )


def _intersection(sets: Iterable[set[str]]) -> set[str]:
    sets = list(sets)
    if not sets:
        return set()
    out = set(sets[0])
    for s in sets[1:]:
        out &= set(s)
    return out


def _format_build_timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")


def _parse_weight_pct(line: str) -> float | None:
    m = re.search(r"\((\d+(?:\.\d+)?)%\)", line)
    if not m:
        return None
    return float(m.group(1))


def _normalize_weights(lines: list[str]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for ln in lines:
        label = re.sub(r"\s*\(.*?\)\s*", "", ln)
        label = re.sub(r"^[-•]\s*", "", label).strip(": ").strip()
        out.append({"label": label, "pct": _parse_weight_pct(ln), "raw": ln})
    return out


def _svg_pie_slices(values: list[float]) -> list[dict[str, float]]:
    total = sum(values) or 1.0
    start = -90.0
    out: list[dict[str, float]] = []
    for v in values:
        sweep = 360.0 * (v / total)
        out.append({"start": start, "sweep": sweep})
        start += sweep
    return out


def _arc_path(start_deg: float, sweep_deg: float, r: float = 46.0) -> str:
    start = math.radians(start_deg)
    end = math.radians(start_deg + sweep_deg)
    x1, y1 = r * math.cos(start), r * math.sin(start)
    x2, y2 = r * math.cos(end), r * math.sin(end)
    large = 1 if sweep_deg > 180 else 0
    return f"M 0 0 L {x1:.4f} {y1:.4f} A {r} {r} 0 {large} 1 {x2:.4f} {y2:.4f} Z"


def _report_sort_key(key: str, order: list[str]) -> tuple[int, str]:
    try:
        idx = order.index(key)
    except ValueError:
        idx = 999
    return (idx, key.lower())


def _rewrite_built_by(html: str, name: str) -> str:
    html = re.sub(
        r'(<div class="label"[^>]*>\s*Built by\s*</div>\s*<div class="value"[^>]*>)([^<]*)(</div>)',
        lambda m: m.group(1) + name + m.group(3),
        html,
        flags=re.I,
    )
    html = re.sub(r"(·\s*Built by\s*)([^<·]+)", r"\1" + name, html, flags=re.I)
    html = re.sub(r"(Built by\s*)(tools/[A-Za-z0-9_./-]+)", r"\1" + name, html, flags=re.I)
    return html


def _rewrite_theme_colors(html: str) -> str:
    html = re.sub(r"--nav:\s*#[0-9a-fA-F]{6};", "--nav: #E74C3C;", html)
    html = re.sub(r"--nav2:\s*#[0-9a-fA-F]{6};", "--nav2:#E74C3C;", html)
    html = re.sub(r"--accent:\s*#[0-9a-fA-F]{6};", "--accent: #E74C3C;", html)
    return html


def _rewrite_titles(html: str, label: str) -> str:
    title = f"{label} Combination Prioritization Report"
    html = re.sub(r"<title>.*?</title>", f"<title>{title}</title>", html, count=1, flags=re.I | re.S)
    html = re.sub(r"<h1>.*?</h1>", f"<h1>{title}</h1>", html, count=1, flags=re.I | re.S)
    return html


def _rewrite_heatmap_scale(html: str) -> str:
    score_hue = r"""
      function scoreHue(score){
        let t = Math.max(0, Math.min(100, score)) / 100;
        if (t <= 0.5) {
          const tt = t / 0.5;
          return 0 + (60 * tt); // red -> yellow
        }
        const tt = (t - 0.5) / 0.5;
        return 60 + (60 * tt); // yellow -> green
      }
    """.strip()

    colorize = r"""
      function colorizeTableByColumn(table){
        const tbody = table.tBodies && table.tBodies[0];
        if (!tbody) return;
        const rows = Array.from(tbody.rows);
        if (!rows.length) return;
        const headers = table.tHead && table.tHead.rows && table.tHead.rows[0]
          ? Array.from(table.tHead.rows[0].cells).map(c => (c.textContent || '').trim().toLowerCase())
          : [];

        const colCount = Math.max(...rows.map(r => r.cells.length));
        for (let col = 0; col < colCount; col++){
          const h = headers[col] || '';
          if (h === 'rank' || h === 'combination' || h === 'genes' || h === 'gene') continue;

          for (const r of rows){
            const cell = r.cells[col];
            if (!cell) continue;
            const v = parseFloat((cell.textContent || '').trim());
            if (!Number.isFinite(v)) continue;

            const hue = scoreHue(v);
            const sat = 88;
            const lightTop = 86;
            const lightBot = 74;
            cell.classList.add('scorecell');
            cell.style.background = `linear-gradient(180deg, hsla(${hue}, ${sat}%, ${lightTop}%, .92), hsla(${hue}, ${sat}%, ${lightBot}%, .78))`;
            cell.style.boxShadow = 'inset 0 0 0 1px rgba(12,19,35,.14), 0 10px 18px rgba(10,16,30,.08)';
            cell.style.borderRadius = '12px';
            cell.style.backgroundClip = 'padding-box';
          }
        }
      }
    """.strip()

    html, _ = re.subn(
        r"function colorizeTableByColumn\(table\)\{.*?\n\s*\}\n\s*\n\s*function colorizeAllScoreTables\(\)\{",
        colorize + "\n\n      function colorizeAllScoreTables(){",
        html,
        count=1,
        flags=re.S,
    )
    html, _ = re.subn(
        r"function scoreHue\(score\)\{.*?\n\s*\}",
        score_hue,
        html,
        count=1,
        flags=re.S,
    )
    return html


def _normalize_report_title(label: str) -> str:
    return f"{label} Combination Prioritization Report"


def _normalize_build_meta(build_meta: str | None, name: str) -> str | None:
    if not build_meta:
        return build_meta
    return re.sub(r"(·\s*Built by\s*)([^<·]+)", r"\1" + name, build_meta, flags=re.I)


def _rewrite_executive_summary(html: str, report: ReportData) -> str:
    if BeautifulSoup is None:
        return html

    m = re.search(r"<section\s+id=\"executive\"[^>]*>.*?</section>", html, flags=re.S | re.I)
    if not m:
        return html

    section_html = m.group(0)
    soup = BeautifulSoup(section_html, "html.parser")
    sec = soup.find("section", attrs={"id": re.compile(r"^executive$", re.I)})
    if not sec:
        return html

    md = sec.find("div", class_="md")
    if not md:
        return html

    def find_block(head_patterns: list[str]) -> tuple[Any | None, Any | None]:
        head = None
        for tag in md.find_all(["h2", "h3", "h4"]):
            t = tag.get_text(" ", strip=True)
            if any(re.search(p, t, flags=re.I) for p in head_patterns):
                head = tag
                break
        if not head:
            return None, None
        ul = None
        for sib in head.next_siblings:
            if getattr(sib, "name", None) in {"h2", "h3", "h4"}:
                break
            if getattr(sib, "name", None) == "ul":
                ul = sib
                break
        return head, ul

    _h_w, ul_w = find_block([r"scoring\s+weights", r"default\s+scoring", r"how\s+combinations\s+are\s+scored"])
    if ul_w is None:
        lines = report.score_strategy_lines or [
            "Clinical Validation (30%)",
            "Disease Association (30%) = DEG (0.40), BioBridge (0.25), ULTRA (0.25), PrimeKG (0.10)",
            "Opportunity (20%)",
            "Safety (10%) (missing → 50.0)",
            "Novelty (10%)",
        ]
        ul_w = soup.new_tag("ul")
        for ln in lines:
            li = soup.new_tag("li")
            li.string = ln
            ul_w.append(li)

    _h_t, ul_t = find_block([r"top\s+combinations", r"top\s+combos", r"top\s+ranked\s+combinations", r"top\s+recommendation"])
    if ul_t is None:
        ul_t = soup.new_tag("ul")
        top = sorted(report.combos, key=lambda c: c.get("rank") or 9999)[:3]
        for c in top:
            name = c.get("name") or "—"
            overall = c.get("overall")
            li = soup.new_tag("li")
            strong = soup.new_tag("strong")
            strong.string = f"#{c.get('rank') or ''} {name}".strip()
            li.append(strong)
            if overall is not None:
                li.append(" — Overall ")
                s2 = soup.new_tag("strong")
                s2.string = f"{overall:.1f}"
                li.append(s2)
            ul_t.append(li)

    md.clear()
    h3a = soup.new_tag("h3")
    h3a.string = "1. Scoring weights"
    md.append(h3a)
    md.append(ul_w)

    h3b = soup.new_tag("h3")
    h3b.string = "2. Top combinations (by Overall score)"
    md.append(h3b)
    md.append(ul_t)

    return html[: m.start()] + str(sec) + html[m.end() :]


def _patch_flat_heatmap_columns(html: str) -> str:
    """
    Embedded reports apply per-column heatmap fills to score tables.
    If a column is constant (max == min), remove the fill to avoid misleading color.
    """

    script = r"""
<script>
(() => {
  function clearCell(cell){
    if (!cell || !cell.style) return;
    cell.style.background = 'transparent';
    cell.style.boxShadow = 'none';
    cell.style.borderRadius = '0';
    cell.style.backgroundClip = 'border-box';
  }

  function fixTable(table){
    const tbody = table.tBodies && table.tBodies[0];
    if (!tbody) return;
    const rows = Array.from(tbody.rows || []);
    if (!rows.length) return;
    const headers = (table.tHead && table.tHead.rows && table.tHead.rows[0])
      ? Array.from(table.tHead.rows[0].cells).map(c => (c.textContent || '').trim().toLowerCase())
      : [];
    const colCount = Math.max(0, ...rows.map(r => r.cells ? r.cells.length : 0));
    for (let col = 0; col < colCount; col++){
      const h = headers[col] || '';
      if (h === 'rank' || h === 'combination' || h === 'genes' || h === 'gene') continue;
      const vals = [];
      for (const r of rows){
        const cell = r.cells[col];
        if (!cell) continue;
        const v = parseFloat((cell.textContent || '').trim());
        if (Number.isFinite(v)) vals.push(v);
      }
      if (!vals.length) continue;
      const min = Math.min(...vals);
      const max = Math.max(...vals);
      if (Math.abs(max - min) < 1e-9){
        for (const r of rows){
          const cell = r.cells[col];
          if (cell) clearCell(cell);
        }
      }
    }
  }

  function run(){
    document.querySelectorAll('table.score-color-table').forEach(fixTable);
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', run, { once: true });
  else run();
})();
</script>
""".strip()

    m = re.search(r"</body\\s*>", html, flags=re.I)
    if m:
        return html[: m.start()] + script + "\\n" + html[m.start() :]
    return html + "\\n" + script + "\\n"


def _glob_paths(root: Path, pattern: str) -> list[Path]:
    if re.match(r"^[A-Za-z]:[\\/]", pattern) or pattern.startswith(("/", "\\")):
        return sorted({Path(p).resolve() for p in glob.glob(pattern, recursive=True)})
    return sorted({p.resolve() for p in root.glob(pattern)})


def build(
    *,
    root: Path,
    output: Path,
    report_glob: str,
    built_by: str,
    order: list[str],
    labels: dict[str, str],
) -> None:
    paths = _glob_paths(root, report_glob)
    if not paths:
        raise SystemExit(f"No reports found for glob: {report_glob}")

    reports = sorted((_parse_report(p) for p in paths), key=lambda r: _report_sort_key(r.key, order))

    embedded_b64: dict[str, str] = {}
    for r in reports:
        raw = r.path.read_bytes()
        try:
            txt = raw.decode("utf-8")
        except UnicodeDecodeError:
            txt = raw.decode("utf-8", errors="replace")
        label = labels.get(r.key, r.key)
        txt = _rewrite_titles(txt, label)
        txt = _rewrite_theme_colors(txt)
        txt = _rewrite_built_by(txt, built_by)
        txt = _rewrite_executive_summary(txt, r)
        txt = _rewrite_heatmap_scale(txt)
        txt = _patch_flat_heatmap_columns(txt)
        embedded_b64[r.key] = b64encode(txt.encode("utf-8")).decode("ascii")

    genes_by_ind = {r.key: set(r.genes) for r in reports}
    combos_by_ind = {r.key: set(c["name"] for c in r.combos if c.get("name")) for r in reports}

    all_genes = set().union(*genes_by_ind.values()) if genes_by_ind else set()
    shared_all = _intersection(list(genes_by_ind.values()))

    unique_genes: dict[str, list[str]] = {}
    for k, gset in genes_by_ind.items():
        others = set().union(*[v for kk, v in genes_by_ind.items() if kk != k])
        unique_genes[k] = sorted(gset - others)

    gene_membership_counts: dict[str, int] = {g: 0 for g in all_genes}
    for g in all_genes:
        gene_membership_counts[g] = sum(1 for s in genes_by_ind.values() if g in s)

    top_genes = sorted(gene_membership_counts.items(), key=lambda kv: (-kv[1], kv[0]))[:20]
    shared_2plus = sorted(
        (g for g, n in gene_membership_counts.items() if n >= 2),
        key=lambda g: (-gene_membership_counts[g], g),
    )

    shared_combos_all = _intersection(list(combos_by_ind.values()))
    unique_combos: dict[str, int] = {}
    for k, cset in combos_by_ind.items():
        others = set().union(*[v for kk, v in combos_by_ind.items() if kk != k])
        unique_combos[k] = len(cset - others)

    strategy_src = next((r for r in reports if r.score_strategy_lines), reports[0])
    strategy_lines = strategy_src.score_strategy_lines or [
        "Clinical Validation (30%)",
        "Disease Association (30%) = DEG (0.40), BioBridge (0.25), ULTRA (0.25), PrimeKG (0.10)",
        "Opportunity (20%)",
        "Safety (10%) (missing → 50.0)",
        "Novelty (10%)",
    ]
    strategy_items = _normalize_weights(strategy_lines)

    pie_labels = [labels.get(r.key, r.key) for r in reports]
    pie_values = [float(len(r.combos)) for r in reports]
    pie_slices = _svg_pie_slices(pie_values)

    summary_data: dict[str, Any] = {
        "built_utc": _format_build_timestamp(),
        "built_by": built_by,
        "reports": [
            {
                "key": r.key,
                "label": labels.get(r.key, r.key),
                "title": _normalize_report_title(labels.get(r.key, r.key)),
                "path": str(r.path),
                "build_meta": _normalize_build_meta(r.build_meta, built_by),
                "n_combos": len(r.combos),
                "n_genes": len(set(r.genes)),
                "top_combo": (r.combos[0] if r.combos else None),
            }
            for r in reports
        ],
        "genes": {
            "total": len(all_genes),
            "shared_all": sorted(shared_all),
            "shared_2plus": [{"gene": g, "n": gene_membership_counts[g]} for g in shared_2plus],
            "unique": unique_genes,
            "by_indication": {k: sorted(v) for k, v in genes_by_ind.items()},
            "top_by_indications": [{"gene": g, "n": n} for g, n in top_genes],
            "distribution": {
                "by_n_indications": {
                    str(i): sum(1 for _g, n in gene_membership_counts.items() if n == i)
                    for i in range(1, len(reports) + 1)
                }
            },
        },
        "combos": {
            "total": sum(len(v) for v in combos_by_ind.values()),
            "shared_all": sorted(shared_combos_all),
            "unique_counts": unique_combos,
            "by_indication": {k: sorted(v) for k, v in combos_by_ind.items()},
        },
        "combo_scores_by_indication": {
            r.key: [
                {"rank": c.get("rank"), "name": c.get("name"), "overall": c.get("overall")}
                for c in r.combos
                if c.get("name") and (c.get("overall") is not None)
            ]
            for r in reports
        },
        "combo_genes_by_name": {
            name: sorted(
                {
                    g
                    for r in reports
                    for c in r.combos
                    if c.get("name") == name
                    for g in (c.get("genes") or [])
                }
            )
            for name in sorted({c.get("name") for r in reports for c in r.combos if c.get("name")})
        },
        "score_strategy": {"source": strategy_src.key, "lines": strategy_lines, "items": strategy_items},
        "pie": {"labels": pie_labels, "values": pie_values, "slices": pie_slices},
    }

    output.write_text(_render_html(summary_data, embedded_b64), encoding="utf-8")


def _render_html(data: dict[str, Any], embedded_reports_b64: dict[str, str]) -> str:
    data_json = json.dumps(data, ensure_ascii=False)
    embedded_json = json.dumps(embedded_reports_b64, ensure_ascii=False)

    tabs = [{"id": "summary", "label": "Summary"}] + [
        {"id": _slug(r["key"]), "label": r.get("label") or r["key"]} for r in data["reports"]
    ]

    # Keep palette consistent with the GI2 per-indication report theme.
    palette = ["#E74C3C", "#1b998b", "#f2b705", "#0b1b3a", "#b21e31", "#137a6d"]

    pie_paths = "".join(
        (
            lambda i, sl, label, value: (
                f'<path class="pie-slice" d="{_arc_path(sl["start"], sl["sweep"])}" '
                f'fill="{palette[i % len(palette)]}" opacity="0.92" '
                f'data-label="{label}" data-value="{int(value)}">'
                f"<title>{label}: {int(value)} lists</title>"
                f"</path>"
            )
        )(i, sl, data["pie"]["labels"][i], data["pie"]["values"][i])
        for i, sl in enumerate(data["pie"]["slices"])
    )
    pie_legend = "".join(
        f"""
                <div class="legend-item">
                  <span class="swatch" style="background:{palette[i % len(palette)]}"></span>
                  <div class="legend-text">
                    <div class="legend-label">{label}</div>
                    <div class="legend-sub">{int(value)} lists</div>
                  </div>
                </div>
        """
        for i, (label, value) in enumerate(zip(data["pie"]["labels"], data["pie"]["values"]))
    )

    report_summary_rows = "".join(
        (
            lambda r, tc: f"""
            <tr>
              <td>{(r.get('label') or r['key'])}</td>
              <td>{r['n_genes']}</td>
              <td>{r['n_combos']}</td>
              <td class="mono">{(tc.get('name') or '—')}</td>
              <td class="mono">{('%.1f' % tc.get('overall')) if tc.get('overall') is not None else '—'}</td>
            </tr>
            """
        )(r, r.get("top_combo") or {})
        for r in data["reports"]
    )

    shared_gene_items = "".join(f'<span class="chip">{g}</span>' for g in data["genes"]["shared_all"][:80])
    if len(data["genes"]["shared_all"]) > 80:
        shared_gene_items += f'<span class="chip subtle">+{len(data["genes"]["shared_all"]) - 80} more</span>'
    shared_gene_block = shared_gene_items or '<span class="faint">None.</span>'

    shared_2plus_items = "".join(
        f'<span class="chip">{x["gene"]}<span class="subtle" style="margin-left:8px;">{x["n"]}/{len(data["reports"])}</span></span>'
        for x in data["genes"]["shared_2plus"][:60]
    )
    if len(data["genes"]["shared_2plus"]) > 60:
        shared_2plus_items += f'<span class="chip subtle">+{len(data["genes"]["shared_2plus"]) - 60} more</span>'
    shared_2plus_block = shared_2plus_items or '<span class="faint">None.</span>'

    uniq_blocks = []
    for r in data["reports"]:
        k = r["key"]
        glist = data["genes"]["unique"].get(k, [])
        chips = "".join(f'<span class="chip">{g}</span>' for g in glist[:36])
        more = f'<span class="chip subtle">+{len(glist) - 36} more</span>' if len(glist) > 36 else ""
        uniq_blocks.append(
            f"""
            <div class="mini-card">
              <div class="mini-head">
                <div class="mini-title">{(r.get('label') or k)} unique genes</div>
                <div class="mini-count">{len(glist)}</div>
              </div>
              <div class="chipwrap">{chips}{more}</div>
            </div>
            """
        )
    uniq_blocks_html = "".join(uniq_blocks)

    ind_cols: list[tuple[str, str]] = [(r["key"], (r.get("label") or r["key"])) for r in data["reports"]]
    bar_options = "".join(
        f'<option value="{k}"{" selected" if i == 0 else ""}>{lbl}</option>'
        for i, (k, lbl) in enumerate(ind_cols)
    )

    genes_by_ind: dict[str, list[str]] = data["genes"].get("by_indication", {})
    lists_by_ind: dict[str, list[str]] = data["combos"].get("by_indication", {})
    gene_sets = {k: set(v) for k, v in genes_by_ind.items()}
    all_gene_names = sorted({g for glist in genes_by_ind.values() for g in glist})
    all_list_names = sorted({name for llist in lists_by_ind.values() for name in llist})
    combo_genes_by_name: dict[str, list[str]] = data.get("combo_genes_by_name", {})

    def mark(present: bool) -> str:
        return '<span class="mark ok" aria-label="yes">✓</span>' if present else '<span class="mark no" aria-label="no">✕</span>'

    gene_matrix_rows = []
    for g in all_gene_names:
        tds = [f'<td class="mono">{g}</td>']
        for key, _lbl in ind_cols:
            tds.append(f"<td class='center'>{mark(g in gene_sets.get(key, set()))}</td>")
        gene_matrix_rows.append("<tr>" + "".join(tds) + "</tr>")
    gene_matrix_table = f"""
      <table class="matrix">
        <thead><tr><th>Gene</th>{''.join(f'<th>{lbl}</th>' for _k, lbl in ind_cols)}</tr></thead>
        <tbody>{''.join(gene_matrix_rows) or ''}</tbody>
      </table>
    """

    list_matrix_rows = []
    for name in all_list_names:
        tds = [f'<td class="mono">{name}</td>']
        for key, _lbl in ind_cols:
            genes = combo_genes_by_name.get(name, [])
            ind_genes = gene_sets.get(key, set())
            ok = bool(genes) and all(g in ind_genes for g in genes)
            missing = [g for g in genes if g not in ind_genes]
            title = ""
            if genes:
                title = f" title='genes: {', '.join(genes)}" + (f" | missing: {', '.join(missing)}'" if missing else "'")
            tds.append(f"<td class='center'{title}>{mark(ok)}</td>")
        list_matrix_rows.append("<tr>" + "".join(tds) + "</tr>")
    list_matrix_table = f"""
      <table class="matrix">
        <thead><tr><th>List</th>{''.join(f'<th>{lbl}</th>' for _k, lbl in ind_cols)}</tr></thead>
        <tbody>{''.join(list_matrix_rows) or ''}</tbody>
      </table>
    """

    dist_rows = []
    for i in range(1, len(data["reports"]) + 1):
        n = data["genes"]["distribution"]["by_n_indications"][str(i)]
        width = 0.0 if data["genes"]["total"] == 0 else (n / data["genes"]["total"]) * 100.0
        dist_rows.append(
            f"""
            <div class="distrow">
              <div class="distlabel">{i} indication{'s' if i != 1 else ''}</div>
              <div class="distbar"><div class="distfill" style="width:{width:.2f}%"></div></div>
              <div class="distval">{n}</div>
            </div>
            """
        )
    dist_rows_html = "".join(dist_rows)

    weights = [x for x in data["score_strategy"]["items"] if x.get("pct") is not None]
    total_pct = sum(float(x["pct"]) for x in weights) or 100.0
    weight_blocks = []
    for i, w in enumerate(weights):
        pct = float(w["pct"])
        width = 100.0 * (pct / total_pct)
        color = palette[i % len(palette)]
        weight_blocks.append(
            f"""
            <div class="wrow">
              <div class="wlabel">{w['label']}</div>
              <div class="wbar"><div class="wfill" style="width:{width:.2f}%; background:{color}"></div></div>
              <div class="wval">{pct:.0f}%</div>
            </div>
            """
        )
    weight_blocks_html = "".join(weight_blocks) or '<div class="faint">No explicit weights found.</div>'

    disease_raw = next(
        (x["raw"] for x in data["score_strategy"]["items"] if "DEG" in x["raw"] and "weights" in x["raw"]),
        "",
    )
    disease_seg = ""
    m = re.search(
        r"weights\s*([0-9.]+)\s*/\s*([0-9.]+)\s*/\s*([0-9.]+)\s*/\s*([0-9.]+)",
        disease_raw,
    )
    if m:
        names = ["DEG", "BioBridge", "ULTRA", "PrimeKG"]
        vals = [float(m.group(i)) for i in range(1, 5)]
        s = sum(vals) or 1.0
        acc = 0.0
        segs = []
        for name, v in zip(names, vals):
            w = 100.0 * (v / s)
            color = {"DEG": "#d7263d", "BioBridge": "#1b998b", "ULTRA": "#f2b705", "PrimeKG": "#0b1b3a"}[name]
            segs.append(
                f'<div class="seg" title="{name} {v:.2f}" style="left:{acc:.4f}%; width:{w:.4f}%; background:{color}"></div>'
            )
            acc += w
        disease_seg = f"""
          <div class="subweights">
            <div class="subhead">Disease component blend</div>
            <div class="segbar" aria-label="Disease component weights">{''.join(segs)}</div>
            <div class="sublegend">
              <span class="tag" style="--c:#d7263d">DEG</span>
              <span class="tag" style="--c:#1b998b">BioBridge</span>
              <span class="tag" style="--c:#f2b705">ULTRA</span>
              <span class="tag" style="--c:#0b1b3a">PrimeKG</span>
            </div>
          </div>
        """

    combos_shared = data["combos"]["shared_all"]
    combo_shared_block = (
        "".join(f'<span class="chip">{c}</span>' for c in combos_shared[:24])
        if combos_shared
        else '<span class="faint">None shared across all indications.</span>'
    )

    unique_combo_rows = "".join(
        f"<tr><td>{(r.get('label') or r['key'])}</td><td>{r['n_combos']}</td><td>{data['combos']['unique_counts'].get(r['key'], 0)}</td></tr>"
        for r in data["reports"]
    )

    report_panels = "".join(
        f"""
      <section class="panel" data-panel="{_slug(r['key'])}">
        <div class="panel-head">
          <div>
            <div class="panel-kicker">Indication</div>
            <h2 class="panel-title">{_normalize_report_title((r.get('label') or r['key']))}</h2>
            <div class="panel-sub">Built by {data.get('built_by','')}</div>
          </div>
          <div class="panel-actions">
            <div class="metric"><div class="metric-k">Lists</div><div class="metric-v">{r['n_combos']}</div></div>
            <div class="metric"><div class="metric-k">Genes</div><div class="metric-v">{r['n_genes']}</div></div>
          </div>
        </div>
        <div class="iframe-wrap">
          <iframe class="reportframe" title="{r['key']} report" data-embed="{r['key']}" loading="lazy"></iframe>
          <div class="frame-status" data-state="loading" aria-live="polite">
            <div class="sbox">
              <div class="stitle">Loading {(r.get('label') or r['key'])} report</div>
              <div class="smsg" data-msg>Waiting to start…</div>
              <div class="serr">Open DevTools Console for details.</div>
            </div>
          </div>
        </div>
      </section>
        """
        for r in data["reports"]
    )

    return f"""<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <title>GI2 Combo Portfolio — Summary + Tabs (Offline)</title>
    <style>
      :root{{
        --paper: #f6f3ee;
        --ink: #0c1323;
        --muted: rgba(12,19,35,.70);
        --faint: rgba(12,19,35,.52);
        --stroke: rgba(12,19,35,.12);
        --stroke2: rgba(12,19,35,.18);
        --shadow: 0 18px 60px rgba(10,16,30,.10);
        --shadow2: 0 8px 22px rgba(10,16,30,.10);

        --nav: #E74C3C;
        --nav2:#E74C3C;
        --accent: #E74C3C;
        --accent2:#1b998b;
        --warn:#f2b705;

        --radius: 18px;
        --radius2: 12px;
      }}
      *{{ box-sizing:border-box; }}
      html,body{{ height:100%; }}
      body{{
        margin:0;
        color:var(--ink);
        background:
          radial-gradient(1200px 900px at 16% 8%, rgba(215,38,61,.10), transparent 55%),
          radial-gradient(900px 700px at 84% 22%, rgba(27,153,139,.12), transparent 55%),
          linear-gradient(180deg, #fff, var(--paper));
        font-family: ui-sans-serif, system-ui, -apple-system, "Segoe UI", Arial, sans-serif;
        font-size: 16px;
        line-height:1.5;
      }}

      .topbar{{
        position:sticky; top:0; z-index:40;
        background: linear-gradient(180deg, rgba(231,76,60,.98), rgba(231,76,60,.94));
        border-bottom: 1px solid rgba(255,255,255,.10);
        box-shadow: 0 18px 44px rgba(0,0,0,.20);
        backdrop-filter: blur(10px);
      }}
      .topbar-inner{{
        max-width: 1760px;
        margin: 0 auto;
        padding: 16px 18px 14px;
        display:grid;
        grid-template-columns: 1fr auto;
        gap: 14px;
        align-items:end;
      }}
      .brand{{
        color: rgba(255,255,255,.92);
        padding: 14px 14px 12px 14px;
        border: 1px solid rgba(255,255,255,.12);
        border-radius: 16px;
        background: rgba(255,255,255,.06);
        box-shadow: 0 14px 30px rgba(0,0,0,.18);
      }}
      .brand .kicker{{
        font-size: 11px;
        text-transform: uppercase;
        letter-spacing: .16em;
        opacity:.78;
      }}
      .brand .title{{
        margin-top: 6px;
        font-family: ui-serif, "Iowan Old Style", "Palatino Linotype", Palatino, serif;
        font-weight: 760;
        letter-spacing: .02em;
        font-size: 22px;
        line-height: 1.1;
      }}
      .brand .meta{{ margin-top: 8px; font-size: 12px; opacity: .72; }}

      .tabs{{ display:flex; gap: 10px; flex-wrap: wrap; justify-content:flex-end; }}
      .tab{{
        appearance:none;
        border: 1px solid rgba(255,255,255,.16);
        background: rgba(255,255,255,.06);
        color: rgba(255,255,255,.92);
        padding: 10px 12px;
        border-radius: 999px;
        font-weight: 780;
        letter-spacing: .01em;
        cursor:pointer;
        transition: transform .16s ease, background .16s ease, border-color .16s ease;
      }}
      .tab:hover{{ transform: translateY(-1px); background: rgba(255,255,255,.10); border-color: rgba(255,255,255,.20); }}
      .tab[data-active="true"]{{
        background: rgba(215,38,61,.16);
        border-color: rgba(215,38,61,.40);
        box-shadow: inset 0 0 0 1px rgba(215,38,61,.20);
      }}

      .wrap{{ max-width: 1760px; margin: 0 auto; padding: 20px 18px 84px; }}
      .panel{{ display:none; }}
      .panel[data-active="true"]{{ display:block; }}

      .summary-stack{{ display:grid; grid-template-columns: 1fr; gap: 16px; }}
      .summary-row{{ display:grid; grid-template-columns: 1fr 1fr; gap: 16px; align-items:start; }}
      .stack{{ display:grid; grid-template-columns: 1fr; gap: 16px; }}
      .hero-split{{ display:grid; grid-template-columns: 1.15fr .85fr; gap: 14px; align-items:start; margin-top: 14px; }}
      .matrix-split{{ display:grid; grid-template-columns: 1fr 1fr; gap: 14px; align-items:start; margin-top: 14px; }}
      .block-head{{ display:flex; justify-content:space-between; align-items:center; gap: 10px; margin-bottom: 10px; }}
      .block-title{{ font-weight: 900; color: var(--muted); }}
      .select{{
        appearance:none;
        border: 1px solid var(--stroke2);
        background: rgba(255,255,255,.78);
        border-radius: 999px;
        padding: 8px 12px;
        font-weight: 850;
        color: var(--ink);
        box-shadow: inset 0 0 0 1px rgba(12,19,35,.05);
      }}
      .chartbox{{
        border-radius: 16px;
        border: 1px solid var(--stroke);
        background: rgba(255,255,255,.62);
        box-shadow: var(--shadow2);
        padding: 12px 12px;
      }}
      .chartbox svg{{ width:100%; height: 320px; display:block; }}
      .chartnote{{ margin-top: 10px; font-size: 12px; color: var(--faint); }}

      .table-scroll{{
        border-radius: 16px;
        border: 1px solid var(--stroke);
        background: rgba(255,255,255,.62);
        box-shadow: var(--shadow2);
        overflow:auto;
        max-height: 380px;
      }}
      table.matrix{{ width:100%; }}
      table.matrix thead th{{ position: sticky; top: 0; background: rgba(246,243,238,.92); backdrop-filter: blur(6px); z-index: 2; }}
      .center{{ text-align:center; }}
      .mark{{
        display:inline-flex;
        width: 22px;
        height: 22px;
        align-items:center;
        justify-content:center;
        border-radius: 8px;
        font-weight: 950;
        line-height: 1;
        user-select:none;
      }}
      .mark.ok{{ background: rgba(27,153,139,.16); color: #137a6d; border: 1px solid rgba(27,153,139,.35); }}
      .mark.no{{ background: rgba(215,38,61,.14); color: #b21e31; border: 1px solid rgba(215,38,61,.30); }}

      @media (max-width: 1080px){{
        .topbar-inner{{ grid-template-columns: 1fr; }}
        .tabs{{ justify-content:flex-start; }}
        .summary-row{{ grid-template-columns: 1fr; }}
        .hero-split{{ grid-template-columns: 1fr; }}
        .matrix-split{{ grid-template-columns: 1fr; }}
      }}

      .card{{
        background: rgba(255,255,255,.70);
        border: 1px solid var(--stroke);
        border-radius: var(--radius);
        box-shadow: var(--shadow2);
        padding: 18px 18px;
        backdrop-filter: blur(8px);
      }}
      .card h2{{ margin: 0 0 10px; font-family: ui-serif, "Iowan Old Style", "Palatino Linotype", Palatino, serif; font-size: 22px; letter-spacing: .01em; }}
      .subtle{{ color: var(--faint); }}
      .faint{{ color: var(--faint); }}
      .mono{{ font-variant-numeric: tabular-nums; font-feature-settings: "tnum" 1; }}

      .kpi-row{{ display:grid; grid-template-columns: repeat(4, 1fr); gap: 12px; margin-top: 10px; }}
      @media (max-width: 980px){{ .kpi-row{{ grid-template-columns: 1fr 1fr; }} }}
      .kpi{{ border-radius: 16px; border: 1px solid var(--stroke); background: rgba(255,255,255,.72); padding: 12px 12px; box-shadow: inset 0 0 0 1px rgba(12,19,35,.06); }}
      .kpi .k{{ font-size: 11px; text-transform: uppercase; letter-spacing: .16em; color: var(--faint); font-weight: 850; }}
      .kpi .v{{ margin-top: 6px; font-size: 26px; font-weight: 850; }}

      table{{ width:100%; border-collapse:separate; border-spacing:0; }}
      th,td{{ padding: 10px 10px; border-bottom: 1px solid var(--stroke); vertical-align: top; }}
      thead th{{ text-align:left; font-size: 11px; text-transform: uppercase; letter-spacing: .16em; color: var(--faint); font-weight: 850; border-bottom: 1px solid var(--stroke2); }}
      .chipwrap{{ display:flex; flex-wrap:wrap; gap: 8px; margin-top: 10px; }}
      .chip{{ display:inline-flex; align-items:center; border: 1px solid var(--stroke); background: rgba(255,255,255,.78); border-radius: 999px; padding: 6px 10px; font-weight: 820; font-size: 13px; box-shadow: inset 0 0 0 1px rgba(12,19,35,.05); }}
      .chip.subtle{{ opacity:.75; font-weight: 780; }}

      .mini-grid{{ display:grid; grid-template-columns: 1fr 1fr; gap: 12px; margin-top: 12px; }}
      @media (max-width: 980px){{ .mini-grid{{ grid-template-columns: 1fr; }} }}
      .mini-card{{ border-radius: 16px; border: 1px solid var(--stroke); background: rgba(255,255,255,.64); padding: 12px 12px; }}
      .mini-head{{ display:flex; justify-content:space-between; align-items:baseline; gap: 10px; }}
      .mini-title{{ font-weight: 860; }}
      .mini-count{{ font-weight: 900; color: var(--accent); }}

      .viz{{ display:grid; grid-template-columns: 180px 1fr; gap: 14px; align-items:center; margin-top: 10px; }}
      @media (max-width: 560px){{ .viz{{ grid-template-columns: 1fr; }} }}
      .legend{{ display:grid; gap: 10px; }}
      .legend-item{{ display:flex; align-items:center; gap: 10px; }}
      .swatch{{ width: 12px; height: 12px; border-radius: 99px; box-shadow: 0 0 0 3px rgba(12,19,35,.08); }}
      .legend-label{{ font-weight: 860; }}
      .legend-sub{{ font-size: 12px; color: var(--faint); }}

      .tooltip{{
        position: fixed;
        z-index: 60;
        pointer-events: none;
        padding: 10px 10px;
        border-radius: 14px;
        border: 1px solid rgba(12,19,35,.16);
        background: rgba(255,255,255,.92);
        box-shadow: 0 18px 44px rgba(0,0,0,.18);
        backdrop-filter: blur(10px);
        color: var(--ink);
        display:none;
        min-width: 160px;
      }}
      .tooltip .t1{{ font-weight: 950; }}
      .tooltip .t2{{ margin-top: 4px; color: var(--muted); font-weight: 850; }}

      .distrow{{ display:grid; grid-template-columns: 120px 1fr 36px; gap: 10px; align-items:center; margin-top: 8px; }}
      .distlabel{{ font-weight: 820; color: var(--muted); font-size: 13px; }}
      .distbar{{ height: 10px; border-radius: 999px; background: rgba(12,19,35,.10); overflow:hidden; border: 1px solid rgba(12,19,35,.10); }}
      .distfill{{ height: 100%; background: linear-gradient(90deg, rgba(215,38,61,.95), rgba(27,153,139,.90)); }}
      .distval{{ text-align:right; font-weight: 900; }}

      .wrow{{ display:grid; grid-template-columns: 1fr 1.6fr 56px; gap: 10px; align-items:center; margin-top: 10px; }}
      .wlabel{{ font-weight: 860; }}
      .wbar{{ height: 12px; border-radius: 999px; background: rgba(12,19,35,.10); overflow:hidden; border: 1px solid rgba(12,19,35,.10); }}
      .wfill{{ height: 100%; border-radius: 999px; }}
      .wval{{ text-align:right; font-weight: 900; }}

      .subweights{{ margin-top: 14px; }}
      .subhead{{ font-size: 11px; text-transform: uppercase; letter-spacing: .16em; color: var(--faint); font-weight: 850; }}
      .segbar{{ position:relative; height: 12px; border-radius: 999px; background: rgba(12,19,35,.10); overflow:hidden; border: 1px solid rgba(12,19,35,.10); margin-top: 8px; }}
      .seg{{ position:absolute; top:0; bottom:0; }}
      .sublegend{{ display:flex; flex-wrap:wrap; gap: 8px; margin-top: 10px; }}
      .tag{{ display:inline-flex; align-items:center; gap: 8px; padding: 6px 10px; border-radius: 999px; border: 1px solid var(--stroke); background: rgba(255,255,255,.70); font-weight: 820; font-size: 12px; }}
      .tag::before{{ content:""; width: 9px; height: 9px; border-radius: 99px; background: var(--c); box-shadow: 0 0 0 3px rgba(12,19,35,.08); }}

      .panel-head{{ display:flex; justify-content:space-between; align-items:flex-end; gap: 16px; margin: 8px 0 14px; }}
      .panel-kicker{{ font-size: 11px; text-transform: uppercase; letter-spacing: .16em; color: var(--faint); font-weight: 850; }}
      .panel-title{{ margin: 6px 0 0; font-family: ui-serif, "Iowan Old Style", "Palatino Linotype", Palatino, serif; font-size: 30px; line-height: 1.08; }}
      .panel-sub{{ margin-top: 8px; color: var(--muted); max-width: 92ch; }}
      .panel-metrics{{ display:flex; gap: 10px; flex-wrap:wrap; justify-content:flex-end; align-items:end; }}
      .panel-actions{{ display:flex; gap: 10px; flex-wrap:wrap; justify-content:flex-end; align-items:end; }}
      .metric{{ padding: 10px 12px; border-radius: 16px; border: 1px solid var(--stroke); background: rgba(255,255,255,.62); box-shadow: var(--shadow2); }}
      .metric-k{{ font-size: 11px; text-transform: uppercase; letter-spacing: .16em; color: var(--faint); font-weight: 850; }}
      .metric-v{{ margin-top: 6px; font-size: 20px; font-weight: 900; }}

      .iframe-wrap{{ border-radius: var(--radius); overflow:hidden; border: 1px solid var(--stroke); box-shadow: var(--shadow); background: rgba(255,255,255,.62); position:relative; }}
      .frame-status{{
        position:absolute;
        inset: 0;
        display:flex;
        align-items:center;
        justify-content:center;
        padding: 18px;
        text-align:center;
        background: rgba(255,255,255,.78);
        backdrop-filter: blur(10px);
        border-radius: var(--radius);
      }}
      .frame-status[data-state="hidden"]{{ display:none; }}
      .frame-status .sbox{{
        max-width: 64ch;
        border: 1px solid var(--stroke);
        background: rgba(255,255,255,.86);
        border-radius: 16px;
        padding: 14px 14px;
        box-shadow: var(--shadow2);
      }}
      .frame-status .stitle{{ font-weight: 900; }}
      .frame-status .smsg{{ margin-top: 8px; color: var(--muted); }}
      .frame-status .serr{{ margin-top: 10px; color: rgba(215,38,61,.92); font-weight: 860; display:none; }}
      .frame-status[data-state="error"] .serr{{ display:block; }}
      /* Default: auto-size iframes to content height (page-scroll, not iframe-scroll). */
      .reportframe{{ width:100%; height: 960px; min-height: 720px; border: 0; background: transparent; }}
    </style>
  </head>
  <body>
    <div class="topbar">
      <div class="topbar-inner">
        <div class="brand">
          <div class="kicker">Offline Portfolio</div>
          <div class="title">GI2 Combo Prioritization — Summary + Indication Tabs</div>
          <div class="meta">{data['built_utc']} · Built by {data.get('built_by','')} · {len(data['reports'])} indications</div>
        </div>
        <div class="tabs" role="tablist" aria-label="Report tabs">
          {''.join(f'<button class="tab" type="button" role="tab" data-tab="{t["id"]}">{t["label"]}</button>' for t in tabs)}
        </div>
      </div>
    </div>

    <div class="wrap">
      <section class="panel" data-panel="summary">
        <div class="panel-head">
          <div>
            <div class="panel-kicker">Portfolio Summary</div>
            <h2 class="panel-title">Cross‑Indication Overview</h2>
            <div class="panel-sub">High-level summary across indications, with lightweight visuals and set overlap views. Built by {data.get('built_by','')}.</div>
          </div>
        </div>

        <div class="summary-stack">
          <div class="card">
            <h2>Executive Summary</h2>
            <div class="subtle">Cross-indication roll-up with score visualization and coverage matrices.</div>

            <div class="kpi-row">
              <div class="kpi"><div class="k">Indications</div><div class="v">{len(data['reports'])}</div></div>
              <div class="kpi"><div class="k">Total Genes</div><div class="v">{data['genes']['total']}</div></div>
              <div class="kpi"><div class="k">Total Lists</div><div class="v">{sum(r['n_combos'] for r in data['reports'])}</div></div>
              <div class="kpi"><div class="k">Shared ≥2</div><div class="v">{len(data['genes']['shared_2plus'])}</div></div>
            </div>

            <div class="hero-split">
              <div>
                <div class="block-head">
                  <div class="block-title">Indication overview</div>
                  <div class="subtle">Top list and top overall from each report’s rankings table.</div>
                </div>
                <div class="table-scroll">
                  <table>
                    <thead><tr><th>Indication</th><th>Genes</th><th>Lists</th><th>Top List</th><th>Top Overall</th></tr></thead>
                    <tbody>
                      {report_summary_rows}
                    </tbody>
                  </table>
                </div>
              </div>

              <div>
                <div class="block-head">
                  <div class="block-title">Overall scores (ranked lists)</div>
                  <select class="select" id="barSelect" aria-label="Select indication">
                    {bar_options}
                  </select>
                </div>
                <div class="chartbox">
                  <svg id="barSvg" viewBox="0 0 720 320" preserveAspectRatio="none" aria-label="Overall scores bar chart"></svg>
                  <div class="chartnote">Bars are sorted by rank; color is single red.</div>
                </div>
              </div>
            </div>

            <div class="matrix-split">
              <div>
                <div class="block-head">
                  <div class="block-title">Single-gene targets by indication</div>
                  <div class="subtle">Green ✓ indicates the gene appears in that indication report.</div>
                </div>
                <div class="table-scroll">
                  {gene_matrix_table}
                </div>
              </div>
              <div>
                <div class="block-head">
                  <div class="block-title">Combination gene coverage by indication</div>
                  <div class="subtle">Green ✓ means all genes for that list are present in the indication’s gene set.</div>
                </div>
                <div class="table-scroll">
                  {list_matrix_table}
                </div>
              </div>
            </div>
          </div>

          <div class="summary-row">
            <div class="card">
              <h2>Shared / Unique</h2>
              <div class="subtle">Shared and unique gene/list signals across the provided indications.</div>

              <div style="margin-top: 14px;">
                <div class="block-title">Genes shared across ≥2 indications</div>
                <div class="chipwrap">{shared_2plus_block}</div>
              </div>

              <div style="margin-top: 14px;">
                <div class="block-title">Genes shared across all indications</div>
                <div class="chipwrap">{shared_gene_block}</div>
              </div>

              <div style="margin-top: 14px;">
                <div class="block-title">Unique genes by indication</div>
                <div class="mini-grid">{uniq_blocks_html}</div>
              </div>

              <div style="margin-top: 14px;">
                <div class="block-title">Gene overlap distribution</div>
                <div class="subtle">Counts of genes present in 1..N indications.</div>
                {dist_rows_html}
              </div>

              <div style="margin-top: 14px;">
                <div class="block-title">Combination names shared across all indications</div>
                <div class="chipwrap">{combo_shared_block}</div>
              </div>

              <div style="margin-top: 14px;">
                <div class="block-title">Unique combination counts</div>
                <div class="table-scroll" style="max-height: 240px;">
                  <table>
                    <thead><tr><th>Indication</th><th>Lists</th><th>Unique Lists</th></tr></thead>
                    <tbody>{unique_combo_rows}</tbody>
                  </table>
                </div>
              </div>
            </div>

            <div class="stack">
              <div class="card">
                <h2>Indication Distribution</h2>
                <div class="subtle"># lists per indication (hover to inspect).</div>
                <div class="viz">
                  <div class="chartbox" style="padding: 14px;">
                    <svg viewBox="-50 -50 100 100" role="img" aria-label="Pie chart">
                      {pie_paths}
                    </svg>
                  </div>
                  <div class="legend">{pie_legend}</div>
                </div>
              </div>

              <div class="card">
                <h2>Overall Score Strategy</h2>
                <div class="subtle">Weights are extracted from an input report when possible; fallback is canonical defaults.</div>
                {weight_blocks_html}
                {disease_seg}
                <details style="margin-top: 14px;">
                  <summary style="font-weight: 900; cursor: pointer;">Raw strategy text</summary>
                  <div class="subtle" style="margin-top: 10px;">{data['score_strategy']['lines']}</div>
                </details>
              </div>
            </div>
          </div>
        </div>
      </section>

      {report_panels}
    </div>

    <div class="tooltip" id="tooltip"><div class="t1"></div><div class="t2"></div></div>

    <script id="portfolio-data" type="application/json">{data_json}</script>
    <script id="embedded-reports" type="application/json">{embedded_json}</script>
    <script>
      const DATA = JSON.parse(document.getElementById('portfolio-data').textContent);
      const EMBED = JSON.parse(document.getElementById('embedded-reports').textContent);

      function qs(sel, root=document) {{ return root.querySelector(sel); }}
      function qsa(sel, root=document) {{ return Array.from(root.querySelectorAll(sel)); }}

      function setActiveTab(tabId) {{
        qsa('.tab').forEach(b => b.dataset.active = (b.dataset.tab === tabId) ? 'true' : 'false');
        qsa('.panel').forEach(p => p.dataset.active = (p.dataset.panel === tabId) ? 'true' : 'false');
        if (tabId !== 'summary') {{
          const panel = qs(`.panel[data-panel="${{tabId}}"]`);
          if (panel) loadPanelIframe(panel);
        }}
      }}

      function b64ToBytes(b64, onProgress) {{
        // Chunked decode to avoid atob() limits on large strings.
        const chunkSize = 256 * 1024; // divisible by 4
        const padding = b64.endsWith('==') ? 2 : (b64.endsWith('=') ? 1 : 0);
        const outLen = Math.floor((b64.length * 3) / 4) - padding;
        const out = new Uint8Array(outLen);
        let outPos = 0;
        for (let offset = 0; offset < b64.length; offset += chunkSize) {{
          const slice = b64.slice(offset, offset + chunkSize);
          const bin = atob(slice);
          for (let i = 0; i < bin.length; i++) out[outPos++] = bin.charCodeAt(i);
          if (onProgress) onProgress(Math.min(1, (offset + chunkSize) / b64.length));
        }}
        return out;
      }}

      function writeDocFromBytes(doc, bytes) {{
        const td = new TextDecoder('utf-8');
        const step = 1024 * 1024;
        doc.open();
        for (let i = 0; i < bytes.length; i += step) {{
          const part = bytes.subarray(i, i + step);
          doc.write(td.decode(part, {{ stream: true }}));
        }}
        doc.write(td.decode());
        doc.close();
      }}

      function measureIframeHeight(iframe) {{
        try {{
          const doc = iframe.contentDocument;
          if (!doc) return 0;
          const de = doc.documentElement;
          const b = doc.body;
          const h = Math.max(
            de ? de.scrollHeight : 0,
            de ? de.offsetHeight : 0,
            b ? b.scrollHeight : 0,
            b ? b.offsetHeight : 0
          );
          return h || 0;
        }} catch (_e) {{
          return 0;
        }}
      }}

      function startAutoHeight(iframe) {{
        if (!iframe) return;
        let ticks = 0;
        let best = 0;
        const set = () => {{
          const h = measureIframeHeight(iframe);
          if (h > best) best = h;
          if (best > 0) {{
            iframe.style.height = (best + 24) + 'px';
          }}
        }};
        set();
        const id = setInterval(() => {{
          ticks += 1;
          set();
          if (ticks >= 40) clearInterval(id);
        }}, 250);
        iframe.dataset.autoh = String(id);
      }}

      function loadPanelIframe(panel) {{
        const iframe = qs('iframe.reportframe', panel);
        const status = qs('.frame-status', panel);
        if (!iframe || !status) return;
        if (iframe.dataset.loaded === 'true') {{
          status.dataset.state = 'hidden';
          return;
        }}
        const key = iframe.dataset.embed;
        const b64 = EMBED[key];
        if (!b64) {{
          status.dataset.state = 'error';
          qs('[data-msg]', status).textContent = 'Missing embedded report payload.';
          return;
        }}
        const msg = qs('[data-msg]', status);
        if (msg) msg.textContent = 'Decoding embedded report…';
        try {{
          const bytes = b64ToBytes(b64, (p) => {{
            if (msg) msg.textContent = 'Decoding embedded report… ' + Math.round(p * 100) + '%';
          }});
          if (msg) msg.textContent = 'Rendering…';
          const doc = iframe.contentDocument;
          if (!doc) throw new Error('iframe.contentDocument unavailable');
          writeDocFromBytes(doc, bytes);
          iframe.dataset.loaded = 'true';
          status.dataset.state = 'hidden';
          startAutoHeight(iframe);
        }} catch (e) {{
          status.dataset.state = 'error';
          if (msg) msg.textContent = String(e && e.message ? e.message : e);
          console.error(e);
        }}
      }}

      // Tabs
      qsa('.tab').forEach(btn => {{
        btn.addEventListener('click', () => setActiveTab(btn.dataset.tab));
      }});
      setActiveTab('summary');

      // Pie hover tooltip
      const tooltip = qs('#tooltip');
      qsa('.pie-slice').forEach(p => {{
        p.addEventListener('mousemove', (ev) => {{
          const label = p.dataset.label || '';
          const value = p.dataset.value || '';
          tooltip.style.display = 'block';
          qs('.t1', tooltip).textContent = label;
          qs('.t2', tooltip).textContent = `${{value}} lists`;
          const pad = 14;
          const w = tooltip.offsetWidth;
          const h = tooltip.offsetHeight;
          let x = ev.clientX + 14;
          let y = ev.clientY + 14;
          if (x + w + pad > window.innerWidth) x = ev.clientX - w - 14;
          if (y + h + pad > window.innerHeight) y = ev.clientY - h - 14;
          tooltip.style.left = x + 'px';
          tooltip.style.top = y + 'px';
        }});
        p.addEventListener('mouseleave', () => {{
          tooltip.style.display = 'none';
        }});
      }});

      window.addEventListener('resize', () => {{
        const active = qs('.panel[data-active="true"]');
        const iframe = active ? qs('iframe.reportframe', active) : null;
        if (active && iframe && iframe.dataset.loaded === 'true') startAutoHeight(iframe);
      }});

      // Bar chart (pure SVG, single color red)
      const barSvg = qs('#barSvg');
      const barSelect = qs('#barSelect');
      function renderBars(key) {{
        const rows = (DATA.combo_scores_by_indication[key] || []).slice();
        rows.sort((a,b) => (a.rank||9999) - (b.rank||9999));
        const top = rows.slice(0, 12);
        const w = 720, h = 320;
        const padL = 72, padR = 16, padT = 20, padB = 88;
        const innerW = w - padL - padR;
        const innerH = h - padT - padB;
        const maxV = Math.max(1, ...top.map(r => r.overall || 0));
        const bw = innerW / Math.max(1, top.length);
        const fill = '#d7263d';
        let out = '';
        out += `<rect x="0" y="0" width="${{w}}" height="${{h}}" fill="transparent"></rect>`;
        out += `<line x1="${{padL}}" y1="${{padT+innerH}}" x2="${{w-padR}}" y2="${{padT+innerH}}" stroke="rgba(12,19,35,.18)" stroke-width="1"/>`;
        for (let i=0;i<top.length;i++) {{
          const r = top[i];
          const v = (r.overall || 0);
          const bh = (v / maxV) * innerH;
          const x = padL + i*bw + 6;
          const y = padT + (innerH - bh);
          const rw = Math.max(6, bw - 12);
          const name = String(r.name || '').replace(/&/g,'&amp;').replace(/</g,'&lt;');
          out += `<rect x="${{x}}" y="${{y}}" width="${{rw}}" height="${{bh}}" rx="8" fill="${{fill}}" opacity="0.92"><title>${{name}}: ${{v.toFixed(1)}}</title></rect>`;
          const lbl = String(i+1);
          out += `<text x="${{x+rw/2}}" y="${{padT+innerH+18}}" text-anchor="middle" font-size="10" fill="rgba(12,19,35,.70)" font-weight="800">${{lbl}}</text>`;
        }}
        out += `<text x="${{padL}}" y="${{padT+innerH+44}}" font-size="11" fill="rgba(12,19,35,.62)" font-weight="900">Rank</text>`;
        out += `<text x="${{padL}}" y="${{padT-4}}" font-size="11" fill="rgba(12,19,35,.62)" font-weight="900">Overall</text>`;
        barSvg.innerHTML = out;
      }}
      barSelect?.addEventListener('change', () => renderBars(barSelect.value));
      if (barSelect) renderBars(barSelect.value);
    </script>
  </body>
</html>
"""


def _parse_label_pairs(pairs: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for p in pairs:
        if "=" not in p:
            raise SystemExit(f"Invalid --label value: {p!r} (expected KEY=Label)")
        k, v = p.split("=", 1)
        k = k.strip()
        v = v.strip()
        if not k or not v:
            raise SystemExit(f"Invalid --label value: {p!r} (empty KEY or Label)")
        out[k] = v
    return out


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Final stitch step: combine multiple per-indication *_Combo_Prioritization_Report_offline.html files into one offline portfolio HTML."
    )
    ap.add_argument("--root", default=".", help="Working directory to resolve relative globs (default: current dir).")
    ap.add_argument("--glob", default=DEFAULT_GLOB, help=f"Input glob (default: {DEFAULT_GLOB}).")
    ap.add_argument("--out", default="GI2_Combo_Portfolio_Report_offline.html", help="Output HTML file path.")
    ap.add_argument("--built-by", default=DEFAULT_BUILT_BY, help=f'Built-by rewrite text (default: "{DEFAULT_BUILT_BY}").')
    ap.add_argument(
        "--order",
        default=",".join(DEFAULT_ORDER),
        help=f"Tab order as comma-separated keys (default: {','.join(DEFAULT_ORDER)}).",
    )
    ap.add_argument(
        "--label",
        action="append",
        default=[],
        help="Label override as KEY=Label (repeatable). Defaults rename SSc/ATD/HS.",
    )
    args = ap.parse_args()

    root = Path(args.root).resolve()
    out = Path(args.out)
    if not out.is_absolute():
        out = (root / out).resolve()

    order = [x.strip() for x in str(args.order).split(",") if x.strip()]
    labels = dict(DEFAULT_LABELS)
    labels.update(_parse_label_pairs(args.label))

    build(
        root=root,
        output=out,
        report_glob=str(args.glob),
        built_by=str(args.built_by),
        order=order,
        labels=labels,
    )

    print(f"Wrote: {out}")


if __name__ == "__main__":
    main()
