#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, UTC
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

PHASE_ORDER_DEFAULT = [
    'Unknown',
    'Preclinical',
    'IND filed',
    'Phase I',
    'Phase I/II',
    'Phase II',
    'Phase III',
    'Pre-registered',
    'Registered',
    'Marketed',
]


def split_multi(val: Any) -> list[str]:
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return []
    parts = re.split(r"\s*\n\s*|\s*;\s*", str(val))
    return [p.strip() for p in parts if p and p.strip()]


def clean_name(val: Any) -> str:
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return ''
    s = str(val).strip()
    if not s:
        return ''
    s = s.split('\n')[0].strip()
    s = re.sub(r"\s*\([^)]*\)\s*$", "", s).strip()
    return s


def one_line(val: Any, limit: int = 360) -> str:
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return ''
    s = re.sub(r"\s+", " ", str(val)).strip()
    return s if len(s) <= limit else (s[: limit - 1] + '…')


def find_col(cols: list[str], candidates: list[str]) -> str | None:
    cols_l = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in cols_l:
            return cols_l[cand.lower()]
    for cand in candidates:
        cl = cand.lower()
        for c in cols:
            if cl in c.lower():
                return c
    return None


def norm_phase(val: Any) -> str:
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return 'Unknown'
    s = str(val).strip()
    if not s:
        return 'Unknown'
    sl = s.lower()
    if sl.startswith('launched') or 'on market' in sl or 'marketed' in sl:
        return 'Marketed'
    if sl.startswith('registered'):
        return 'Registered'
    if 'pre-registered' in sl:
        return 'Pre-registered'
    if 'phase iii' in sl:
        return 'Phase III'
    if 'phase i/ii' in sl:
        return 'Phase I/II'
    if 'phase ii' in sl:
        return 'Phase II'
    if re.search(r"\bphase\s+i\b", sl):
        return 'Phase I'
    if 'preclinical' in sl:
        return 'Preclinical'
    if 'ind filed' in sl:
        return 'IND filed'
    return 'Unknown'


@dataclass
class Rule:
    name: str
    patterns: list[str]

    def matches(self, blob: str) -> bool:
        if not self.patterns:
            return False
        for p in self.patterns:
            if re.search(p, blob, flags=re.IGNORECASE):
                return True
        return False


def load_taxonomy(path: Path | None) -> dict[str, Any]:
    if not path:
        return {}
    obj = yaml.safe_load(path.read_text(encoding='utf-8'))
    return obj or {}


def compile_rules(items: list[dict[str, Any]] | None) -> list[Rule]:
    rules: list[Rule] = []
    for it in items or []:
        name = str(it.get('name') or it.get('tag') or '').strip()
        pats = it.get('patterns') or []
        if not name:
            continue
        rules.append(Rule(name=name, patterns=[str(p) for p in pats]))
    return rules


def classify_first(blob: str, rules: list[Rule], fallback: str = 'Other') -> str:
    for r in rules:
        if r.matches(blob):
            return r.name
    return fallback


def classify_many(blob: str, rules: list[Rule]) -> list[str]:
    out: list[str] = []
    for r in rules:
        if r.matches(blob):
            out.append(r.name)
    return out


def modality_from(drug_type: str, product_category: str) -> str:
    dt = (drug_type or '').lower()
    pc = (product_category or '').lower()
    if 'small molecules' in dt or 'small molecule' in dt:
        return 'Small molecule'
    if 'peptide' in dt:
        return 'Peptide'
    if 'herbal' in dt:
        return 'Herbal'
    if any(k in pc for k in ['probiotic', 'live bacterial', 'microbiome']):
        return 'Microbiome'
    if any(k in pc for k in ['car t', 'treg', 'mesenchymal', 'stem cell', 'cell']):
        return 'Cell therapy'
    if 'biologics' in dt or 'biologic' in dt:
        if 'antibod' in pc:
            return 'Antibody'
        if 'fusion protein' in pc or 'recombinant' in pc:
            return 'Protein'
        return 'Biologic'
    if 'antibod' in pc:
        return 'Antibody'
    return 'Other'


def main() -> None:
    ap = argparse.ArgumentParser(description='Generate a CI interactive HTML dashboard from an Excel export.')
    ap.add_argument('--input', required=True, help='Input .xlsx file')
    ap.add_argument('--output', required=True, help='Output .html file')
    ap.add_argument(
        '--taxonomy',
        default=str(Path(__file__).resolve().parents[1] / 'references' / 'taxonomy_default.yaml'),
        help='Taxonomy YAML (rules for target type/family + optional focus tag rules)',
    )
    ap.add_argument(
        '--template',
        default=str(Path(__file__).resolve().parents[1] / 'assets' / 'dashboard_template.html'),
        help='Dashboard HTML template (must contain __PAYLOAD__ placeholder)',
    )
    ap.add_argument(
        '--context-keywords',
        default='',
        help='Regex/keywords for context phase filtering (comma-separated; matches dev-status Condition/Indication)',
    )
    ap.add_argument('--sheet', default='', help='Optional main sheet name override')
    args = ap.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output)
    tpl_path = Path(args.template)

    taxonomy = load_taxonomy(Path(args.taxonomy) if args.taxonomy else None)
    phase_order = taxonomy.get('phase_order') or PHASE_ORDER_DEFAULT
    phase_rank = {p: i for i, p in enumerate(phase_order)}

    target_type_rules = compile_rules(taxonomy.get('target_type_rules'))
    target_family_rules = compile_rules(taxonomy.get('target_family_rules'))
    focus_tag_rules = compile_rules(taxonomy.get('focus_tag_rules'))

    xl = pd.ExcelFile(in_path)
    sheets = xl.sheet_names

    # Select main sheet
    main_sheet = args.sheet.strip() if args.sheet.strip() else None
    if not main_sheet:
        for cand in ['Product List', 'Products', 'Drugs & Biologics']:
            if cand in sheets:
                main_sheet = cand
                break
    if not main_sheet:
        best = None
        for s in sheets:
            df = pd.read_excel(in_path, sheet_name=s)
            if best is None or df.shape[0] > best[1]:
                best = (s, df.shape[0])
        main_sheet = best[0]

    main = pd.read_excel(in_path, sheet_name=main_sheet)

    dev = None
    for cand in ['Development Status', 'Development', 'Clinical Studies']:
        if cand in sheets:
            dev = pd.read_excel(in_path, sheet_name=cand)
            break

    # Column mapping (heuristic)
    cols = list(main.columns)

    col_entry = find_col(cols, ['Entry Number', 'Entry', 'ID'])
    col_name = find_col(cols, ['Display Name', 'Drug Name', 'Drug Name (All)', 'Main Name', 'Name', 'Compound', 'Asset'])
    col_brand = find_col(cols, ['Brand Name', 'Brand'])
    col_generic = find_col(cols, ['Generic Name', 'Generic'])
    col_code = find_col(cols, ['Code Name', 'Code'])
    col_org = find_col(cols, ['Organization', 'Company', 'Sponsor'])
    col_active = find_col(cols, ['Under Active Development', 'Active'])
    col_drug_type = find_col(cols, ['Drug Type', 'Type', 'Modality'])
    col_product_cat = find_col(cols, ['Product Category', 'Category'])
    col_phase = find_col(cols, ['Highest Phase', 'Phase', 'Development Phase'])
    col_target = find_col(cols, ['Target', 'Targets'])
    col_moa = find_col(cols, ['Mechanism of Action', 'MOA'])
    col_mm = find_col(cols, ['Molecular Mechanism'])
    col_cm = find_col(cols, ['Cellular Mechanism'])
    col_cond = find_col(cols, ['Condition', 'Indication', 'Disease'])

    if not col_name and cols:
        col_name = cols[0]

    # Build context phase/routes map from dev sheet if available
    context_re = None
    if args.context_keywords.strip():
        parts = [p.strip() for p in args.context_keywords.split(',') if p.strip()]
        if parts:
            context_re = re.compile('|'.join(parts), flags=re.IGNORECASE)

    context_phase_map: dict[int, tuple[str, str]] = defaultdict(lambda: ('Unknown', ''))
    routes_map: dict[int, set[str]] = defaultdict(set)

    if dev is not None:
        dev_cols = list(dev.columns)
        d_entry = find_col(dev_cols, ['Entry Number', 'Entry', 'ID'])
        d_phase = find_col(dev_cols, ['Phase', 'Development Phase'])
        d_cond = find_col(dev_cols, ['Condition', 'Indication', 'Disease'])
        d_route = find_col(dev_cols, ['Administration Route', 'Route'])

        if d_entry and d_phase:
            for _, row in dev.iterrows():
                try:
                    en = int(row.get(d_entry))
                except Exception:
                    continue

                cond_val = str(row.get(d_cond) or '') if d_cond else ''
                if context_re is not None and not context_re.search(cond_val):
                    continue

                ph_raw = row.get(d_phase)
                ph = norm_phase(ph_raw)
                if phase_rank.get(ph, 0) > phase_rank.get(context_phase_map[en][0], 0):
                    context_phase_map[en] = (ph, one_line(ph_raw, 90))

                if d_route:
                    for rt in split_multi(row.get(d_route)):
                        rt2 = rt.strip()
                        if rt2:
                            routes_map[en].add(rt2)

    entries: list[dict[str, Any]] = []
    for i, row in main.iterrows():
        entry_number: int | None = None
        if col_entry and not pd.isna(row.get(col_entry)):
            try:
                entry_number = int(row.get(col_entry))
            except Exception:
                entry_number = None

        drug_all = row.get(col_name) if col_name else None
        generic = row.get(col_generic) if col_generic else None
        brand = row.get(col_brand) if col_brand else None
        code = row.get(col_code) if col_code else None

        display = clean_name(generic) or clean_name(brand) or clean_name(code) or clean_name(drug_all) or f'Row {i+1}'

        targets = split_multi(row.get(col_target)) if col_target else []
        conditions = split_multi(row.get(col_cond)) if col_cond else []
        moa = one_line(row.get(col_moa), 800) if col_moa else ''
        mm = one_line(row.get(col_mm), 800) if col_mm else ''
        cm = one_line(row.get(col_cm), 800) if col_cm else ''

        drug_type = one_line(row.get(col_drug_type), 90) if col_drug_type else ''
        product_category = one_line(row.get(col_product_cat), 200) if col_product_cat else ''
        modality = modality_from(drug_type, product_category)

        blob = ' '.join([*(targets or []), moa, mm, cm, drug_type, product_category])

        _ = classify_first(blob, target_type_rules, fallback='Other')
        fams = classify_many(blob, target_family_rules)
        if not fams:
            fams = ['Other'] if targets else []
        primary_family = fams[0] if fams else 'Other'

        # Focus tags
        focus_tags: list[str] = []
        if focus_tag_rules and conditions:
            cond_blob = ' '.join(conditions)
            for rr in focus_tag_rules:
                if rr.matches(cond_blob):
                    focus_tags.append(rr.name)
            focus_tags = sorted(set(focus_tags))

        overall_phase_raw = one_line(row.get(col_phase), 90) if col_phase else ''
        overall_phase = norm_phase(row.get(col_phase)) if col_phase else 'Unknown'

        context_phase, context_phase_raw = ('Unknown', '')
        if entry_number is not None:
            context_phase, context_phase_raw = context_phase_map[entry_number]
        if context_phase == 'Unknown':
            context_phase = overall_phase
            context_phase_raw = overall_phase_raw

        routes = sorted(routes_map.get(entry_number or -1, set())) if routes_map else []

        org = one_line(row.get(col_org), 180) if col_org else ''
        active = str(row.get(col_active) or '').strip() if col_active else ''

        entries.append({
            'entryNumber': entry_number if entry_number is not None else (i + 1),
            'displayName': display,
            'brandName': clean_name(brand),
            'genericName': clean_name(generic),
            'organization': org,
            'underActiveDevelopment': active,
            'drugType': drug_type,
            'productCategory': product_category,
            'modality': modality,
            'targets': targets,
            'targetFamilies': fams,
            'targetFamilyPrimary': primary_family,
            'mechanismOfAction': moa,
            'molecularMechanism': mm,
            'cellularMechanism': cm,
            'conditions': conditions,
            'focusTags': focus_tags,
            'overallPhaseRaw': overall_phase_raw,
            'overallPhase': overall_phase,
            'contextPhaseRaw': context_phase_raw,
            'contextPhase': context_phase,
            'routes': routes,
        })

    data = {
        'meta': {
            'sourceFile': str(in_path.name),
            'generatedAt': datetime.now(UTC).replace(microsecond=0).isoformat().replace('+00:00', 'Z'),
            'count': len(entries),
            'notes': [
                "Phase modes: 'Context phase' can be computed from dev-status rows filtered by context keywords; 'Overall phase' uses the primary sheet phase column.",
                f"Main sheet: {main_sheet}",
            ],
        },
        'entries': entries,
        'vocab': {
            'phaseOrder': phase_order,
        },
    }

    payload = json.dumps(data, ensure_ascii=False).replace('</', '<\\/')
    template = tpl_path.read_text(encoding='utf-8')
    if '__PAYLOAD__' not in template:
        raise SystemExit(f'Template missing __PAYLOAD__ placeholder: {tpl_path}')

    out_html = template.replace('__PAYLOAD__', payload)
    out_path.write_text(out_html, encoding='utf-8')
    print(f'Wrote {out_path} ({out_path.stat().st_size/1024:.1f} KB)')


if __name__ == '__main__':
    main()
