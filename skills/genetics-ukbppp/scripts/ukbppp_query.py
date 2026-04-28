#!/usr/bin/env python3
"""UKB-PPP disease signatures query tool.

Query plasma proteomics protein-disease associations from the UK Biobank
Pharma Proteomics Project. Supports downloading data from S3, resolving
disease names to ICD-10 codes, and exporting filtered/enriched results.
"""

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

S3_BUCKET = "s3://tec-dev-usvga-11158-ukb-sumstats-share/UKBPPP_disease_signatures"
AWS_PROFILE = "cmp-dev"

S3_FILES = {
    "signatures": "final.signatures.csv.gz",
    "protein_info": "protein_info.tsv",
    "icd10_info": "ICD10_info.tsv",
}

COLUMN_MAP = {
    "logit": {
        "beta": "logit_beta",
        "se": "logit_se",
        "pval": "logit_pval",
        "Ncontrol": "logit_Ncontrol",
        "Ncase": "logit_Ncase",
        "fdr": "logit_FDR",
    },
    "ARD": {
        "beta": "ARD_beta",
        "se": "ARD_se",
        "pval": "ARD_pval",
        "Ncontrol": "ARD_Ncontrol",
        "Ncase": "ARD_Ncase",
        "fdr": "ARD_FDR",
    },
    "CoxPH": {
        "beta": "CoxPH_beta",
        "se": "CoxPH_se",
        "pval": "CoxPH_pval",
        "Ncontrol": "CoxPH_Ncontrol",
        "Ncase": "CoxPH_Ncase",
        "fdr": "CoxPH_FDR",
    },
}

# Source data uses "ARD_fdr" (lowercase) — rename early in pipeline
SOURCE_RENAME = {"ARD_fdr": "ARD_FDR"}

MODEL_DESCRIPTIONS = {
    "logit": "Logistic regression Y ~ Covariates + Protein X (all patients regardless of diagnosis time)",
    "ARD": "Automatic Relevance Determination, Protein X ~ Covariates + Diseases (within +/-5 years of sample collection)",
    "CoxPH": "Cox Proportional Hazard regression (patients diagnosed after sample collection only)",
}

ALL_MODELS = ["logit", "ARD", "CoxPH"]

NET_SCORE_MAP = {0: 0, 1: 30, 2: 60, 3: 100, -1: -30, -2: -60, -3: -100}


def _fmt(v, sci=False):
    if pd.isna(v):
        return "—"
    if sci:
        return f"{v:.3e}"
    return f"{v:.4f}"


def _fmt_int(v):
    if pd.isna(v):
        return "—"
    return str(int(v))


def check_sso_auth(profile: str) -> bool:
    try:
        result = subprocess.run(
            ["aws", "s3", "ls", S3_BUCKET + "/", "--profile", profile],
            capture_output=True,
            text=True,
            timeout=30,
        )
    except subprocess.TimeoutExpired:
        return False
    return result.returncode == 0


def download_data(cache_dir: Path, profile: str) -> dict:
    cache_dir.mkdir(parents=True, exist_ok=True)
    downloaded = {}

    for key, filename in S3_FILES.items():
        local_path = cache_dir / filename
        if key == "signatures" and local_path.exists():
            print(f"Cached: {local_path}", file=sys.stderr)
            downloaded[key] = str(local_path)
            continue

        s3_path = f"{S3_BUCKET}/{filename}"
        timeout = 600 if key == "signatures" else 60
        try:
            result = subprocess.run(
                ["aws", "s3", "cp", s3_path, str(local_path), "--profile", profile],
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            print(json.dumps({"status": "error", "file": filename, "error": f"Download timed out after {timeout}s"}))
            sys.exit(1)
        if result.returncode != 0:
            print(json.dumps({"status": "error", "file": filename, "error": result.stderr.strip()}))
            sys.exit(1)

        print(f"Downloaded: {local_path}", file=sys.stderr)
        downloaded[key] = str(local_path)

    print(json.dumps({"status": "ok", "files": downloaded}))
    return downloaded


def resolve_disease(query: str, cache_dir: Path) -> None:
    icd10_path = cache_dir / S3_FILES["icd10_info"]
    if not icd10_path.exists():
        print(json.dumps({"status": "error", "error": f"ICD10_info.tsv not found at {icd10_path}. Run --action download first."}))
        sys.exit(1)

    icd10_df = pd.read_csv(icd10_path, sep="\t", dtype=str)
    icd10_3char = icd10_df[icd10_df["coding"].str.match(r"^[A-Z]\d{2}$", na=False)]

    if re.match(r"^[A-Z]\d{2}$", query.strip().upper()):
        code = query.strip().upper()
        match = icd10_3char[icd10_3char["coding"] == code]
        if match.empty:
            print(json.dumps({"status": "ok", "query": query, "matches": [], "message": f"ICD10 code {code} not found in reference data."}))
        else:
            matches = [{"coding": r["coding"], "meaning": r["meaning"]} for _, r in match.iterrows()]
            print(json.dumps({"status": "ok", "query": query, "matches": matches}))
        return

    query_lower = query.strip().lower()
    hits = icd10_df[icd10_df["meaning"].str.lower().str.contains(query_lower, na=False, regex=False)]

    if hits.empty:
        print(json.dumps({"status": "ok", "query": query, "matches": [], "message": f"No ICD10 codes found matching '{query}'."}))
        return

    parent_codes = set()
    for code in hits["coding"]:
        if re.match(r"^[A-Z]\d{2}$", str(code)):
            parent_codes.add(str(code))
        elif re.match(r"^[A-Z]\d{2,}$", str(code)):
            parent_codes.add(str(code)[:3])

    result_rows = icd10_3char[icd10_3char["coding"].isin(parent_codes)].sort_values("coding")
    matches = [{"coding": r["coding"], "meaning": r["meaning"]} for _, r in result_rows.iterrows()]

    print(json.dumps({"status": "ok", "query": query, "matches": matches, "total_hits": len(hits), "parent_codes": len(matches)}))


def load_data(cache_dir: Path):
    sig_path = cache_dir / S3_FILES["signatures"]
    prot_path = cache_dir / S3_FILES["protein_info"]
    icd10_path = cache_dir / S3_FILES["icd10_info"]

    for p in [sig_path, prot_path, icd10_path]:
        if not p.exists():
            print(json.dumps({"status": "error", "error": f"{p.name} not found. Run --action download first."}))
            sys.exit(1)

    print("Loading signatures data (~2M rows)...", file=sys.stderr)
    signatures_df = pd.read_csv(sig_path, compression="gzip")
    print(f"Loaded {len(signatures_df)} rows.", file=sys.stderr)

    protein_info_df = pd.read_csv(prot_path, sep="\t")
    icd10_info_df = pd.read_csv(icd10_path, sep="\t", dtype=str)

    return signatures_df, protein_info_df, icd10_info_df


def validate_genes(genes: list, signatures_df: pd.DataFrame):
    available = set(signatures_df["Assay"].unique())
    found = [g for g in genes if g in available]
    missing = [g for g in genes if g not in available]
    return found, missing


def get_selected_columns(models: list) -> list:
    cols = ["UKBPPP_ProteinID", "ICD10", "Assay"]
    for model in models:
        for stat_col in COLUMN_MAP[model].values():
            cols.append(stat_col)
    return cols


def filter_and_enrich(
    signatures_df: pd.DataFrame,
    protein_info_df: pd.DataFrame,
    icd10_info_df: pd.DataFrame,
    genes: list,
    icd10_codes: list,
    models: list,
) -> pd.DataFrame:
    required_cols = {"Assay", "ICD10", "UKBPPP_ProteinID"}
    missing = required_cols - set(signatures_df.columns)
    if missing:
        print(json.dumps({"status": "error", "error": f"Missing required columns in signatures data: {missing}"}))
        sys.exit(1)

    if "ARD_fdr" in signatures_df.columns:
        signatures_df = signatures_df.rename(columns=SOURCE_RENAME)

    mask = signatures_df["Assay"].isin(genes) & signatures_df["ICD10"].isin(icd10_codes)
    filtered = signatures_df.loc[mask].copy()

    if filtered.empty:
        return filtered

    selected_cols = get_selected_columns(models)
    existing_cols = [c for c in selected_cols if c in filtered.columns]
    filtered = filtered[existing_cols].copy()

    icd10_3char = icd10_info_df[icd10_info_df["coding"].str.match(r"^[A-Z]\d{2}$", na=False)].copy()
    icd10_3char = icd10_3char.rename(columns={"meaning": "disease_name"})
    filtered = filtered.merge(icd10_3char[["coding", "disease_name"]], left_on="ICD10", right_on="coding", how="left")
    filtered = filtered.drop(columns=["coding"], errors="ignore")

    filtered = filtered.merge(
        protein_info_df[["UKBPPP_ProteinID", "olink_target_fullname", "ensembl_id"]],
        on="UKBPPP_ProteinID",
        how="left",
    )

    front_cols = ["Assay", "UKBPPP_ProteinID", "ICD10", "disease_name"]
    back_cols = ["olink_target_fullname", "ensembl_id"]
    model_cols = [c for c in filtered.columns if c not in front_cols + back_cols]
    ordered = front_cols + model_cols + back_cols
    ordered = [c for c in ordered if c in filtered.columns]
    filtered = filtered[ordered]

    return filtered


def _score_gene_indication(row, models: list) -> int:
    sig_pos = 0
    sig_neg = 0
    for model in models:
        mcols = COLUMN_MAP[model]
        beta = row.get(mcols["beta"])
        fdr = row.get(mcols["fdr"])
        if pd.isna(beta) or pd.isna(fdr) or fdr >= 0.05:
            continue
        if beta > 0:
            sig_pos += 1
        else:
            sig_neg += 1
    net = sig_pos - sig_neg
    return NET_SCORE_MAP.get(net, 0)


def compute_scores(
    enriched_df: pd.DataFrame,
    genes: list,
    found_genes: list,
    missing_genes: list,
    icd10_codes: list,
    models: list,
) -> dict:
    scores = {}
    for gene in genes:
        gene_scores = {}
        if gene in missing_genes:
            for icd_code in icd10_codes:
                gene_scores[icd_code] = 0
        else:
            gene_data = enriched_df[enriched_df["Assay"] == gene]
            for icd_code in icd10_codes:
                pair_data = gene_data[gene_data["ICD10"] == icd_code]
                if pair_data.empty:
                    gene_scores[icd_code] = 0
                else:
                    gene_scores[icd_code] = _score_gene_indication(pair_data.iloc[0], models)
        indication_vals = [gene_scores[c] for c in icd10_codes]
        gene_scores["_total"] = sum(indication_vals) / len(indication_vals) if indication_vals else 0
        scores[gene] = gene_scores
    return scores


def compute_combo_scores(
    gene_scores: dict,
    combos: list,
    icd10_codes: list,
) -> dict:
    combo_scores = {}
    for combo in combos:
        label = "+".join(combo)
        cs = {}
        for icd_code in icd10_codes:
            vals = [gene_scores.get(g, {}).get(icd_code, 0) for g in combo]
            cs[icd_code] = sum(vals) / len(vals) if vals else 0
        indication_vals = [cs[c] for c in icd10_codes]
        cs["_total"] = sum(indication_vals) / len(indication_vals) if indication_vals else 0
        combo_scores[label] = cs
    return combo_scores


def generate_summary(
    enriched_df: pd.DataFrame,
    genes: list,
    found_genes: list,
    missing_genes: list,
    icd10_codes: list,
    icd10_info_df: pd.DataFrame,
    models: list,
    gene_scores: dict,
    combo_scores: dict = None,
    combos: list = None,
) -> str:
    icd10_3char = icd10_info_df[icd10_info_df["coding"].str.match(r"^[A-Z]\d{2}$", na=False)]
    code_to_name = dict(zip(icd10_3char["coding"], icd10_3char["meaning"]))

    icd_display = ", ".join(f"{c} ({code_to_name.get(c, 'Unknown')})" for c in icd10_codes)
    model_display = ", ".join(models)

    lines = [
        "# UKB-PPP Query Summary",
        "",
        f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"**Models:** {model_display}",
        f"**ICD10 codes queried:** {icd_display}",
        "**FDR summary threshold:** < 0.05",
        "",
        "## Scoring Method",
        "",
        "For each gene-indication pair, only FDR-significant models (FDR<0.05) are counted.",
        "Net direction score: net +1→30, +2→60, +3→100, -1→-30, -2→-60, -3→-100, 0→0.",
        "Opposing directions cancel (e.g., 1 pos + 1 neg = net 0 → score 0).",
        "Gene total = average across all indication scores. Missing genes score 0.",
        "",
    ]

    # --- Individual Gene Rankings ---
    ranked = sorted(genes, key=lambda g: gene_scores.get(g, {}).get("_total", 0), reverse=True)
    icd_headers = " | ".join(f"{c} Score" for c in icd10_codes)
    icd_dashes = " | ".join("---:" for _ in icd10_codes)
    def _gene_status(gene):
        if gene in missing_genes:
            return "Not on panel"
        gs = gene_scores.get(gene, {})
        if any(gs.get(c, 0) != 0 for c in icd10_codes):
            return "Significant"
        return "Not significant"

    lines.append("## Individual Gene Rankings")
    lines.append("")
    lines.append(f"| Rank | Gene | {icd_headers} | Total Score | Status |")
    lines.append(f"|------|------|{icd_dashes}|------------:|--------|")
    for rank, gene in enumerate(ranked, 1):
        gs = gene_scores.get(gene, {})
        icd_vals = " | ".join(str(gs.get(c, 0)) for c in icd10_codes)
        total = gs.get("_total", 0)
        status = _gene_status(gene)
        lines.append(f"| {rank} | {gene} | {icd_vals} | {total:.1f} | {status} |")
    lines.append("")

    # --- Individual Gene Details ---
    for gene in ranked:
        lines.append(f"## {gene}")
        lines.append("")

        gs = gene_scores.get(gene, {})
        total = gs.get("_total", 0)
        lines.append(f"**Total Score: {total:.1f}**")
        lines.append("")

        if gene in missing_genes:
            lines.append("**Gene not found in UKB-PPP data.** Score = 0 for all indications.")
            lines.append("")
            continue

        gene_data = enriched_df[enriched_df["Assay"] == gene]

        if gene_data.empty:
            lines.append("No data found for this gene in the queried diseases.")
            lines.append("")
            continue

        for icd_code in icd10_codes:
            disease_name = code_to_name.get(icd_code, "Unknown")
            pair_data = gene_data[gene_data["ICD10"] == icd_code]
            indication_score = gs.get(icd_code, 0)

            lines.append(f"### {gene} x {icd_code} ({disease_name}) — Score: {indication_score}")
            lines.append("")

            if pair_data.empty:
                lines.append("No data for this gene-disease pair.")
                lines.append("")
                continue

            for _, row in pair_data.iterrows():
                if len(pair_data) > 1:
                    protein_id = row.get("UKBPPP_ProteinID", "")
                    lines.append(f"**Assay: {protein_id}**")
                    lines.append("")

                lines.append("| Model | Beta | SE | P-value | FDR | Ncase | Ncontrol | Sig (FDR<0.05) |")
                lines.append("|-------|------|----|---------|-----|-------|----------|----------------|")

                sig_count = 0
                directions = {}
                available_models = 0

                for model in models:
                    mcols = COLUMN_MAP[model]
                    beta = row.get(mcols["beta"])
                    se = row.get(mcols["se"])
                    pval = row.get(mcols["pval"])
                    fdr = row.get(mcols["fdr"])
                    ncase = row.get(mcols["Ncase"])
                    ncontrol = row.get(mcols["Ncontrol"])

                    if pd.isna(beta):
                        lines.append(f"| {model} | — | — | — | — | — | — | N/A |")
                        continue

                    available_models += 1
                    is_sig = not pd.isna(fdr) and fdr < 0.05
                    if is_sig:
                        sig_count += 1
                    sig_str = "**Yes**" if is_sig else "No"
                    directions[model] = "positive" if beta > 0 else "negative"

                    lines.append(
                        f"| {model} | {_fmt(beta)} | {_fmt(se)} | {_fmt(pval, True)} | {_fmt(fdr, True)} | {_fmt_int(ncase)} | {_fmt_int(ncontrol)} | {sig_str} |"
                    )

                lines.append("")

                if available_models > 0:
                    dir_summary = ", ".join(f"{m} ({d})" for m, d in directions.items())
                    lines.append(f"**Concordance:** Significant in {sig_count}/{available_models} models. Directions: {dir_summary}.")
                else:
                    lines.append("**Concordance:** No model data available for this pair.")

                lines.append("")

    # --- Combo Section ---
    if combos and combo_scores:
        combo_ranked = sorted(combo_scores.keys(), key=lambda k: combo_scores[k].get("_total", 0), reverse=True)

        lines.append("## Combo Rankings")
        lines.append("")
        lines.append(f"| Rank | Combo | {icd_headers} | Total Score |")
        lines.append(f"|------|-------|{icd_dashes}|------------:|")
        for rank, label in enumerate(combo_ranked, 1):
            cs = combo_scores[label]
            icd_vals = " | ".join(f"{cs.get(c, 0):.1f}" for c in icd10_codes)
            total = cs.get("_total", 0)
            lines.append(f"| {rank} | {label} | {icd_vals} | {total:.1f} |")
        lines.append("")

        for label in combo_ranked:
            combo_genes = label.split("+")
            cs = combo_scores[label]
            total = cs.get("_total", 0)

            lines.append(f"## Combo: {label}")
            lines.append("")
            lines.append(f"**Total Combo Score: {total:.1f}**")
            lines.append("")

            gene_score_headers = " | ".join(f"{g} Score" for g in combo_genes)
            gene_score_dashes = " | ".join("---:" for _ in combo_genes)
            lines.append(f"| ICD10 | Disease | {gene_score_headers} | Combo Score |")
            lines.append(f"|-------|---------|{gene_score_dashes}|------------:|")

            for icd_code in icd10_codes:
                disease_name = code_to_name.get(icd_code, "Unknown")
                gene_vals = " | ".join(str(gene_scores.get(g, {}).get(icd_code, 0)) for g in combo_genes)
                combo_val = cs.get(icd_code, 0)
                lines.append(f"| {icd_code} | {disease_name} | {gene_vals} | {combo_val:.1f} |")

            lines.append("")

    return "\n".join(lines)


def export_results(
    enriched_df: pd.DataFrame,
    genes: list,
    found_genes: list,
    missing_genes: list,
    output_dir: Path,
    models: list,
    icd10_codes: list,
    icd10_info_df: pd.DataFrame,
    combos: list = None,
) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)

    icd10_3char = icd10_info_df[icd10_info_df["coding"].str.match(r"^[A-Z]\d{2}$", na=False)]
    code_to_name = dict(zip(icd10_3char["coding"], icd10_3char["meaning"]))
    icd_display = ", ".join(f"{c} ({code_to_name.get(c, 'Unknown')})" for c in icd10_codes)

    icd10_suffix = "_".join(icd10_codes)
    file_summary = {}

    for gene in found_genes:
        gene_df = enriched_df[enriched_df["Assay"] == gene]
        out_path = output_dir / f"{gene}_{icd10_suffix}_ukbppp.tsv"
        gene_df.to_csv(out_path, sep="\t", index=False)
        file_summary[gene] = {"file": str(out_path), "rows": len(gene_df)}
        print(f"Exported: {out_path} ({len(gene_df)} rows)", file=sys.stderr)

    for gene in missing_genes:
        file_summary[gene] = {"rows": 0, "note": "Gene not found in UKB-PPP data (no TSV exported)"}
        print(f"Skipped (not on Olink panel): {gene}", file=sys.stderr)

    gene_scores = compute_scores(enriched_df, genes, found_genes, missing_genes, icd10_codes, models)

    combo_scores_dict = None
    if combos:
        combo_scores_dict = compute_combo_scores(gene_scores, combos, icd10_codes)

    summary_md = generate_summary(
        enriched_df, genes, found_genes, missing_genes, icd10_codes, icd10_info_df, models,
        gene_scores, combo_scores_dict, combos,
    )
    summary_path = output_dir / "ukbppp_summary.md"
    with open(summary_path, "w") as f:
        f.write(summary_md)
    print(f"Summary: {summary_path}", file=sys.stderr)

    log_content = _generate_log(genes, found_genes, missing_genes, models, icd10_codes, code_to_name, file_summary)
    log_path = output_dir / "ukbppp_query.log"
    with open(log_path, "w") as f:
        f.write(log_content)
    print(f"Log: {log_path}", file=sys.stderr)

    result = {
        "status": "ok",
        "files": file_summary,
        "summary_file": str(summary_path),
        "log_file": str(log_path),
        "genes_found": found_genes,
        "genes_missing": missing_genes,
    }
    print(json.dumps(result))
    return result


def _generate_log(
    genes: list,
    found_genes: list,
    missing_genes: list,
    models: list,
    icd10_codes: list,
    code_to_name: dict,
    file_summary: dict,
) -> str:
    icd_display = ", ".join(f"{c} ({code_to_name.get(c, 'Unknown')})" for c in icd10_codes)
    model_display = ", ".join(models)

    lines = [
        "=== UKB-PPP Query Log ===",
        f"Timestamp: {datetime.now().isoformat()}",
        f"Genes queried: {', '.join(genes)}",
        f"Genes found: {', '.join(found_genes) if found_genes else 'None'}",
        f"Genes NOT found: {', '.join(missing_genes) if missing_genes else 'None'}",
        f"ICD10 codes: {icd_display}",
        f"Models exported: {model_display}",
        "Summary FDR threshold: 0.05",
        "",
        "--- Column Descriptions ---",
    ]

    for model in models:
        prefix = model
        lines.append(f"{prefix}_beta: Effect size")
        lines.append(f"{prefix}_se: Standard error of the estimated effect size")
        lines.append(f"{prefix}_pval: Raw p-value")
        lines.append(f"{prefix}_Ncontrol: Number of control subjects")
        lines.append(f"{prefix}_Ncase: Number of case subjects")
        lines.append(f"{prefix}_FDR: False discovery rate (Benjamini-Hochberg)")
    lines.append("")

    lines.append("--- Model Notes ---")
    for model in models:
        lines.append(f"{model}: {MODEL_DESCRIPTIONS[model]}")
    lines.append("")

    lines.append("--- Results Summary ---")
    for gene in genes:
        info = file_summary.get(gene, {})
        if "note" in info and "not found" in info["note"].lower():
            lines.append(f"{gene}: NOT FOUND in UKB-PPP data (no TSV exported)")
        else:
            rows = info.get("rows", 0)
            lines.append(f"{gene}: {rows} rows exported to {info.get('file', gene + '_ukbppp.tsv')}")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="UKB-PPP disease signatures query tool")
    parser.add_argument("--action", required=True, choices=["download", "resolve-disease", "query"])
    parser.add_argument("--cache-dir", default=str(Path.home() / "tmp"), help="Local cache directory (default: ~/tmp)")
    parser.add_argument("--profile", default=AWS_PROFILE, help=f"AWS profile (default: {AWS_PROFILE})")
    parser.add_argument("--disease", help="Disease name or ICD10 code (for resolve-disease action)")
    parser.add_argument("--genes", help="Comma-separated gene symbols (for query action)")
    parser.add_argument("--icd10", help="Comma-separated ICD10 3-char codes (for query action)")
    parser.add_argument("--models", default="all", help="Model selection: all, logit, ARD, CoxPH (default: all)")
    parser.add_argument("--combos", help="Combo groups, e.g. 'IL11+OSM,OSM+GREM1' (optional)")
    parser.add_argument("--output-dir", default=".", help="Output directory (default: current directory)")
    args = parser.parse_args()

    cache_dir = Path(args.cache_dir)

    if args.action == "download":
        download_data(cache_dir, args.profile)

    elif args.action == "resolve-disease":
        if not args.disease:
            print(json.dumps({"status": "error", "error": "--disease is required for resolve-disease action"}))
            sys.exit(1)
        resolve_disease(args.disease, cache_dir)

    elif args.action == "query":
        if not args.genes or not args.icd10:
            print(json.dumps({"status": "error", "error": "--genes and --icd10 are required for query action"}))
            sys.exit(1)

        genes = [g.strip().upper() for g in args.genes.split(",") if g.strip()]
        icd10_codes = [c.strip().upper() for c in args.icd10.split(",") if c.strip()]

        if args.models == "all":
            models = ALL_MODELS
        elif args.models in COLUMN_MAP:
            models = [args.models]
        else:
            print(json.dumps({"status": "error", "error": f"Invalid model: {args.models}. Choose from: all, logit, ARD, CoxPH"}))
            sys.exit(1)

        combos = None
        if args.combos:
            combos = [
                [g.strip().upper() for g in group.split("+")]
                for group in args.combos.split(",")
                if group.strip()
            ]

        output_dir = Path(args.output_dir)

        signatures_df, protein_info_df, icd10_info_df = load_data(cache_dir)
        found_genes, missing_genes = validate_genes(genes, signatures_df)

        if found_genes:
            enriched_df = filter_and_enrich(signatures_df, protein_info_df, icd10_info_df, found_genes, icd10_codes, models)
        else:
            enriched_df = pd.DataFrame()

        export_results(enriched_df, genes, found_genes, missing_genes, output_dir, models, icd10_codes, icd10_info_df, combos)


if __name__ == "__main__":
    main()
