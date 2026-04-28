#!/usr/bin/env python3
"""DrugnomeAI Individual + Combo Druggability Scoring.

Reads the DrugnomeAI master table and produces:
  1. drugnome_individual_score.tsv  — per-gene druggability scores
  2. drugnome_combo_score.tsv       — per-combo mean scores (only if combos provided)

Scores are numeric only. Interpretation is done dynamically by the
calling agent (Claude Code / Codex) after TSVs are generated.

If output files already exist, only missing genes/combos are appended.
"""

import argparse
import os
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DEFAULT_MASTER = Path(
    os.path.expanduser(
        "~/.claude/skills/drugnomeai/references/"
        "druggability_master_table.csv"
    )
)
DEFAULT_OUTPUT_DIR = Path.cwd()

# Columns to extract for individual scoring
INDIVIDUAL_COLS = [
    "Gene_Name",
    "composite_score",
    "pharos_mean_proba",
    "tier_mean_proba",
    "best_modality",
    "best_modality_proba",
    "modality_specificity",
    "small_mol_proba",
    "antibody_proba",
    "protac_proba",
    "consensus_tier",
    "n_runs_top_decile",
    "n_runs_above_75perc",
    "is_novel_everywhere",
]

# Per-run proba columns for detailed individual output
PER_RUN_PROBA_COLS = [
    "tclin_proba", "tchem_proba", "tclin_tchem_proba",
    "tier1_proba", "tier12_proba", "tier123A_proba", "tier123AB_proba",
    "tclin_tier1_proba",
    "small_mol_proba", "antibody_proba", "protac_proba",
]

PER_RUN_KNOWN_COLS = [
    "tclin_known", "tchem_known", "tclin_tchem_known",
    "tier1_known", "tier12_known", "tier123A_known", "tier123AB_known",
    "tclin_tier1_known",
    "small_mol_known", "antibody_known", "protac_known",
]

# Numeric columns to average for combo scoring
NUMERIC_MEAN_COLS = [
    "composite_score",
    "pharos_mean_proba",
    "tier_mean_proba",
    "best_modality_proba",
    "modality_specificity",
    "n_runs_top_decile",
    "n_runs_above_75perc",
    "small_mol_proba",
    "antibody_proba",
    "protac_proba",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def get_combo_consensus_tier(score: float) -> str:
    if score >= 0.7:
        return "High"
    if score >= 0.3:
        return "Moderate"
    return "Low"


def get_majority_modality(modalities: list) -> str:
    counts = Counter(modalities)
    max_count = counts.most_common(1)[0][1]
    tied = sorted(m for m, c in counts.items() if c == max_count)
    return tied[0]


def parse_gene_list(text: str) -> list:
    """Parse comma/space/semicolon-separated gene names."""
    import re
    return [g.strip() for g in re.split(r"[,;\s]+", text) if g.strip()]


def parse_combos(combo_strs: list) -> list:
    """Parse combo strings like 'GENE1+GENE2+GENE3'."""
    combos = []
    for s in combo_strs:
        genes = [g.strip() for g in s.split("+") if g.strip()]
        if len(genes) >= 2:
            combos.append(genes)
    return combos


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------
def load_master_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Master table not found at: {path}\n"
            f"Expected location: ~/.claude/skills/drugnomeai/references/druggability_master_table.csv\n"
            f"Copy from DrugnomeAI output or regenerate via the skill's Section F."
        )
    return pd.read_csv(path)


def compute_individual_rows(master: pd.DataFrame, genes: list) -> pd.DataFrame:
    missing = set(genes) - set(master["Gene_Name"])
    if missing:
        print(f"  WARNING: genes not found in master table: {sorted(missing)}")
        genes = [g for g in genes if g not in missing]

    all_cols = INDIVIDUAL_COLS + [
        c for c in PER_RUN_PROBA_COLS + PER_RUN_KNOWN_COLS
        if c not in INDIVIDUAL_COLS
    ]
    available_cols = [c for c in all_cols if c in master.columns]
    subset = master[master["Gene_Name"].isin(genes)][available_cols].copy()
    return subset.sort_values("composite_score", ascending=False).reset_index(drop=True)


def compute_combo_row(genes: list, master: pd.DataFrame) -> dict:
    rows = master[master["Gene_Name"].isin(genes)]
    combo_id = "+".join(genes)

    means = {}
    for col in NUMERIC_MEAN_COLS:
        means[col] = float(np.nanmean(rows[col].values))

    member_scores = [
        float(rows.loc[rows["Gene_Name"] == g, "composite_score"].iloc[0])
        for g in genes if g in rows["Gene_Name"].values
    ]
    member_tiers = [
        str(rows.loc[rows["Gene_Name"] == g, "consensus_tier"].iloc[0])
        for g in genes if g in rows["Gene_Name"].values
    ]
    member_mods = [
        str(rows.loc[rows["Gene_Name"] == g, "best_modality"].iloc[0])
        for g in genes if g in rows["Gene_Name"].values
    ]

    majority_mod = get_majority_modality(member_mods)
    combo_tier = get_combo_consensus_tier(means["composite_score"])

    return {
        "Combo_ID": combo_id,
        "Genes": ",".join(genes),
        "Size": len(genes),
        "Combo_Composite_Score": round(means["composite_score"], 6),
        "Combo_Pharos_Mean": round(means["pharos_mean_proba"], 6),
        "Combo_Tier_Mean": round(means["tier_mean_proba"], 6),
        "Combo_Best_Modality_Proba": round(means["best_modality_proba"], 6),
        "Combo_Modality_Specificity": round(means["modality_specificity"], 6),
        "Combo_Avg_Runs_Top_Decile": round(means["n_runs_top_decile"], 2),
        "Combo_Avg_Runs_Above_75perc": round(means["n_runs_above_75perc"], 2),
        "Combo_Small_Mol_Mean": round(means["small_mol_proba"], 6),
        "Combo_Antibody_Mean": round(means["antibody_proba"], 6),
        "Combo_PROTAC_Mean": round(means["protac_proba"], 6),
        "Combo_Best_Modality": majority_mod,
        "Combo_Consensus_Tier": combo_tier,
        "Member_Scores": "|".join(f"{s:.4f}" for s in member_scores),
        "Member_Tiers": "|".join(member_tiers),
        "Member_Modalities": "|".join(member_mods),
    }


def load_existing_ids(path: Path, id_col: str) -> set:
    if not path.exists():
        return set()
    try:
        df = pd.read_csv(path, sep="\t", usecols=[id_col])
        return set(df[id_col].dropna().astype(str))
    except Exception:
        return set()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="DrugnomeAI individual + combo druggability scoring"
    )
    parser.add_argument(
        "--genes", type=str, required=True,
        help="Comma-separated gene names (e.g., 'TYK2,JAK1,IL17A')"
    )
    parser.add_argument(
        "--combos", type=str, nargs="*", default=None,
        help="Combo definitions as GENE1+GENE2+GENE3 (e.g., 'TYK2+JAK1' 'IL17A+IL17F+IL4R')"
    )
    parser.add_argument(
        "--master-table", type=str, default=str(DEFAULT_MASTER),
        help="Path to master druggability table CSV"
    )
    parser.add_argument(
        "--output-dir", type=str, default=str(DEFAULT_OUTPUT_DIR),
        help="Output directory for TSV files"
    )
    parser.add_argument(
        "--report", action="store_true", default=False,
        help="Flag placeholder — report generation is handled by the calling "
             "agent (Claude Code), not this script"
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    individual_tsv = output_dir / "drugnome_individual_score.tsv"
    combo_tsv = output_dir / "drugnome_combo_score.tsv"

    genes = parse_gene_list(args.genes)
    combos = parse_combos(args.combos) if args.combos else []

    # Also collect genes from combos
    combo_genes = {g for combo in combos for g in combo}
    all_genes = sorted(set(genes) | combo_genes)

    print("Loading master table ...")
    master = load_master_table(Path(args.master_table))
    found = set(master["Gene_Name"]) & set(all_genes)
    print(f"  {len(master)} genes in table, {len(found)}/{len(all_genes)} target genes found.")

    # --- Individual scores ---
    existing_genes = load_existing_ids(individual_tsv, "Gene_Name")
    missing_genes = sorted(set(all_genes) - existing_genes)

    if not missing_genes:
        print(f"\nIndividual: all {len(all_genes)} genes already in {individual_tsv.name}")
    else:
        print(f"\nIndividual: scoring {len(missing_genes)} genes ...")
        indiv_df = compute_individual_rows(master, missing_genes)

        if not existing_genes:
            indiv_df.to_csv(individual_tsv, sep="\t", index=False)
            print(f"  Created {individual_tsv.name} with {len(indiv_df)} genes")
        else:
            indiv_df.to_csv(
                individual_tsv, sep="\t", index=False, mode="a", header=False
            )
            print(f"  Appended {len(indiv_df)} genes to {individual_tsv.name}")

    # --- Combo scores (only if combos provided) ---
    if combos:
        existing_combos = load_existing_ids(combo_tsv, "Combo_ID")
        all_combo_ids = ["+".join(c) for c in combos]
        missing_combos = [
            (cid, cg) for cid, cg in zip(all_combo_ids, combos)
            if cid not in existing_combos
        ]

        if not missing_combos:
            print(f"\nCombo: all {len(combos)} combos already in {combo_tsv.name}")
        else:
            print(f"\nCombo: scoring {len(missing_combos)} combinations ...")
            combo_rows = [
                compute_combo_row(cg, master) for _, cg in missing_combos
            ]
            combo_df = pd.DataFrame(combo_rows)

            if not existing_combos:
                combo_df.to_csv(combo_tsv, sep="\t", index=False)
                print(f"  Created {combo_tsv.name} with {len(combo_df)} combos")
            else:
                combo_df.to_csv(
                    combo_tsv, sep="\t", index=False, mode="a", header=False
                )
                print(f"  Appended {len(combo_df)} combos to {combo_tsv.name}")

    # --- Summary ---
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    if individual_tsv.exists():
        idf = pd.read_csv(individual_tsv, sep="\t")
        print(f"\nIndividual: {len(idf)} genes -> {individual_tsv}")
        top = idf.nlargest(min(5, len(idf)), "composite_score")[
            ["Gene_Name", "composite_score", "consensus_tier", "best_modality"]
        ]
        for _, r in top.iterrows():
            print(
                f"  {r['Gene_Name']:12s}  composite={r['composite_score']:.4f}  "
                f"tier={r['consensus_tier']:8s}  modality={r['best_modality']}"
            )

    if combos and combo_tsv.exists():
        cdf = pd.read_csv(combo_tsv, sep="\t")
        print(f"\nCombo: {len(cdf)} combos -> {combo_tsv}")
        top = cdf.nlargest(min(5, len(cdf)), "Combo_Composite_Score")[
            ["Combo_ID", "Combo_Composite_Score", "Combo_Consensus_Tier", "Combo_Best_Modality"]
        ]
        for _, r in top.iterrows():
            print(
                f"  {r['Combo_ID']:40s}  composite={r['Combo_Composite_Score']:.4f}  "
                f"tier={r['Combo_Consensus_Tier']:8s}  modality={r['Combo_Best_Modality']}"
            )

    if args.report:
        print("\n[--report flag set: interpretation + markdown report should be "
              "generated by the calling agent (Claude Code / Codex)]")

    print("\nDone.")


if __name__ == "__main__":
    main()
