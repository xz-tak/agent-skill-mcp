#!/usr/bin/env python3
"""
Target Evaluation Prompt Generator - Enhanced CLI

Generates target evaluation documents by combining:
- Cortellis API data for targets
- Pathway analysis (KEGG, Reactome, MSigDB) via pathwaydb-query
- PPI data (STRING, IntAct, BioGRID) via interactdb-query (combo only)
- GPT-powered placeholder population

Output: Markdown prompt file (not DOCX)

Required User Inputs:
- targets: Target gene(s)
- indication: Primary indication
- related_diseases: Related diseases for context (user-defined, NOT auto-populated)

Usage:
    python target_evaluation_gen.py --targets "TYK2,JAK1" --indication "Crohn's disease" --related-diseases "ulcerative colitis, MASH, PSC"
    python target_evaluation_gen.py --targets "GLP2R" --indication "IBD" --related-diseases "Crohn's disease, ulcerative colitis, short bowel syndrome"
"""

import argparse
import asyncio
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import pandas as pd

# Paths
SKILL_ROOT = Path(__file__).parent.parent
# Use v2 template by default (comprehensive format from template_Target_Evaluation_Prompt_v2.docx)
MARKDOWN_TEMPLATE_PATH = SKILL_ROOT / "templates" / "prompt_template_v2.md"
# Fallback to original template if v2 doesn't exist
MARKDOWN_TEMPLATE_PATH_LEGACY = SKILL_ROOT / "templates" / "prompt_template.md"
GPT_PROMPTS_PATH = SKILL_ROOT / "templates" / "gpt_prompts.json"
CORTELLIS_SCRIPT = Path.home() / ".claude/skills/cortellis/cortellis_api/scripts/cortellis_gene_query.py"

# Pathway and interaction skill paths
PATHWAYDB_SCRIPT = Path.home() / ".claude/skills/pathwaydb-query/scripts"
INTERACTDB_SCRIPT = Path.home() / ".claude/skills/interactdb-query/scripts"

# Core placeholders that form info.txt
CORE_PLACEHOLDERS = [
    "target",
    "target_full_name",
    "indication",
    "modality_primary",
    "related_diseases",
    "function_summary",
    "canonical_pathway",
    "expression_cells",
    "effector_cells",
    "moa_rationale",
    "combo_opportunities",
    "provided_sources",
    # New placeholders for pathway/interaction data
    "shared_pathways",
    "unique_pathways_T1",
    "unique_pathways_T2",
    "interaction_path",
    "bridge_proteins",
    "combo_effect_type",
    "pathway_count_shared",
]


def sanitize_name(name: str) -> str:
    """Sanitize name for filesystem use."""
    return re.sub(r'[^\w\-]', '_', name.replace("'", "").replace(" ", "_"))


def parse_targets(targets_str: str) -> List[str]:
    """Parse comma-separated target string into list."""
    return [t.strip().upper() for t in targets_str.split(",") if t.strip()]


def get_output_dir_name(targets: List[str], indication: str) -> str:
    """Generate output directory name."""
    targets_part = "_".join(targets)
    indication_part = sanitize_name(indication)
    return f"{targets_part}_{indication_part}"


def run_cortellis_query(target: str, output_dir: Path) -> Dict:
    """Run Cortellis API query for a single target."""
    print(f"[Cortellis] Querying {target}...")

    cmd = [
        sys.executable,
        str(CORTELLIS_SCRIPT),
        target,
        "--all",
        "--excel",
        "--output-dir", str(output_dir)
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )

        if result.returncode != 0:
            print(f"[Cortellis] Warning: Query for {target} returned non-zero: {result.stderr}")

        # Load the generated JSON
        json_path = output_dir / f"{target}_cortellis_data.json"
        if json_path.exists():
            with open(json_path) as f:
                return json.load(f)
        else:
            print(f"[Cortellis] Warning: No JSON output for {target}")
            return {}

    except subprocess.TimeoutExpired:
        print(f"[Cortellis] Error: Query for {target} timed out")
        return {}
    except Exception as e:
        print(f"[Cortellis] Error querying {target}: {e}")
        return {}


def run_codex_gpt(prompt: str, timeout: int = 180, max_retries: int = 3, reasoning: str = "medium") -> str:
    """Run codex CLI to generate content via GPT with retry logic.

    Args:
        prompt: The prompt to send to GPT
        timeout: Timeout in seconds for each attempt
        max_retries: Maximum number of retry attempts
        reasoning: Reasoning effort level (low/medium/high/xhigh)

    Returns:
        GPT response string or "[TO BE PROVIDED]" if all retries fail
    """
    cmd = [
        "codex", "exec",
        "--skip-git-repo-check",
        "--sandbox", "read-only",
        "-c", f"model_reasoning_effort={reasoning}"
    ]

    for attempt in range(max_retries):
        try:
            result = subprocess.run(
                cmd,
                input=prompt,
                capture_output=True,
                text=True,
                timeout=timeout
            )

            if result.returncode == 0 and result.stdout.strip():
                return result.stdout.strip()

            # Empty response or non-zero return code - retry
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt  # Exponential backoff: 1s, 2s, 4s
                print(f"[Codex] Attempt {attempt + 1} failed, retrying in {wait_time}s...")
                time.sleep(wait_time)
                continue
            else:
                print(f"[Codex] Warning after {max_retries} attempts: {result.stderr}")

        except subprocess.TimeoutExpired:
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
                print(f"[Codex] Timeout on attempt {attempt + 1}, retrying in {wait_time}s...")
                time.sleep(wait_time)
                continue
            else:
                print(f"[Codex] Error: GPT request timed out after {max_retries} attempts")

        except FileNotFoundError:
            print("[Codex] Error: codex CLI not found, falling back to placeholder")
            return "[TO BE PROVIDED]"

        except Exception as e:
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
                print(f"[Codex] Error on attempt {attempt + 1}: {e}, retrying in {wait_time}s...")
                time.sleep(wait_time)
                continue
            else:
                print(f"[Codex] Error after {max_retries} attempts: {e}")

    return "[TO BE PROVIDED]"  # All retries exhausted


def load_gpt_prompts() -> Dict[str, str]:
    """Load GPT prompt templates."""
    if GPT_PROMPTS_PATH.exists():
        with open(GPT_PROMPTS_PATH) as f:
            return json.load(f)
    return {}


def extract_from_cortellis(data: Dict, field: str) -> str:
    """Extract specific field from Cortellis data."""
    if not data:
        return ""

    # Navigation paths for common fields
    # Note: related_diseases is now user-provided (required input), not extracted from Cortellis
    paths = {
        "target_full_name": ["annotation", "target_name"],
        "function_summary": ["annotation", "function_description"],
        "canonical_pathway": ["annotation", "pathway_info"],
        "expression_cells": ["annotation", "expression_pattern"],
    }

    if field in paths:
        current = data
        for key in paths[field]:
            if isinstance(current, dict) and key in current:
                current = current[key]
            else:
                return ""
        return str(current) if current else ""

    return ""


def query_pathways(targets: List[str], output_dir: Path) -> Dict:
    """
    Query pathway databases for all targets using pathwaydb-query skill.

    Returns pathway results and saves CSV files.
    """
    print("\n[Pathways] Querying KEGG, Reactome, MSigDB...")

    # Add pathwaydb-query scripts to path
    if str(PATHWAYDB_SCRIPT) not in sys.path:
        sys.path.insert(0, str(PATHWAYDB_SCRIPT))

    try:
        from multi_gene_analysis_async import query_multiple_genes_async, export_results_to_tables

        # Query all databases
        pathway_results = asyncio.run(query_multiple_genes_async(
            genes=targets,
            kegg_organism="hsa",
            msigdb_collections=["H", "C2"],
            reactome_species=9606,
            network_compatible_names=True
        ))

        # Export results - pass working_dir to keep outputs in main directory
        output_prefix = str(output_dir / "pathways")
        export_results_to_tables(
            gene_results=pathway_results,
            output_prefix=output_prefix,
            genes=targets,
            working_dir=str(output_dir)  # Ensure pathwaydb outputs go to main output dir
        )

        print(f"[Pathways] Saved to {output_dir}/")
        return pathway_results

    except ImportError as e:
        print(f"[Pathways] Warning: Could not import pathwaydb-query: {e}")
        return {}
    except Exception as e:
        print(f"[Pathways] Error: {e}")
        return {}


def analyze_pathway_overlap(pathway_results: Dict, targets: List[str]) -> Dict:
    """
    Analyze shared and unique pathways for each target.

    Args:
        pathway_results: {gene: {database: set(pathways)}}
        targets: List of target genes

    Returns:
        Dict with shared_pathways, unique_pathways, pathway_counts
    """
    if not pathway_results:
        return {
            "shared_pathways": [],
            "unique_pathways": {t: [] for t in targets},
            "pathway_counts": {t: 0 for t in targets}
        }

    # Combine pathways across databases for each gene
    gene_pathways = {}
    for gene in targets:
        if gene in pathway_results:
            all_pathways = set()
            for db_name, pathways in pathway_results[gene].items():
                all_pathways.update(pathways)
            gene_pathways[gene] = all_pathways
        else:
            gene_pathways[gene] = set()

    # Find shared pathways (present in ALL targets)
    if len(targets) >= 2 and all(gene_pathways.get(t) for t in targets):
        shared = set.intersection(*[gene_pathways[t] for t in targets])
    else:
        shared = set()

    # Find unique pathways per target
    unique = {}
    for target in targets:
        target_pathways = gene_pathways.get(target, set())
        other_pathways = set()
        for other_target in targets:
            if other_target != target:
                other_pathways.update(gene_pathways.get(other_target, set()))
        unique[target] = list(target_pathways - other_pathways)[:20]

    return {
        "shared_pathways": list(shared)[:30],
        "unique_pathways": unique,
        "pathway_counts": {t: len(gene_pathways.get(t, [])) for t in targets}
    }


def query_interactions(targets: List[str], output_dir: Path) -> Dict:
    """
    Query PPI databases for shortest paths between combo targets.

    Only called for multi-target (combo) evaluations.
    """
    if len(targets) < 2:
        return {}

    print("\n[Interactions] Finding shortest paths in STRING, IntAct, BioGRID...")

    # Add interactdb-query scripts to path
    if str(INTERACTDB_SCRIPT) not in sys.path:
        sys.path.insert(0, str(INTERACTDB_SCRIPT))

    try:
        from unified_query import query_shortest_paths_all_databases

        # Find shortest paths between targets
        interaction_results = query_shortest_paths_all_databases(
            gene_list=targets,
            species=9606,
            max_distance=50,
            export_results=True,
            output_dir=str(output_dir)
        )

        print(f"[Interactions] Completed path queries")
        return interaction_results

    except ImportError as e:
        print(f"[Interactions] Warning: Could not import interactdb-query: {e}")
        return {}
    except Exception as e:
        print(f"[Interactions] Error: {e}")
        return {}


def extract_interaction_summary(interaction_results: Dict) -> Dict:
    """
    Extract summary information from interaction query results.

    Returns:
        Dict with direct_interaction, shortest_path, interaction_scores, bridge_proteins
    """
    summary = {
        "direct_interaction": False,
        "shortest_path": {
            "path": [],
            "hops": None,
            "bridge_proteins": []
        },
        "interaction_scores": {}
    }

    if not interaction_results:
        return summary

    # Check each database for paths
    best_path = None
    best_hops = float('inf')

    for db_name in ['string', 'intact', 'biogrid']:
        db_result = interaction_results.get(db_name)
        if isinstance(db_result, str):  # Error message
            continue
        if not db_result:
            continue

        # db_result is Dict of paths: {(geneA, geneB): {path, hops, distance, scores}}
        for pair_key, path_info in db_result.items():
            if isinstance(path_info, dict):
                hops = path_info.get('hops', float('inf'))
                path = path_info.get('path', [])

                # Check for direct interaction
                if hops == 1:
                    summary["direct_interaction"] = True

                # Track best (shortest) path
                if hops < best_hops:
                    best_hops = hops
                    best_path = path

                # Extract scores
                scores = path_info.get('scores', [])
                valid_scores = [s for s in scores if s is not None]
                if valid_scores:
                    avg_score = sum(valid_scores) / len(valid_scores)
                    summary["interaction_scores"][db_name.upper()] = round(avg_score, 2)

    # Set shortest path info
    if best_path:
        summary["shortest_path"]["path"] = best_path
        summary["shortest_path"]["hops"] = best_hops

        # Extract bridge proteins (intermediate nodes)
        if len(best_path) > 2:
            summary["shortest_path"]["bridge_proteins"] = best_path[1:-1]

    return summary


def interpret_combo_effect(
    pathway_analysis: Dict,
    interaction_summary: Dict,
    targets: List[str],
    indication: str
) -> Dict:
    """
    Determine combo effect type combining real data + GPT interpretation.

    Args:
        pathway_analysis: Result from analyze_pathway_overlap()
        interaction_summary: Result from extract_interaction_summary()
        targets: List of target genes
        indication: Disease indication

    Returns:
        Dict with data_effect_type, gpt_interpretation, final_effect_type, evidence
    """
    shared_pathways = pathway_analysis.get("shared_pathways", [])
    has_direct_interaction = interaction_summary.get("direct_interaction", False)
    bridge_proteins = interaction_summary.get("shortest_path", {}).get("bridge_proteins", [])

    # Initial classification from real data
    if len(shared_pathways) > 5 and has_direct_interaction:
        data_effect_type = "synergistic"
    elif len(shared_pathways) > 2 or bridge_proteins:
        data_effect_type = "complementary"
    else:
        data_effect_type = "additive"

    # Build evidence summary
    evidence = {
        "shared_pathway_count": len(shared_pathways),
        "has_direct_interaction": has_direct_interaction,
        "bridge_protein_count": len(bridge_proteins)
    }

    # GPT synthesis to interpret and validate the effect type
    gpt_prompt = f"""Analyze combo therapy potential for {', '.join(targets)} in {indication}.

Real data findings:
- Shared pathways ({len(shared_pathways)}): {', '.join(shared_pathways[:5]) if shared_pathways else 'None found'}
- Direct interaction: {has_direct_interaction}
- Bridge proteins: {', '.join(bridge_proteins[:5]) if bridge_proteins else 'None'}
- Data-based effect type: {data_effect_type}

Based on this data, provide:
1. Confirmed effect type (synergistic/complementary/additive)
2. Scientific rationale (2-3 sentences)
3. Key mechanistic insight

Format as:
EFFECT TYPE: <type>
RATIONALE: <rationale>
MECHANISTIC INSIGHT: <insight>"""

    print("[GPT] Generating combo effect interpretation...")
    gpt_response = run_codex_gpt(gpt_prompt)

    # Parse GPT response
    gpt_interpretation = {
        "raw_response": gpt_response,
        "rationale": "",
        "mechanistic_insight": ""
    }

    if gpt_response and gpt_response != "[TO BE PROVIDED]":
        lines = gpt_response.split('\n')
        for line in lines:
            if line.startswith("RATIONALE:"):
                gpt_interpretation["rationale"] = line.replace("RATIONALE:", "").strip()
            elif line.startswith("MECHANISTIC INSIGHT:"):
                gpt_interpretation["mechanistic_insight"] = line.replace("MECHANISTIC INSIGHT:", "").strip()

    return {
        "data_effect_type": data_effect_type,
        "gpt_interpretation": gpt_interpretation,
        "final_effect_type": data_effect_type,  # Can be overridden by GPT
        "evidence": evidence
    }


def populate_placeholders(
    targets: List[str],
    indication: str,
    cortellis_data: Dict[str, Dict],
    pathway_analysis: Dict = None,
    interaction_summary: Dict = None,
    combo_effect: Dict = None,
    related_diseases: str = None
) -> Dict[str, str]:
    """Populate all placeholders with values."""

    replacements = {}
    gpt_prompts = load_gpt_prompts()

    # Set core values from input (user-provided required fields)
    replacements["target"] = ", ".join(targets) if len(targets) > 1 else targets[0]
    replacements["indication"] = indication
    replacements["related_diseases"] = related_diseases if related_diseases else "[TO BE PROVIDED]"
    replacements["provided_sources"] = "provided files in local directory: *.json,xlsx,csv"

    # Merge Cortellis data for all targets
    merged_context = []
    for target, data in cortellis_data.items():
        if data:
            merged_context.append(f"Target {target}: {json.dumps(data.get('annotation', {}), indent=2)[:2000]}")

    context_str = "\n".join(merged_context)

    # Add pathway analysis placeholders
    if pathway_analysis:
        shared = pathway_analysis.get("shared_pathways", [])
        unique = pathway_analysis.get("unique_pathways", {})
        counts = pathway_analysis.get("pathway_counts", {})

        replacements["shared_pathways"] = ", ".join(shared[:10]) if shared else "No shared pathways identified"
        replacements["pathway_count_shared"] = str(len(shared))

        # Unique pathways for each target
        if len(targets) >= 1:
            replacements["unique_pathways_T1"] = ", ".join(unique.get(targets[0], [])[:5]) or "None identified"
        if len(targets) >= 2:
            replacements["unique_pathways_T2"] = ", ".join(unique.get(targets[1], [])[:5]) or "None identified"
        else:
            # Single target - no T2
            replacements["unique_pathways_T2"] = "N/A (single target)"
    else:
        replacements["shared_pathways"] = "[Pathway query failed]"
        replacements["pathway_count_shared"] = "0"
        replacements["unique_pathways_T1"] = "[Pathway query failed]"
        replacements["unique_pathways_T2"] = "N/A (single target)" if len(targets) < 2 else "[Pathway query failed]"

    # Add interaction analysis placeholders
    if interaction_summary:
        path = interaction_summary.get("shortest_path", {}).get("path", [])
        bridge = interaction_summary.get("shortest_path", {}).get("bridge_proteins", [])

        replacements["interaction_path"] = " → ".join(path) if path else "No direct path found"
        replacements["bridge_proteins"] = ", ".join(bridge) if bridge else "None"
    else:
        # Single target - combo fields not applicable
        replacements["interaction_path"] = "N/A (single target)"
        replacements["bridge_proteins"] = "N/A (single target)"

    # Add combo effect placeholders
    if combo_effect:
        replacements["combo_effect_type"] = combo_effect.get("final_effect_type", "unknown")
    else:
        # Single target - combo effect not applicable
        replacements["combo_effect_type"] = "N/A (single target)"

    # Populate remaining core placeholders
    for placeholder in CORE_PLACEHOLDERS:
        if placeholder in replacements:
            continue

        # Try to extract from Cortellis data first
        for target, data in cortellis_data.items():
            extracted = extract_from_cortellis(data, placeholder)
            if extracted:
                replacements[placeholder] = extracted
                break

        # If not found, generate via codex
        if placeholder not in replacements:
            prompt_template = gpt_prompts.get(placeholder,
                f"For target(s) {', '.join(targets)} in {indication}, provide a concise value for {placeholder}."
            )
            prompt = prompt_template.format(
                targets=", ".join(targets),
                indication=indication,
                context=context_str[:3000]
            )

            print(f"[GPT] Generating {placeholder}...")
            replacements[placeholder] = run_codex_gpt(prompt)

        # Fallback
        if placeholder not in replacements:
            replacements[placeholder] = "[TO BE PROVIDED]"

    return replacements


def analyze_combo_effects(
    targets: List[str],
    indication: str,
    cortellis_data: Dict[str, Dict],
    pathway_analysis: Dict,
    interaction_summary: Dict
) -> Dict:
    """
    Generate enhanced combo_analysis.json with real data + GPT interpretation.
    """
    if len(targets) < 2:
        return {}

    # Extract pathways and indications for each target from Cortellis
    target_pathways = {}
    target_indications = {}

    for target, data in cortellis_data.items():
        pathways = []
        indications = []

        if data:
            annotation = data.get("annotation", {})
            if "pathway_info" in annotation:
                pathways.extend(annotation["pathway_info"] if isinstance(annotation["pathway_info"], list) else [annotation["pathway_info"]])

            drugs = data.get("drug", {}).get("Drug", [])
            for drug in drugs[:10]:
                if "indication" in drug:
                    indications.append(drug["indication"])

        target_pathways[target] = pathways
        target_indications[target] = indications

    # Find shared pathways from Cortellis
    if len(target_pathways) >= 2:
        pathway_sets = [set(p) for p in target_pathways.values() if p]
        cortellis_shared_pathways = list(set.intersection(*pathway_sets)) if pathway_sets else []
    else:
        cortellis_shared_pathways = []

    # Find shared indications
    if len(target_indications) >= 2:
        indication_sets = [set(i) for i in target_indications.values() if i]
        shared_indications = list(set.intersection(*indication_sets)) if indication_sets else []
    else:
        shared_indications = []

    # Get combo effect interpretation
    combo_effect = interpret_combo_effect(
        pathway_analysis,
        interaction_summary,
        targets,
        indication
    )

    # Build enhanced analysis structure
    analysis = {
        "targets": targets,
        "indication": indication,
        "effect_classification": {
            "data_effect_type": combo_effect["data_effect_type"],
            "final_effect_type": combo_effect["final_effect_type"],
            "evidence": combo_effect["evidence"]
        },
        "pathway_analysis": {
            "shared_pathways": pathway_analysis.get("shared_pathways", [])[:20],
            "unique_pathways": pathway_analysis.get("unique_pathways", {}),
            "pathway_counts": pathway_analysis.get("pathway_counts", {}),
            "cortellis_shared_pathways": cortellis_shared_pathways[:10]
        },
        "interaction_analysis": {
            "direct_interaction": interaction_summary.get("direct_interaction", False),
            "shortest_path": interaction_summary.get("shortest_path", {}),
            "interaction_scores": interaction_summary.get("interaction_scores", {})
        },
        "gpt_interpretation": {
            "rationale": combo_effect["gpt_interpretation"].get("rationale", ""),
            "mechanistic_insight": combo_effect["gpt_interpretation"].get("mechanistic_insight", "")
        },
        "shared_indications": shared_indications[:10],
        "unique_mechanisms": {
            target: list(set(target_pathways.get(target, [])) - set(cortellis_shared_pathways))[:5]
            for target in targets
        }
    }

    return analysis


def save_info_txt(replacements: Dict[str, str], output_path: Path) -> None:
    """Save replacements to info.txt in TSV format."""
    with open(output_path, "w") as f:
        for key, value in replacements.items():
            # Escape newlines in values
            escaped_value = str(value).replace("\n", " ").replace("\t", " ")
            f.write(f"${{{key}}}\t{escaped_value}\n")

    print(f"[Info] Saved to {output_path}")


def validate_placeholders(info_path: Path) -> Dict[str, Any]:
    """Validate that all placeholders are populated in info.txt.

    Returns:
        Dict with:
        - passed: bool
        - total_fields: int
        - populated_fields: int
        - failed_fields: List[str] - fields still showing [TO BE PROVIDED]
        - completion_rate: float (0.0 to 1.0)
    """
    failed_fields = []
    total_fields = 0
    populated_fields = 0

    with open(info_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or not line.startswith('${'):
                continue

            parts = line.split('\t', 1)
            if len(parts) != 2:
                continue

            placeholder_name = parts[0]  # e.g., ${target}
            value = parts[1]
            total_fields += 1

            # Check for unpopulated values
            if value in ['[TO BE PROVIDED]', '']:
                failed_fields.append(placeholder_name)
            elif 'N/A (single target)' in value:
                # N/A is acceptable for combo-only fields on single targets
                populated_fields += 1
            else:
                populated_fields += 1

    completion_rate = populated_fields / total_fields if total_fields > 0 else 0.0

    return {
        'passed': len(failed_fields) == 0,
        'total_fields': total_fields,
        'populated_fields': populated_fields,
        'failed_fields': failed_fields,
        'completion_rate': completion_rate
    }


def retry_failed_placeholders(
    info_path: Path,
    failed_fields: List[str],
    gpt_prompts: Dict[str, str],
    context: Dict[str, str],
    max_qc_rounds: int = 5
) -> Dict[str, Any]:
    """Retry GPT calls for failed placeholders until 100% completion.

    Args:
        info_path: Path to info.txt
        failed_fields: List of placeholder names that failed (e.g., "${expression_cells}")
        gpt_prompts: Dict of placeholder_key -> prompt template
        context: Dict with targets, indication, context for prompt formatting
        max_qc_rounds: Maximum retry rounds (default: 5 for 100% target)

    Returns:
        Dict with retry results
    """
    results = {
        'rounds_run': 0,
        'fields_recovered': [],
        'fields_still_failed': list(failed_fields)
    }

    for round_num in range(max_qc_rounds):
        if not results['fields_still_failed']:
            print(f"[QC] ✓ 100% completion achieved after {round_num} rounds")
            break

        results['rounds_run'] = round_num + 1
        print(f"[QC] Retry round {round_num + 1}/{max_qc_rounds} for {len(results['fields_still_failed'])} failed fields")

        # Read current info.txt
        with open(info_path, 'r') as f:
            lines = f.readlines()

        # Retry each failed field
        newly_recovered = []
        for field in results['fields_still_failed']:
            # Get prompt template for this field
            field_key = field.replace('${', '').replace('}', '')
            if field_key not in gpt_prompts:
                print(f"[QC] No prompt template for {field_key}, skipping")
                continue

            prompt_template = gpt_prompts[field_key]

            # Format prompt with context
            try:
                prompt = prompt_template.format(**context)
            except KeyError as e:
                print(f"[QC] Missing context key for {field_key}: {e}")
                continue

            # Try GPT with fresh call (increased timeout for retry)
            print(f"[QC] Retrying {field_key}...")
            response = run_codex_gpt(prompt, timeout=240, max_retries=2)

            if response and response != '[TO BE PROVIDED]':
                # Update the line in info.txt
                for i, line in enumerate(lines):
                    if line.startswith(field + '\t'):
                        lines[i] = f"{field}\t{response.replace(chr(10), ' ').replace(chr(9), ' ')}\n"
                        newly_recovered.append(field)
                        print(f"[QC] ✓ Recovered {field_key}")
                        break

        # Write updated info.txt
        with open(info_path, 'w') as f:
            f.writelines(lines)

        # Update tracking
        results['fields_recovered'].extend(newly_recovered)
        results['fields_still_failed'] = [
            f for f in results['fields_still_failed']
            if f not in newly_recovered
        ]

        # Add delay between rounds to let API recover (skip delay on last round or if all recovered)
        if results['fields_still_failed'] and round_num < max_qc_rounds - 1:
            print(f"[QC] Waiting 30s before next retry round to let API recover...")
            time.sleep(30)

    return results


def generate_qc_report(
    output_dir: Path,
    validation_result: Dict[str, Any],
    retry_result: Dict[str, Any],
    targets: List[str],
    indication: str
) -> Path:
    """Generate QC report as JSON.

    Returns:
        Path to QC report file
    """
    total_fields = validation_result['total_fields']
    final_failed = len(retry_result.get('fields_still_failed', []))
    final_populated = total_fields - final_failed

    report = {
        'timestamp': datetime.now().isoformat(),
        'targets': targets,
        'indication': indication,
        'initial_validation': {
            'total_fields': validation_result['total_fields'],
            'populated_fields': validation_result['populated_fields'],
            'completion_rate': f"{validation_result['completion_rate']:.1%}",
            'failed_fields': validation_result['failed_fields']
        },
        'retry_pass': {
            'rounds_run': retry_result.get('rounds_run', 0),
            'fields_recovered': retry_result.get('fields_recovered', []),
            'fields_still_failed': retry_result.get('fields_still_failed', [])
        },
        'final_status': {
            'passed': final_failed == 0,
            'final_completion_rate': f"{final_populated / total_fields:.1%}" if total_fields > 0 else "N/A",
            'final_populated_fields': final_populated,
            'final_failed_fields': final_failed
        }
    }

    report_path = output_dir / 'qc_report.json'
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2)

    return report_path


def load_info_txt_as_replacements(info_path: Path) -> Dict[str, str]:
    """Load info.txt back into a replacements dict for markdown regeneration."""
    replacements = {}
    with open(info_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or not line.startswith('${'):
                continue
            parts = line.split('\t', 1)
            if len(parts) == 2:
                key = parts[0].replace('${', '').replace('}', '')
                replacements[key] = parts[1]
    return replacements


def generate_markdown_prompt(replacements: Dict[str, str], output_path: Path, is_combo: bool = False) -> Path:
    """
    Generate markdown prompt file from template and replacements.

    Args:
        replacements: Dict of placeholder -> value
        output_path: Path for output markdown file
        is_combo: Whether this is a combo target evaluation

    Returns:
        Path to generated markdown file
    """
    # Load template - prefer v2, fallback to legacy, then default
    if MARKDOWN_TEMPLATE_PATH.exists():
        template = MARKDOWN_TEMPLATE_PATH.read_text()
        print(f"[Markdown] Using v2 template: {MARKDOWN_TEMPLATE_PATH}")
    elif MARKDOWN_TEMPLATE_PATH_LEGACY.exists():
        template = MARKDOWN_TEMPLATE_PATH_LEGACY.read_text()
        print(f"[Markdown] Using legacy template: {MARKDOWN_TEMPLATE_PATH_LEGACY}")
    else:
        # Default template
        template = get_default_markdown_template(is_combo)
        print("[Markdown] Using built-in default template")

    # Replace all placeholders
    pattern = re.compile(r'\$\{([^}]+)\}')

    def replace_match(match):
        key = match.group(1)
        return str(replacements.get(key, f"[{key}]"))

    content = pattern.sub(replace_match, template)

    # Write output
    output_path.write_text(content)
    print(f"[Markdown] Saved to {output_path}")

    return output_path


def get_default_markdown_template(is_combo: bool = False) -> str:
    """Return default markdown template."""
    combo_section = """
## Combo Analysis

### Effect Classification
- **Effect Type:** ${combo_effect_type}
- **Evidence:**
  - Shared pathways: ${pathway_count_shared}
  - Direct interaction: ${interaction_path}
  - Bridge proteins: ${bridge_proteins}

### Interaction Path
${interaction_path}

### Bridge Proteins
${bridge_proteins}
""" if is_combo else ""

    return f"""# Target Evaluation: ${{target}}

## Overview
- **Target(s):** ${{target}}
- **Full Name:** ${{target_full_name}}
- **Indication:** ${{indication}}
- **Primary Modality:** ${{modality_primary}}

## Biological Context

### Function Summary
${{function_summary}}

### Canonical Pathway
${{canonical_pathway}}

### Expression Pattern
${{expression_cells}}

### Effector Cells
${{effector_cells}}

## Pathway Analysis

### Shared Pathways (${{pathway_count_shared}})
${{shared_pathways}}

### Unique Pathways
**${{target}} (Target 1):** ${{unique_pathways_T1}}
**Target 2:** ${{unique_pathways_T2}}
{combo_section}
## MOA Rationale
${{moa_rationale}}

## Combination Opportunities
${{combo_opportunities}}

## Decision Framework

### Claims
1. ${{claim_1}}
2. ${{claim_2}}

### Evidence
- ${{evidence_1a}}
- ${{evidence_1b}}
- ${{evidence_2a}}
- ${{evidence_2b}}

### Counterarguments
- ${{counter_1}}
- ${{counter_2}}

### Reconciliation
${{how_we_reconcile}}

### Confidence Assessment
- Claim 1: ${{confidence_1}}
- Claim 2: ${{confidence_2}}

### Key Assumptions
${{assump_1}}

## Risk Assessment

| Failure Mode | Impact | Mitigation |
|--------------|--------|------------|
| ${{failure_mode_1}} | ${{impact_1}} | ${{mitigation_1}} |
| ${{failure_mode_2}} | ${{impact_2}} | ${{mitigation_2}} |

## Options Analysis

### Option A (Base Case)
${{option_A}}

### Option B (Alternative)
${{option_B}}

### Option C (De-risked)
${{option_C}}

## Logic Chain
1. ${{logic_1}}
2. ${{logic_2}}

## Validation
- **Key Test:** ${{test_1}}
- **What Would Change Our Mind:** ${{flip_1}}

## Executive Summary
${{one_paragraph_why}}

## Decision

**Statement:** ${{decision_statement}}

**Top 3 Risks:**
${{top3_risks}}

**Next Actions:**
${{next_actions}}

---
*Generated by /agenticboost-prompt*
*Sources: ${{provided_sources}}*
*Review ID: ${{dr_id}} | Date: ${{when}} | Reviewer: ${{who}}*
"""


def main():
    parser = argparse.ArgumentParser(
        description="Generate target evaluation documents with pathway and interaction data"
    )
    parser.add_argument(
        "--targets", "-t",
        required=True,
        help="Comma-separated list of targets (e.g., 'TYK2,JAK1')"
    )
    parser.add_argument(
        "--indication", "-i",
        required=True,
        help="Primary indication (e.g., 'Crohn\\'s disease')"
    )
    parser.add_argument(
        "--output-dir", "-o",
        help="Output directory (default: auto-generated)"
    )
    parser.add_argument(
        "--template",
        help="Custom markdown template path"
    )
    parser.add_argument(
        "--related-diseases", "-r",
        required=True,
        help="Related diseases, comma-separated (e.g., 'Crohn\\'s disease, ulcerative colitis, MASH')"
    )

    args = parser.parse_args()

    # Parse inputs
    targets = parse_targets(args.targets)
    indication = args.indication
    related_diseases = args.related_diseases
    is_combo = len(targets) > 1

    print(f"\n{'='*60}")
    print(f"Target Evaluation Generator (Enhanced)")
    print(f"{'='*60}")
    print(f"Targets: {targets}")
    print(f"Indication: {indication}")
    print(f"Related diseases: {related_diseases}")
    print(f"Mode: {'Combo' if is_combo else 'Single Target'}")
    print(f"{'='*60}\n")

    # Setup output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = Path.cwd() / get_output_dir_name(targets, indication)

    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"[Setup] Output directory: {output_dir}")

    # Step 1: Run Cortellis queries in parallel
    print("\n[Phase 1] Running Cortellis queries...")
    cortellis_data = {}

    with ThreadPoolExecutor(max_workers=len(targets)) as executor:
        futures = {
            executor.submit(run_cortellis_query, target, output_dir): target
            for target in targets
        }

        for future in as_completed(futures):
            target = futures[future]
            try:
                cortellis_data[target] = future.result()
            except Exception as e:
                print(f"[Error] Cortellis query failed for {target}: {e}")
                cortellis_data[target] = {}

    # Step 2: Run pathway queries
    print("\n[Phase 2] Running pathway queries...")
    pathway_results = query_pathways(targets, output_dir)
    pathway_analysis = analyze_pathway_overlap(pathway_results, targets)

    # Step 3: Run interaction queries (combo only)
    interaction_summary = {}
    if is_combo:
        print("\n[Phase 3] Running interaction queries...")
        interaction_results = query_interactions(targets, output_dir)
        interaction_summary = extract_interaction_summary(interaction_results)
    else:
        print("\n[Phase 3] Skipping interaction queries (single target)")

    # Step 4: Analyze combo effects (combo only)
    combo_analysis = {}
    combo_effect = {}
    if is_combo:
        print("\n[Phase 4] Analyzing combo effects...")
        combo_effect = interpret_combo_effect(
            pathway_analysis, interaction_summary, targets, indication
        )
        combo_analysis = analyze_combo_effects(
            targets, indication, cortellis_data, pathway_analysis, interaction_summary
        )

        # Save combo_analysis.json
        combo_path = output_dir / "combo_analysis.json"
        with open(combo_path, "w") as f:
            json.dump(combo_analysis, f, indent=2)
        print(f"[Combo] Saved to {combo_path}")
    else:
        print("\n[Phase 4] Skipping combo analysis (single target)")

    # Step 5: Populate placeholders
    print("\n[Phase 5] Populating placeholders...")
    replacements = populate_placeholders(
        targets, indication, cortellis_data,
        pathway_analysis, interaction_summary, combo_effect if is_combo else None,
        related_diseases=related_diseases
    )

    # Step 6: Save info.txt
    targets_str = "_".join(targets)
    info_path = output_dir / f"{targets_str}_info.txt"
    save_info_txt(replacements, info_path)

    # Step 7: Generate markdown prompt
    print("\n[Phase 6] Generating markdown prompt...")
    indication_sanitized = sanitize_name(indication)

    # Use custom template if provided
    if args.template:
        custom_template = Path(args.template)
        if custom_template.exists():
            template_content = custom_template.read_text()
            # Save as the active template
            MARKDOWN_TEMPLATE_PATH.parent.mkdir(parents=True, exist_ok=True)
            MARKDOWN_TEMPLATE_PATH.write_text(template_content)

    output_md = output_dir / f"{targets_str}_({indication_sanitized})_Target_Evaluation_Prompt.md"
    generate_markdown_prompt(replacements, output_md, is_combo)

    # === QC PASS ===
    print(f"\n{'='*60}")
    print("[QC] Validating placeholder population...")
    print(f"{'='*60}")

    validation = validate_placeholders(info_path)

    print(f"[QC] Initial completion: {validation['completion_rate']:.1%} ({validation['populated_fields']}/{validation['total_fields']} fields)")

    if not validation['passed']:
        print(f"[QC] {len(validation['failed_fields'])} fields need retry:")
        for field in validation['failed_fields']:
            print(f"      - {field}")

        # Load GPT prompts for retry
        gpt_prompts = load_gpt_prompts()

        # Merge Cortellis context for retry prompts
        merged_context = []
        for target, data in cortellis_data.items():
            if data:
                merged_context.append(f"Target {target}: {json.dumps(data.get('annotation', {}), indent=2)[:2000]}")
        context_str = "\n".join(merged_context)

        # Build context for prompt formatting
        prompt_context = {
            'targets': ", ".join(targets),
            'target': targets[0],
            'target2': targets[1] if len(targets) > 1 else '',
            'indication': indication,
            'context': context_str[:3000]
        }

        # Retry failed fields until 100% or max 5 rounds
        retry_result = retry_failed_placeholders(
            info_path=info_path,
            failed_fields=validation['failed_fields'],
            gpt_prompts=gpt_prompts,
            context=prompt_context,
            max_qc_rounds=5
        )

        # Re-validate after retry
        final_validation = validate_placeholders(info_path)
        print(f"\n[QC] Final completion: {final_validation['completion_rate']:.1%} ({final_validation['populated_fields']}/{final_validation['total_fields']} fields)")

        # Regenerate markdown with updated info if fields were recovered
        if retry_result['fields_recovered']:
            print("[QC] Regenerating markdown with recovered fields...")
            updated_replacements = load_info_txt_as_replacements(info_path)
            generate_markdown_prompt(updated_replacements, output_md, is_combo)
    else:
        retry_result = {'rounds_run': 0, 'fields_recovered': [], 'fields_still_failed': []}
        final_validation = validation

    # Generate QC report
    qc_report_path = generate_qc_report(
        output_dir=output_dir,
        validation_result=validation,
        retry_result=retry_result,
        targets=targets,
        indication=indication
    )
    print(f"[QC] Report saved: {qc_report_path}")

    # Final status
    print(f"\n{'='*60}")
    if not retry_result.get('fields_still_failed'):
        print("✓ QC PASSED: 100% completion achieved")
    else:
        print(f"⚠️  QC INCOMPLETE: {len(retry_result['fields_still_failed'])} fields could not be populated after {retry_result['rounds_run']} rounds:")
        for field in retry_result['fields_still_failed']:
            print(f"   - {field}")
    print(f"{'='*60}")

    # Summary
    print(f"\n{'='*60}")
    print("Generation Complete!")
    print(f"{'='*60}")
    print(f"Output directory: {output_dir}")
    print(f"Files generated:")
    for f in sorted(output_dir.iterdir()):
        print(f"  - {f.name}")
    print(f"{'='*60}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
