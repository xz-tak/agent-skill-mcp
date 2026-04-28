#!/usr/bin/env python3
"""
Placeholder Population Helper

Extracts placeholders from DOCX template and assists with population.
Can be used standalone or as a module by target_evaluation_gen.py.

Usage:
    python populate_info.py --template path/to/template.docx --extract
    python populate_info.py --template path/to/template.docx --populate --target TYK2 --indication "Crohn's"
"""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple
from zipfile import ZipFile

# Skill paths
SKILL_ROOT = Path(__file__).parent.parent
GPT_PROMPTS_PATH = SKILL_ROOT / "templates" / "gpt_prompts.json"


def extract_placeholders_from_docx(docx_path: Path) -> Set[str]:
    """Extract all ${placeholder} patterns from DOCX file."""
    placeholders = set()
    pattern = re.compile(r'\$\{([^}]+)\}')

    # DOCX is a ZIP file containing XML
    with ZipFile(docx_path, 'r') as zf:
        # Check document.xml (main content)
        for xml_file in ['word/document.xml', 'word/header1.xml', 'word/footer1.xml']:
            try:
                content = zf.read(xml_file).decode('utf-8')
                matches = pattern.findall(content)
                placeholders.update(matches)
            except KeyError:
                continue

    return placeholders


def categorize_placeholders(placeholders: Set[str]) -> Dict[str, List[str]]:
    """Categorize placeholders by type."""
    categories = {
        "core": [],
        "pathway_interaction": [],
        "decision_framework": [],
        "risk_assessment": [],
        "scoring": [],
        "analysis": [],
        "metadata": [],
        "other": []
    }

    # Core placeholders (info.txt standard)
    core_keys = {
        "target", "target_full_name", "indication", "modality_primary",
        "related_diseases", "function_summary", "canonical_pathway",
        "expression_cells", "effector_cells", "moa_rationale",
        "combo_opportunities", "provided_sources"
    }

    # Pathway and interaction placeholders (new in v2.0)
    pathway_interaction_keys = {
        "shared_pathways", "unique_pathways_T1", "unique_pathways_T2",
        "interaction_path", "bridge_proteins", "combo_effect_type",
        "pathway_count_shared"
    }

    # Decision framework
    decision_keys = {
        "claim_1", "claim_2", "evidence_1a", "evidence_1b",
        "evidence_2a", "evidence_2b", "counter_1", "counter_2",
        "how_we_reconcile", "confidence_1", "confidence_2", "assump_1"
    }

    # Risk assessment
    risk_keys = {
        "failure_mode_1", "failure_mode_2", "mitigation_1", "mitigation_2",
        "impact_1", "impact_2"
    }

    # Scoring
    scoring_pattern = re.compile(r'^[ws]\d+$')

    # Analysis
    analysis_keys = {
        "option_A", "option_B", "option_C", "logic_1", "logic_2",
        "test_1", "flip_1", "one_paragraph_why", "decision_statement",
        "top3_risks", "next_actions"
    }

    # Metadata
    metadata_keys = {"dr_id", "when", "who", "base_1"}

    for p in placeholders:
        if p in core_keys:
            categories["core"].append(p)
        elif p in pathway_interaction_keys:
            categories["pathway_interaction"].append(p)
        elif p in decision_keys:
            categories["decision_framework"].append(p)
        elif p in risk_keys:
            categories["risk_assessment"].append(p)
        elif scoring_pattern.match(p) or p in {"total", "criteria_list_with_weights"}:
            categories["scoring"].append(p)
        elif p in analysis_keys:
            categories["analysis"].append(p)
        elif p in metadata_keys:
            categories["metadata"].append(p)
        else:
            categories["other"].append(p)

    # Sort each category
    for cat in categories:
        categories[cat] = sorted(categories[cat])

    return categories


def load_gpt_prompts() -> Dict[str, str]:
    """Load GPT prompt templates."""
    if GPT_PROMPTS_PATH.exists():
        with open(GPT_PROMPTS_PATH) as f:
            return json.load(f)
    return {}


def run_codex_gpt(prompt: str, timeout: int = 60) -> str:
    """Run codex CLI to generate content."""
    cmd = ["codex", "exec", "--skip-git-repo-check", "--sandbox", "read-only"]

    try:
        result = subprocess.run(
            cmd,
            input=prompt,
            capture_output=True,
            text=True,
            timeout=timeout
        )
        return result.stdout.strip() if result.returncode == 0 else "[TO BE PROVIDED]"
    except Exception:
        return "[TO BE PROVIDED]"


def populate_placeholder(
    placeholder: str,
    targets: List[str],
    indication: str,
    context: str = "",
    gpt_prompts: Dict[str, str] = None
) -> str:
    """Populate a single placeholder using GPT."""
    if gpt_prompts is None:
        gpt_prompts = load_gpt_prompts()

    # Get template or use default
    template = gpt_prompts.get(
        placeholder,
        f"For target(s) {{targets}} in {{indication}}, provide a value for {placeholder}."
    )

    # Format prompt
    prompt = template.format(
        targets=", ".join(targets),
        indication=indication,
        context=context[:2000] if context else ""
    )

    return run_codex_gpt(prompt)


def generate_info_txt(
    placeholders: Set[str],
    targets: List[str],
    indication: str,
    context: str = "",
    output_path: Path = None,
    interactive: bool = False
) -> Dict[str, str]:
    """Generate complete info.txt content."""

    replacements = {}
    gpt_prompts = load_gpt_prompts()

    # Set core values
    replacements["target"] = ", ".join(targets) if len(targets) > 1 else targets[0]
    replacements["indication"] = indication
    replacements["provided_sources"] = "provided files in local directory: *.json,xlsx,csv"

    # Categorize remaining placeholders
    remaining = placeholders - set(replacements.keys())
    categories = categorize_placeholders(remaining)

    # Process by category
    for category, keys in categories.items():
        if not keys:
            continue

        print(f"\n[{category.upper()}] Processing {len(keys)} placeholders...")

        for key in keys:
            if interactive:
                # Ask user first
                user_value = input(f"  {key} (press Enter to auto-generate): ").strip()
                if user_value:
                    replacements[key] = user_value
                    continue

            # Auto-generate via GPT
            print(f"  Generating {key}...")
            replacements[key] = populate_placeholder(
                key, targets, indication, context, gpt_prompts
            )

    # Save if output path provided
    if output_path:
        with open(output_path, "w") as f:
            for key, value in sorted(replacements.items()):
                escaped = value.replace("\n", " ").replace("\t", " ")
                f.write(f"${{{key}}}\t{escaped}\n")
        print(f"\nSaved to {output_path}")

    return replacements


def main():
    parser = argparse.ArgumentParser(description="DOCX Placeholder Population Helper")

    parser.add_argument(
        "--template", "-t",
        required=True,
        help="Path to DOCX template"
    )
    parser.add_argument(
        "--extract", "-e",
        action="store_true",
        help="Extract and list all placeholders"
    )
    parser.add_argument(
        "--populate", "-p",
        action="store_true",
        help="Populate placeholders"
    )
    parser.add_argument(
        "--target",
        help="Target gene(s), comma-separated"
    )
    parser.add_argument(
        "--indication",
        help="Primary indication"
    )
    parser.add_argument(
        "--output", "-o",
        help="Output path for info.txt"
    )
    parser.add_argument(
        "--interactive", "-i",
        action="store_true",
        help="Interactive mode - prompt for user values"
    )
    parser.add_argument(
        "--context-file",
        help="Path to context file (e.g., Cortellis JSON)"
    )

    args = parser.parse_args()

    template_path = Path(args.template)
    if not template_path.exists():
        print(f"Error: Template not found: {template_path}")
        sys.exit(1)

    # Extract placeholders
    placeholders = extract_placeholders_from_docx(template_path)

    if args.extract:
        categories = categorize_placeholders(placeholders)
        print(f"\n{'='*60}")
        print(f"Placeholders in {template_path.name}")
        print(f"{'='*60}")
        print(f"Total: {len(placeholders)}\n")

        for category, keys in categories.items():
            if keys:
                print(f"[{category.upper()}] ({len(keys)} items)")
                for key in keys:
                    print(f"  ${{{key}}}")
                print()

    if args.populate:
        if not args.target or not args.indication:
            print("Error: --target and --indication required for --populate")
            sys.exit(1)

        targets = [t.strip().upper() for t in args.target.split(",")]

        # Load context if provided
        context = ""
        if args.context_file:
            context_path = Path(args.context_file)
            if context_path.exists():
                context = context_path.read_text()[:5000]

        output_path = Path(args.output) if args.output else None

        generate_info_txt(
            placeholders,
            targets,
            args.indication,
            context,
            output_path,
            args.interactive
        )


if __name__ == "__main__":
    main()
