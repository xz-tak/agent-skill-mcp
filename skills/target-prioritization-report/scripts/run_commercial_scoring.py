#!/usr/bin/env python3
"""
Commercial Scoring Preparation Script for Target Prioritization Report.

Extracts CI context from HTML dashboard, builds research prompts per combo,
and outputs a structured JSON ready for deep-research scoring.

Usage:
    python scripts/run_commercial_scoring.py \
        --workdir /path/to/project \
        --indication ibd \
        --config pipeline_config.json

Outputs:
    results/commercial_<indication>/commercial_prompts.json  (research prompts)
    results/commercial_<indication>/commercial_scores.json   (scores, after research)

The script has two modes:
    --prepare : Generate prompts + CI context (run BEFORE deep research)
    --finalize : Validate and format scores (run AFTER deep research fills scores)
"""

import argparse
import json
import re
import sys
from pathlib import Path


SKILL_DIR = Path(__file__).parent.parent
TEMPLATE_FILE = SKILL_DIR / "references" / "commercial-research-prompt-template.md"

COMMERCIAL_WEIGHTS = {
    "market_opportunity": 0.40,
    "competitive_profile": 0.40,
    "strategic_fit": 0.20,
}


def load_research_template():
    """Load the standard prompt template from references/commercial-research-prompt-template.md.

    The template is indication/combo-agnostic. Variables like {combo_name},
    {indication_full}, {ci_context} are filled per-combo at runtime.
    """
    if TEMPLATE_FILE.exists():
        content = TEMPLATE_FILE.read_text()
        m = re.search(r"```\n(.*?)```", content, re.DOTALL)
        if m:
            return m.group(1)
    return _FALLBACK_TEMPLATE


_FALLBACK_TEMPLATE = """Research the commercial viability of a {combo_name} combination therapy for {indication_full}.

## Assessment Criteria

Score each dimension on a continuous 0-100 scale (higher = better opportunity / lower risk).

### 1. MARKET OPPORTUNITY (weight: 40%)
- Total addressable market size for {indication_full}
- Patient population (prevalent + incident), growth trajectory
- Unmet need: what current SOC fails to address that this combo could
- Reimbursement/payer landscape considerations

### 2. COMPETITIVE PROFILE (weight: 40%)
- Existing marketed therapies targeting {gene_list} (monotherapy or combo)
- Pipeline programs with similar mechanisms (phase, timeline)
- Patent/exclusivity landscape
- Differentiation: what makes this combo unique vs existing options
- Biosimilar/generic risk for component monotherapies

### 3. STRATEGIC FIT (weight: 20%)
- Alignment with GI/immunology therapeutic area focus
- Leverageability of existing programs and capabilities
- (Score modality-agnostic, only on TA relevance)

## Local Pipeline Context (from CI Dashboard)
{ci_context}

## Cortellis Summary
{cortellis_context}

## Required Output Format
Return a JSON object with this exact schema:
{{
  "combo": "{combo_name}",
  "indication": "{indication}",
  "commercial_score": <float 0-100>,
  "confidence": "<High|Medium|Low>",
  "components": {{
    "market_opportunity": {{
      "score": <float 0-100>,
      "rationale": "<1-2 sentence justification>",
      "key_evidence": ["<evidence point 1>", "<evidence point 2>"]
    }},
    "competitive_profile": {{
      "score": <float 0-100>,
      "rationale": "<1-2 sentence justification>",
      "key_evidence": ["<evidence point 1>", "<evidence point 2>"]
    }},
    "strategic_fit": {{
      "score": <float 0-100>,
      "rationale": "<1-2 sentence justification>",
      "key_evidence": ["<evidence point 1>", "<evidence point 2>"]
    }}
  }},
  "citations": [
    {{"source": "<url or reference>", "relevance": "<what it supports>"}}
  ],
  "narrative_summary": "<one paragraph executive summary>"
}}

Scoring guidance:
- Market: {indication_full} is a {market_hint} market. Calibrate thresholds accordingly.
- Competition: score LOW (0-30) if highly crowded with little differentiation, HIGH (70-100) if first-in-class with clear advantages.
- Strategic: GI/immunology core = 80-100, adjacent autoimmune = 60-79, distant TA = 20-59.
- commercial_score = 0.40*market + 0.40*competitive + 0.20*strategic

CRITICAL: Scores must DIFFERENTIATE between combos within the same indication.
- All combos share the same market/indication context, so Market scores may be similar.
- Competitive Profile should be the PRIMARY differentiator between combos:
  * A combo targeting a novel mechanism with no direct competitors scores HIGH (75-95).
  * A combo where BOTH targets have marketed drugs (e.g., JAK1+TYK2) scores LOWER (30-55)
    because the competitive landscape is already crowded.
  * A combo with ONE validated target + ONE novel target scores MODERATE (50-75).
- Use the CI pipeline data above to ground your competitive assessment in REAL pipeline counts.
- Combos with targets that have many Phase III+ programs are MORE competitive = LOWER score.
- Combos with targets that have NO pipeline programs offer HIGHER novelty = HIGHER score.
"""

RESEARCH_TEMPLATE = None  # Loaded lazily from template file


def extract_ci_data(ci_html_path):
    """Extract entries from CI HTML dashboard."""
    content = Path(ci_html_path).read_text(errors="ignore")
    m = re.search(r'<script[^>]*id="data"[^>]*>(.*?)</script>', content, re.DOTALL)
    if not m:
        return []
    return json.loads(m.group(1)).get("entries", [])


def build_gene_ci_map(entries, genes):
    """Map each gene to its CI pipeline entries."""
    target_map = {}
    for entry in entries:
        for t in entry.get("targets", []):
            tl = t.lower()
            for gene in genes:
                if gene.lower() in tl:
                    if gene not in target_map:
                        target_map[gene] = []
                    target_map[gene].append({
                        "drug": entry.get("displayName", ""),
                        "phase": entry.get("ibdPhase", entry.get("overallPhase", "")),
                        "org": entry.get("organization", "")[:60],
                        "active": entry.get("underActiveDevelopment", ""),
                        "modality": entry.get("modality", ""),
                    })
                    break
    return target_map


def build_combo_ci_context(combo, target_map):
    """Build CI summary text for a combo."""
    lines = []
    for g in combo:
        drugs = target_map.get(g, [])
        active = [d for d in drugs if d["active"] == "Yes"]
        if active:
            drug_lines = [f"  - {d['drug']} ({d['phase']}, {d['org']})" for d in active[:5]]
            lines.append(f"{g}: {len(active)} active drugs\n" + "\n".join(drug_lines))
        else:
            lines.append(f"{g}: No active drugs in pipeline")
    return "\n".join(lines)


def build_cortellis_context(combo, cortellis_dir):
    """Build Cortellis summary text for a combo."""
    lines = []
    for g in combo:
        summary_file = cortellis_dir / f"{g}_summary.md"
        if summary_file.exists():
            content = summary_file.read_text(errors="ignore")
            first_lines = "\n".join(content.split("\n")[:10])
            lines.append(f"{g}: {first_lines}")
        else:
            lines.append(f"{g}: No Cortellis summary available")
    return "\n".join(lines)


def estimate_market_hint(indication):
    """Provide market-size hint for dynamic calibration."""
    hints = {
        "ibd": "large ($25B+ global)",
        "uc": "large ($15B+ global, IBD subset)",
        "cd": "large ($15B+ global, IBD subset)",
        "ssc": "moderate-to-small ($2-5B global)",
        "ipf": "moderate ($5-8B global)",
        "nash": "large (projected $20B+)",
        "ra": "large ($20B+ global)",
    }
    return hints.get(indication.lower(), "size to be determined by research")


def prepare_prompts(workdir, indication, config):
    """Generate research prompts for all combos."""
    combos = config["combos"]
    genes = config["genes_unique"]
    diseases = config.get("diseases", {})
    indication_full = diseases.get("primary", indication.upper())
    subtypes = diseases.get("subtypes", [])
    if subtypes:
        indication_full = f"{indication_full} ({', '.join(subtypes[:2])})"

    ci_dir = workdir / "results" / f"ci_{indication}"
    ci_entries = []
    for html_file in ci_dir.glob("*.html"):
        ci_entries = extract_ci_data(html_file)
        if ci_entries:
            break

    target_map = build_gene_ci_map(ci_entries, genes)

    cortellis_dir = workdir / "results" / f"cortellis_{indication}"
    market_hint = estimate_market_hint(indication)

    prompts = []
    for combo in combos:
        combo_name = " + ".join(combo)
        ci_context = build_combo_ci_context(combo, target_map)
        cortellis_context = build_cortellis_context(combo, cortellis_dir)

        template = load_research_template()
        prompt = (template
            .replace("{combo_name}", combo_name)
            .replace("{indication}", indication.upper())
            .replace("{indication_full}", indication_full)
            .replace("{gene_list}", ", ".join(combo))
            .replace("{ci_context}", ci_context)
            .replace("{cortellis_context}", cortellis_context)
            .replace("{market_hint}", market_hint)
        )

        prompts.append({
            "combo": combo,
            "combo_name": combo_name,
            "prompt": prompt,
            "ci_context": ci_context,
            "cortellis_context": cortellis_context,
        })

    output_dir = workdir / "results" / f"commercial_{indication}"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "commercial_prompts.json"
    output_file.write_text(json.dumps(prompts, indent=2))
    print(f"Prepared {len(prompts)} research prompts → {output_file}")
    return prompts


def validate_scores(scores_file):
    """Validate commercial_scores.json structure."""
    scores = json.loads(Path(scores_file).read_text())
    required_keys = {"combo", "commercial_score", "confidence", "components", "narrative_summary"}
    component_keys = {"market_opportunity", "competitive_profile", "strategic_fit"}

    errors = []
    for i, entry in enumerate(scores):
        missing = required_keys - set(entry.keys())
        if missing:
            errors.append(f"Entry {i} ({entry.get('combo', '?')}): missing keys {missing}")
        if "components" in entry:
            comp_missing = component_keys - set(entry["components"].keys())
            if comp_missing:
                errors.append(f"Entry {i}: missing components {comp_missing}")
            for key in component_keys & set(entry["components"].keys()):
                score = entry["components"][key].get("score")
                if score is None or not (0 <= score <= 100):
                    errors.append(f"Entry {i}: {key}.score invalid ({score})")

        computed = (
            COMMERCIAL_WEIGHTS["market_opportunity"] * entry.get("components", {}).get("market_opportunity", {}).get("score", 0)
            + COMMERCIAL_WEIGHTS["competitive_profile"] * entry.get("components", {}).get("competitive_profile", {}).get("score", 0)
            + COMMERCIAL_WEIGHTS["strategic_fit"] * entry.get("components", {}).get("strategic_fit", {}).get("score", 0)
        )
        declared = entry.get("commercial_score", 0)
        if abs(computed - declared) > 1.0:
            errors.append(f"Entry {i} ({entry.get('combo')}): declared score {declared:.1f} != computed {computed:.1f}")

    if errors:
        print(f"Validation FAILED ({len(errors)} errors):")
        for e in errors:
            print(f"  - {e}")
        return False
    print(f"Validation PASSED: {len(scores)} entries, all valid.")
    return True


def main():
    parser = argparse.ArgumentParser(description="Commercial Scoring for Target Prioritization")
    parser.add_argument("--workdir", type=str, required=True, help="Project working directory")
    parser.add_argument("--indication", type=str, required=True, help="Indication label (e.g., ibd)")
    parser.add_argument("--config", type=str, default="pipeline_config.json", help="Pipeline config file")
    parser.add_argument("--mode", choices=["prepare", "validate"], default="prepare",
                       help="prepare=generate prompts, validate=check scores file")
    parser.add_argument("--scores-file", type=str, help="Path to scores JSON (for validate mode)")
    args = parser.parse_args()

    workdir = Path(args.workdir)
    config_path = workdir / args.config
    if not config_path.exists():
        print(f"Error: config not found at {config_path}")
        sys.exit(1)

    config = json.loads(config_path.read_text())

    if args.mode == "prepare":
        prepare_prompts(workdir, args.indication, config)
    elif args.mode == "validate":
        scores_file = args.scores_file or str(workdir / "results" / f"commercial_{args.indication}" / "commercial_scores.json")
        if not Path(scores_file).exists():
            print(f"Error: scores file not found at {scores_file}")
            sys.exit(1)
        valid = validate_scores(scores_file)
        sys.exit(0 if valid else 1)


if __name__ == "__main__":
    main()
