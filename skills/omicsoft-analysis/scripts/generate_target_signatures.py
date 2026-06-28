#!/usr/bin/env python3
"""Generate functional gene signatures for target genes via Codex GPT-5.5 xhigh.

For each target gene, queries Codex to produce a signature of 10-50 co-regulated
genes based on perturbation evidence (CRISPR screens, Perturb-seq, OE studies,
pathway databases). Results are cached to the output JSON file: genes already
present are skipped on re-runs.
"""

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


CODEX_MODEL = "gpt-5.5"
CODEX_REASONING = "xhigh"
CODEX_TIMEOUT_SEC = 180


def build_prompt(gene: str, disease: str, tissue: str) -> str:
    return f"""You are a functional genomics expert specializing in perturbation biology.

Generate a gene signature for the target gene: {gene}
Disease context: {disease} in {tissue}

Requirements:
1. List 10-50 genes that are CO-REGULATED with {gene} under perturbation
2. These genes should go DOWN when {gene} is knocked out/knocked down/inhibited, OR go UP when {gene} is overexpressed/activated
3. Prioritize genes with HIGH-CONFIDENCE functional genomics evidence:
   - CRISPR knockout/knockdown screens (DepMap, Perturb-seq)
   - Overexpression studies
   - Pathway databases (KEGG, Reactome) showing direct functional interaction
4. Include genes from key pathways that {gene} participates in
5. ALL gene symbols must be HUMAN (ALL UPPERCASE)
6. Prioritize {tissue}-relevant genes when possible

Output ONLY valid JSON (no markdown, no code blocks, no explanation outside JSON):
{{
  "gene": "{gene}",
  "signature_genes": ["GENE1", "GENE2", ...],
  "perturbation_direction": "KO/KD downregulated" or "OE/activation upregulated",
  "key_pathways": ["pathway1", "pathway2", ...],
  "rationale": [
    {{"gene": "GENE1", "confidence": "high|medium", "evidence": "specific study or database", "pathway": "relevant pathway"}},
    ...
  ]
}}"""


def build_simpler_prompt(gene: str, disease: str, tissue: str) -> str:
    return f"""Generate a functional gene signature for {gene} in the context of {disease} ({tissue}).

Return ONLY a JSON object with this exact schema (no markdown, no prose):
{{
  "gene": "{gene}",
  "signature_genes": ["GENE1", "GENE2", ...],
  "perturbation_direction": "KO/KD downregulated",
  "key_pathways": ["pathway1"],
  "rationale": [{{"gene": "GENE1", "confidence": "high", "evidence": "source", "pathway": "name"}}]
}}

Include 10-50 human gene symbols (UPPERCASE) that are co-regulated with {gene} based on CRISPR screens, Perturb-seq, or pathway databases."""


def extract_json(stdout: str) -> dict:
    """Pull the first balanced JSON object out of Codex stdout.

    Codex sometimes prepends/appends status lines or wraps output in fences,
    so we scan for the first '{' and bracket-match to its closing '}'.
    """
    text = stdout.strip()
    if text.startswith("```"):
        first_nl = text.find("\n")
        if first_nl != -1:
            text = text[first_nl + 1 :]
        if text.endswith("```"):
            text = text[:-3]
        text = text.strip()

    start = text.find("{")
    if start == -1:
        raise ValueError("No '{' found in Codex output")

    depth = 0
    in_str = False
    escape = False
    for i in range(start, len(text)):
        ch = text[i]
        if in_str:
            if escape:
                escape = False
            elif ch == "\\":
                escape = True
            elif ch == '"':
                in_str = False
            continue
        if ch == '"':
            in_str = True
        elif ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return json.loads(text[start : i + 1])

    raise ValueError("Unbalanced JSON in Codex output")


def call_codex(prompt: str) -> str:
    result = subprocess.run(
        ["codex", "exec", "--skip-git-repo-check", prompt],
        capture_output=True,
        text=True,
        timeout=CODEX_TIMEOUT_SEC,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"codex exec exited with code {result.returncode}: {result.stderr[:500]}"
        )
    return result.stdout


def generate_signature(gene: str, disease: str, tissue: str) -> dict:
    """Call Codex for one gene; retry once with a simpler prompt on parse failure."""
    last_error: Exception | None = None
    last_stdout: str = ""

    for attempt, prompt_builder in enumerate(
        (build_prompt, build_simpler_prompt), start=1
    ):
        prompt = prompt_builder(gene, disease, tissue)
        try:
            stdout = call_codex(prompt)
            last_stdout = stdout
            print(
                f"  [attempt {attempt}] codex returned {len(stdout)} chars",
                file=sys.stderr,
            )
            parsed = extract_json(stdout)
            parsed["_raw_response_len"] = len(stdout)
            return parsed
        except (json.JSONDecodeError, ValueError, RuntimeError, subprocess.TimeoutExpired) as exc:
            last_error = exc
            print(
                f"  [attempt {attempt}] failed: {type(exc).__name__}: {exc}",
                file=sys.stderr,
            )
            if last_stdout:
                preview = last_stdout[:500].replace("\n", " ")
                print(f"  [attempt {attempt}] stdout preview: {preview!r}", file=sys.stderr)

    raise RuntimeError(f"Failed to generate signature for {gene}: {last_error}")


def load_existing(output_path: Path) -> dict:
    if not output_path.exists():
        return {}
    try:
        with output_path.open("r") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError) as exc:
        print(
            f"warning: could not load existing output {output_path}: {exc}; starting fresh",
            file=sys.stderr,
        )
        return {}


def save_output(output_path: Path, payload: dict) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = output_path.with_suffix(output_path.suffix + ".tmp")
    with tmp.open("w") as f:
        json.dump(payload, f, indent=2)
    tmp.replace(output_path)


def signatures_to_addon_format(sig_json: dict, combo_genes: list[str]) -> str:
    """Convert the full signatures JSON to the --addon-signatures CLI format.

    Format: "GENE1_sig:G1,G2,...;GENE2_sig:GA,GB,..."
    Only genes in combo_genes that have a successful signature are included.
    """
    signatures = sig_json.get("signatures", {})
    parts: list[str] = []
    for gene in combo_genes:
        entry = signatures.get(gene)
        if not entry:
            continue
        genes = entry.get("signature_genes") or []
        if not genes:
            continue
        parts.append(f"{gene}_sig:{','.join(genes)}")
    return ";".join(parts)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate functional gene signatures via Codex GPT-5.5 xhigh."
    )
    parser.add_argument(
        "--genes",
        required=True,
        help="Comma-separated target genes (e.g., TYK2,JAK1,STAT3)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output JSON path (e.g., signatures_auto.json)",
    )
    parser.add_argument(
        "--disease",
        default="IBD",
        help="Disease context for signature generation (default: IBD)",
    )
    parser.add_argument(
        "--tissue",
        default="intestine",
        help="Tissue context (default: intestine)",
    )
    args = parser.parse_args()

    genes = [g.strip().upper() for g in args.genes.split(",") if g.strip()]
    if not genes:
        print("error: --genes produced an empty list", file=sys.stderr)
        return 2

    output_path = Path(args.output)
    existing = load_existing(output_path)

    payload = {
        "generated_at": existing.get("generated_at")
        or datetime.now(timezone.utc).isoformat(),
        "model": CODEX_MODEL,
        "reasoning": CODEX_REASONING,
        "disease_context": args.disease,
        "tissue_context": args.tissue,
        "genes_requested": genes,
        "signatures": dict(existing.get("signatures") or {}),
        "failed_genes": list(existing.get("failed_genes") or []),
    }
    # Reset failed_genes for genes we're about to retry this run.
    payload["failed_genes"] = [g for g in payload["failed_genes"] if g not in genes]

    total = len(genes)
    for idx, gene in enumerate(genes, start=1):
        if gene in payload["signatures"]:
            print(
                f"[{idx}/{total}] Skipping {gene} (already in cache)",
                file=sys.stderr,
            )
            continue

        print(f"[{idx}/{total}] Generating signature for {gene}...", file=sys.stderr)
        try:
            sig = generate_signature(gene, args.disease, args.tissue)
        except Exception as exc:
            print(f"  -> FAILED: {exc}", file=sys.stderr)
            if gene not in payload["failed_genes"]:
                payload["failed_genes"].append(gene)
            save_output(output_path, payload)
            continue

        entry = {
            "signature_genes": sig.get("signature_genes", []),
            "perturbation_direction": sig.get("perturbation_direction", ""),
            "key_pathways": sig.get("key_pathways", []),
            "rationale": sig.get("rationale", []),
        }
        payload["signatures"][gene] = entry
        if gene in payload["failed_genes"]:
            payload["failed_genes"].remove(gene)
        n = len(entry["signature_genes"])
        print(f"  -> OK: {n} signature genes", file=sys.stderr)
        save_output(output_path, payload)

    save_output(output_path, payload)
    print(
        f"\nDone. signatures: {len(payload['signatures'])}, failed: {len(payload['failed_genes'])}",
        file=sys.stderr,
    )
    print(f"Output: {output_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
