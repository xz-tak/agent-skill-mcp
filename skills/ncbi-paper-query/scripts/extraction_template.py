#!/usr/bin/env python3
"""
Paper Extraction Template for Codex

This template is used by Codex to iteratively extract structured information
from downloaded papers. Codex will modify this template based on user requirements.

Usage:
    # With study name (recommended) - uses ./output/{study_name}/ structure
    python extraction_template.py --study-name ibd_crohn \\
        --fields "study_design,key_findings,biomarkers,genes_mentioned"

    # With explicit paths (advanced)
    python extraction_template.py \\
        --input-dir ./output/my_study/downloads \\
        --output ./output/my_study/extraction_results.csv \\
        --fields "study_design,key_findings"
"""

import os
import re
import csv
import json
import argparse
from pathlib import Path
from typing import List, Dict, Optional, Any
from dataclasses import dataclass, field, asdict

# Skill directory - defaults to ~/.claude/skills/ncbi-paper-query
SKILL_DIR = Path(os.environ.get("SKILL_DIR",
    Path.home() / ".claude" / "skills" / "ncbi-paper-query"))

# PDF parsing
try:
    from PyPDF2 import PdfReader
    PDF_AVAILABLE = True
except ImportError:
    PDF_AVAILABLE = False
    print("Warning: PyPDF2 not installed. PDF extraction disabled.")

# HTML parsing
try:
    from bs4 import BeautifulSoup
    HTML_AVAILABLE = True
except ImportError:
    HTML_AVAILABLE = False
    print("Warning: BeautifulSoup not installed. HTML extraction disabled.")


@dataclass
class ExtractionResult:
    """Container for extracted information from a single paper."""
    pmid: str
    title: str
    extraction_status: str = "pending"
    error_message: Optional[str] = None

    # Extraction fields (populated dynamically based on user requirements)
    study_design: Optional[str] = None
    sample_size: Optional[str] = None
    methodology: Optional[str] = None
    key_findings: Optional[str] = None
    statistical_results: Optional[str] = None
    biomarkers: Optional[str] = None
    genes_mentioned: Optional[str] = None
    pathways: Optional[str] = None
    dataset_types: Optional[str] = None
    omics_platforms: Optional[str] = None
    sample_characteristics: Optional[str] = None
    therapeutic_targets: Optional[str] = None
    drug_candidates: Optional[str] = None
    cell_types_analyzed: Optional[str] = None
    conclusions: Optional[str] = None


class PaperExtractor:
    """Extract structured information from papers."""

    # Field extraction patterns (rule-based fallback)
    PATTERNS = {
        "sample_size": [
            r"n\s*=\s*(\d+)",
            r"(\d+)\s+patients",
            r"(\d+)\s+subjects",
            r"sample\s+size\s*[:\s]+(\d+)",
            r"(\d+)\s+participants",
            r"(\d+)\s+samples",
        ],
        "study_design": [
            r"(randomized controlled trial|RCT)",
            r"(cohort study|prospective study|retrospective study)",
            r"(case-control|case control)",
            r"(cross-sectional|meta-analysis|systematic review)",
            r"(observational study|clinical trial|pilot study)",
            r"(single-cell RNA-seq|scRNA-seq|bulk RNA-seq)",
        ],
        "genes_mentioned": [
            # Common gene name patterns
            r"\b([A-Z][A-Z0-9]{2,})\b",  # Uppercase gene symbols
        ],
        "biomarkers": [
            r"(biomarker|marker|indicator|signature)",
            r"(CRP|IL-\d+|TNF|IFN|CD\d+)",
        ],
    }

    def __init__(self, fields: List[str]):
        """
        Initialize extractor with list of fields to extract.

        Args:
            fields: List of field names to extract (e.g., ["study_design", "key_findings"])
        """
        self.fields = fields

    def extract_from_pdf(self, pdf_path: str) -> str:
        """Extract text content from PDF file."""
        if not PDF_AVAILABLE:
            return ""

        try:
            reader = PdfReader(pdf_path)
            text = ""
            for page in reader.pages:
                text += page.extract_text() + "\n"
            return text
        except Exception as e:
            print(f"Error extracting PDF {pdf_path}: {e}")
            return ""

    def extract_from_html(self, html_path: str) -> str:
        """Extract text content from HTML file."""
        if not HTML_AVAILABLE:
            return ""

        try:
            with open(html_path, 'r', encoding='utf-8') as f:
                soup = BeautifulSoup(f.read(), 'html.parser')

            # Remove script and style elements
            for element in soup(['script', 'style', 'nav', 'footer', 'header']):
                element.decompose()

            return soup.get_text(separator='\n', strip=True)
        except Exception as e:
            print(f"Error extracting HTML {html_path}: {e}")
            return ""

    def extract_field(self, text: str, field_name: str) -> Optional[str]:
        """
        Extract a specific field from text using patterns.

        This is a rule-based fallback. Codex will typically use LLM-based
        extraction for more accurate results.
        """
        patterns = self.PATTERNS.get(field_name, [])

        for pattern in patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                return match.group(0)

        return None

    def extract_all(self, text: str, pmid: str, title: str) -> ExtractionResult:
        """
        Extract all requested fields from paper text.

        Args:
            text: Full text content of the paper
            pmid: PubMed ID
            title: Paper title

        Returns:
            ExtractionResult with populated fields
        """
        result = ExtractionResult(pmid=pmid, title=title)

        try:
            for field in self.fields:
                value = self.extract_field(text, field)
                if hasattr(result, field):
                    setattr(result, field, value)

            result.extraction_status = "success"
        except Exception as e:
            result.extraction_status = "error"
            result.error_message = str(e)

        return result


def find_paper_files(input_dir: str) -> List[Dict[str, str]]:
    """
    Find all paper files in the input directory.

    Looks for PDFs in both flat and nested (PMID_*/) directory structures.

    Returns list of dicts with 'path', 'pmid', 'type' keys.
    """
    papers = []
    input_path = Path(input_dir)

    # Look for PDFs in flat structure
    for pdf_file in input_path.glob("*.pdf"):
        # Extract PMID from filename (assumed format: PMID_title.pdf or PMID.pdf)
        pmid = pdf_file.stem.split("_")[0]
        papers.append({
            "path": str(pdf_file),
            "pmid": pmid,
            "type": "pdf"
        })

    # Look for PDFs in nested PMID_*/ structure (from ncbi_paper_query downloads)
    for pmid_dir in input_path.glob("PMID_*"):
        if pmid_dir.is_dir():
            pmid = pmid_dir.name.replace("PMID_", "")
            for pdf_file in pmid_dir.glob("*.pdf"):
                papers.append({
                    "path": str(pdf_file),
                    "pmid": pmid,
                    "type": "pdf"
                })

    # Look for HTML files
    for html_file in input_path.glob("*.html"):
        pmid = html_file.stem.split("_")[0]
        papers.append({
            "path": str(html_file),
            "pmid": pmid,
            "type": "html"
        })

    return papers


def main():
    parser = argparse.ArgumentParser(
        description="Extract structured information from downloaded papers"
    )
    parser.add_argument(
        "--study-name", "-s",
        help="Study name (creates ./output/{study_name}/ structure). "
             "Input defaults to ./output/{study_name}/downloads, "
             "output defaults to ./output/{study_name}/extraction_results.csv"
    )
    parser.add_argument(
        "--input-dir", "-i",
        help="Directory containing downloaded papers (PDF/HTML). "
             "Defaults to ./output/{study_name}/downloads if --study-name provided"
    )
    parser.add_argument(
        "--output", "-o",
        help="Output CSV file path. "
             "Defaults to ./output/{study_name}/extraction_results.csv if --study-name provided"
    )
    parser.add_argument(
        "--fields", "-f",
        default="study_design,key_findings,biomarkers,genes_mentioned,sample_size,methodology",
        help="Comma-separated list of fields to extract. Available fields: "
             "study_design, sample_size, methodology, key_findings, statistical_results, "
             "biomarkers, genes_mentioned, pathways, dataset_types, omics_platforms, "
             "sample_characteristics, therapeutic_targets, drug_candidates, "
             "cell_types_analyzed, conclusions"
    )
    parser.add_argument(
        "--metadata", "-m",
        help="Path to results CSV with PMID/Title metadata (e.g., results.csv from ncbi_paper_query)"
    )

    args = parser.parse_args()

    # Resolve paths based on study_name or explicit arguments
    if args.study_name:
        base_dir = Path("./output") / args.study_name
        input_dir = args.input_dir or str(base_dir / "downloads")
        output_path = args.output or str(base_dir / "extraction_results.csv")
        metadata_path = args.metadata or str(base_dir / "results.csv")
    else:
        input_dir = args.input_dir or "./output/downloads"
        output_path = args.output or "./output/extraction_results.csv"
        metadata_path = args.metadata

    # Parse fields
    fields = [f.strip() for f in args.fields.split(",")]
    print(f"Extracting fields: {fields}")
    print(f"Input directory: {input_dir}")
    print(f"Output file: {output_path}")

    # Find paper files
    papers = find_paper_files(input_dir)
    print(f"Found {len(papers)} papers to process")

    if not papers:
        print("No papers found in input directory")
        return

    # Load metadata if provided
    metadata = {}
    if metadata_path and os.path.exists(metadata_path):
        import pandas as pd
        df = pd.read_csv(metadata_path)
        for _, row in df.iterrows():
            pmid = str(row.get("PMID", ""))
            metadata[pmid] = {
                "title": row.get("Title", ""),
                "journal": row.get("Journal", ""),
                "year": row.get("Publication_Year", "")
            }
        print(f"Loaded metadata for {len(metadata)} papers from {metadata_path}")

    # Initialize extractor
    extractor = PaperExtractor(fields)

    # Process each paper
    results = []
    for i, paper in enumerate(papers, 1):
        print(f"Processing {i}/{len(papers)}: PMID {paper['pmid']}")

        # Extract text based on file type
        if paper["type"] == "pdf":
            text = extractor.extract_from_pdf(paper["path"])
        else:
            text = extractor.extract_from_html(paper["path"])

        # Get title from metadata or use placeholder
        title = metadata.get(paper["pmid"], {}).get("title", "Unknown")

        # Extract fields
        result = extractor.extract_all(text, paper["pmid"], title)
        results.append(result)

    # Write results to CSV
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        # Determine columns based on fields
        columns = ["pmid", "title", "extraction_status"] + fields + ["error_message"]
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()

        for result in results:
            row = asdict(result)
            writer.writerow({k: row.get(k, "") for k in columns})

    print(f"\nResults saved to: {output_file}")
    print(f"Processed {len(results)} papers")

    # Summary
    success = sum(1 for r in results if r.extraction_status == "success")
    errors = sum(1 for r in results if r.extraction_status == "error")
    print(f"Success: {success}, Errors: {errors}")


if __name__ == "__main__":
    main()
