#!/usr/bin/env python3
"""
NCBI Paper Query - CLI Entry Point

Run with: python -m ncbi_paper_query [args]

Example:
    python -m ncbi_paper_query --disease "Crohn's disease" --tissue intestine \
        --organism human --if-cutoff 7.0 --year-cutoff 2020 --max-results 50 \
        --abstract-only --output my_study
"""

import argparse
import sys

from . import retrieve_publications, EntityMatcher


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="NCBI Publication Retrieval System",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples (run from your project directory with PYTHONPATH set):
  export PYTHONPATH=~/.claude/skills/ncbi-paper-query/scripts

  # Basic search (abstract only, fast):
  python -m ncbi_paper_query --disease "Crohn's disease" --tissue intestine \\
      --organism human --if-cutoff 7.0 --year-cutoff 2020 --max-results 50 \\
      --abstract-only --output crohn_intestine

  # Web scraping mode (extracts accessions from article HTML):
  python -m ncbi_paper_query --disease "IBD" --tissue intestine \\
      --organism human --if-cutoff 5.0 --max-results 20 \\
      --web-scrape --output ibd_webscrape

  # Full download mode (downloads PDFs):
  python -m ncbi_paper_query --disease "ulcerative colitis" --tissue colon \\
      --organism human --max-results 10 --output uc_full

Output Structure (created in your current working directory):
  output/{study_name}/
  ├── results.csv          # Main results with all metadata
  ├── accessions.csv       # Omics accessions only
  ├── failed_studies.csv   # Papers that couldn't be accessed
  └── downloads/           # Downloaded PDFs (if applicable)
        """
    )

    parser.add_argument("--disease", "-d", nargs="+", required=True,
                        help="Target disease(s) - multiple values use OR logic. Supports: space-separated (--disease Crohn UC), comma-separated (--disease 'Crohn,UC'), or both")
    parser.add_argument("--interpret-diseases", action="store_true",
                        help="Use Claude to interpret and expand disease terms (e.g., 'IBD' → ['Crohn disease', 'ulcerative colitis'])")
    parser.add_argument("--tissue", "-t", nargs="+", required=True,
                        help="Target tissue(s)/organ(s) - multiple values use OR logic")
    parser.add_argument("--exact", "-e", action="store_true",
                        help="Use exact matching for disease and tissue terms (no Claude or fuzzy expansion). "
                             "OR logic still applies for comma-separated or multiple terms.")
    parser.add_argument("--organism", "-o", nargs="+", default=["human", "mouse", "rat"],
                        help="Organisms to include")
    parser.add_argument("--cell-type", "-c", help="Specific cell type")
    parser.add_argument("--keywords", "-k", nargs="+",
                        default=["single-cell", "single cell", "scRNA-seq", "CITE-seq", "spatial",
                                "omics", "transcriptomics", "proteomics", "CyTOF", "RNA-seq",
                                "gene expression", "sequencing"],
                        help="Technology/method keywords to filter results (OR logic). "
                             "Default: single-cell, CITE-seq, spatial, omics, CyTOF, etc. "
                             "Use --no-keywords to disable keyword filtering.")
    parser.add_argument("--no-keywords", action="store_true",
                        help="Disable default keyword filtering (broader search, more results)")
    parser.add_argument("--if-cutoff", type=float, default=7.0,
                        help="Impact factor cutoff")
    parser.add_argument("--year-cutoff", type=int, default=2015,
                        help="Publication year cutoff")
    parser.add_argument("--no-preprints", action="store_true",
                        help="Exclude preprints")
    parser.add_argument("--include-unknown-if", action="store_true",
                        help="Include papers with unknown impact factors. By default, papers without a known IF are excluded.")
    parser.add_argument("--abstract-only", action="store_true",
                        help="Extract omics from abstract/title only (no full-text access)")
    parser.add_argument("--web-scrape", action="store_true",
                        help="Extract omics from web HTML without downloading PDFs. Uses Playwright to fetch article pages (PMC first, then publisher). Auto-fallback to PDF download if no accessions found.")
    parser.add_argument("--subscription-download", action="store_true",
                        help="Hybrid mode: web-scrape free/open access papers (no download), but download PDFs for subscription papers. Saves bandwidth while still getting paywalled content.")
    parser.add_argument("--subscription-only", action="store_true",
                        help="Only download subscription papers (skip free PMC papers). Requires VPN for institutional access.")
    parser.add_argument("--no-institutional-auth", action="store_true",
                        help="Disable institutional authentication. By default, Takeda credentials from .env are used for subscription papers.")
    parser.add_argument("--max-results", type=int, default=1000,
                        help="Maximum results")
    parser.add_argument("--rounds", "-r", type=int, default=3,
                        help="Number of validation rounds (default 3). Round 1: initial retrieval, Round 2+: validation and proofreading")
    parser.add_argument("--output", "-O",
                        help="Study name for output directory. Creates output/{name}/ with: "
                             "results.csv, accessions.csv, downloads/. "
                             "If not provided, uses timestamped filenames in output/ directory.")
    parser.add_argument("--format", choices=["csv", "excel"], default="csv",
                        help="Output format")

    args = parser.parse_args()

    # Process disease arguments: handle comma-separated values
    diseases = []
    for d in args.disease:
        # Split by comma if present
        if ',' in d:
            diseases.extend([term.strip() for term in d.split(',') if term.strip()])
        else:
            diseases.append(d.strip())

    # Optionally use Claude to interpret/expand disease terms
    if args.interpret_diseases:
        matcher = EntityMatcher()
        # Join all diseases and interpret as a single input
        disease_input = ", ".join(diseases)
        diseases = matcher.interpret_diseases(disease_input)
        print(f"Interpreted diseases: {diseases}")

    # Handle keywords - use None if --no-keywords is set
    keywords = None if args.no_keywords else args.keywords

    df = retrieve_publications(
        disease=diseases,
        tissues=args.tissue,
        organism=args.organism,
        cell_type=args.cell_type,
        keywords=keywords,
        impact_factor_cutoff=args.if_cutoff,
        publication_date_cutoff=args.year_cutoff,
        max_results=args.max_results,
        download_papers=not (args.abstract_only or args.web_scrape or args.subscription_download),
        scrape_web_content=args.web_scrape,
        subscription_download_only=args.subscription_download,
        subscription_only=args.subscription_only,
        use_institutional_auth=not args.no_institutional_auth,
        include_preprints=not args.no_preprints,
        include_unknown_if=args.include_unknown_if,
        output_format=args.format,
        study_name=args.output,
        validation_rounds=args.rounds,
        exact_match=args.exact,
    )

    print(f"\n{'='*60}")
    print(f"Retrieved {len(df)} publications")

    if args.output:
        print(f"Results saved to: output/{args.output}/")

    return 0


if __name__ == "__main__":
    sys.exit(main())
