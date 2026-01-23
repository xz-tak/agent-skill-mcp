"""
NCBI Paper Query - Retrieval Module

Main entry point for publication retrieval with multi-round validation.
"""

import logging
from typing import List, Union

import pandas as pd

from ..config import logger
from ..models import Publication
from ..retrieval import PubMedSearcher
from ..matching import EntityMatcher
from ..metadata import ImpactFactorLookup
from ..download import PaperDownloader, OmicsExtractor
from ..export import ResultsExporter
from .validation import validate_results, compare_rounds


def retrieve_publications(
    disease: Union[str, List[str]],
    tissues: List[str],
    organism: List[str] = None,
    cell_type: str = None,
    keywords: List[str] = None,
    impact_factor_cutoff: float = 7.0,
    publication_date_cutoff: int = 2015,
    include_preprints: bool = True,
    include_unknown_if: bool = False,
    download_papers: bool = True,
    scrape_web_content: bool = False,
    subscription_only: bool = False,
    use_institutional_auth: bool = True,
    max_results: int = 1000,
    study_name: str = None,
    output_format: str = "csv",
    validation_rounds: int = 1,
    exact_match: bool = False
) -> pd.DataFrame:
    """
    Main entry point for publication retrieval with multi-round validation.

    Args:
        disease: Target disease (fuzzy matched)
        tissues: List of target tissues/organs (fuzzy matched, OR logic)
        organism: List of organisms to filter ["human", "mouse", "rat"]
        cell_type: Specific cell type to match (optional)
        keywords: Technology/method keywords to filter (OR logic, e.g., ["single-cell", "omics"])
        impact_factor_cutoff: Minimum journal IF (ignored for preprints)
        publication_date_cutoff: Papers from this year onward
        include_preprints: Include preprint servers
        include_unknown_if: If True, include papers with unknown impact factors.
                           If False (default), exclude papers without a known IF.
        download_papers: Download PDFs when possible
        scrape_web_content: If True, extract omics from web content without downloading PDFs.
                           Uses Playwright to fetch HTML (PMC first, then publisher DOI).
                           Auto-fallback to PDF download if web scraping finds no accessions.
        subscription_only: If True, skip free papers and only download subscription papers
        use_institutional_auth: If True (default), use institutional access
                               credentials from .env for subscription papers
        max_results: Maximum papers to retrieve
        study_name: Study name for output directory structure. If provided, creates:
                   - output/{study_name}/results.csv
                   - output/{study_name}/downloads/
                   - output/{study_name}/accessions.csv
                   If None, uses timestamped filenames in output/ directory.
        output_format: "csv" or "excel"
        validation_rounds: Number of rounds to run (default 1, 3+ for validation)
                          Round 1: Initial retrieval
                          Round 2+: Validation and proofreading
        exact_match: If True, use disease and tissue terms exactly as provided (no expansion).
                    If False (default), expand terms using Claude-based interpretation.

    Returns:
        DataFrame with results
    """
    if organism is None:
        organism = ["human", "mouse", "rat"]

    tissues_str = ", ".join(tissues) if isinstance(tissues, list) else tissues
    logger.info(f"Starting retrieval: disease='{disease}', tissues='{tissues_str}'")
    if scrape_web_content:
        logger.info("Web scraping mode: Will extract omics from web content (PMC first, then publisher)")
    if subscription_only:
        logger.info("Subscription-only mode: Will skip free access papers")
    if use_institutional_auth:
        logger.info("Institutional auth enabled: Will use credentials for subscription papers")

    # Initialize components
    searcher = PubMedSearcher()
    matcher = EntityMatcher()
    if_lookup = ImpactFactorLookup()

    # Initialize exporter with study-specific directory structure
    exporter = ResultsExporter(study_name=study_name)

    # Initialize downloader with study-specific downloads directory
    downloader = PaperDownloader(
        download_dir=str(exporter.downloads_dir),
        subscription_only=subscription_only
    )
    omics_extractor = OmicsExtractor()

    if study_name:
        logger.info(f"Output directory: {exporter.output_dir}")

    # Handle disease and tissue terms based on exact_match flag
    if exact_match:
        # Use terms exactly as provided (no expansion)
        if isinstance(disease, str):
            # Split comma-separated terms for OR logic
            disease_terms = [d.strip() for d in disease.split(',') if d.strip()]
        else:
            disease_terms = disease
        # Split comma-separated tissues for OR logic
        if isinstance(tissues, list):
            all_tissue_terms = []
            for t in tissues:
                all_tissue_terms.extend([x.strip() for x in t.split(',') if x.strip()])
        else:
            all_tissue_terms = [t.strip() for t in tissues.split(',') if t.strip()]
        logger.info(f"Exact match mode - disease terms: {disease_terms}")
        logger.info(f"Exact match mode - tissue terms: {all_tissue_terms}")
    else:
        # Expand search terms using Claude-based interpretation
        disease_terms = matcher.interpret_diseases(disease) if isinstance(disease, str) else disease
        logger.info(f"Expanded disease terms: {disease_terms}")

        # Expand tissue terms using Claude-based interpretation
        tissue_input = ", ".join(tissues) if isinstance(tissues, list) else tissues
        all_tissue_terms = matcher.interpret_tissues(tissue_input)
        logger.info(f"Expanded tissue terms: {all_tissue_terms}")

    # Log keywords if provided
    if keywords:
        logger.info(f"Keyword filters: {keywords}")

    # Build and execute query with expanded terms
    query = searcher.build_query(
        disease=disease_terms,
        tissues=all_tissue_terms,
        organisms=organism,
        cell_type=cell_type,
        date_cutoff=publication_date_cutoff,
        keywords=keywords
    )

    pmids = searcher.search(query, max_results=max_results)

    if not pmids:
        logger.warning("No publications found")
        return pd.DataFrame()

    # Fetch details
    records = searcher.fetch_details(pmids)
    logger.info(f"Fetched details for {len(records)} publications")

    # Process publications
    publications = []

    for record in records:
        # Score and filter
        match_result = matcher.score_publication(
            record, disease, tissues, organism, cell_type
        )

        if not match_result["matched"]:
            continue

        # Check impact factor
        journal = record.get("journal", "")
        is_preprint = if_lookup.is_preprint(journal)
        impact_factor = if_lookup.get_impact_factor(journal) if not is_preprint else None

        # Apply IF filter with proper handling of preprints and unknown IF
        if is_preprint:
            # Preprints: include only if include_preprints=True
            if not include_preprints:
                continue
        else:
            # Non-preprints: check IF
            if impact_factor is not None:
                # Known IF: must meet cutoff
                if impact_factor < impact_factor_cutoff:
                    continue
            else:
                # Unknown IF: include only if include_unknown_if=True
                if not include_unknown_if:
                    continue

        # Create publication object
        pub = Publication(
            pmid=record["pmid"],
            title=record["title"],
            abstract=record["abstract"],
            authors=record["authors"],
            affiliations=record["affiliations"],
            journal=journal,
            publication_year=record["publication_year"],
            doi=record.get("doi"),
            pmc_id=record.get("pmc_id"),
            matched_disease=match_result["matched_disease"],
            matched_tissue=match_result["matched_tissue"],
            matched_organism=match_result["matched_organism"],
            matched_cell_type=match_result["matched_cell_type"],
            relevance_score=match_result["relevance_score"],
            impact_factor=impact_factor,
            is_preprint=is_preprint,
        )

        # Download paper OR scrape web content
        if download_papers:
            # Standard PDF download mode
            success, download_path, access_url, requires_sub = downloader.download_paper(
                record,
                use_institutional_auth=use_institutional_auth
            )
            pub.free_access = success and not requires_sub
            pub.requires_subscription = requires_sub
            pub.download_path = download_path
            pub.access_url = access_url

            # Extract omics from PDF
            if success and download_path:
                datasets = omics_extractor.extract_from_pdf(download_path)
                pub.omics_datasets = omics_extractor.enrich_datasets(datasets)

        elif scrape_web_content:
            # Web scraping mode - extract from HTML without downloading PDF
            # Decoupled logic: web_scraped = True when full access achieved, regardless of accessions
            pmc_id = record.get("pmc_id")
            doi = pub.doi

            logger.debug(f"PMID {pub.pmid}: Attempting web scraping (PMC: {pmc_id}, DOI: {doi})")

            # Step 1: Try to fetch article HTML to check web access
            article_text = downloader.fetch_article_html(pmc_id, doi)

            if article_text and len(article_text) > 1000:
                # Full web access achieved - mark as web_scraped regardless of accessions
                pub.web_scraped = True
                pub.free_access = True  # We got full text without download

                # Set access URL based on source
                if pmc_id:
                    pmc_clean = pmc_id.replace("PMC", "") if pmc_id.startswith("PMC") else pmc_id
                    pub.access_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_clean}/"
                elif doi:
                    pub.access_url = downloader.get_doi_url(doi) or f"https://doi.org/{doi}"
                else:
                    pub.access_url = f"https://pubmed.ncbi.nlm.nih.gov/{pub.pmid}/"

                # Extract omics from the fetched text (without re-fetching)
                datasets = omics_extractor.extract_from_fetched_text(article_text)
                if datasets:
                    pub.omics_datasets = omics_extractor.enrich_datasets(datasets)
                    logger.info(f"PMID {pub.pmid}: Web access successful, extracted {len(datasets)} accessions")
                else:
                    logger.info(f"PMID {pub.pmid}: Web access successful, no accessions found (no PDF fallback)")

                # DO NOT fall back to PDF when web access succeeded
            else:
                # NO web access - fall back to PDF download
                logger.info(f"PMID {pub.pmid}: No web access, falling back to PDF download")
                success, download_path, access_url, requires_sub = downloader.download_paper(
                    record,
                    use_institutional_auth=use_institutional_auth
                )
                pub.free_access = success and not requires_sub
                pub.requires_subscription = requires_sub
                pub.download_path = download_path
                pub.access_url = access_url

                if success and download_path:
                    datasets = omics_extractor.extract_from_pdf(download_path)
                    pub.omics_datasets = omics_extractor.enrich_datasets(datasets)

            # Check subscription status
            pub.requires_subscription = not downloader.has_free_access(pmc_id)

        else:
            # No full-text extraction mode (abstract only)
            pub.access_url = f"https://pubmed.ncbi.nlm.nih.gov/{pub.pmid}/"
            # Still check if requires subscription
            pub.requires_subscription = not downloader.has_free_access(record.get("pmc_id"))

        # Extract omics from abstract (always, as supplementary source)
        text_datasets = omics_extractor.extract_from_text(f"{pub.title} {pub.abstract}", source="abstract")
        for ds in text_datasets:
            if not any(d.accession == ds.accession for d in pub.omics_datasets):
                pub.omics_datasets.append(ds)

        publications.append(pub)
        logger.info(f"Processed PMID {pub.pmid}: {pub.title[:50]}...")

    logger.info(f"Found {len(publications)} relevant publications")

    # Sort by relevance
    publications.sort(key=lambda p: p.relevance_score, reverse=True)

    # Convert to DataFrame for validation
    df = exporter.publications_to_dataframe(publications)

    # Multi-round validation
    all_validation_results = []
    all_round_dfs = [df]

    if validation_rounds >= 1:
        # Round 1 validation
        logger.info(f"=== Round 1: Initial retrieval complete ({len(df)} papers) ===")
        validation = validate_results(df, 1, if_lookup)
        all_validation_results.append(validation)

        if validation["issues"]:
            logger.warning(f"Round 1 issues: {validation['issues']}")
        if validation["corrections"]:
            logger.info(f"Round 1 corrections: {validation['corrections']}")

    # Additional validation rounds (2+)
    for round_num in range(2, validation_rounds + 1):
        logger.info(f"=== Round {round_num}: Validation and proofreading ===")

        # Re-validate IF values with fresh lookup
        if_lookup.get_impact_factor.cache_clear()

        # Validate this round
        validation = validate_results(df, round_num, if_lookup)
        all_validation_results.append(validation)

        # Compare with previous round
        if len(all_round_dfs) > 0:
            comparison = compare_rounds(all_round_dfs[-1], df, round_num - 1, round_num)
            if not comparison["consistent"]:
                logger.warning(f"Round {round_num} inconsistencies: {comparison['differences']}")

        all_round_dfs.append(df.copy())

        # Report issues and corrections
        if validation["issues"]:
            logger.warning(f"Round {round_num} issues: {validation['issues']}")
        if validation["corrections"]:
            logger.info(f"Round {round_num} corrections: {validation['corrections']}")
        if validation["warnings"]:
            logger.info(f"Round {round_num} warnings: {validation['warnings'][:5]}")  # Limit warnings

        # If all rounds pass validation, break early
        if validation["valid"] and len(validation["issues"]) == 0:
            logger.info(f"Round {round_num}: All validations passed")

    # Final summary
    if validation_rounds > 1:
        total_issues = sum(len(v["issues"]) for v in all_validation_results)
        total_warnings = sum(len(v["warnings"]) for v in all_validation_results)
        total_corrections = sum(len(v.get("corrections", [])) for v in all_validation_results)
        logger.info(f"=== Validation Summary: {validation_rounds} rounds, "
                    f"{total_issues} issues, {total_warnings} warnings, {total_corrections} corrections ===")

    # Export results
    if output_format == "excel":
        output_file = exporter.export_excel(publications)
    else:
        output_file = exporter.export_csv(publications)

    # Export accessions to separate file
    accessions_file = exporter.export_accessions(publications)

    # Export failed studies (papers that couldn't be accessed)
    failed_studies_file = exporter.export_failed_studies(publications)

    # Log output summary
    if study_name:
        logger.info(f"=== Output Summary ===")
        logger.info(f"  Study directory: {exporter.output_dir}")
        logger.info(f"  Results: {output_file}")
        if accessions_file:
            logger.info(f"  Accessions: {accessions_file}")
        if failed_studies_file:
            logger.info(f"  Failed studies: {failed_studies_file}")
        logger.info(f"  Downloads: {exporter.downloads_dir}")

    return df
