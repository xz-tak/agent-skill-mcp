"""
NCBI Paper Query - PubMed Searcher

Handles PubMed/NCBI Entrez API queries for publication retrieval.
"""

import time
import logging
from typing import List, Dict, Optional

from Bio import Entrez

from ..config import NCBI_EMAIL, NCBI_API_KEY

logger = logging.getLogger(__name__)


class PubMedSearcher:
    """Handles PubMed/NCBI Entrez API queries."""

    def __init__(self, email: str = NCBI_EMAIL, api_key: str = None):
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.rate_limit_delay = 0.34 if api_key else 1.0  # 3/sec with key, 1/sec without

    def build_query(
        self,
        disease: str,
        tissues: List[str],
        organisms: List[str] = None,
        cell_type: str = None,
        date_cutoff: int = 2015,
        keywords: List[str] = None
    ) -> str:
        """Build an NCBI Entrez query string."""
        query_parts = []

        # Disease (main criteria - supports multiple diseases with OR)
        if isinstance(disease, list):
            if len(disease) == 1:
                query_parts.append(f'("{disease[0]}"[Title/Abstract])')
            else:
                disease_query = " OR ".join([f'"{d}"[Title/Abstract]' for d in disease])
                query_parts.append(f"({disease_query})")
        else:
            query_parts.append(f'("{disease}"[Title/Abstract])')

        # Tissue filter (supports multiple tissues with OR)
        if tissues:
            if len(tissues) == 1:
                query_parts.append(f'("{tissues[0]}"[Title/Abstract])')
            else:
                tissue_query = " OR ".join([f'"{t}"[Title/Abstract]' for t in tissues])
                query_parts.append(f"({tissue_query})")

        # Organism filter - use BOTH [Organism] AND [MeSH] with OR logic
        # [Organism] = sequence taxonomy; [MeSH] = indexed study subjects
        # Either match fulfills the requirement
        if organisms:
            mesh_map = {
                "human": "Humans",
                "mouse": "Mice",
                "rat": "Rats",
                "zebrafish": "Zebrafish",
                "pig": "Swine",
                "monkey": "Primates",
                "dog": "Dogs",
            }
            org_parts = []
            for org in organisms:
                org_lower = org.lower()
                mesh_term = mesh_map.get(org_lower, org.capitalize() + "s")
                # Either [Organism] OR [MeSH] match passes
                org_parts.append(f'("{org_lower}"[Organism] OR "{mesh_term}"[MeSH])')
            org_query = " OR ".join(org_parts)
            query_parts.append(f"({org_query})")

        # Cell type if specified
        if cell_type:
            query_parts.append(f'("{cell_type}"[Title/Abstract])')

        # Keywords filter (OR logic) - focuses on omics/technology terms
        if keywords:
            if len(keywords) == 1:
                query_parts.append(f'("{keywords[0]}"[Title/Abstract])')
            else:
                keywords_query = " OR ".join([f'"{kw}"[Title/Abstract]' for kw in keywords])
                query_parts.append(f"({keywords_query})")

        # Date filter
        query_parts.append(f'("{date_cutoff}"[Date - Publication] : "3000"[Date - Publication])')

        return " AND ".join(query_parts)

    def search(self, query: str, max_results: int = 1000) -> List[str]:
        """Search PubMed and return list of PMIDs."""
        logger.info(f"Searching PubMed with query: {query[:100]}...")

        try:
            handle = Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=max_results,
                sort="relevance"
            )
            record = Entrez.read(handle)
            handle.close()

            pmids = record.get("IdList", [])
            logger.info(f"Found {len(pmids)} publications")
            return pmids

        except Exception as e:
            logger.error(f"PubMed search failed: {e}")
            return []

    def fetch_details(self, pmids: List[str], batch_size: int = 100) -> List[Dict]:
        """Fetch detailed metadata for a list of PMIDs."""
        all_records = []

        for i in range(0, len(pmids), batch_size):
            batch = pmids[i:i + batch_size]
            logger.info(f"Fetching details for PMIDs {i+1}-{i+len(batch)}...")

            try:
                handle = Entrez.efetch(
                    db="pubmed",
                    id=",".join(batch),
                    rettype="xml",
                    retmode="xml"
                )
                records = Entrez.read(handle)
                handle.close()

                for article in records.get("PubmedArticle", []):
                    parsed = self._parse_article(article)
                    if parsed:
                        all_records.append(parsed)

                time.sleep(self.rate_limit_delay)

            except Exception as e:
                logger.error(f"Failed to fetch batch: {e}")
                continue

        return all_records

    def _parse_article(self, article: Dict) -> Optional[Dict]:
        """Parse a PubMed article XML record into a dictionary."""
        try:
            medline = article.get("MedlineCitation", {})
            article_data = medline.get("Article", {})

            # Basic info
            pmid = str(medline.get("PMID", ""))
            title = str(article_data.get("ArticleTitle", ""))

            # Abstract
            abstract_parts = article_data.get("Abstract", {}).get("AbstractText", [])
            if isinstance(abstract_parts, list):
                abstract = " ".join([str(p) for p in abstract_parts])
            else:
                abstract = str(abstract_parts)

            # Journal
            journal_info = article_data.get("Journal", {})
            journal = journal_info.get("Title", "")

            # Publication year
            pub_date = article_data.get("ArticleDate", [{}])
            if pub_date:
                year = int(pub_date[0].get("Year", 0))
            else:
                journal_date = journal_info.get("JournalIssue", {}).get("PubDate", {})
                year = int(journal_date.get("Year", 0))

            # Authors
            author_list = article_data.get("AuthorList", [])
            authors = []
            affiliations = []
            for author in author_list:
                last = author.get("LastName", "")
                first = author.get("ForeName", "")
                if last:
                    authors.append(f"{last} {first}".strip())

                aff_info = author.get("AffiliationInfo", [])
                for aff in aff_info:
                    aff_str = aff.get("Affiliation", "")
                    if aff_str and aff_str not in affiliations:
                        affiliations.append(aff_str)

            # DOI and PMC
            doi = None
            pmc_id = None
            id_list = article.get("PubmedData", {}).get("ArticleIdList", [])
            for id_item in id_list:
                id_type = id_item.attributes.get("IdType", "")
                if id_type == "doi":
                    doi = str(id_item)
                elif id_type == "pmc":
                    pmc_id = str(id_item)

            return {
                "pmid": pmid,
                "title": title,
                "abstract": abstract,
                "journal": journal,
                "publication_year": year,
                "authors": authors,
                "affiliations": affiliations,
                "doi": doi,
                "pmc_id": pmc_id
            }

        except Exception as e:
            logger.error(f"Failed to parse article: {e}")
            return None
