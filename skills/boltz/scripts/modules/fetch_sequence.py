import requests

from modules.utils import STANDARD_AA, retry_request, setup_logger

logger = setup_logger("fetch_sequence")

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


def fetch_uniprot_sequence(
    gene_name: str,
    organism_id: int = 9606,
) -> dict:
    params = {
        "query": f"gene_exact:{gene_name} AND organism_id:{organism_id} AND reviewed:true",
        "format": "json",
        "size": 1,
        "fields": "accession,gene_names,sequence,organism_name",
    }

    def _request():
        resp = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=30)
        resp.raise_for_status()
        return resp.json()

    data = retry_request(_request)
    results = data.get("results", [])

    if not results:
        params["query"] = f"gene_exact:{gene_name} AND organism_id:{organism_id}"

        data = retry_request(_request)
        results = data.get("results", [])

    if not results:
        raise ValueError(
            f"Gene '{gene_name}' not found in UniProt for organism_id={organism_id}. "
            "Check spelling or try gene synonyms."
        )

    entry = results[0]
    accession = entry.get("primaryAccession", "")
    sequence_info = entry.get("sequence", {})
    sequence = sequence_info.get("value", "")
    length = sequence_info.get("length", len(sequence))

    genes = entry.get("genes", [])
    resolved_gene = gene_name
    if genes:
        primary = genes[0].get("geneName", {})
        resolved_gene = primary.get("value", gene_name)

    organism = entry.get("organism", {}).get("scientificName", "")

    logger.info(
        "Found %s (%s) for gene %s — %d amino acids",
        accession, organism, resolved_gene, length,
    )

    return {
        "accession": accession,
        "gene_name": resolved_gene,
        "sequence": sequence,
        "organism": organism,
        "length": length,
    }


def validate_sequence(sequence: str) -> bool:
    return bool(sequence) and all(c in STANDARD_AA for c in sequence.upper())
