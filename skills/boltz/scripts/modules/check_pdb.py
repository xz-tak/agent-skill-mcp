import json
from pathlib import Path
from typing import Optional

import requests

from modules.utils import retry_request, setup_logger

logger = setup_logger("check_pdb")

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA_URL = "https://data.rcsb.org/rest/v1/core/entry"
RCSB_ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity"
RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download"


def search_pdb_by_gene(
    gene_name: str,
    organism_id: int = 9606,
) -> list[dict]:
    query_body = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                        "operator": "exact_match",
                        "value": gene_name,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.taxonomy_lineage.id",
                        "operator": "exact_match",
                        "value": str(organism_id),
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "sort": [
                {
                    "sort_by": "rcsb_entry_info.resolution_combined",
                    "direction": "asc",
                }
            ],
            "paginate": {"start": 0, "rows": 10},
        },
    }

    def _request():
        resp = requests.post(
            RCSB_SEARCH_URL,
            json=query_body,
            timeout=30,
        )
        if resp.status_code == 204:
            return {"result_set": []}
        resp.raise_for_status()
        return resp.json()

    data = retry_request(_request)
    result_set = data.get("result_set", [])

    entries = []
    for entry in result_set:
        pdb_id = entry.get("identifier", "")
        score = entry.get("score", 0)
        entries.append({"pdb_id": pdb_id, "score": score})

    if entries:
        entry_meta = _get_entry_metadata(entries[0]["pdb_id"])
        entity_meta = _get_entity_metadata(entries[0]["pdb_id"], gene_name)
        entries[0]["resolution"] = entry_meta["resolution"]
        entries[0]["title"] = entry_meta["title"]
        entries[0]["method"] = entry_meta["method"]
        entries[0]["deposition_date"] = entry_meta["deposition_date"]
        entries[0].update(entity_meta)

    logger.info(
        "Found %d PDB entries for gene %s (organism_id=%d)",
        len(entries), gene_name, organism_id,
    )
    return entries


def _get_entry_metadata(pdb_id: str) -> dict:
    metadata = {
        "resolution": None,
        "title": None,
        "method": None,
        "deposition_date": None,
    }
    try:
        def _request():
            resp = requests.get(f"{RCSB_DATA_URL}/{pdb_id}", timeout=15)
            resp.raise_for_status()
            return resp.json()

        data = retry_request(_request, max_retries=2)
        resolutions = (
            data.get("rcsb_entry_info", {})
            .get("resolution_combined", [])
        )
        if resolutions:
            metadata["resolution"] = resolutions[0]
        metadata["title"] = data.get("struct", {}).get("title")
        exptl = data.get("exptl", [])
        if exptl:
            metadata["method"] = exptl[0].get("method")
        metadata["deposition_date"] = (
            data.get("pdbx_database_status", {})
            .get("recvd_initial_deposition_date")
        )
    except Exception:
        pass
    return metadata


def _get_entity_metadata(pdb_id: str, gene_name: str) -> dict:
    metadata = {
        "entity_sequence_length": None,
        "entity_description": None,
        "entity_gene_name": None,
        "entity_organism": None,
        "entity_weight_kda": None,
    }
    try:
        for entity_id in range(1, 10):
            def _request(eid=entity_id):
                resp = requests.get(
                    f"{RCSB_ENTITY_URL}/{pdb_id}/{eid}", timeout=15,
                )
                resp.raise_for_status()
                return resp.json()

            try:
                data = retry_request(_request, max_retries=1)
            except Exception:
                break

            src_list = data.get("entity_src_gen", [])
            entity_gene = None
            for src in src_list:
                entity_gene = src.get("pdbx_gene_src_gene")
                if entity_gene:
                    break

            if entity_gene and gene_name.upper() in entity_gene.upper():
                poly = data.get("entity_poly", {})
                metadata["entity_sequence_length"] = poly.get(
                    "rcsb_sample_sequence_length",
                )
                rcsb_entity = data.get("rcsb_polymer_entity", {})
                metadata["entity_description"] = rcsb_entity.get(
                    "pdbx_description",
                )
                weight = rcsb_entity.get("formula_weight")
                if weight is not None:
                    metadata["entity_weight_kda"] = round(weight, 3)
                metadata["entity_gene_name"] = entity_gene
                for src in src_list:
                    org = src.get("pdbx_gene_src_scientific_name")
                    if org:
                        metadata["entity_organism"] = org
                        break
                break
    except Exception:
        pass
    return metadata


def download_pdb(
    pdb_id: str,
    output_dir: Path,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{pdb_id}.pdb"

    def _request():
        resp = requests.get(
            f"{RCSB_DOWNLOAD_URL}/{pdb_id}.pdb",
            timeout=60,
        )
        resp.raise_for_status()
        return resp.content

    content = retry_request(_request)
    output_path.write_bytes(content)
    logger.info("Downloaded PDB %s to %s", pdb_id, output_path)
    return output_path


def check_and_download_best_pdb(
    gene_name: str,
    base_dir: Path,
    organism_id: int = 9606,
    pdb_dir_suffix: str = "_pdb",
) -> Optional[dict]:
    entries = search_pdb_by_gene(gene_name, organism_id)
    if not entries:
        logger.info("No PDB structures found for %s", gene_name)
        return None

    best = entries[0]
    pdb_dir = base_dir / f"{gene_name}{pdb_dir_suffix}"

    try:
        file_path = download_pdb(best["pdb_id"], pdb_dir)
    except Exception as e:
        logger.error("Failed to download PDB %s: %s", best["pdb_id"], e)
        return None

    return {
        "pdb_id": best["pdb_id"],
        "resolution": best.get("resolution"),
        "file_path": str(file_path),
        "title": best.get("title"),
        "method": best.get("method"),
        "deposition_date": best.get("deposition_date"),
        "entity_sequence_length": best.get("entity_sequence_length"),
        "entity_description": best.get("entity_description"),
        "entity_gene_name": best.get("entity_gene_name"),
        "entity_organism": best.get("entity_organism"),
        "entity_weight_kda": best.get("entity_weight_kda"),
    }
