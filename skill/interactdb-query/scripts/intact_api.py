"""
Query IntAct (PSICQUIC) for neighbor proteins of a gene.

Features
--------
- Data source: IntAct molecular interaction database via PSICQUIC REST,
  retrieving PSI-MITAB 2.7 with IntAct-specific confidence scores
  (e.g. 'intact-miscore'). :contentReference[oaicite:0]{index=0}
- Input: gene symbol (MIQL field `geneName`) and species name keyword
  (e.g. 'human' → taxid:9606). :contentReference[oaicite:1]{index=1}
- Query: MIQL expression of the form
      geneName:<gene> AND species:<species>
  sent to the IntAct PSICQUIC endpoint:
      https://www.ebi.ac.uk/Tools/webservices/psicquic/
      intact/webservices/current/search/query/<MIQL> :contentReference[oaicite:2]{index=2}
- Network: binary protein–protein interactions expanded and distributed
  as MITAB records; multiple confidence scores may be present, including
  'intact-miscore:<value>' in the confidence column. :contentReference[oaicite:3]{index=3}
- Neighbor definition: direct (1-hop) interactors of the queried gene in
  IntAct; each row corresponds to a unique neighbor protein.
- Ranking:
    1) hop distance (all = 1 for this implementation),
    2) best IntAct MIscore per neighbor (descending),
    3) number of supporting interactions (descending).
- Output: pandas DataFrame with one row per neighbor, including:
    * neighbor_id       – primary interactor identifier (e.g. UniProt AC)
    * hop               – shortest hop distance from the seed gene (1 here)
    * neighbor_name     – preferred display / gene name (from MITAB aliases)
    * neighbor_taxid    – NCBI taxon ID of the neighbor
    * neighbor_organism – organism name (e.g. 'Homo sapiens')
    * best_miscore      – maximum 'intact-miscore' across supporting edges
    * n_interactions    – number of IntAct interactions with the seed gene
- Default: return top_n = 100 neighbors after ranking.

Notes
-----
- MITAB 2.7 parsing follows the PSICQUIC / PSI-MI specifications:
  interactor IDs (cols 1–2), aliases (cols 5–6), organisms (cols 10–11),
  and confidence scores (col 15). :contentReference[oaicite:4]{index=4}
- This implementation focuses on 1-hop neighbors; n-hop expansion would
  require additional PSICQUIC calls (e.g. /interactor/<id>) and a BFS
  over the interaction graph.
"""

import requests
import pandas as pd
from urllib.parse import quote
from typing import Dict, List, Tuple, Any, Optional, Set
from collections import defaultdict
import heapq


# IntAct PSICQUIC REST base URL (current spec version)
BASE_URL = (
    "https://www.ebi.ac.uk/Tools/webservices/psicquic/"
    "intact/webservices/current/search"
)


def psicquic_query_by_gene(gene, species="human",
                           max_results=20000,
                           fmt="tab27",
                           timeout=60):
    """
    Query IntAct via PSICQUIC using MIQL:
        geneName:<gene> AND species:<species>

    Returns raw MITAB text.
    """
    miql = f"geneName:{gene} AND species:{species}"
    # Properly URL-encode the MIQL query for use in the path
    encoded_miql = quote(miql)

    url = f"{BASE_URL}/query/{encoded_miql}"
    params = {
        "format": fmt,           # MITAB 2.7
        "firstResult": 0,
        "maxResults": max_results,
    }

    resp = requests.get(url, params=params, timeout=timeout)
    resp.raise_for_status()
    return resp.text


def extract_miscore(conf_str):
    """
    Parse the 'confidence score' column (MITAB col 15) and
    pull out intact-miscore:<float> if present.
    """
    if not conf_str or conf_str == "-":
        return None
    for token in conf_str.split("|"):
        token = token.strip().strip('"')
        if token.startswith("intact-miscore:"):
            val = token.split(":", 1)[1]
            try:
                return float(val)
            except ValueError:
                return None
    return None


def parse_taxon(tax_str):
    """
    Parse 'taxid:9606(Homo sapiens)' → ('9606', 'Homo sapiens')
    from MITAB cols 10 / 11.
    """
    if not tax_str or tax_str == "-":
        return None, None

    token = tax_str.split("|")[0].strip()
    if not token.startswith("taxid:"):
        return None, None

    rest = token[len("taxid:"):]
    if "(" in rest:
        tax_id, org = rest.split("(", 1)
        return tax_id.strip(), org.rstrip(")").strip()
    else:
        return rest.strip(), None


def pick_name_from_alias(alias_str):
    """
    MITAB cols 5 / 6 = Aliases:
        uniprotkb:TP53(gene name)|uniprotkb:P53(display_short)|...

    Try to pick a nice display name:
    - Prefer 'display_short' or 'gene name' or 'recommended name'
    - Otherwise, return the first alias's 'name' part.
    """
    if not alias_str or alias_str == "-":
        return None

    best_fallback = None

    for token in alias_str.split("|"):
        token = token.strip()
        if not token:
            continue

        alias_type = ""
        if "(" in token and token.endswith(")"):
            base, alias_type = token.rsplit("(", 1)
            alias_type = alias_type[:-1]
        else:
            base = token

        if ":" in base:
            _, name = base.split(":", 1)
        else:
            name = base

        name = name.strip()
        if alias_type in ("display_short", "gene name", "recommended name"):
            return name

        if best_fallback is None:
            best_fallback = name

    return best_fallback


def pick_primary_id(id_str):
    """
    MITAB col 1 / 2 = Unique identifier, e.g.
        'uniprotkb:P04637' or 'uniprotkb:P04637|refseq:NP_...'

    We pick the first database:identifier pair, and return just the accession,
    e.g. 'P04637'.
    """
    if not id_str or id_str == "-":
        return None

    token = id_str.split("|")[0].strip()
    if ":" in token:
        _, acc = token.split(":", 1)
        return acc.strip()
    return token.strip()


def extract_readable_term(term_str):
    """
    Extract readable term from PSI-MI ontology format.
    E.g. 'psi-mi:"MI:0018"(two hybrid)' -> 'two hybrid'
         'psi-mi:"MI:0915"(physical association)' -> 'physical association'
    """
    if not term_str or term_str == "-":
        return None

    # Take first term if multiple
    token = term_str.split("|")[0].strip()

    # Extract text in parentheses
    if "(" in token and token.endswith(")"):
        return token.split("(", 1)[1].rstrip(")")

    return token


def extract_pmids(pub_str):
    """
    Extract PubMed IDs from publication column.
    E.g. 'pubmed:12345|pubmed:67890' -> '12345|67890'
    """
    if not pub_str or pub_str == "-":
        return None

    pmids = []
    for token in pub_str.split("|"):
        token = token.strip()
        if token.startswith("pubmed:"):
            pmid = token.split(":", 1)[1]
            if pmid and pmid not in pmids:
                pmids.append(pmid)

    return "|".join(pmids) if pmids else None


def parse_mitab(mitab_text):
    """
    Parse a MITAB 2.7 text block into a list of interaction dicts
    with the fields we care about.

    Key MITAB 2.7 columns:
    0-1: Interactor IDs
    4-5: Aliases (gene names)
    6: Detection method
    8: Publication IDs
    9-10: Taxon IDs
    11: Interaction type
    12: Source database
    14: Confidence scores
    """
    interactions = []

    for line in mitab_text.splitlines():
        if not line or line.startswith("#"):
            continue

        cols = line.split("\t")
        # Need at least up to confidence score (col 15, index 14)
        if len(cols) < 15:
            continue

        idA_raw, idB_raw = cols[0], cols[1]
        aliasA, aliasB = cols[4], cols[5]
        detection_method = cols[6] if len(cols) > 6 else "-"
        publication_ids = cols[8] if len(cols) > 8 else "-"
        taxA, taxB = cols[9], cols[10]
        interaction_type = cols[11] if len(cols) > 11 else "-"
        source_db = cols[12] if len(cols) > 12 else "-"
        conf = cols[14]

        miscore = extract_miscore(conf)
        taxidA, orgA = parse_taxon(taxA)
        taxidB, orgB = parse_taxon(taxB)
        nameA = pick_name_from_alias(aliasA)
        nameB = pick_name_from_alias(aliasB)

        interactions.append(
            {
                "idA_raw": idA_raw,
                "idB_raw": idB_raw,
                "idA": pick_primary_id(idA_raw),
                "idB": pick_primary_id(idB_raw),
                "nameA": nameA,
                "nameB": nameB,
                "taxidA": taxidA,
                "organismA": orgA,
                "taxidB": taxidB,
                "organismB": orgB,
                "miscore": miscore,
                "detection_method": detection_method,
                "publication_ids": publication_ids,
                "interaction_type": interaction_type,
                "source_database": source_db,
            }
        )

    return interactions


def get_direct_neighbors(
    gene,
    species="human",
    top_n=100,
    miql_max_results=20000,
    organism_filter=None,
):
    """
    Get direct (1-hop) IntAct neighbors for a human gene.

    Ranking:
      1) all are hop = 1
      2) sort by best IntAct MIscore (descending)
      3) tie-breaker: number of supporting interactions

    Returns a pandas DataFrame with columns:
      ['neighbor_id', 'hop', 'neighbor_name',
       'neighbor_taxid', 'neighbor_organism',
       'best_miscore', 'n_interactions']
    """
    mitab = psicquic_query_by_gene(
        gene, species=species, max_results=miql_max_results
    )
    interactions = parse_mitab(mitab)

    gene_upper = gene.upper()
    neighbors = {}

    # Parse organism filter if provided
    filter_organisms = None
    if organism_filter:
        if isinstance(organism_filter, str):
            filter_organisms = {o.strip().lower() for o in organism_filter.split(",")}
        else:
            filter_organisms = {str(o).lower() for o in organism_filter}

    for it in interactions:
        miscore = it["miscore"] if it["miscore"] is not None else 0.0

        # Extract edge annotations
        detection = extract_readable_term(it.get("detection_method", "-"))
        interaction_type = extract_readable_term(it.get("interaction_type", "-"))
        pmids = extract_pmids(it.get("publication_ids", "-"))
        source_db = it.get("source_database", "IntAct")

        # Check if gene is interactor A → neighbor = B
        if it["nameA"] and it["nameA"].upper() == gene_upper:
            nid = it["idB"]

            # Apply organism filter
            if filter_organisms and it["organismB"]:
                if it["organismB"].lower() not in filter_organisms:
                    continue

            if nid:
                rec = neighbors.get(nid)
                if rec is None:
                    neighbors[nid] = {
                        "neighbor_id": nid,
                        "hop": 1,
                        "neighbor_name": it["nameB"],
                        "neighbor_taxid": it["taxidB"],
                        "neighbor_organism": it["organismB"],
                        "best_miscore": miscore,
                        "n_interactions": 1,
                        "detection_methods": detection,
                        "interaction_types": interaction_type,
                        "publications": pmids,
                        "source_database": source_db,
                        "path": f"{gene}-{it['nameB']}" if it["nameB"] else None,
                    }
                else:
                    rec["n_interactions"] += 1
                    if miscore > rec["best_miscore"]:
                        rec["best_miscore"] = miscore
                    if not rec["neighbor_name"] and it["nameB"]:
                        rec["neighbor_name"] = it["nameB"]
                    if not rec["neighbor_organism"] and it["organismB"]:
                        rec["neighbor_organism"] = it["organismB"]
                        rec["neighbor_taxid"] = it["taxidB"]

                    # Accumulate unique edge annotations
                    if detection and detection not in (rec["detection_methods"] or ""):
                        rec["detection_methods"] = f"{rec['detection_methods']}|{detection}" if rec["detection_methods"] else detection
                    if interaction_type and interaction_type not in (rec["interaction_types"] or ""):
                        rec["interaction_types"] = f"{rec['interaction_types']}|{interaction_type}" if rec["interaction_types"] else interaction_type
                    if pmids:
                        existing_pmids = set((rec["publications"] or "").split("|")) if rec["publications"] else set()
                        new_pmids = set(pmids.split("|"))
                        all_pmids = existing_pmids | new_pmids
                        rec["publications"] = "|".join(sorted(all_pmids))

        # Check if gene is interactor B → neighbor = A
        if it["nameB"] and it["nameB"].upper() == gene_upper:
            nid = it["idA"]

            # Apply organism filter
            if filter_organisms and it["organismA"]:
                if it["organismA"].lower() not in filter_organisms:
                    continue

            if nid:
                rec = neighbors.get(nid)
                if rec is None:
                    neighbors[nid] = {
                        "neighbor_id": nid,
                        "hop": 1,
                        "neighbor_name": it["nameA"],
                        "neighbor_taxid": it["taxidA"],
                        "neighbor_organism": it["organismA"],
                        "best_miscore": miscore,
                        "n_interactions": 1,
                        "detection_methods": detection,
                        "interaction_types": interaction_type,
                        "publications": pmids,
                        "source_database": source_db,
                        "path": f"{gene}-{it['nameA']}" if it["nameA"] else None,
                    }
                else:
                    rec["n_interactions"] += 1
                    if miscore > rec["best_miscore"]:
                        rec["best_miscore"] = miscore
                    if not rec["neighbor_name"] and it["nameA"]:
                        rec["neighbor_name"] = it["nameA"]
                    if not rec["neighbor_organism"] and it["organismA"]:
                        rec["neighbor_organism"] = it["organismA"]
                        rec["neighbor_taxid"] = it["taxidA"]

                    # Accumulate unique edge annotations
                    if detection and detection not in (rec["detection_methods"] or ""):
                        rec["detection_methods"] = f"{rec['detection_methods']}|{detection}" if rec["detection_methods"] else detection
                    if interaction_type and interaction_type not in (rec["interaction_types"] or ""):
                        rec["interaction_types"] = f"{rec['interaction_types']}|{interaction_type}" if rec["interaction_types"] else interaction_type
                    if pmids:
                        existing_pmids = set((rec["publications"] or "").split("|")) if rec["publications"] else set()
                        new_pmids = set(pmids.split("|"))
                        all_pmids = existing_pmids | new_pmids
                        rec["publications"] = "|".join(sorted(all_pmids))

    if not neighbors:
        return pd.DataFrame(
            columns=[
                "neighbor_id",
                "hop",
                "neighbor_name",
                "neighbor_taxid",
                "neighbor_organism",
                "best_miscore",
                "n_interactions",
                "detection_methods",
                "interaction_types",
                "publications",
                "source_database",
                "path",
            ]
        )

    df = pd.DataFrame(neighbors.values())

    # Rank: hop (all 1) then MIscore desc, then number of interactions desc
    df = df.sort_values(
        by=["hop", "best_miscore", "n_interactions"],
        ascending=[True, False, False],
    )

    return df.head(top_n)


def get_neighbors_multihop(
    gene,
    species="human",
    top_n=100,
    max_hops=3,
    min_miscore=0.0,
    miql_max_results=20000,
    organism_filter=None,
):
    """
    Get IntAct neighbors with automatic multi-hop BFS expansion.

    This function iteratively expands from 1-hop to 2-hop to 3-hop neighbors
    until top_n neighbors are found or max_hops is reached.

    **Algorithm**: Breadth-First Search (BFS) with iterative expansion
    **Principle**: If 1-hop neighbors < top_n, expand to 2-hop, then 3-hop, etc.

    Parameters
    ----------
    gene : str
        Gene name to query
    species : str
        Species keyword for MIQL query (default: "human")
    top_n : int
        Target number of neighbors to return (default: 100)
    max_hops : int
        Maximum number of hops to expand (default: 3)
    min_miscore : float
        Minimum IntAct MI-score threshold (0.0-1.0, default: 0.0)
    miql_max_results : int
        Maximum results per PSICQUIC query (default: 20000)
    organism_filter : str, optional
        Filter by organism (e.g., "homo sapiens" or "homo sapiens,mus musculus")

    Returns
    -------
    pd.DataFrame
        Neighbors sorted by (hop, best_miscore desc, n_interactions desc)
        Each neighbor includes hop distance and path from seed gene

    Examples
    --------
    >>> neighbors = get_neighbors_multihop("GREM1", top_n=100, max_hops=3, min_miscore=0.4)
    >>> print(f"Found {len(neighbors)} neighbors")
    >>> print(f"1-hop: {sum(neighbors['hop'] == 1)}")
    >>> print(f"2-hop: {sum(neighbors['hop'] == 2)}")
    """
    gene = gene.strip()
    gene_upper = gene.upper()

    # Track all neighbors across hops
    all_neighbors: Dict[str, Dict] = {}
    visited: Set[str] = {gene_upper}

    # BFS frontier (genes to expand from)
    frontier: Set[str] = {gene_upper}

    # Track paths: gene_name -> [seed, intermediate, ..., gene]
    paths: Dict[str, List[str]] = {gene_upper: [gene]}

    # Parse organism filter once
    filter_organisms = None
    if organism_filter:
        if isinstance(organism_filter, str):
            filter_organisms = {o.strip().lower() for o in organism_filter.split(",")}
        else:
            filter_organisms = {str(o).lower() for o in organism_filter}

    for hop in range(1, max_hops + 1):
        if not frontier:
            break

        # Query neighbors of all genes in current frontier
        new_frontier: Set[str] = set()

        for current_gene in frontier:
            # Get 1-hop neighbors of current gene
            try:
                mitab = psicquic_query_by_gene(
                    current_gene, species=species, max_results=miql_max_results
                )
                interactions = parse_mitab(mitab)
            except Exception:
                continue

            for it in interactions:
                miscore = it["miscore"] if it["miscore"] is not None else 0.0

                # Apply MI-score filter
                if miscore < min_miscore:
                    continue

                # Extract edge annotations
                detection = extract_readable_term(it.get("detection_method", "-"))
                interaction_type = extract_readable_term(it.get("interaction_type", "-"))
                pmids = extract_pmids(it.get("publication_ids", "-"))
                source_db = it.get("source_database", "IntAct")

                # Process both sides of the interaction
                for idx, (name_key, id_key, taxid_key, org_key) in enumerate([
                    ("nameA", "idA", "taxidA", "organismA"),
                    ("nameB", "idB", "taxidB", "organismB")
                ]):
                    neighbor_name = it[name_key]
                    neighbor_id = it[id_key]
                    neighbor_taxid = it[taxid_key]
                    neighbor_organism = it[org_key]

                    if not neighbor_name or not neighbor_id:
                        continue

                    neighbor_name_upper = neighbor_name.upper()

                    # Skip seed gene
                    if neighbor_name_upper == gene_upper:
                        continue

                    # Skip if already visited
                    if neighbor_name_upper in visited:
                        continue

                    # Apply organism filter
                    if filter_organisms and neighbor_organism:
                        if neighbor_organism.lower() not in filter_organisms:
                            continue

                    # Build path for this neighbor
                    if current_gene in paths:
                        neighbor_path = paths[current_gene] + [neighbor_name]
                    else:
                        neighbor_path = [gene, neighbor_name]

                    paths[neighbor_name_upper] = neighbor_path

                    # Add or update neighbor record
                    if neighbor_name_upper not in all_neighbors:
                        all_neighbors[neighbor_name_upper] = {
                            "neighbor_id": neighbor_id,
                            "hop": hop,
                            "neighbor_name": neighbor_name,
                            "neighbor_taxid": neighbor_taxid,
                            "neighbor_organism": neighbor_organism,
                            "best_miscore": miscore,
                            "n_interactions": 1,
                            "detection_methods": detection,
                            "interaction_types": interaction_type,
                            "publications": pmids,
                            "source_database": source_db,
                            "path": "-".join(neighbor_path),
                        }
                    else:
                        rec = all_neighbors[neighbor_name_upper]
                        rec["n_interactions"] += 1
                        if miscore > rec["best_miscore"]:
                            rec["best_miscore"] = miscore

                        # Accumulate unique edge annotations
                        if detection and detection not in (rec["detection_methods"] or ""):
                            rec["detection_methods"] = f"{rec['detection_methods']}|{detection}" if rec["detection_methods"] else detection
                        if interaction_type and interaction_type not in (rec["interaction_types"] or ""):
                            rec["interaction_types"] = f"{rec['interaction_types']}|{interaction_type}" if rec["interaction_types"] else interaction_type
                        if pmids:
                            existing_pmids = set((rec["publications"] or "").split("|")) if rec["publications"] else set()
                            new_pmids = set(pmids.split("|"))
                            all_pmids = existing_pmids | new_pmids
                            rec["publications"] = "|".join(sorted(all_pmids))

                    # Add to new frontier for next hop expansion
                    new_frontier.add(neighbor_name_upper)
                    visited.add(neighbor_name_upper)

        frontier = new_frontier

        # Stop if we have enough neighbors
        if len(all_neighbors) >= top_n:
            break

    if not all_neighbors:
        return pd.DataFrame(
            columns=[
                "neighbor_id",
                "hop",
                "neighbor_name",
                "neighbor_taxid",
                "neighbor_organism",
                "best_miscore",
                "n_interactions",
                "detection_methods",
                "interaction_types",
                "publications",
                "source_database",
                "path",
            ]
        )

    # Sort by hop (ascending), then by best_miscore (descending), then by n_interactions (descending)
    df = pd.DataFrame(all_neighbors.values())
    df = df.sort_values(
        by=["hop", "best_miscore", "n_interactions"],
        ascending=[True, False, False],
    )

    return df.head(top_n)


def find_shortest_paths_intact(
    gene_list: List[str],
    species: str = "human",
    max_distance: int = 3,
    min_miscore: float = 0.4,
    organism_filter: Optional[str] = None,
    miql_max_results: int = 50000,
) -> Dict[Tuple[str, str], Dict[str, Any]]:
    """
    Find shortest paths between all pairs of query genes using IntAct network.

    **Algorithm**: Dijkstra's shortest path with score-based edge weighting
    **Edge weight formula**: weight = 1.0 - miscore
    **Principle**: Higher MI-score = lower weight = shorter distance (closer nodes)

    Parameters
    ----------
    gene_list : List[str]
        List of gene names to find paths between
    species : str
        Species keyword for MIQL query (default: "human")
    max_distance : int
        Maximum number of intermediate nodes (default: 3)
    min_miscore : float
        Minimum IntAct MI-score for edges (0.0-1.0, default: 0.4)
    organism_filter : str, optional
        Filter by organism (e.g., "human" or "human,mouse")
    miql_max_results : int
        Maximum results from PSICQUIC query (default: 50000)

    Returns
    -------
    Dict[Tuple[str, str], Dict]
        Dictionary mapping (geneA, geneB) tuples to path info:
        {
            'path': ['A', 'D', 'B'],
            'distance': 1.25,
            'hops': 2,
            'scores': [0.85, 0.92],  # MI-scores along path
            'algorithm': 'Dijkstra'
        }
    """
    if len(gene_list) < 2:
        return {}

    # Query IntAct for all genes together
    miql_genes = " OR ".join([f"geneName:{g}" for g in gene_list])
    miql = f"({miql_genes}) AND species:{species}"
    encoded_miql = quote(miql)

    url = f"{BASE_URL}/query/{encoded_miql}"
    params = {
        "format": "tab27",
        "firstResult": 0,
        "maxResults": miql_max_results,
    }

    try:
        resp = requests.get(url, params=params, timeout=120)
        resp.raise_for_status()
        mitab = resp.text
    except Exception as e:
        print(f"Error querying IntAct: {e}")
        return {}

    interactions = parse_mitab(mitab)

    if not interactions:
        return {}

    # Build weighted graph: edge weight = 1 - miscore (lower score = higher cost)
    graph: Dict[str, Dict[str, float]] = defaultdict(dict)
    edge_scores: Dict[Tuple[str, str], float] = {}
    gene_set = {g.upper() for g in gene_list}

    # Parse organism filter
    filter_organisms = None
    if organism_filter:
        filter_organisms = {o.strip().lower() for o in organism_filter.split(",")}

    for it in interactions:
        name_a = it["nameA"]
        name_b = it["nameB"]
        organism_a = it["organismA"]
        organism_b = it["organismB"]
        miscore = it["miscore"] if it["miscore"] is not None else 0.0

        # Apply filters
        if miscore < min_miscore:
            continue

        if not name_a or not name_b:
            continue

        # Apply organism filter
        if filter_organisms:
            if organism_a and organism_a.lower() not in filter_organisms:
                continue
            if organism_b and organism_b.lower() not in filter_organisms:
                continue

        # Edge weight: higher MI-score = lower distance
        weight = 1.0 - miscore

        # Add edge (undirected)
        name_a_upper = name_a.upper()
        name_b_upper = name_b.upper()

        graph[name_a_upper][name_b_upper] = weight
        graph[name_b_upper][name_a_upper] = weight

        # Store score for reporting
        edge_key = tuple(sorted([name_a_upper, name_b_upper]))
        # Keep best (highest) score if duplicate edge
        if edge_key not in edge_scores or miscore > edge_scores[edge_key]:
            edge_scores[edge_key] = miscore

    # Find shortest paths between all pairs using Dijkstra
    results = {}
    query_genes_upper = [g.upper() for g in gene_list]

    for i, gene_a in enumerate(query_genes_upper):
        for gene_b in query_genes_upper[i+1:]:
            path, distance, scores = _dijkstra_path_intact(
                graph, gene_a, gene_b, edge_scores, max_distance
            )

            if path:
                results[(gene_a, gene_b)] = {
                    'path': path,
                    'distance': distance,
                    'hops': len(path) - 1,
                    'scores': scores,
                    'algorithm': 'Dijkstra',
                    'weight_formula': '1.0 - miscore'
                }

    return results


def _dijkstra_path_intact(
    graph: Dict[str, Dict[str, float]],
    start: str,
    end: str,
    edge_scores: Dict[Tuple[str, str], float],
    max_hops: int
) -> Tuple[Optional[List[str]], float, List[float]]:
    """
    Dijkstra's algorithm for IntAct network.

    Returns (path, total_distance, mi_scores_along_path)
    """
    if start not in graph or end not in graph:
        return None, float('inf'), []

    # Priority queue: (distance, node, path, hop_count)
    pq = [(0.0, start, [start], 0)]
    visited = set()

    while pq:
        dist, node, path, hops = heapq.heappop(pq)

        if node in visited:
            continue

        visited.add(node)

        if node == end:
            # Extract MI-scores along path
            scores = []
            for i in range(len(path) - 1):
                edge_key = tuple(sorted([path[i], path[i+1]]))
                score = edge_scores.get(edge_key, 0.0)
                scores.append(score)

            return path, dist, scores

        if hops >= max_hops:
            continue

        for neighbor, weight in graph.get(node, {}).items():
            if neighbor not in visited:
                heapq.heappush(
                    pq,
                    (dist + weight, neighbor, path + [neighbor], hops + 1)
                )

    return None, float('inf'), []


if __name__ == "__main__":
    # Example: get top 100 neighbors of TP53 in human
    gene = "TP53"
    neighbors_df = get_direct_neighbors(gene, species="human", top_n=100)
    print(neighbors_df.head(20))
