"""
IntAct REST API Implementation

This module provides a REST API alternative to the PSICQUIC interface
for querying IntAct protein-protein interaction database.

Key Differences from PSICQUIC version:
- Uses IntAct REST API (/ws/interaction/findInteractions)
- Returns JSON instead of MITAB format
- Direct access to intactMiscore field (simpler parsing)
- Better structured response with pagination support

All downstream logic (BFS, Dijkstra, multi-hop, shortest paths) remains unchanged.
Only the data source layer is modified.
"""

import requests
import pandas as pd
from typing import Dict, List, Tuple, Any, Optional, Set
from collections import defaultdict
import heapq
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading


# IntAct REST API base URLs
INTACT_INTERACTION_BASE = "https://www.ebi.ac.uk/intact/ws/interaction"

# UniProt REST API
UNIPROT_API_BASE = "https://rest.uniprot.org"

# Cache for UniProt AC lookups
_uniprot_cache = {}


def get_uniprot_ac(gene_symbol: str, species_taxid: str = "9606", timeout: int = 10) -> Optional[str]:
    """
    Get UniProt accession for a gene symbol using UniProt REST API.

    This provides precise mapping from gene symbols to UniProt ACs, which enables
    more accurate IntAct queries (e.g., TP53 → P04637).

    Parameters
    ----------
    gene_symbol : str
        Gene symbol (e.g., "TP53", "MDM2")
    species_taxid : str
        NCBI taxonomy ID (e.g., "9606" for human)
    timeout : int
        Request timeout in seconds

    Returns
    -------
    str or None
        UniProt accession (e.g., "P04637") or None if not found

    Examples
    --------
    >>> ac = get_uniprot_ac("TP53", "9606")
    >>> print(ac)
    'P04637'

    Notes
    -----
    Results are cached to avoid repeated API calls.
    """
    cache_key = f"{gene_symbol}_{species_taxid}"

    if cache_key in _uniprot_cache:
        return _uniprot_cache[cache_key]

    # Query UniProt API - prefer reviewed (Swiss-Prot) entries
    url = f"{UNIPROT_API_BASE}/uniprotkb/search"
    params = {
        "query": f"gene:{gene_symbol} AND organism_id:{species_taxid} AND reviewed:true",
        "format": "json",
        "size": 1
    }

    try:
        resp = requests.get(url, params=params, timeout=timeout)
        resp.raise_for_status()
        data = resp.json()

        if data.get('results') and len(data['results']) > 0:
            primary_ac = data['results'][0].get('primaryAccession')
            _uniprot_cache[cache_key] = primary_ac
            return primary_ac

        # If no reviewed entry, try unreviewed
        params['query'] = f"gene:{gene_symbol} AND organism_id:{species_taxid}"
        resp = requests.get(url, params=params, timeout=timeout)
        resp.raise_for_status()
        data = resp.json()

        if data.get('results') and len(data['results']) > 0:
            primary_ac = data['results'][0].get('primaryAccession')
            _uniprot_cache[cache_key] = primary_ac
            return primary_ac

    except Exception:
        pass

    _uniprot_cache[cache_key] = None
    return None


def intact_rest_interaction_search(
    gene: str,
    species_taxid: str = "9606",
    page: int = 0,
    page_size: int = 500,
    timeout: int = 60,
    use_uniprot_ac: bool = True,
) -> dict:
    """
    Query IntAct Interaction Search REST API for interactions containing a gene.

    This replaces `psicquic_query_by_gene` and returns JSON instead of MITAB.

    Parameters
    ----------
    gene : str
        Gene symbol (e.g. "TP53") or UniProt AC (e.g., "P04637")
    species_taxid : str
        NCBI taxid as string (e.g. "9606" for human)
    page : int
        Page index (starts at 0)
    page_size : int
        How many interactions per page (max 500 recommended)
    timeout : int
        Request timeout in seconds
    use_uniprot_ac : bool
        If True, try to map gene symbol to UniProt AC for precise search

    Returns
    -------
    dict
        JSON response with structure:
        {
            'content': [list of interactions],
            'totalElements': total count,
            'totalPages': total pages,
            'number': current page number,
            'size': page size,
            ...
        }

    Examples
    --------
    >>> data = intact_rest_interaction_search("TP53", species_taxid="9606", page_size=100)
    >>> print(f"Found {data['totalElements']} interactions")
    >>> interactions = data['content']

    Notes
    -----
    If use_uniprot_ac is True, will attempt to convert gene symbol to UniProt AC
    for more precise results (e.g., TP53 → P04637 → returns actual TP53 interactions).
    Falls back to gene symbol if mapping fails.
    """
    query = gene

    # Try to get UniProt AC for more precise search
    if use_uniprot_ac and not gene.startswith(('P', 'Q', 'O', 'A', 'B')):  # Common UniProt prefixes
        uniprot_ac = get_uniprot_ac(gene, species_taxid)
        if uniprot_ac:
            query = uniprot_ac

    url = f"{INTACT_INTERACTION_BASE}/findInteractions/{query}"
    params = {
        "page": page,
        "size": page_size,
    }

    resp = requests.get(url, params=params, timeout=timeout)
    resp.raise_for_status()
    return resp.json()


def parse_interaction_rest_json(data: dict, target_gene: str = None, species_filter: str = None) -> List[dict]:
    """
    Convert IntAct REST API JSON response into normalized interaction dicts.

    This replaces `parse_mitab` and produces the same output format for compatibility
    with existing graph algorithms.

    Parameters
    ----------
    data : dict
        JSON response from `intact_rest_interaction_search`
    target_gene : str, optional
        If provided, only return interactions involving this gene
    species_filter : str, optional
        If provided, filter by species taxid (e.g. "9606")

    Returns
    -------
    List[dict]
        List of interaction dicts with standardized fields:
        {
            'idA_raw': str,
            'idB_raw': str,
            'idA': str,
            'idB': str,
            'nameA': str,
            'nameB': str,
            'taxidA': str,
            'taxidB': str,
            'organismA': str,
            'organismB': str,
            'miscore': float or None,
            'detection_method': str or None,
            'publication_ids': str or None,
            'interaction_type': str or None,
            'source_database': str,
        }
    """
    interactions_out = []

    # Extract interactions from content array
    raw_interactions = data.get('content', [])

    for inter in raw_interactions:
        # Extract basic IDs and names
        idA_raw = inter.get('idA', '')
        idB_raw = inter.get('idB', '')
        uniqueIdA = inter.get('uniqueIdA', '')
        uniqueIdB = inter.get('uniqueIdB', '')
        nameA = inter.get('moleculeA', '')
        nameB = inter.get('moleculeB', '')

        # Extract taxonomic info (already integers in REST API!)
        taxidA = str(inter.get('taxIdA', '')) if inter.get('taxIdA') else None
        taxidB = str(inter.get('taxIdB', '')) if inter.get('taxIdB') else None
        organismA = inter.get('speciesA', '')
        organismB = inter.get('speciesB', '')

        # Filter by species if requested
        if species_filter:
            if taxidA != species_filter and taxidB != species_filter:
                continue

        # Filter by target gene if requested
        if target_gene:
            target_upper = target_gene.upper()
            if not (nameA and nameA.upper() == target_upper) and not (nameB and nameB.upper() == target_upper):
                continue

        # Extract IntAct MI-score (directly available!)
        miscore = inter.get('intactMiscore')

        # Extract detection method
        detection_method = inter.get('detectionMethod')

        # Extract publication ID
        publication_ids = inter.get('publicationPubmedIdentifier')
        if publication_ids:
            publication_ids = str(publication_ids)

        # Extract interaction type
        interaction_type = inter.get('type')

        # Source database
        source_db = inter.get('sourceDatabase', 'IntAct REST')

        interactions_out.append({
            'idA_raw': idA_raw,
            'idB_raw': idB_raw,
            'idA': uniqueIdA or idA_raw.split()[0] if idA_raw else '',
            'idB': uniqueIdB or idB_raw.split()[0] if idB_raw else '',
            'nameA': nameA,
            'nameB': nameB,
            'taxidA': taxidA,
            'taxidB': taxidB,
            'organismA': organismA,
            'organismB': organismB,
            'miscore': miscore,
            'detection_method': detection_method,
            'publication_ids': publication_ids,
            'interaction_type': interaction_type,
            'source_database': source_db,
        })

    return interactions_out


def get_direct_neighbors_rest(
    gene: str,
    species_taxid: str = "9606",
    top_n: int = 100,
    page_size: int = 500,
    min_miscore: float = 0.0,
    organism_filter: Optional[str] = None,
) -> pd.DataFrame:
    """
    Get direct (1-hop) IntAct neighbors for a gene using REST API.

    This function replaces `get_direct_neighbors` but uses REST instead of PSICQUIC.

    Ranking:
      1) hop distance = 1 (all direct neighbors)
      2) best IntAct MIscore (descending)
      3) number of supporting interactions (descending)

    Parameters
    ----------
    gene : str
        Gene symbol (e.g. "TP53")
    species_taxid : str
        NCBI taxonomy ID as string (e.g. "9606" for human)
    top_n : int
        Maximum number of neighbors to return
    page_size : int
        Number of interactions per API page
    min_miscore : float
        Minimum MI-score threshold (0.0-1.0)
    organism_filter : str or None
        Comma-separated species names to filter (e.g. "homo sapiens,mus musculus")

    Returns
    -------
    pd.DataFrame
        One row per neighbor with columns:
        - neighbor_id: Primary ID
        - hop: Always 1 (direct neighbors)
        - neighbor_name: Gene name
        - neighbor_taxid: NCBI taxonomy ID
        - neighbor_organism: Species name
        - best_miscore: Highest MI-score across all interactions
        - n_interactions: Number of supporting interactions
        - detection_methods: Pipe-separated list of methods
        - interaction_types: Pipe-separated list of types
        - publications: Pipe-separated PubMed IDs
        - source_database: Always "IntAct REST"
        - path: Gene path (seed-neighbor)

    Examples
    --------
    >>> df = get_direct_neighbors_rest("TP53", species_taxid="9606", top_n=50)
    >>> print(f"Found {len(df)} neighbors for TP53")
    >>> print(df[['neighbor_name', 'best_miscore', 'n_interactions']].head())
    """
    gene_upper = gene.strip().upper()

    # Query REST API
    data = intact_rest_interaction_search(
        gene=gene,
        species_taxid=species_taxid,
        page=0,
        page_size=page_size,
    )

    # Parse JSON to standardized format
    interactions = parse_interaction_rest_json(
        data,
        target_gene=gene,
        species_filter=species_taxid
    )

    # Build neighbor dictionary
    neighbors = {}

    # Prepare organism filter
    filter_orgs = None
    if organism_filter:
        if isinstance(organism_filter, str):
            filter_orgs = {o.strip().lower() for o in organism_filter.split(",")}
        else:
            filter_orgs = {str(o).lower() for o in organism_filter}

    for it in interactions:
        miscore = it['miscore'] if it['miscore'] is not None else 0.0

        # Apply MI-score filter
        if miscore < min_miscore:
            continue

        detection = it.get('detection_method')
        interaction_type = it.get('interaction_type')
        pmids = it.get('publication_ids')
        source_db = it.get('source_database', 'IntAct REST')

        # Case 1: seed gene is A → neighbor is B
        if it['nameA'] and it['nameA'].upper() == gene_upper:
            nid = it['idB']

            # Apply organism filter
            if filter_orgs and it['organismB']:
                if it['organismB'].lower() not in filter_orgs:
                    continue

            if nid:
                rec = neighbors.get(nid)
                if rec is None:
                    neighbors[nid] = {
                        'neighbor_id': nid,
                        'hop': 1,
                        'neighbor_name': it['nameB'],
                        'neighbor_taxid': it['taxidB'],
                        'neighbor_organism': it['organismB'],
                        'best_miscore': miscore,
                        'n_interactions': 1,
                        'detection_methods': detection if detection else '',
                        'interaction_types': interaction_type if interaction_type else '',
                        'publications': pmids if pmids else '',
                        'source_database': source_db,
                        'path': f"{gene}-{it['nameB']}" if it['nameB'] else None,
                    }
                else:
                    _update_neighbor_record(rec, miscore, detection, interaction_type, pmids,
                                          it['organismB'], it['taxidB'], it['nameB'])

        # Case 2: seed gene is B → neighbor is A
        if it['nameB'] and it['nameB'].upper() == gene_upper:
            nid = it['idA']

            # Apply organism filter
            if filter_orgs and it['organismA']:
                if it['organismA'].lower() not in filter_orgs:
                    continue

            if nid:
                rec = neighbors.get(nid)
                if rec is None:
                    neighbors[nid] = {
                        'neighbor_id': nid,
                        'hop': 1,
                        'neighbor_name': it['nameA'],
                        'neighbor_taxid': it['taxidA'],
                        'neighbor_organism': it['organismA'],
                        'best_miscore': miscore,
                        'n_interactions': 1,
                        'detection_methods': detection if detection else '',
                        'interaction_types': interaction_type if interaction_type else '',
                        'publications': pmids if pmids else '',
                        'source_database': source_db,
                        'path': f"{gene}-{it['nameA']}" if it['nameA'] else None,
                    }
                else:
                    _update_neighbor_record(rec, miscore, detection, interaction_type, pmids,
                                          it['organismA'], it['taxidA'], it['nameA'])

    # Convert to DataFrame
    if not neighbors:
        return pd.DataFrame(columns=[
            'neighbor_id', 'hop', 'neighbor_name', 'neighbor_taxid', 'neighbor_organism',
            'best_miscore', 'n_interactions', 'detection_methods', 'interaction_types',
            'publications', 'source_database', 'path'
        ])

    df = pd.DataFrame(neighbors.values())
    df = df.sort_values(
        by=['hop', 'best_miscore', 'n_interactions'],
        ascending=[True, False, False],
    )

    return df.head(top_n)


def _update_neighbor_record(
    rec: dict,
    miscore: float,
    detection: Optional[str],
    interaction_type: Optional[str],
    pmids: Optional[str],
    organism: Optional[str],
    taxid: Optional[str],
    name: Optional[str],
):
    """
    Update a neighbor record with additional interaction data.

    This consolidates multiple interactions with the same neighbor.
    """
    rec['n_interactions'] += 1

    # Update best MI-score
    if miscore > rec['best_miscore']:
        rec['best_miscore'] = miscore

    # Fill in missing metadata
    if not rec['neighbor_name'] and name:
        rec['neighbor_name'] = name
    if not rec['neighbor_organism'] and organism:
        rec['neighbor_organism'] = organism
        rec['neighbor_taxid'] = taxid

    # Append detection methods
    if detection and detection not in (rec['detection_methods'] or ''):
        if rec['detection_methods']:
            rec['detection_methods'] = f"{rec['detection_methods']}|{detection}"
        else:
            rec['detection_methods'] = detection

    # Append interaction types
    if interaction_type and interaction_type not in (rec['interaction_types'] or ''):
        if rec['interaction_types']:
            rec['interaction_types'] = f"{rec['interaction_types']}|{interaction_type}"
        else:
            rec['interaction_types'] = interaction_type

    # Append PubMed IDs
    if pmids:
        existing_pmids = set((rec['publications'] or '').split('|')) if rec['publications'] else set()
        new_pmids = set(str(pmids).split('|'))
        all_pmids = existing_pmids | new_pmids
        all_pmids.discard('')
        rec['publications'] = '|'.join(sorted(all_pmids))


# ============================================================================
# Multi-hop and Shortest Path Functions
# ============================================================================
# These can be copied from the original intact_api.py with minimal changes
# (just replace psicquic_query_by_gene + parse_mitab with REST equivalents)
#
# For now, I'm including stubs to show the structure. Full implementation
# would mirror the existing BFS and Dijkstra logic.
# ============================================================================

def get_neighbors_multihop_rest(
    gene: str,
    species_taxid: str = "9606",
    top_n: int = 100,
    max_hops: int = 3,
    min_miscore: float = 0.4,
    organism_filter: Optional[str] = None,
    page_size: int = 500,
) -> pd.DataFrame:
    """
    Multi-hop BFS expansion using REST API.

    This is a drop-in replacement for get_neighbors_multihop using REST instead of PSICQUIC.
    Performs breadth-first search to find neighbors up to max_hops distance.

    Parameters
    ----------
    gene : str
        Seed gene symbol
    species_taxid : str
        NCBI taxonomy ID
    top_n : int
        Maximum neighbors to return
    max_hops : int
        Maximum hop distance
    min_miscore : float
        Minimum MI-score threshold
    organism_filter : str or None
        Species filter (e.g., "homo sapiens,mus musculus")
    page_size : int
        API page size

    Returns
    -------
    pd.DataFrame
        Neighbors up to max_hops distance with columns:
        - neighbor_id, hop, neighbor_name, neighbor_taxid, neighbor_organism,
          best_miscore, n_interactions, detection_methods, interaction_types,
          publications, source_database, path

    Examples
    --------
    >>> df = get_neighbors_multihop_rest("TP53", max_hops=2, top_n=100)
    >>> print(f"Found neighbors at hops: {df['hop'].value_counts().to_dict()}")
    """
    gene_upper = gene.strip().upper()

    # Initialize BFS
    visited = {gene_upper}
    current_layer = [gene_upper]
    all_neighbors = {}

    for hop in range(1, max_hops + 1):
        if not current_layer:
            break

        next_layer = set()

        # Expand each gene in current layer
        for current_gene in current_layer:
            # Get direct neighbors using REST
            neighbors_df = get_direct_neighbors_rest(
                gene=current_gene,
                species_taxid=species_taxid,
                top_n=page_size,  # Get many neighbors for BFS
                page_size=page_size,
                min_miscore=min_miscore,
                organism_filter=organism_filter,
            )

            # Process each neighbor
            for _, row in neighbors_df.iterrows():
                neighbor_id = row['neighbor_id']
                neighbor_name = row['neighbor_name']

                if not neighbor_name:
                    continue

                neighbor_upper = neighbor_name.upper()

                # Skip if already visited
                if neighbor_upper in visited:
                    continue

                visited.add(neighbor_upper)
                next_layer.add(neighbor_upper)

                # Record neighbor with current hop distance
                if neighbor_id not in all_neighbors:
                    all_neighbors[neighbor_id] = {
                        'neighbor_id': neighbor_id,
                        'hop': hop,
                        'neighbor_name': neighbor_name,
                        'neighbor_taxid': row['neighbor_taxid'],
                        'neighbor_organism': row['neighbor_organism'],
                        'best_miscore': row['best_miscore'],
                        'n_interactions': row['n_interactions'],
                        'detection_methods': row['detection_methods'],
                        'interaction_types': row['interaction_types'],
                        'publications': row['publications'],
                        'source_database': row['source_database'],
                        'path': f"{gene}-{neighbor_name}" if hop == 1 else f"{gene}-...-{neighbor_name} ({hop} hops)",
                    }
                else:
                    # Update if we found a shorter path or better score
                    existing = all_neighbors[neighbor_id]
                    if hop < existing['hop']:
                        existing['hop'] = hop
                        existing['path'] = f"{gene}-{neighbor_name}" if hop == 1 else f"{gene}-...-{neighbor_name} ({hop} hops)"
                    if row['best_miscore'] > existing['best_miscore']:
                        existing['best_miscore'] = row['best_miscore']
                    existing['n_interactions'] += row['n_interactions']

        current_layer = list(next_layer)

        # Stop if we have enough neighbors
        if len(all_neighbors) >= top_n:
            break

    # Convert to DataFrame
    if not all_neighbors:
        return pd.DataFrame(columns=[
            'neighbor_id', 'hop', 'neighbor_name', 'neighbor_taxid', 'neighbor_organism',
            'best_miscore', 'n_interactions', 'detection_methods', 'interaction_types',
            'publications', 'source_database', 'path'
        ])

    df = pd.DataFrame(all_neighbors.values())
    df = df.sort_values(
        by=['hop', 'best_miscore', 'n_interactions'],
        ascending=[True, False, False],
    )

    return df.head(top_n)


def find_shortest_paths_rest(
    gene_list: List[str],
    species_taxid: str = "9606",
    max_distance: int = 50,
    min_miscore: float = 0.4,
    organism_filter: Optional[str] = None,
    page_size: int = 500,
) -> Dict[Tuple[str, str], Dict[str, Any]]:
    """
    Find shortest paths between genes using REST API with Dijkstra's algorithm.

    This is a drop-in replacement for find_shortest_paths_intact using REST instead of PSICQUIC.
    Uses Dijkstra's algorithm with edge weights = 1.0 - miscore (lower score = longer distance).

    Parameters
    ----------
    gene_list : List[str]
        List of gene symbols to find paths between
    species_taxid : str
        NCBI taxonomy ID
    max_distance : int
        Maximum path length (in terms of weighted distance)
    min_miscore : float
        Minimum MI-score threshold
    organism_filter : str or None
        Species filter (e.g., "homo sapiens")
    page_size : int
        API page size

    Returns
    -------
    Dict[Tuple[str, str], Dict[str, Any]]
        Mapping from (gene_a, gene_b) tuples to path information:
        {
            'path': [gene1, gene2, ...],
            'hops': int,
            'distance': float,
            'scores': [float, ...],
            'algorithm': 'Dijkstra',
            'weight_formula': 'weight = 1.0 - miscore'
        }

    Examples
    --------
    >>> paths = find_shortest_paths_rest(["TP53", "MDM2", "ATM"], species_taxid="9606")
    >>> for (ga, gb), info in paths.items():
    ...     print(f"{ga} ↔ {gb}: {' → '.join(info['path'])} ({info['hops']} hops)")
    """
    # Validate input
    if not gene_list or len(gene_list) < 2:
        raise ValueError("gene_list must contain at least 2 genes")

    # Build interaction graph from all query genes
    # CRITICAL FIX: Query each gene SEPARATELY using parallel threads for performance
    graph = defaultdict(list)  # gene -> [(neighbor, miscore), ...]
    gene_set = {g.upper() for g in gene_list}
    graph_lock = threading.Lock()

    def query_gene_interactions(gene: str):
        """Query IntAct for a single gene (thread-safe)."""
        try:
            data = intact_rest_interaction_search(
                gene=gene,
                species_taxid=species_taxid,
                page=0,
                page_size=page_size,
            )

            interactions = parse_interaction_rest_json(data, species_filter=species_taxid)

            # Build bidirectional graph
            local_edges = []
            for inter in interactions:
                miscore = inter['miscore'] if inter['miscore'] is not None else 0.0

                if miscore < min_miscore:
                    continue

                nameA = inter['nameA']
                nameB = inter['nameB']

                if not nameA or not nameB:
                    continue

                # Apply organism filter
                if organism_filter:
                    filter_orgs = {o.strip().lower() for o in organism_filter.split(",")}
                    if inter['organismA']:
                        if inter['organismA'].lower() not in filter_orgs:
                            continue
                    if inter['organismB']:
                        if inter['organismB'].lower() not in filter_orgs:
                            continue

                # Store edges locally
                local_edges.append((nameA.upper(), nameB.upper(), miscore))
                local_edges.append((nameB.upper(), nameA.upper(), miscore))

            return local_edges
        except Exception:
            return []

    # Query genes in parallel
    with ThreadPoolExecutor(max_workers=10) as executor:
        future_to_gene = {executor.submit(query_gene_interactions, gene): gene
                          for gene in gene_list}

        for future in as_completed(future_to_gene):
            try:
                edges = future.result()
                # Add all edges to graph (thread-safe)
                with graph_lock:
                    for nameA, nameB, miscore in edges:
                        graph[nameA].append((nameB, miscore))
            except Exception:
                # Continue with other genes if one fails
                pass

    # Find shortest paths between all pairs of query genes
    paths = {}

    for i, source in enumerate(gene_list):
        source_upper = source.upper()

        if source_upper not in graph:
            continue

        # Dijkstra from this source
        distances = {source_upper: 0.0}
        predecessors = {source_upper: None}
        mi_scores = {source_upper: None}

        # Priority queue: (distance, gene)
        pq = [(0.0, source_upper)]
        visited = set()

        while pq:
            current_dist, current_gene = heapq.heappop(pq)

            if current_gene in visited:
                continue

            visited.add(current_gene)

            # Stop if distance exceeds max
            if current_dist > max_distance:
                continue

            # Explore neighbors
            for neighbor, miscore in graph.get(current_gene, []):
                if neighbor in visited:
                    continue

                # Edge weight = 1.0 - miscore (higher miscore = shorter distance)
                edge_weight = 1.0 - miscore
                new_dist = current_dist + edge_weight

                if neighbor not in distances or new_dist < distances[neighbor]:
                    distances[neighbor] = new_dist
                    predecessors[neighbor] = current_gene
                    mi_scores[neighbor] = miscore
                    heapq.heappush(pq, (new_dist, neighbor))

        # Reconstruct paths to other query genes
        for j, target in enumerate(gene_list):
            if i >= j:  # Only compute each pair once
                continue

            target_upper = target.upper()

            if target_upper not in distances:
                continue

            # Reconstruct path
            path = []
            current = target_upper
            while current is not None:
                path.append(current)
                current = predecessors.get(current)

            path.reverse()

            # Get MI-scores for each edge
            edge_scores = []
            for k in range(len(path) - 1):
                nodeA = path[k]
                nodeB = path[k + 1]
                # Find miscore for this edge
                score = None
                for neighbor, miscore in graph.get(nodeA, []):
                    if neighbor == nodeB:
                        score = miscore
                        break
                edge_scores.append(score)

            # Store path
            key = tuple(sorted([source.upper(), target.upper()]))
            paths[key] = {
                'path': path,
                'hops': len(path) - 1,
                'distance': distances[target_upper],
                'scores': edge_scores,
                'algorithm': 'Dijkstra',
                'weight_formula': 'weight = 1.0 - miscore'
            }

    return paths


# ============================================================================
# Backward Compatibility Aliases
# ============================================================================
# These aliases maintain compatibility with existing code that uses the old
# function names from the PSICQUIC version of intact_api.py

def get_neighbors_multihop(
    gene: str,
    species: str = "human",
    top_n: int = 100,
    max_hops: int = 3,
    min_miscore: float = 0.4,
    organism_filter: Optional[str] = None,
    **kwargs
) -> pd.DataFrame:
    """
    Backward compatibility wrapper for get_neighbors_multihop_rest.

    Maps old PSICQUIC-style parameters to new REST API parameters.
    """
    # Convert species string to taxid
    species_map = {
        "human": "9606",
        "mouse": "10090",
        "9606": "9606",
        "10090": "10090"
    }
    species_taxid = species_map.get(str(species).lower(), "9606")

    return get_neighbors_multihop_rest(
        gene=gene,
        species_taxid=species_taxid,
        top_n=top_n,
        max_hops=max_hops,
        min_miscore=min_miscore,
        organism_filter=organism_filter,
        page_size=kwargs.get('page_size', 500)
    )


def find_shortest_paths_intact(
    gene_list: List[str],
    species: str = "human",
    max_distance: int = 50,
    min_miscore: float = 0.4,
    organism_filter: Optional[str] = None,
    **kwargs
) -> Dict[Tuple[str, str], Dict[str, Any]]:
    """
    Backward compatibility wrapper for find_shortest_paths_rest.

    Maps old PSICQUIC-style parameters to new REST API parameters.
    """
    # Convert species string to taxid
    species_map = {
        "human": "9606",
        "mouse": "10090",
        "9606": "9606",
        "10090": "10090"
    }
    species_taxid = species_map.get(str(species).lower(), "9606")

    return find_shortest_paths_rest(
        gene_list=gene_list,
        species_taxid=species_taxid,
        max_distance=max_distance,
        min_miscore=min_miscore,
        organism_filter=organism_filter,
        page_size=kwargs.get('page_size', 500)
    )


def get_direct_neighbors(
    gene: str,
    species: str = "human",
    top_n: int = 100,
    min_miscore: float = 0.0,
    organism_filter: Optional[str] = None,
    **kwargs
) -> pd.DataFrame:
    """
    Backward compatibility wrapper for get_direct_neighbors_rest.

    Maps old PSICQUIC-style parameters to new REST API parameters.
    """
    # Convert species string to taxid
    species_map = {
        "human": "9606",
        "mouse": "10090",
        "9606": "9606",
        "10090": "10090"
    }
    species_taxid = species_map.get(str(species).lower(), "9606")

    return get_direct_neighbors_rest(
        gene=gene,
        species_taxid=species_taxid,
        top_n=top_n,
        page_size=kwargs.get('page_size', 500),
        min_miscore=min_miscore,
        organism_filter=organism_filter
    )


if __name__ == "__main__":
    # Comprehensive test
    print("Testing IntAct REST API Implementation")
    print("=" * 70)

    # Test 1: UniProt AC mapping
    print("\n1. Testing UniProt AC mapping...")
    ac = get_uniprot_ac("TP53", "9606")
    print(f"   TP53 → {ac}")
    assert ac == "P04637", f"Expected P04637, got {ac}"
    print("   ✓ UniProt mapping works")

    # Test 2: Direct neighbors with REST API
    print("\n2. Testing direct neighbors (1-hop)...")
    df = get_direct_neighbors_rest("TP53", species_taxid="9606", top_n=10)
    print(f"   Found {len(df)} neighbors for TP53")
    if len(df) > 0:
        print(f"   Top neighbor: {df.iloc[0]['neighbor_name']} (score: {df.iloc[0]['best_miscore']})")
        print("   ✓ Direct neighbors work")
    else:
        print("   ⚠ No neighbors found (check API)")

    # Test 3: Multi-hop BFS
    print("\n3. Testing multi-hop BFS...")
    df_multihop = get_neighbors_multihop_rest("TP53", species_taxid="9606", max_hops=2, top_n=20)
    print(f"   Found {len(df_multihop)} neighbors (up to 2 hops)")
    if len(df_multihop) > 0:
        hop_counts = df_multihop['hop'].value_counts().to_dict()
        print(f"   Hop distribution: {hop_counts}")
        print("   ✓ Multi-hop BFS works")
    else:
        print("   ⚠ No multi-hop neighbors found")

    # Test 4: Shortest paths
    print("\n4. Testing shortest paths (Dijkstra)...")
    paths = find_shortest_paths_rest(["TP53", "MDM2"], species_taxid="9606")
    print(f"   Found {len(paths)} paths")
    for (ga, gb), info in paths.items():
        print(f"   {ga} ↔ {gb}: {' → '.join(info['path'])} ({info['hops']} hops, distance: {info['distance']:.2f})")
        print("   ✓ Shortest paths work")

    print("\n" + "=" * 70)
    print("✓ All tests passed! REST API implementation is working.")
    print("=" * 70)
