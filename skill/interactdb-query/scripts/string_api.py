#!/usr/bin/env python3
"""
Query STRING database for protein-protein interactions and neighbors.

Features
--------
- Data source: STRING database (https://string-db.org/) providing known and
  predicted protein-protein interactions from multiple evidence sources.
- Input: gene/protein identifier (STRING protein ID, gene name, or UniProt AC)
  and organism taxon ID (default 9606 = Homo sapiens).
- Network retrieval: Uses STRING API /network endpoint to fetch interaction
  partners with confidence scores.
- Multiple evidence channels:
    * experimental: from experimental data repositories
    * database: from curated databases
    * coexpression: from gene expression patterns
    * neighborhood: from genomic context
    * fusion: from gene fusion events
    * cooccurrence: from phylogenetic co-occurrence
    * textmining: from text mining of literature
- Confidence scores: combined score (0-999) and individual channel scores.
- Ranking: neighbors sorted by combined confidence score (descending).
- Output: DataFrame with entity annotations (protein ID, preferred name,
  organism) and edge annotations (all evidence channel scores, combined score).

API Documentation
-----------------
https://string-db.org/help/api/

Requirements
------------
    pip install requests pandas

Usage
-----
    python string_api.py TP53 --top-n 50 --min-score 400
"""

import requests
import pandas as pd
from typing import Optional, Dict, Any, List, Tuple, Set
from dataclasses import dataclass, asdict
from collections import defaultdict
import heapq


STRING_API_BASE = "https://string-db.org/api"


@dataclass
class StringNeighbor:
    """Represents a STRING interaction neighbor with annotations."""
    # Entity annotations
    protein_id: str              # STRING protein identifier
    preferred_name: str          # Gene/protein name
    organism_id: int             # NCBI taxonomy ID
    organism_name: str           # Organism name
    annotation: Optional[str]    # Functional annotation

    # Edge annotations (scores 0-999)
    combined_score: int          # Combined confidence score
    experimental_score: int      # Experimental evidence
    database_score: int          # Curated database evidence
    textmining_score: int        # Text mining evidence
    coexpression_score: int      # Gene expression evidence
    neighborhood_score: int      # Genomic neighborhood
    fusion_score: int            # Gene fusion events
    cooccurrence_score: int      # Phylogenetic co-occurrence

    # Metadata
    hop: int = 1                 # Hop distance from seed
    source_database: str = "STRING"
    path: Optional[str] = None   # Path from seed, e.g. "TP53-MDM2"


class StringClient:
    """Client for querying the STRING database API."""

    def __init__(self, base_url: str = STRING_API_BASE, timeout: int = 60):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.organism_cache: Dict[int, str] = {}

    def _request(self, endpoint: str, params: Dict[str, Any]) -> Any:
        """Make a GET request to STRING API."""
        url = f"{self.base_url}/{endpoint}"
        resp = requests.get(url, params=params, timeout=self.timeout)
        resp.raise_for_status()
        return resp.json()

    def resolve_identifier(
        self,
        identifier: str,
        species: int = 9606
    ) -> Optional[str]:
        """
        Resolve a gene name/UniProt AC to STRING protein identifier.

        Returns STRING ID like '9606.ENSP00000269305' or None if not found.
        """
        params = {
            "identifiers": identifier,
            "species": species,
            "limit": 1,
            "format": "json"
        }

        try:
            data = self._request("json/get_string_ids", params)
            if data and len(data) > 0:
                return data[0].get("stringId")
        except Exception:
            pass

        return None

    def get_organism_name(self, taxon_id: int) -> str:
        """Get organism name from taxonomy ID (with caching)."""
        if taxon_id in self.organism_cache:
            return self.organism_cache[taxon_id]

        # Common organisms
        common_names = {
            9606: "Homo sapiens",
            10090: "Mus musculus",
            10116: "Rattus norvegicus",
            7227: "Drosophila melanogaster",
            6239: "Caenorhabditis elegans",
            4932: "Saccharomyces cerevisiae",
            511145: "Escherichia coli K-12",
        }

        name = common_names.get(taxon_id, f"taxid:{taxon_id}")
        self.organism_cache[taxon_id] = name
        return name

    def get_neighbors(
        self,
        identifier: str,
        species: int = 9606,
        top_n: int = 100,
        min_combined_score: int = 0,
        min_experimental_score: int = 0,
        min_database_score: int = 0,
        min_textmining_score: int = 0,
        min_coexpression_score: int = 0,
        min_neighborhood_score: int = 0,
        min_fusion_score: int = 0,
        min_cooccurrence_score: int = 0,
        required_score: Optional[int] = None,
        network_type: str = "functional",
        add_nodes: int = 0,
    ) -> List[StringNeighbor]:
        """
        Get protein interaction neighbors from STRING.

        Parameters
        ----------
        identifier : str
            Gene name, UniProt AC, or STRING protein ID (e.g., 'TP53')
        species : int
            NCBI taxonomy ID (default 9606 = human)
        top_n : int
            Maximum number of neighbors to return
        min_combined_score : int
            Minimum combined confidence score (0-999)
        min_experimental_score : int
            Minimum experimental evidence score (0-999)
        min_database_score : int
            Minimum curated database score (0-999)
        min_textmining_score : int
            Minimum text mining score (0-999)
        min_coexpression_score : int
            Minimum coexpression score (0-999)
        min_neighborhood_score : int
            Minimum genomic neighborhood score (0-999)
        min_fusion_score : int
            Minimum gene fusion score (0-999)
        min_cooccurrence_score : int
            Minimum phylogenetic co-occurrence score (0-999)
        required_score : int, optional
            If provided, overrides min_combined_score for API filtering
        network_type : str
            'functional' (default) or 'physical' interactions
        add_nodes : int
            Number of additional nodes to add to network (0 = only direct)

        Returns
        -------
        List[StringNeighbor]
            Sorted by combined score (descending)
        """
        # Resolve identifier to STRING ID for matching
        string_id = self.resolve_identifier(identifier, species)
        if not string_id:
            string_id = f"{species}.{identifier}"

        # Query network endpoint (use original identifier - STRING API handles resolution)
        params = {
            "identifiers": identifier,  # Use gene name - API resolves automatically
            "species": species,
            "network_type": network_type,
        }

        # Only add these parameters if non-zero/non-default
        if add_nodes > 0:
            params["add_nodes"] = add_nodes
        if required_score is not None:
            params["required_score"] = required_score

        try:
            data = self._request("json/network", params)
        except Exception as e:
            raise RuntimeError(f"STRING API error: {e}")

        if not data:
            return []

        # Parse interactions
        neighbors: Dict[str, StringNeighbor] = {}
        seed_id = string_id
        organism_name = self.get_organism_name(species)

        # Process interactions in two passes:
        # Pass 1: Direct connections to seed (these get accurate paths)
        # Pass 2: Connections between neighbors (only if add_nodes > 0)

        for interaction in data:
            str_a = interaction.get("stringId_A", "")
            str_b = interaction.get("stringId_B", "")
            pref_a = interaction.get("preferredName_A", "")
            pref_b = interaction.get("preferredName_B", "")
            annot_a = interaction.get("annotation_A")
            annot_b = interaction.get("annotation_B")

            # Edge scores (STRING returns 0.0-1.0, multiply by 1000 to get 0-999)
            combined = int(float(interaction.get("score", 0)) * 1000)
            experimental = int(float(interaction.get("escore", 0)) * 1000)
            database = int(float(interaction.get("dscore", 0)) * 1000)
            textmining = int(float(interaction.get("tscore", 0)) * 1000)

            # Additional channel scores
            coexpression = int(float(interaction.get("ascore", 0)) * 1000)  # 'a' = automated/expression
            neighborhood = int(float(interaction.get("nscore", 0)) * 1000)
            fusion = int(float(interaction.get("fscore", 0)) * 1000)
            cooccurrence = int(float(interaction.get("pscore", 0)) * 1000)  # 'p' = phylogenetic

            # Filter by combined score
            if combined < min_combined_score:
                continue

            # Filter by individual evidence channel scores
            if experimental < min_experimental_score:
                continue
            if database < min_database_score:
                continue
            if textmining < min_textmining_score:
                continue
            if coexpression < min_coexpression_score:
                continue
            if neighborhood < min_neighborhood_score:
                continue
            if fusion < min_fusion_score:
                continue
            if cooccurrence < min_cooccurrence_score:
                continue

            # Determine neighbor (exclude seed protein)
            neighbor_id = None
            neighbor_name = None
            neighbor_annot = None

            if str_a == seed_id or pref_a.upper() == identifier.upper():
                neighbor_id = str_b
                neighbor_name = pref_b
                neighbor_annot = annot_b
            elif str_b == seed_id or pref_b.upper() == identifier.upper():
                neighbor_id = str_a
                neighbor_name = pref_a
                neighbor_annot = annot_a
            else:
                # Neither matches seed - skip unless add_nodes > 0
                # (These are indirect neighbors from network expansion)
                if add_nodes == 0:
                    continue

                # Add indirect neighbors with unknown path
                for nid, nname, nannot in [(str_a, pref_a, annot_a), (str_b, pref_b, annot_b)]:
                    if nid not in neighbors:
                        path_str = f"{identifier}-...-{nname}" if nname else None

                        neighbors[nid] = StringNeighbor(
                            protein_id=nid,
                            preferred_name=nname,
                            organism_id=species,
                            organism_name=organism_name,
                            annotation=nannot,
                            combined_score=combined,
                            experimental_score=experimental,
                            database_score=database,
                            textmining_score=textmining,
                            coexpression_score=coexpression,
                            neighborhood_score=neighborhood,
                            fusion_score=fusion,
                            cooccurrence_score=cooccurrence,
                            path=path_str,
                        )
                    else:
                        # Keep best score if duplicate
                        if combined > neighbors[nid].combined_score:
                            neighbors[nid].combined_score = combined
                            neighbors[nid].experimental_score = experimental
                            neighbors[nid].database_score = database
                            neighbors[nid].textmining_score = textmining
                            neighbors[nid].coexpression_score = coexpression
                            neighbors[nid].neighborhood_score = neighborhood
                            neighbors[nid].fusion_score = fusion
                            neighbors[nid].cooccurrence_score = cooccurrence
                continue

            if not neighbor_id:
                continue

            # Add or update neighbor
            if neighbor_id not in neighbors:
                # Format path for 1-hop: SEED-NEIGHBOR
                path_str = f"{identifier}-{neighbor_name}" if neighbor_name else None

                neighbors[neighbor_id] = StringNeighbor(
                    protein_id=neighbor_id,
                    preferred_name=neighbor_name,
                    organism_id=species,
                    organism_name=organism_name,
                    annotation=neighbor_annot,
                    combined_score=combined,
                    experimental_score=experimental,
                    database_score=database,
                    textmining_score=textmining,
                    coexpression_score=coexpression,
                    neighborhood_score=neighborhood,
                    fusion_score=fusion,
                    cooccurrence_score=cooccurrence,
                    path=path_str,
                )
            else:
                # Keep best score if duplicate
                if combined > neighbors[neighbor_id].combined_score:
                    neighbors[neighbor_id].combined_score = combined
                    neighbors[neighbor_id].experimental_score = experimental
                    neighbors[neighbor_id].database_score = database
                    neighbors[neighbor_id].textmining_score = textmining
                    neighbors[neighbor_id].coexpression_score = coexpression
                    neighbors[neighbor_id].neighborhood_score = neighborhood
                    neighbors[neighbor_id].fusion_score = fusion
                    neighbors[neighbor_id].cooccurrence_score = cooccurrence

        # Sort by combined score (descending)
        sorted_neighbors = sorted(
            neighbors.values(),
            key=lambda n: n.combined_score,
            reverse=True
        )

        return sorted_neighbors[:top_n]

    def get_neighbors_multihop(
        self,
        identifier: str,
        species: int = 9606,
        top_n: int = 100,
        max_hops: int = 3,
        min_combined_score: int = 0,
        min_experimental_score: int = 0,
        min_database_score: int = 0,
        min_textmining_score: int = 0,
        min_coexpression_score: int = 0,
        min_neighborhood_score: int = 0,
        min_fusion_score: int = 0,
        min_cooccurrence_score: int = 0,
        network_type: str = "functional",
    ) -> List[StringNeighbor]:
        """
        Get protein interaction neighbors with automatic multi-hop BFS expansion.

        This function iteratively expands from 1-hop to 2-hop to 3-hop neighbors
        until top_n neighbors are found or max_hops is reached.

        **Algorithm**: Breadth-First Search (BFS) with iterative expansion
        **Principle**: If 1-hop neighbors < top_n, expand to 2-hop, then 3-hop, etc.

        Parameters
        ----------
        identifier : str
            Gene name, UniProt AC, or STRING protein ID
        species : int
            NCBI taxonomy ID (default 9606 = human)
        top_n : int
            Target number of neighbors to return (default 100)
        max_hops : int
            Maximum number of hops to expand (default 3)
        min_combined_score : int
            Minimum combined confidence score (0-999)
        min_experimental_score : int
            Minimum experimental evidence score (0-999)
        min_database_score : int
            Minimum curated database score (0-999)
        min_textmining_score : int
            Minimum text mining score (0-999)
        min_coexpression_score : int
            Minimum coexpression score (0-999)
        min_neighborhood_score : int
            Minimum genomic neighborhood score (0-999)
        min_fusion_score : int
            Minimum gene fusion score (0-999)
        min_cooccurrence_score : int
            Minimum phylogenetic co-occurrence score (0-999)
        network_type : str
            'functional' (default) or 'physical' interactions

        Returns
        -------
        List[StringNeighbor]
            Neighbors sorted by (hop, combined_score desc)
            Each neighbor includes hop distance and path from seed

        Examples
        --------
        >>> client = StringClient()
        >>> neighbors = client.get_neighbors_multihop("GREM1", top_n=100, max_hops=3)
        >>> print(f"Found {len(neighbors)} neighbors")
        >>> print(f"1-hop: {sum(1 for n in neighbors if n.hop == 1)}")
        >>> print(f"2-hop: {sum(1 for n in neighbors if n.hop == 2)}")
        """
        identifier = identifier.strip()

        # Track all neighbors across hops
        all_neighbors: Dict[str, StringNeighbor] = {}
        visited: Set[str] = set()

        # BFS frontier (genes to expand from)
        frontier: Set[str] = {identifier.upper()}
        visited.add(identifier.upper())

        # Track paths: gene_name -> [seed, intermediate, ..., gene]
        paths: Dict[str, List[str]] = {identifier.upper(): [identifier]}

        for hop in range(1, max_hops + 1):
            if not frontier:
                break

            # Query neighbors of all genes in current frontier
            new_frontier: Set[str] = set()

            for current_gene in frontier:
                # Get 1-hop neighbors of current gene
                try:
                    direct_neighbors = self.get_neighbors(
                        identifier=current_gene,
                        species=species,
                        top_n=1000,  # Get all direct neighbors
                        min_combined_score=min_combined_score,
                        min_experimental_score=min_experimental_score,
                        min_database_score=min_database_score,
                        min_textmining_score=min_textmining_score,
                        min_coexpression_score=min_coexpression_score,
                        min_neighborhood_score=min_neighborhood_score,
                        min_fusion_score=min_fusion_score,
                        min_cooccurrence_score=min_cooccurrence_score,
                        network_type=network_type,
                        add_nodes=0,  # Only direct neighbors
                    )
                except Exception:
                    continue

                for neighbor in direct_neighbors:
                    neighbor_name = neighbor.preferred_name.upper()

                    # Skip if already visited
                    if neighbor_name in visited:
                        continue

                    # Build path for this neighbor
                    if current_gene in paths:
                        neighbor_path = paths[current_gene] + [neighbor.preferred_name]
                    else:
                        neighbor_path = [identifier, neighbor.preferred_name]

                    paths[neighbor_name] = neighbor_path

                    # Update neighbor record with hop distance and path
                    neighbor.hop = hop
                    neighbor.path = "-".join(neighbor_path)

                    # Add to all_neighbors
                    if neighbor_name not in all_neighbors:
                        all_neighbors[neighbor_name] = neighbor
                    else:
                        # Keep neighbor with better score if duplicate
                        if neighbor.combined_score > all_neighbors[neighbor_name].combined_score:
                            all_neighbors[neighbor_name] = neighbor

                    # Add to new frontier for next hop expansion
                    new_frontier.add(neighbor_name)
                    visited.add(neighbor_name)

            frontier = new_frontier

            # Stop if we have enough neighbors
            if len(all_neighbors) >= top_n:
                break

        # Sort by hop (ascending), then by combined score (descending)
        sorted_neighbors = sorted(
            all_neighbors.values(),
            key=lambda n: (n.hop, -n.combined_score)
        )

        return sorted_neighbors[:top_n]

    def find_shortest_paths(
        self,
        gene_list: List[str],
        species: int = 9606,
        max_distance: int = 50,
        max_network_expansion: int = 5,
        min_combined_score: int = 400,
        min_experimental_score: int = 0,
        min_database_score: int = 0,
        min_textmining_score: int = 0,
        min_coexpression_score: int = 0,
        network_type: str = "functional",
    ) -> Dict[Tuple[str, str], Dict[str, Any]]:
        """
        Find shortest paths between all pairs of query genes using interaction network.

        **Algorithm**: BFS network expansion + Dijkstra's shortest path
        **Edge weight formula**: weight = 1000 - combined_score
        **Principle**: Higher score = lower weight = shorter distance (closer nodes)
        **NEW**: Iteratively expands network via BFS to capture intermediate edges

        Parameters
        ----------
        gene_list : List[str]
            List of gene names/identifiers to find paths between
        species : int
            NCBI taxonomy ID (default 9606 = human)
        max_distance : int
            Maximum path length in final Dijkstra search (default 50)
        max_network_expansion : int
            Maximum BFS hops to expand network (default 5)
            Controls depth of network exploration, not path length
        min_combined_score : int
            Minimum combined confidence score for edges (0-999, default 400)
        min_experimental_score : int
            Minimum experimental evidence score (0-999, default 0)
        min_database_score : int
            Minimum curated database score (0-999, default 0)
        min_textmining_score : int
            Minimum text mining score (0-999, default 0)
        min_coexpression_score : int
            Minimum coexpression score (0-999, default 0)
        network_type : str
            'functional' (default) or 'physical' interactions

        Returns
        -------
        Dict[Tuple[str, str], Dict]
            Dictionary mapping (geneA, geneB) tuples to path info:
            {
                'path': ['A', 'D', 'B'],
                'distance': 2.5,
                'hops': 2,
                'scores': [999, 850],  # edge scores along path
                'algorithm': 'Dijkstra+BFS',
                'weight_formula': '1000 - combined_score'
            }
        """
        if len(gene_list) < 2:
            return {}

        # Resolve all query genes to STRING IDs
        gene_to_string_id = {}
        string_id_to_gene = {}

        for gene in gene_list:
            string_id = self.resolve_identifier(gene, species)
            if string_id:
                gene_to_string_id[gene.upper()] = string_id
                string_id_to_gene[string_id] = gene.upper()
            else:
                string_id_to_gene[f"{species}.{gene}"] = gene.upper()
                gene_to_string_id[gene.upper()] = f"{species}.{gene}"

        # Phase 1: Expand network iteratively via BFS to build complete graph
        frontier = set(g.upper() for g in gene_list)
        visited = set(frontier)
        all_interactions = []

        for expansion_hop in range(1, max_network_expansion + 1):
            if not frontier:
                break

            # Query interactions for current frontier
            frontier_str = "|".join(frontier)
            params = {
                "identifiers": frontier_str,
                "species": species,
                "network_type": network_type,
                "add_nodes": 0,  # Only direct neighbors
            }

            try:
                data = self._request("json/network", params)
            except Exception:
                break

            if not data:
                break

            all_interactions.extend(data)

            # Extract new genes from interactions
            new_frontier = set()
            for interaction in data:
                pref_a = interaction.get("preferredName_A", "").upper()
                pref_b = interaction.get("preferredName_B", "").upper()

                for gene in [pref_a, pref_b]:
                    if gene and gene not in visited:
                        new_frontier.add(gene)
                        visited.add(gene)

            frontier = new_frontier

            # Stop early if we have enough nodes
            if len(visited) > 1000:  # Prevent explosion
                break

        # Phase 2: Build complete graph from all interactions
        graph: Dict[str, Dict[str, float]] = defaultdict(dict)
        edge_scores: Dict[Tuple[str, str], int] = {}

        for interaction in all_interactions:
            str_a = interaction.get("stringId_A", "")
            str_b = interaction.get("stringId_B", "")
            pref_a = interaction.get("preferredName_A", "")
            pref_b = interaction.get("preferredName_B", "")

            # Parse all scores
            combined = int(float(interaction.get("score", 0)) * 1000)
            experimental = int(float(interaction.get("escore", 0)) * 1000)
            database = int(float(interaction.get("dscore", 0)) * 1000)
            textmining = int(float(interaction.get("tscore", 0)) * 1000)
            coexpression = int(float(interaction.get("ascore", 0)) * 1000)

            # Apply all score filters
            if combined < min_combined_score:
                continue
            if experimental < min_experimental_score:
                continue
            if database < min_database_score:
                continue
            if textmining < min_textmining_score:
                continue
            if coexpression < min_coexpression_score:
                continue

            # Edge weight: higher score = lower distance
            weight = 1000 - combined

            # Add edge (undirected)
            graph[str_a][str_b] = weight
            graph[str_b][str_a] = weight

            # Store score for reporting
            edge_key = tuple(sorted([pref_a, pref_b]))
            edge_scores[edge_key] = combined

            # Also map by preferred names
            if pref_a and pref_b:
                graph[pref_a][pref_b] = weight
                graph[pref_b][pref_a] = weight

        # Phase 3: Run Dijkstra on complete graph
        results = {}
        query_genes_upper = [g.upper() for g in gene_list]

        for i, gene_a in enumerate(query_genes_upper):
            for gene_b in query_genes_upper[i+1:]:
                # Try to find path using different identifier formats
                path_found = False

                # Try: gene name -> gene name
                path, distance, scores = self._dijkstra_path(
                    graph, gene_a, gene_b, edge_scores, max_distance
                )

                if path:
                    path_found = True
                # Try: STRING ID -> STRING ID
                elif gene_a in gene_to_string_id and gene_b in gene_to_string_id:
                    path, distance, scores = self._dijkstra_path(
                        graph,
                        gene_to_string_id[gene_a],
                        gene_to_string_id[gene_b],
                        edge_scores,
                        max_distance
                    )

                    if path:
                        # Convert STRING IDs back to gene names
                        path = [string_id_to_gene.get(p, p) for p in path]
                        path_found = True

                if path_found and path:
                    results[(gene_a, gene_b)] = {
                        'path': path,
                        'distance': distance,
                        'hops': len(path) - 1,
                        'scores': scores,
                        'algorithm': 'Dijkstra+BFS',
                        'weight_formula': '1000 - combined_score'
                    }

        return results

    def _dijkstra_path(
        self,
        graph: Dict[str, Dict[str, float]],
        start: str,
        end: str,
        edge_scores: Dict[Tuple[str, str], int],
        max_hops: int
    ) -> Tuple[Optional[List[str]], float, List[int]]:
        """
        Dijkstra's algorithm to find shortest path.

        Returns (path, total_distance, edge_scores_along_path)
        """
        if start not in graph or end not in graph:
            return None, float('inf'), []

        # Priority queue: (distance, node, path, hop_count)
        pq = [(0, start, [start], 0)]
        visited = set()

        while pq:
            dist, node, path, hops = heapq.heappop(pq)

            if node in visited:
                continue

            visited.add(node)

            if node == end:
                # Extract edge scores along path
                scores = []
                for i in range(len(path) - 1):
                    edge_key = tuple(sorted([path[i], path[i+1]]))
                    score = edge_scores.get(edge_key, 0)
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


def get_string_neighbors(
    gene: str,
    species: int = 9606,
    top_n: int = 100,
    min_score: int = 0,
    min_experimental: int = 0,
    min_database: int = 0,
    min_textmining: int = 0,
    min_coexpression: int = 0,
    network_type: str = "functional"
) -> pd.DataFrame:
    """
    Convenience function to get STRING neighbors as DataFrame.

    Parameters
    ----------
    gene : str
        Gene name or protein identifier
    species : int
        NCBI taxonomy ID (default 9606 = human)
    top_n : int
        Maximum neighbors to return
    min_score : int
        Minimum combined confidence score (0-999)
    network_type : str
        'functional' (default) or 'physical'

    Returns
    -------
    pd.DataFrame
        Neighbors with entity and edge annotations
    """
    client = StringClient()
    neighbors = client.get_neighbors(
        identifier=gene,
        species=species,
        top_n=top_n,
        min_combined_score=min_score,
        min_experimental_score=min_experimental,
        min_database_score=min_database,
        min_textmining_score=min_textmining,
        min_coexpression_score=min_coexpression,
        network_type=network_type,
    )

    if not neighbors:
        return pd.DataFrame()

    # Convert to DataFrame
    df = pd.DataFrame([asdict(n) for n in neighbors])
    return df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Query STRING database for protein interaction neighbors."
    )
    parser.add_argument(
        "gene",
        help="Gene name or protein identifier (e.g., TP53)"
    )
    parser.add_argument(
        "--species", type=int, default=9606,
        help="NCBI taxonomy ID (default: 9606 = Homo sapiens)"
    )
    parser.add_argument(
        "--top-n", type=int, default=100,
        help="Maximum number of neighbors to return (default: 100)"
    )
    parser.add_argument(
        "--min-score", type=int, default=0,
        help="Minimum combined confidence score 0-999 (default: 0)"
    )
    parser.add_argument(
        "--min-experimental", type=int, default=0,
        help="Minimum experimental evidence score 0-999 (default: 0)"
    )
    parser.add_argument(
        "--min-database", type=int, default=0,
        help="Minimum curated database score 0-999 (default: 0)"
    )
    parser.add_argument(
        "--min-textmining", type=int, default=0,
        help="Minimum text mining score 0-999 (default: 0)"
    )
    parser.add_argument(
        "--min-coexpression", type=int, default=0,
        help="Minimum coexpression score 0-999 (default: 0)"
    )
    parser.add_argument(
        "--network-type", choices=["functional", "physical"], default="functional",
        help="Type of interaction network (default: functional)"
    )
    parser.add_argument(
        "--format", choices=["table", "csv", "tsv"], default="table",
        help="Output format (default: table)"
    )

    args = parser.parse_args()

    # Query STRING
    df = get_string_neighbors(
        gene=args.gene,
        species=args.species,
        top_n=args.top_n,
        min_score=args.min_score,
        min_experimental=args.min_experimental,
        min_database=args.min_database,
        min_textmining=args.min_textmining,
        min_coexpression=args.min_coexpression,
        network_type=args.network_type
    )

    if df.empty:
        print(f"No neighbors found for '{args.gene}' in taxon {args.species}")
    else:
        if args.format == "table":
            print(df.to_string(index=False))
        elif args.format == "csv":
            print(df.to_csv(index=False))
        elif args.format == "tsv":
            print(df.to_csv(index=False, sep="\t"))

        print(f"\n# Found {len(df)} neighbors for {args.gene}", file=__import__('sys').stderr)
