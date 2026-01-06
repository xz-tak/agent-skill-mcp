#!/usr/bin/env python3
"""
Query neighbor proteins from BioGRID for a given gene in human (or other taxa).

Features
--------
- Uses the BioGRID REST /interactions endpoint with format=json/jsonExtended.
- Given a seed gene (e.g. TP53) and a taxId (default 9606 = Homo sapiens),
  finds neighbor genes up to N hops in the interaction graph.
- Ranks neighbors:
    1) by hop distance (1, 2, 3, ...)
    2) within each hop by best QUANTITATION score (descending) if available.
- Returns a table including neighbor symbol, Entrez ID, organism ID, organism
  name (if /organisms is queried), hop distance, best score and count of
  supporting interactions.
- Supports filters:
    * include / exclude Experimental System names (evidenceList/includeEvidence)
    * min / max QUANTITATION (confidence/score)
    * EXPERIMENTAL_SYSTEM_TYPE (e.g. "physical" or "genetic")
    * throughputTag (any/low/high)
    * interSpeciesExcluded / selfInteractionsExcluded.

Requirements
------------
    pip install requests pandas

Set your key:
    export BIOGRID_API_KEY="your_32_char_key"
"""

import os
from dataclasses import dataclass, asdict
from typing import Dict, Any, List, Optional, Iterable, Set, Tuple
import requests
from collections import defaultdict
import heapq
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

BIOGRID_BASE_URL = "https://webservice.thebiogrid.org"


@dataclass
class NeighborRecord:
    gene_symbol: str
    entrez_id: Optional[int]
    organism_id: int
    organism_name: Optional[str]
    hop: int
    best_score: Optional[float]
    interaction_count: int
    # Edge annotations
    interaction_types: Optional[str] = None  # e.g. "physical|genetic"
    detection_methods: Optional[str] = None  # Experimental System names
    evidence_types: Optional[str] = None     # Experimental System Types
    source_database: str = "BioGRID"
    publications: Optional[str] = None       # PMIDs
    path: Optional[str] = None               # Path from seed, e.g. "TP53-MDM2-ATM"


class BioGRIDClient:
    def __init__(self, access_key: str, base_url: str = BIOGRID_BASE_URL, timeout: int = 60):
        self.access_key = access_key
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout

    # ---------- low-level REST helper ----------

    def _request(self, path: str, params: Dict[str, Any]) -> Any:
        """
        Generic GET wrapper around the BioGRID REST API.

        Automatically injects accesskey and returns parsed JSON.
        Raises RuntimeError if BioGRID returns an error JSON.
        """
        url = f"{self.base_url}/{path.lstrip('/')}"
        params = dict(params)
        params["accesskey"] = self.access_key

        resp = requests.get(url, params=params, timeout=self.timeout)
        resp.raise_for_status()
        data = resp.json()

        # REST errors are returned as JSON with STATUS="Error"
        if isinstance(data, dict) and data.get("STATUS") == "Error":
            raise RuntimeError(f"BioGRID error: {data}")
        return data

    # ---------- helper: organisms ----------

    def get_organism_map(self) -> Dict[int, str]:
        """
        Call /organisms?format=json and build TAX_ID -> name mapping.

        See /organisms access point docs.  
        """
        data = self._request("organisms/", {"format": "json"})
        mapping: Dict[int, str] = {}
        for _, info in data.items():
            try:
                tid = int(info["TAX_ID"])
            except (KeyError, ValueError):
                continue
            name = (
                info.get("OFFICIAL_NAME")
                or info.get("SYSTEMATIC_NAME")
                or info.get("OFFICIAL_SYMBOL")
                or str(tid)
            )
            mapping[tid] = name
        return mapping

    # ---------- helper: fetch interactions for a set of genes ----------

    def get_interactions_for_genes(
        self,
        genes: Iterable[str],
        tax_id: str = "9606",
        evidence_include: Optional[Iterable[str]] = None,
        evidence_exclude: Optional[Iterable[str]] = None,
        experimental_system_types: Optional[Iterable[str]] = None,
        throughput_tag: str = "any",
        inter_species_excluded: bool = True,
        self_interactions_excluded: bool = True,
        start: int = 0,
        max_results: int = 10000,
        format: str = "json",
    ) -> Dict[str, Any]:
        """
        Fetch interactions for one or more genes using the BioGRID /interactions
        endpoint (format=json or jsonExtended).  

        geneList is constructed as a pipe-separated list, searched in OFFICIAL_SYMBOL
        (searchNames=true) with the given taxId constraint.  
        """
        genes = [g for g in genes if g]
        if not genes:
            return {}

        if evidence_include and evidence_exclude:
            raise ValueError("Use either evidence_include or evidence_exclude, not both")

        params: Dict[str, Any] = {
            "format": format,
            "start": start,
            "max": max_results,
            "geneList": "|".join(sorted(set(genes))),
            "searchNames": "true",
            "includeInteractors": "true",
            "includeInteractorInteractions": "false",  # first-order only for given geneList
            "taxId": tax_id,
            "interSpeciesExcluded": str(inter_species_excluded).lower(),
            "selfInteractionsExcluded": str(self_interactions_excluded).lower(),
            "throughputTag": throughput_tag,
        }

        # Evidence filtering via evidenceList / includeEvidence
        # evidenceList is a | separated list of Experimental System names. 
        if evidence_include:
            params["evidenceList"] = "|".join(evidence_include)
            params["includeEvidence"] = "true"
        elif evidence_exclude:
            params["evidenceList"] = "|".join(evidence_exclude)
            params["includeEvidence"] = "false"

        data = self._request("interactions/", params)

        # Optional local filter on EXPERIMENTAL_SYSTEM_TYPE (physical / genetic)
        if experimental_system_types:
            keep_types = {t.lower() for t in experimental_system_types}
            data = {
                iid: rec
                for iid, rec in data.items()
                if rec.get("EXPERIMENTAL_SYSTEM_TYPE", "").lower() in keep_types
            }

        return data

    # ---------- main: neighbors / n-hop graph expansion ----------

    def get_neighbors(
        self,
        seed_gene: str,
        tax_id: str = "9606",
        max_hops: int = 1,
        max_neighbors: int = 100,
        min_score: Optional[float] = None,
        max_score: Optional[float] = None,
        evidence_include: Optional[Iterable[str]] = None,
        evidence_exclude: Optional[Iterable[str]] = None,
        experimental_system_types: Optional[Iterable[str]] = None,
        throughput_tag: str = "any",
        inter_species_excluded: bool = True,
        self_interactions_excluded: bool = True,
        organism_map: Optional[Dict[int, str]] = None,
    ) -> List[NeighborRecord]:
        """
        BFS-style expansion starting from seed_gene.

        hop = 1: genes that interact directly with seed_gene.
        hop = 2: genes discovered by expanding from hop=1 frontier, etc.

        Ranking:
            - neighbors are grouped by hop
            - within each hop they are sorted by best QUANTITATION score (desc)
              across all interactions that touch that neighbor.  
        """
        seed_gene = seed_gene.strip()
        visited: Set[str] = {seed_gene}
        frontier: Set[str] = {seed_gene}

        neighbors: Dict[str, NeighborRecord] = {}
        # Track paths: gene -> shortest path from seed
        paths: Dict[str, List[str]] = {seed_gene: [seed_gene]}

        hop = 1
        while frontier and hop <= max_hops:
            interactions = self.get_interactions_for_genes(
                genes=frontier,
                tax_id=tax_id,
                evidence_include=evidence_include,
                evidence_exclude=evidence_exclude,
                experimental_system_types=experimental_system_types,
                throughput_tag=throughput_tag,
                inter_species_excluded=inter_species_excluded,
                self_interactions_excluded=self_interactions_excluded,
                format="json",
            )

            new_frontier: Set[str] = set()

            for _, rec in interactions.items():
                a = rec.get("OFFICIAL_SYMBOL_A")
                b = rec.get("OFFICIAL_SYMBOL_B")
                entrez_a = rec.get("ENTREZ_GENE_A")
                entrez_b = rec.get("ENTREZ_GENE_B")
                org_a = rec.get("ORGANISM_A")
                org_b = rec.get("ORGANISM_B")
                quant = rec.get("QUANTITATION", "-")

                # Edge annotation fields
                exp_system = rec.get("EXPERIMENTAL_SYSTEM", "")
                exp_system_type = rec.get("EXPERIMENTAL_SYSTEM_TYPE", "")
                pubmed_id = rec.get("PUBMED_ID", "")

                # QUANTITATION is the numeric score field used for P-values,
                # confidence scores, etc.; "-" if none reported.
                score: Optional[float] = None
                if quant not in ("-", "", None):
                    try:
                        score = float(quant)
                    except ValueError:
                        score = None

                # Optionally drop interactions based on score thresholds
                if score is not None:
                    if min_score is not None and score < min_score:
                        continue
                    if max_score is not None and score > max_score:
                        continue

                # Determine which gene in the pair is in the current frontier (source)
                # and which is the neighbor (target)
                source_gene = None
                for check_gene in (a, b):
                    if check_gene in frontier:
                        source_gene = check_gene
                        break

                for gene, entrez, org in ((a, entrez_a, org_a), (b, entrez_b, org_b)):
                    if gene is None:
                        continue
                    if gene == seed_gene:
                        # We only want *neighbors*, not the seed itself
                        continue
                    if gene == source_gene:
                        # Skip the source gene itself in this pair
                        continue

                    try:
                        org_id = int(org) if org is not None else 0
                    except ValueError:
                        org_id = 0

                    # Build path for this neighbor
                    if gene not in paths and source_gene and source_gene in paths:
                        paths[gene] = paths[source_gene] + [gene]

                    if gene not in neighbors:
                        organism_name = None
                        if organism_map and org_id in organism_map:
                            organism_name = organism_map[org_id]

                        # Format path as "A-B-C"
                        path_str = "-".join(paths.get(gene, [seed_gene, gene]))

                        neighbors[gene] = NeighborRecord(
                            gene_symbol=gene,
                            entrez_id=int(entrez) if entrez not in (None, "-", "") else None,
                            organism_id=org_id,
                            organism_name=organism_name,
                            hop=hop,
                            best_score=score,
                            interaction_count=1,
                            detection_methods=exp_system if exp_system else None,
                            evidence_types=exp_system_type if exp_system_type else None,
                            publications=pubmed_id if pubmed_id else None,
                            path=path_str,
                        )
                    else:
                        n = neighbors[gene]
                        n.interaction_count += 1
                        if score is not None:
                            if n.best_score is None or score > n.best_score:
                                n.best_score = score

                        # Accumulate unique edge annotations
                        if exp_system and exp_system not in (n.detection_methods or ""):
                            n.detection_methods = f"{n.detection_methods}|{exp_system}" if n.detection_methods else exp_system
                        if exp_system_type and exp_system_type not in (n.evidence_types or ""):
                            n.evidence_types = f"{n.evidence_types}|{exp_system_type}" if n.evidence_types else exp_system_type
                        if pubmed_id and pubmed_id not in (n.publications or ""):
                            n.publications = f"{n.publications}|{pubmed_id}" if n.publications else pubmed_id

                    # Expand BFS frontier only for genes in the chosen tax_id
                    if gene not in visited and (not tax_id or str(org_id) == str(tax_id)):
                        new_frontier.add(gene)
                        visited.add(gene)

            frontier = new_frontier
            hop += 1

        # Sort neighbors: primary key = hop, secondary = best_score (descending, None last)
        def sort_key(n: NeighborRecord):
            score_key = n.best_score if n.best_score is not None else float("-inf")
            return (n.hop, -score_key)

        sorted_neighbors = sorted(neighbors.values(), key=sort_key)

        return sorted_neighbors[:max_neighbors]

    def find_shortest_paths(
        self,
        gene_list: List[str],
        tax_id: str = "9606",
        max_distance: int = 50,
        min_score: Optional[float] = 0.4,
        experimental_system_types: Optional[Iterable[str]] = None,
        throughput_tag: str = "any",
        inter_species_excluded: bool = True,
        self_interactions_excluded: bool = True,
    ) -> Dict[Tuple[str, str], Dict[str, Any]]:
        """
        Find shortest paths between all pairs of query genes using BioGRID network.

        **Algorithm**: Dijkstra's shortest path with score-based edge weighting
        **Edge weight formula**: weight = 1.0 / (score + 0.01)
        **Principle**: Higher QUANTITATION score = lower weight = shorter distance (closer nodes)

        Parameters
        ----------
        gene_list : List[str]
            List of gene symbols to find paths between
        tax_id : str
            NCBI taxonomy ID (default "9606" = human)
        max_distance : int
            Maximum number of intermediate nodes (default 3)
        min_score : float, optional
            Minimum QUANTITATION score to include edges
        experimental_system_types : Iterable[str], optional
            Filter by system type (e.g., ["physical"])
        throughput_tag : str
            Filter by throughput (any/low/high)
        inter_species_excluded : bool
            Exclude inter-species interactions
        self_interactions_excluded : bool
            Exclude self-interactions

        Returns
        -------
        Dict[Tuple[str, str], Dict]
            Dictionary mapping (geneA, geneB) tuples to path info:
            {
                'path': ['A', 'D', 'B'],
                'distance': 2.5,
                'hops': 2,
                'scores': [0.8, 0.9]  # QUANTITATION scores along path
            }
        """
        if len(gene_list) < 2:
            return {}

        # CRITICAL FIX: Query each gene SEPARATELY to ensure all direct edges are captured
        # (BioGRID may return different results when querying genes together vs separately)
        all_interactions = {}
        interactions_lock = threading.Lock()

        def query_gene_interactions(gene: str):
            """Query BioGRID for a single gene (thread-safe)."""
            try:
                return self.get_interactions_for_genes(
                    genes=[gene],  # Query one gene at a time
                    tax_id=tax_id,
                    experimental_system_types=experimental_system_types,
                    throughput_tag=throughput_tag,
                    inter_species_excluded=inter_species_excluded,
                    self_interactions_excluded=self_interactions_excluded,
                    max_results=50000,
                    format="json",
                )
            except Exception:
                return {}

        # Query genes in parallel for performance
        with ThreadPoolExecutor(max_workers=10) as executor:
            future_to_gene = {executor.submit(query_gene_interactions, gene): gene
                              for gene in gene_list}

            for future in as_completed(future_to_gene):
                try:
                    interactions = future.result()
                    # Merge interactions (thread-safe)
                    with interactions_lock:
                        all_interactions.update(interactions)
                except Exception:
                    # Continue with other genes if one fails
                    pass

        if not all_interactions:
            return {}

        # Build weighted graph
        graph: Dict[str, Dict[str, float]] = defaultdict(dict)
        edge_scores: Dict[Tuple[str, str], Optional[float]] = {}

        for _, rec in all_interactions.items():
            a = rec.get("OFFICIAL_SYMBOL_A")
            b = rec.get("OFFICIAL_SYMBOL_B")
            quant = rec.get("QUANTITATION", "-")

            if not a or not b:
                continue

            # Parse score
            score: Optional[float] = None
            if quant not in ("-", "", None):
                try:
                    score = float(quant)
                except ValueError:
                    score = None

            # Apply score filter
            if min_score is not None and score is not None and score < min_score:
                continue

            # Edge weight: higher score = lower distance
            # Use 1/(score+0.01) to convert score to distance
            if score is not None and score > 0:
                weight = 1.0 / (score + 0.01)
            else:
                weight = 1.0  # Default weight for edges without scores

            # Add edge (undirected)
            a_upper = a.upper()
            b_upper = b.upper()

            # Keep edge with best (lowest) weight if duplicate
            edge_key1 = (a_upper, b_upper)
            edge_key2 = (b_upper, a_upper)

            if edge_key1 not in graph.get(a_upper, {}) or weight < graph[a_upper][b_upper]:
                graph[a_upper][b_upper] = weight
                graph[b_upper][a_upper] = weight

            # Store score for reporting (keep best/highest score)
            edge_tuple = tuple(sorted([a_upper, b_upper]))
            if edge_tuple not in edge_scores or (score is not None and (edge_scores[edge_tuple] is None or score > edge_scores[edge_tuple])):
                edge_scores[edge_tuple] = score

        # Find shortest paths between all pairs using Dijkstra
        results = {}
        query_genes_upper = [g.upper() for g in gene_list]

        for i, gene_a in enumerate(query_genes_upper):
            for gene_b in query_genes_upper[i+1:]:
                path, distance, scores = self._dijkstra_path(
                    graph, gene_a, gene_b, edge_scores, max_distance
                )

                if path:
                    results[(gene_a, gene_b)] = {
                        'path': path,
                        'distance': distance,
                        'hops': len(path) - 1,
                        'scores': scores,
                        'algorithm': 'Dijkstra',
                        'weight_formula': '1.0 / (score + 0.01)'
                    }

        return results

    def _dijkstra_path(
        self,
        graph: Dict[str, Dict[str, float]],
        start: str,
        end: str,
        edge_scores: Dict[Tuple[str, str], Optional[float]],
        max_hops: int
    ) -> Tuple[Optional[List[str]], float, List[Optional[float]]]:
        """
        Dijkstra's algorithm for BioGRID network.

        Returns (path, total_distance, quantitation_scores_along_path)
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
                # Extract QUANTITATION scores along path
                scores = []
                for i in range(len(path) - 1):
                    edge_key = tuple(sorted([path[i], path[i+1]]))
                    score = edge_scores.get(edge_key)
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


# ---------- CLI usage ----------

if __name__ == "__main__":
    import argparse
    try:
        import pandas as pd
    except ImportError:
        pd = None

    parser = argparse.ArgumentParser(
        description="Query BioGRID neighbor proteins for a human gene."
    )
    parser.add_argument("gene", help="Seed gene symbol, e.g. TP53")
    parser.add_argument(
        "--api-key", dest="api_key", default="d3367ed24eeea8fe8718f4993ed63ec9",
        help="BioGRID access key (or set BIOGRID_API_KEY env var)."
    )
    parser.add_argument(
        "--tax-id", default="9606",
        help="NCBI taxonomy ID for seed organism (default: 9606 = human)."
    )
    parser.add_argument(
        "--max-hops", type=int, default=1,
        help="Maximum hop distance (graph radius)."
    )
    parser.add_argument(
        "--max-neighbors", type=int, default=100,
        help="Maximum number of neighbors to return."
    )
    parser.add_argument(
        "--min-score", type=float, default=None,
        help="Minimum QUANTITATION score to keep."
    )
    parser.add_argument(
        "--max-score", type=float, default=None,
        help="Maximum QUANTITATION score to keep."
    )
    parser.add_argument(
        "--include-evidence", type=str, default=None,
        help="Comma-separated Experimental System names to INCLUDE "
             "(evidenceList, includeEvidence=true)."
    )
    parser.add_argument(
        "--exclude-evidence", type=str, default=None,
        help="Comma-separated Experimental System names to EXCLUDE "
             "(evidenceList, includeEvidence=false)."
    )
    parser.add_argument(
        "--system-types", type=str, default=None,
        help="Comma-separated EXPERIMENTAL_SYSTEM_TYPE values to keep, "
             "e.g. physical,genetic."
    )
    parser.add_argument(
        "--throughput", type=str, default="any", choices=["any", "low", "high"],
        help="Filter by throughputTag (any/low/high)."
    )
    parser.add_argument(
        "--allow-interspecies", action="store_true",
        help="Include inter-species interactions (interSpeciesExcluded=false)."
    )
    parser.add_argument(
        "--keep-self-interactions", action="store_true",
        help="Keep self-interactions (selfInteractionsExcluded=false)."
    )
    parser.add_argument(
        "--resolve-organisms", action="store_true",
        help="Resolve organism IDs to names via /organisms endpoint."
    )

    args = parser.parse_args()

    api_key = args.api_key or os.environ.get("BIOGRID_API_KEY")
    if not api_key:
        raise SystemExit("ERROR: Provide a BioGRID API key via --api-key or BIOGRID_API_KEY.")

    client = BioGRIDClient(api_key)
    organism_map = client.get_organism_map() if args.resolve_organisms else None

    evidence_include = args.include_evidence.split(",") if args.include_evidence else None
    evidence_exclude = args.exclude_evidence.split(",") if args.exclude_evidence else None
    system_types = args.system_types.split(",") if args.system_types else None

    neighbors = client.get_neighbors(
        seed_gene=args.gene,
        tax_id=args.tax_id,
        max_hops=args.max_hops,
        max_neighbors=args.max_neighbors,
        min_score=args.min_score,
        max_score=args.max_score,
        evidence_include=evidence_include,
        evidence_exclude=evidence_exclude,
        experimental_system_types=system_types,
        throughput_tag=args.throughput,
        inter_species_excluded=not args.allow_interspecies,
        self_interactions_excluded=not args.keep_self_interactions,
        organism_map=organism_map,
    )

    rows = [asdict(n) for n in neighbors]

    if not rows:
        print("No neighbors found.")
    else:
        if pd:
            df = pd.DataFrame(rows)
            print(df.to_string(index=False))
        else:
            # Plain TSV output
            headers = list(rows[0].keys())
            print("\t".join(headers))
            for row in rows:
                print("\t".join(str(row[h]) for h in headers))
