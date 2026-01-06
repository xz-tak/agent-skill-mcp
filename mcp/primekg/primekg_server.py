#!/usr/bin/env python3
"""
PrimeKG MCP Server

A Model Context Protocol server for querying the PrimeKG biomedical knowledge graph.
Provides tools for exploring entity types, relation types, resources, and entity connections.
"""

import os
import pickle
import logging
from pathlib import Path
from typing import Optional, List, Dict, Any
import polars as pl
from fastmcp import FastMCP
import primekg_setup

# Initialize logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Initialize MCP server
mcp = FastMCP("primekg")

# Global data storage (loaded on first use)
_data_cache = {
    "nodes": None,
    "edges": None,
    "entity_index": None,  # For fast entity lookup
    "ent2name_dict": None,
}

# MCP root directory
MCP_ROOT = Path(__file__).parent

# Default data paths
DATA_DIR = None  # Will be set by find_primekg_data()


def find_primekg_data() -> Optional[Path]:
    """
    Find existing PrimeKG data by checking multiple locations.

    Search order:
    1. ultra-inference/data/
    2. ultraquery-inference/data/
    3. Local primekg/data/ (after auto-download)

    Returns:
        Path to valid PrimeKG data directory, or None if not found
    """
    # Define search locations (relative to MCP root's parent directory)
    mcp_services_root = MCP_ROOT.parent
    search_locations = [
        mcp_services_root / "ultra-inference" / "data",
        mcp_services_root / "ultraquery-inference" / "data",
        MCP_ROOT / "data",  # Local data directory
    ]

    logger.info("Searching for existing PrimeKG data...")

    for location in search_locations:
        logger.info(f"  Checking: {location}")

        # Check if directory exists
        if not location.exists():
            logger.info(f"    Directory does not exist")
            continue

        # Validate required files exist
        required_files = [
            location / "primekg.csv",
            location / "primekg1" / "raw" / "nodes.txt",
        ]

        all_exist = all(f.exists() for f in required_files)

        if all_exist:
            logger.info(f"  ✓ Found valid PrimeKG data at: {location}")
            return location
        else:
            missing = [f.name for f in required_files if not f.exists()]
            logger.info(f"    Missing files: {missing}")

    logger.info("  No existing PrimeKG data found")
    return None


def load_data():
    """Load PrimeKG data files into memory (lazy loading)."""
    global DATA_DIR

    if _data_cache["nodes"] is None:
        # Find or download PrimeKG data
        if DATA_DIR is None:
            DATA_DIR = find_primekg_data()

            # If no data found, auto-download
            if DATA_DIR is None:
                logger.info("No existing PrimeKG data found. Starting auto-download...")
                local_data_dir = MCP_ROOT / "data"
                dataset_path = local_data_dir / "primekg1"

                try:
                    # Run setup
                    result = primekg_setup.setup_primekg(
                        dataset_path=str(dataset_path),
                        force_redownload=False,
                    )

                    if result["status"] in ["completed", "already_exists"]:
                        DATA_DIR = local_data_dir
                        logger.info(f"✓ PrimeKG data ready at: {DATA_DIR}")
                    else:
                        raise RuntimeError("Failed to setup PrimeKG data")

                except Exception as e:
                    logger.error(f"Failed to auto-download PrimeKG data: {e}")
                    raise RuntimeError(
                        f"Could not find or download PrimeKG data. Error: {e}"
                    )

        logger.info(f"Loading PrimeKG data from {DATA_DIR}...")

        # Load nodes
        _data_cache["nodes"] = pl.read_csv(
            os.path.join(DATA_DIR, "primekg1", "raw", "nodes.txt"),
            separator="\t",
            has_header=True,
            schema={
                "source_id": pl.Utf8,
                "name": pl.Utf8,
                "type": pl.Categorical,
                "source": pl.Categorical,
                "source_label": pl.Utf8,
            },
        )

        # Load edges from the full CSV (more detailed than graph.txt)
        _data_cache["edges"] = pl.read_csv(
            os.path.join(DATA_DIR, "primekg.csv"),
            has_header=True,
            schema={
                "relation": pl.Categorical,
                "display_relation": pl.Categorical,
                "x_index": pl.Int64,
                "x_id": pl.Utf8,
                "x_type": pl.Categorical,
                "x_name": pl.Utf8,
                "x_source": pl.Categorical,
                "y_index": pl.Int64,
                "y_id": pl.Utf8,
                "y_type": pl.Categorical,
                "y_name": pl.Utf8,
                "y_source": pl.Categorical,
            },
        ).with_columns(
            (pl.col("x_source") + ":" + pl.col("x_id")).alias("x_source_label"),
            (pl.col("y_source") + ":" + pl.col("y_id")).alias("y_source_label"),
        )

        # Create entity index for fast lookup
        _data_cache["entity_index"] = {
            row["source_label"]: idx
            for idx, row in enumerate(_data_cache["nodes"].iter_rows(named=True))
        }

        # Load pickle dictionaries if they exist
        try:
            with open(os.path.join(DATA_DIR, "ent2name_dict.pkl"), "rb") as f:
                _data_cache["ent2name_dict"] = pickle.load(f)
        except FileNotFoundError:
            logger.warning("Warning: ent2name_dict.pkl not found")
            _data_cache["ent2name_dict"] = {}

        logger.info(
            f"✓ Loaded {len(_data_cache['nodes'])} nodes and {len(_data_cache['edges'])} edges"
        )


def get_entity_types() -> Dict[str, Any]:
    """
    Get all unique entity types in PrimeKG.

    Returns a dictionary containing:
    - entity_types: List of unique entity types
    - count_by_type: Count of entities for each type
    """
    load_data()

    type_counts = (
        _data_cache["nodes"]
        .group_by("type")
        .agg(pl.len().alias("count"))
        .sort("count", descending=True)
    )

    return {
        "entity_types": type_counts["type"].to_list(),
        "count_by_type": {
            row["type"]: row["count"] for row in type_counts.iter_rows(named=True)
        },
        "total_entities": len(_data_cache["nodes"]),
    }


def get_relation_types() -> Dict[str, Any]:
    """
    Get all unique relation types in PrimeKG.

    Returns a dictionary containing:
    - relation_types: List of unique relation types
    - display_relations: List of human-readable relation names
    - count_by_relation: Count of edges for each relation type
    """
    load_data()

    relation_counts = (
        _data_cache["edges"]
        .group_by(["display_relation", "relation"])
        .agg(pl.len().alias("count"))
        .sort("count", descending=True)
    )

    return {
        "relation_types": relation_counts["display_relation"].to_list(),
        "count_by_relation": {
            row["relation"]: {
                "display_name": row["display_relation"],
                "count": row["count"],
            }
            for row in relation_counts.iter_rows(named=True)
        },
        "total_edges": len(_data_cache["edges"]),
    }


def get_data_sources() -> Dict[str, Any]:
    """
    Get all data sources/resources included in PrimeKG.

    Returns a dictionary containing:
    - sources: List of unique data sources
    - count_by_source: Count of entities from each source
    - source_by_type: Breakdown of sources by entity type
    """
    load_data()

    source_counts = (
        _data_cache["nodes"]
        .group_by("source")
        .agg(pl.count().alias("count"))
        .sort("count", descending=True)
    )

    source_type_counts = (
        _data_cache["nodes"]
        .group_by(["source", "type"])
        .agg(pl.count().alias("count"))
        .sort(["source", "count"], descending=[False, True])
    )

    # Organize by source -> types
    source_by_type = {}
    for row in source_type_counts.iter_rows(named=True):
        source = row["source"]
        if source not in source_by_type:
            source_by_type[source] = {}
        source_by_type[source][row["type"]] = row["count"]

    return {
        "sources": source_counts["source"].to_list(),
        "count_by_source": {
            row["source"]: row["count"] for row in source_counts.iter_rows(named=True)
        },
        "source_by_type": source_by_type,
        "total_sources": len(source_counts),
    }


def check_entity_exists(entity_id: str) -> Dict[str, Any]:
    """
    Check if an entity exists in PrimeKG and return its details.

    Args:
        entity_id: Entity ID in format "SOURCE:ID" (e.g., "NCBI:9796", "DRUGBANK:DB00001")

    Returns:
        Dictionary with entity information if found, or error message if not found.
    """
    load_data()

    # Check in entity index
    if entity_id in _data_cache["entity_index"]:
        idx = _data_cache["entity_index"][entity_id]
        entity_row = _data_cache["nodes"].row(idx, named=True)

        # Get human-readable name if available
        human_name = _data_cache["ent2name_dict"].get(entity_id, entity_row.get("name"))

        return {
            "exists": True,
            "entity_id": entity_id,
            "name": entity_row.get("name"),
            "human_name": human_name,
            "type": entity_row.get("type"),
            "source": entity_row.get("source"),
            "source_id": entity_row.get("source_id"),
        }
    else:
        return {
            "exists": False,
            "entity_id": entity_id,
            "message": f"Entity {entity_id} not found in PrimeKG",
        }


def check_relation_exists(relation_type: str) -> Dict[str, Any]:
    """
    Check if a relation type exists in PrimeKG.

    Args:
        relation_type: Relation type (e.g., "protein_protein", "drug_drug", "disease_protein")

    Returns:
        Dictionary with relation information if found, or error message if not found.
    """
    load_data()

    relation_info = (
        _data_cache["edges"]
        .filter(pl.col("relation") == relation_type)
        .select(["relation", "display_relation"])
        .unique()
    )

    if len(relation_info) > 0:
        row = relation_info.row(0, named=True)
        count = _data_cache["edges"].filter(pl.col("relation") == relation_type).height

        return {
            "exists": True,
            "relation_type": relation_type,
            "display_name": row["display_relation"],
            "edge_count": count,
        }
    else:
        return {
            "exists": False,
            "relation_type": relation_type,
            "message": f"Relation type {relation_type} not found in PrimeKG",
        }


def get_entity_connections(
    entity_id: str,
    relation_type: Optional[str] = None,
    limit: int = 100,
    dataframe_output: bool = False,
    save_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Get all connections (edges) for a specific entity.

    Args:
        entity_id: Entity ID in format "SOURCE:ID" (e.g., "NCBI:9796")
        relation_type: Optional filter by display relation type (e.g., "ppi", "indication")
        limit: Maximum number of connections to return (default: 100)
        dataframe_output: If True, return connections as a DataFrame and save to file
        save_dir: Directory to save the connections file if dataframe_output is True

    Returns:
        Dictionary containing:
        - entity_id: The queried entity
        - entity_info: Details about the entity
        - connections: List of connected entities with relation information
        - total_connections: Total number of connections (before limit)
    """
    load_data()

    # First check if entity exists
    entity_check = check_entity_exists(entity_id)
    if not entity_check["exists"]:
        return entity_check

    # Find connections where entity is either source (x) or target (y)
    edges_as_source = _data_cache["edges"].filter(pl.col("x_source_label") == entity_id)
    edges_as_target = _data_cache["edges"].filter(pl.col("y_source_label") == entity_id)

    # Apply relation type filter if provided (use display_relation)
    if relation_type:
        edges_as_source = edges_as_source.filter(
            pl.col("display_relation") == relation_type
        )
        edges_as_target = edges_as_target.filter(
            pl.col("display_relation") == relation_type
        )

    total_connections = len(edges_as_source) + len(edges_as_target)
    if dataframe_output:
        save_dir = save_dir if save_dir is not None else DATA_DIR
        os.makedirs(save_dir, exist_ok=True)
        combined_edges = pl.concat([edges_as_source, edges_as_target])
        combined_edges.drop(
            ["relation", "x_index", "x_source", "x_id", "y_index", "y_source", "y_id"]
        ).rename(
            {
                "x_name": "head",
                "x_type": "head_type",
                "x_source_label": "head_label",
                "display_relation": "relation",
                "y_name": "tail",
                "y_type": "tail_type",
                "y_source_label": "tail_label",
            }
        ).write_csv(
            os.path.join(save_dir or DATA_DIR, f"{entity_id}_primekg_connections.tsv"),
            separator="\t",
            include_header=True,
        )

    # Prepare outgoing connections (entity is source)
    outgoing = []
    for row in edges_as_source.head(limit // 2).iter_rows(named=True):
        outgoing.append(
            {
                "direction": "outgoing",
                "relation": row["display_relation"],
                "relation_type": row["relation"],
                "target_entity_id": f"{row['y_source']}:{row['y_id']}",
                "target_name": row["y_name"],
                "target_type": row["y_type"],
            }
        )

    # Prepare incoming connections (entity is target)
    incoming = []
    for row in edges_as_target.head(limit // 2).iter_rows(named=True):
        incoming.append(
            {
                "direction": "incoming",
                "relation": row["display_relation"],
                "relation_type": row["relation"],
                "source_entity_id": f"{row['x_source']}:{row['x_id']}",
                "source_name": row["x_name"],
                "source_type": row["x_type"],
            }
        )

    return {
        "entity_id": entity_id,
        "entity_info": entity_check,
        "connections": outgoing + incoming,
        "total_connections": total_connections,
        "returned_connections": len(outgoing) + len(incoming),
        "limit": limit,
    }


@mcp.tool()
def search_entities(
    search_term: str,
    entity_type: Optional[str] = None,
    source: Optional[str] = None,
    limit: int = 20,
) -> Dict[str, Any]:
    """
    Search for entities by name or ID.

    Args:
        search_term: Text to search for in entity names
        entity_type: Optional filter by entity type (e.g., "gene/protein", "drug", "disease")
        source: Optional filter by data source (e.g., "NCBI", "DRUGBANK")
        limit: Maximum number of results to return (default: 20)

    Returns:
        Dictionary containing:
        - search_term: The query term
        - filters: Applied filters
        - results: List of matching entities
        - total_matches: Total number of matches (before limit)
    """
    load_data()

    # Start with all nodes
    filtered = _data_cache["nodes"]

    # Apply search term (case-insensitive partial match on name)
    if search_term:
        filtered = filtered.filter(
            pl.col("name").str.to_lowercase().str.contains(search_term.lower())
        )

    # Apply type filter
    if entity_type:
        filtered = filtered.filter(pl.col("type") == entity_type)

    # Apply source filter
    if source:
        filtered = filtered.filter(pl.col("source") == source)

    total_matches = len(filtered)
    results = []

    for row in filtered.head(limit).iter_rows(named=True):
        entity_id = row["source_label"]
        human_name = _data_cache["ent2name_dict"].get(entity_id, row["name"])

        results.append(
            {
                "entity_id": entity_id,
                "name": row["name"],
                "human_name": human_name,
                "type": row["type"],
                "source": row["source"],
                "source_id": row["source_id"],
            }
        )

    return {
        "search_term": search_term,
        "filters": {
            "entity_type": entity_type,
            "source": source,
        },
        "results": results,
        "total_matches": total_matches,
        "returned_results": len(results),
        "limit": limit,
    }


@mcp.tool()
def get_entity_neighborhoods(
    entity_id: str,
    max_depth: int = 2,
    relation_types: Optional[List[str]] = None,
    max_neighbors_per_level: int = 10,
) -> Dict[str, Any]:
    """
    Get the neighborhood of an entity up to a certain depth (multi-hop connections).

    Args:
        entity_id: Entity ID in format "SOURCE:ID" (e.g., "NCBI:9796")
        max_depth: Maximum depth to traverse (1 or 2, default: 2)
        relation_types: Optional list of display relation types to follow (e.g., ["ppi"])
        max_neighbors_per_level: Max neighbors to explore at each level (default: 10)

    Returns:
        Dictionary containing the entity's neighborhood graph structure.
    """
    load_data()

    # Check entity exists
    entity_check = check_entity_exists(entity_id)
    if not entity_check["exists"]:
        return entity_check

    if max_depth > 2:
        return {
            "error": "Maximum depth of 2 is supported to avoid excessive computation"
        }

    # Get 1-hop neighbors
    first_hop = get_entity_connections(
        entity_id,
        relation_type=relation_types[0] if relation_types else None,
        limit=max_neighbors_per_level * 2,
    )

    neighborhood = {
        "center_entity": entity_check,
        "depth": max_depth,
        "first_hop_neighbors": first_hop["connections"][:max_neighbors_per_level],
        "first_hop_count": first_hop["total_connections"],
    }

    # Get 2-hop neighbors if requested
    if max_depth >= 2:
        second_hop_entities = {}
        for neighbor in neighborhood["first_hop_neighbors"][
            :5
        ]:  # Limit to avoid explosion
            neighbor_id = neighbor.get("target_entity_id") or neighbor.get(
                "source_entity_id"
            )
            if neighbor_id:
                neighbor_connections = get_entity_connections(
                    neighbor_id,
                    relation_type=relation_types[0] if relation_types else None,
                    limit=5,
                )
                second_hop_entities[neighbor_id] = {
                    "entity": neighbor,
                    "connections": neighbor_connections["connections"],
                }

        neighborhood["second_hop_neighbors"] = second_hop_entities

    return neighborhood


@mcp.tool()
def get_statistics() -> Dict[str, Any]:
    """
    Get overall statistics about the PrimeKG dataset.

    Returns:
        Dictionary with comprehensive dataset statistics.
    """
    load_data()

    entity_stats = get_entity_types()
    relation_stats = get_relation_types()
    source_stats = get_data_sources()

    return {
        "dataset": "PrimeKG",
        "total_entities": entity_stats["total_entities"],
        "total_edges": relation_stats["total_edges"],
        "total_entity_types": len(entity_stats["entity_types"]),
        "total_relation_types": len(relation_stats["relation_types"]),
        "total_data_sources": source_stats["total_sources"],
        "entity_type_distribution": entity_stats["count_by_type"],
        "relation_type_distribution": relation_stats["count_by_relation"],
        "data_source_distribution": source_stats["count_by_source"],
    }


mcp.tool(get_entity_types)
mcp.tool(get_relation_types)
mcp.tool(get_data_sources)
mcp.tool(check_entity_exists)
mcp.tool(check_relation_exists)
mcp.tool(get_entity_connections)


@mcp.tool()
def setup_primekg_data(
    force_redownload: bool = False,
    train_frac: float = 0.8,
    test_frac: float = 0.1,
    valid_frac: float = 0.1,
    seed: int = 42,
) -> Dict[str, Any]:
    """
    Download and setup PrimeKG data with train/test/valid splits.

    This function allows manual control over data setup. By default, data is
    automatically downloaded on first query, but this tool can be used to:
    - Force re-download of data
    - Customize train/test/valid split ratios
    - Pre-download data before running queries

    Args:
        force_redownload: If True, re-download even if files exist (default: False)
        train_frac: Training set fraction (default: 0.8)
        test_frac: Test set fraction (default: 0.1)
        valid_frac: Validation set fraction (default: 0.1)
        seed: Random seed for reproducibility (default: 42)

    Returns:
        Dictionary with setup results:
        - success: bool
        - status: str ("already_exists" or "completed")
        - nodes_count: int (number of unique entities)
        - train_edges: int (number of training edges)
        - test_edges: int (number of test edges)
        - valid_edges: int (number of validation edges)
        - data_path: str (path to data directory)
    """
    try:
        # Setup data in local directory
        local_data_dir = MCP_ROOT / "data"
        dataset_path = local_data_dir / "primekg1"

        logger.info(f"Setting up PrimeKG data at {dataset_path}...")

        # Run setup
        result = primekg_setup.setup_primekg(
            dataset_path=str(dataset_path),
            force_redownload=force_redownload,
            train_frac=train_frac,
            test_frac=test_frac,
            valid_frac=valid_frac,
            seed=seed,
        )

        # Update global DATA_DIR if setup succeeded
        global DATA_DIR
        if result["status"] in ["completed", "already_exists"]:
            DATA_DIR = local_data_dir
            logger.info(f"✓ DATA_DIR set to: {DATA_DIR}")

        return {"success": True, "data_path": str(local_data_dir), **result}

    except Exception as e:
        logger.error(f"Failed to setup PrimeKG: {e}")
        import traceback

        traceback.print_exc()
        return {"success": False, "error": str(e), "error_type": type(e).__name__}


if __name__ == "__main__":
    # Run the MCP server
    mcp.run()
