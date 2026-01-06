#!/usr/bin/env python3
"""
ULTRA Inference MCP Server

On-demand inference for ULTRA foundation model on PrimeKG biomedical knowledge graph.
"""
import sys
import os
import logging
import json

# Add src directory to path to import ultra package
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from fastmcp import FastMCP
import torch
from torch_geometric.data import Data
from easydict import EasyDict
import pickle
from typing import Optional, Any
import polars as pl
import time
from difflib import get_close_matches
from ultra.models import Ultra
from ultra import util, tasks, data_util, datasets
import primekg_setup

# Initialize logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Initialize MCP server
mcp = FastMCP("ultra_inference", stateless_http=True)

# Base paths - MCP server is self-contained
MCP_ROOT = os.path.dirname(os.path.abspath(__file__))
ULTRA_ROOT = os.path.dirname(MCP_ROOT)

# Default configuration
DEFAULT_CONFIG = {
    "dataset": {"class": "PrimeKG1", "root": os.path.join(MCP_ROOT, "data")},
    "model": {
        "relation_model": {
            "class": "RelNBFNet",
            "input_dim": 64,
            "hidden_dims": [64, 64, 64, 64, 64, 64],
            "message_func": "distmult",
            "aggregate_func": "sum",
            "short_cut": True,
            "layer_norm": True,
        },
        "entity_model": {
            "class": "EntityNBFNet",
            "input_dim": 64,
            "hidden_dims": [64, 64, 64, 64, 64, 64],
            "message_func": "distmult",
            "aggregate_func": "sum",
            "short_cut": True,
            "layer_norm": True,
        },
    },
    "checkpoint": os.getenv(
        "ULTRA_CHECKPOINT_PATH",
        os.path.join(MCP_ROOT, "ckpts/ultra_primekg_50g_ft_epoch_1.pth"),
    ),
    "output_dir": os.getenv("ULTRA_OUTPUT_DIR", os.path.join(MCP_ROOT, "output")),
    "batch_size": 4,
    "device": "auto",
}


def load_pickle(file_path: str) -> Any:
    """
    Load a pickle file. Helper function to avoid code duplication.

    Args:
        file_path: Path to pickle file

    Returns:
        Deserialized object from pickle file
    """
    with open(file_path, "rb") as f:
        return pickle.load(f)


def ensure_dataset_dictionaries(
    dataset_root: str, dataset_name: str, dataset=None
) -> None:
    """
    Ensure all required pickle files exist for dataset translation.

    This function creates:
    - id2ent_dict.pkl, ent2id_dict.pkl, id2rel_dict.pkl, rel2id_dict.pkl in dataset_root/
    - ent2name_dict.pkl in dataset_root/

    These files are critical for maintaining consistency between the dataset embeddings
    and the translation of entity/relation IDs to human-readable names.

    Args:
        dataset_root: Path to data directory (e.g., /path/to/data)
        dataset_name: Name of dataset (e.g., "primekg1" - lowercase dataset.name)
        dataset: Optional pre-loaded dataset object to avoid reloading
    """

    dataset_dir = os.path.join(dataset_root, dataset_name.lower())

    # Check if entity/relation dictionaries exist
    dict_files = [
        "id2ent_dict.pkl",
        "ent2id_dict.pkl",
        "id2rel_dict.pkl",
        "rel2id_dict.pkl",
    ]
    all_exist = all(os.path.exists(os.path.join(dataset_root, f)) for f in dict_files)

    if not all_exist:
        logger.info("Creating entity/relation dictionaries from dataset...")

        # Load dataset if not provided
        if dataset is None:
            # dataset_name comes from dataset.name which is lowercase (e.g., "primekg1")
            # We need to find the matching class (e.g., "PrimeKG1")
            logger.info(f"Loading dataset {dataset_name} to extract vocabularies...")

            # Find the dataset class by searching for matching name attribute
            dataset_class = None
            for attr_name in dir(datasets):
                attr = getattr(datasets, attr_name)
                if (
                    isinstance(attr, type)
                    and hasattr(attr, "name")
                    and attr.name == dataset_name.lower()
                ):
                    dataset_class = attr
                    break

            if dataset_class is None:
                raise ValueError(
                    f"Could not find dataset class for name '{dataset_name}'. "
                    f"Available datasets: {[cls for cls in dir(datasets) if 'PrimeKG' in cls]}"
                )

            dataset = dataset_class(root=dataset_root)
            logger.info(
                f"Dataset loaded: {len(dataset[0].edge_index[0])} edges in training graph"
            )

        # Create dictionaries directly in dataset_root
        logger.info(f"Creating dictionary files in {dataset_root}")

        # Create id2ent_dict
        id2ent_dict_path = os.path.join(dataset_root, "id2ent_dict.pkl")
        if not os.path.exists(id2ent_dict_path):
            id2ent_dict = {v: k for k, v in dataset.entity_vocab.items()}
            with open(id2ent_dict_path, "wb") as f:
                pickle.dump(id2ent_dict, f)
            logger.info(f"Created id2ent_dict.pkl with {len(id2ent_dict)} entities")

        # Create ent2id_dict
        ent2id_dict_path = os.path.join(dataset_root, "ent2id_dict.pkl")
        if not os.path.exists(ent2id_dict_path):
            with open(ent2id_dict_path, "wb") as f:
                pickle.dump(dataset.entity_vocab, f)
            logger.info(
                f"Created ent2id_dict.pkl with {len(dataset.entity_vocab)} entities"
            )

        # Create id2rel_dict
        id2rel_dict_path = os.path.join(dataset_root, "id2rel_dict.pkl")
        if not os.path.exists(id2rel_dict_path):
            id2rel_dict = {v: k for k, v in dataset.relation_vocab.items()}
            with open(id2rel_dict_path, "wb") as f:
                pickle.dump(id2rel_dict, f)
            logger.info(f"Created id2rel_dict.pkl with {len(id2rel_dict)} relations")

        # Create rel2id_dict
        rel2id_dict_path = os.path.join(dataset_root, "rel2id_dict.pkl")
        if not os.path.exists(rel2id_dict_path):
            with open(rel2id_dict_path, "wb") as f:
                pickle.dump(dataset.relation_vocab, f)
            logger.info(
                f"Created rel2id_dict.pkl with {len(dataset.relation_vocab)} relations"
            )

        logger.info("Entity/relation dictionaries created successfully")
    else:
        logger.info("Entity/relation dictionaries already exist")

    # Check if ent2name_dict exists
    ent2name_path = os.path.join(dataset_root, "ent2name_dict.pkl")
    if not os.path.exists(ent2name_path):
        logger.info("Creating ent2name_dict from nodes file...")
        nodes_path = os.path.join(dataset_dir, "raw", "nodes.txt")

        if not os.path.exists(nodes_path):
            raise FileNotFoundError(
                f"Nodes file not found at {nodes_path}. "
                f"Please ensure the dataset is properly downloaded."
            )

        # Load nodes file and create mapping
        nodes_df = pl.read_csv(nodes_path, separator="\t", has_header=True)

        # Create source_label -> name mapping
        ent2name_dict = dict(zip(nodes_df["source_label"], nodes_df["name"]))

        # Save to pickle
        with open(ent2name_path, "wb") as f:
            pickle.dump(ent2name_dict, f)

        logger.info(f"Created ent2name_dict with {len(ent2name_dict)} entries")
    else:
        logger.info("ent2name_dict already exists")


class ModelManager:
    """Singleton manager for lazy loading and caching ULTRA model."""

    _instance = None
    _initialized = False

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(ModelManager, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        if not ModelManager._initialized:
            self.model = None
            self.dataset = None
            self.device = None
            self.id2ent_dict = None
            self.ent2id_dict = None
            self.id2rel_dict = None
            self.rel2id_dict = None
            self.ent2name_dict = None
            self.nodes_df = None
            self.filtered_data = None
            self.graph_schema = None  # Cache for graph schema
            # Performance optimization caches (built during model loading)
            self.entity_type_map = None  # source_label -> type (O(1) lookup)
            self.entity_id_list = None  # List of all entity IDs for suggestions
            self.entity_name_list = None  # List of all entity names for suggestions
            self.entity_name_to_id_map = None  # name -> source_label (O(1) lookup)
            self.relation_list = None  # List of all relations for validation
            ModelManager._initialized = True

    def load_model(self):
        """Load model, dataset, and dictionaries on first inference request."""
        if self.model is not None:
            logger.info("Model already loaded, using cached instance")
            return

        logger.info("Loading ULTRA model and dataset for the first time...")

        # Track initialization state for cleanup on error
        initialization_stage = "starting"

        try:
            initialization_stage = "building dataset"
            # Build dataset
            cfg = EasyDict(DEFAULT_CONFIG)
            self.dataset = util.build_dataset(cfg)
            logger.info(f"Loaded dataset: {self.dataset.name}")

            initialization_stage = "setting device"
            # Get device
            if DEFAULT_CONFIG["device"] == "auto":
                self.device = torch.device(
                    "cuda" if torch.cuda.is_available() else "cpu"
                )
            else:
                self.device = torch.device(DEFAULT_CONFIG["device"])
            logger.info(f"Using device: {self.device}")

            initialization_stage = "creating model"
            # Load model
            self.model = Ultra(
                rel_model_cfg=cfg.model.relation_model,
                entity_model_cfg=cfg.model.entity_model,
            )

            initialization_stage = "loading checkpoint"
            # Load checkpoint
            checkpoint_path = DEFAULT_CONFIG["checkpoint"]

            # Check if checkpoint exists
            if not os.path.exists(checkpoint_path):
                raise FileNotFoundError(
                    f"Checkpoint not found: {checkpoint_path}. "
                    f"Please ensure the checkpoint file is present or set ULTRA_CHECKPOINT_PATH."
                )

            # Load checkpoint with validation
            try:
                state = torch.load(
                    checkpoint_path, map_location="cpu", weights_only=False
                )
                if "model" not in state:
                    raise ValueError(
                        f"Invalid checkpoint format: 'model' key not found in {checkpoint_path}"
                    )
                self.model.load_state_dict(state["model"])
            except Exception as e:
                raise RuntimeError(
                    f"Failed to load checkpoint from {checkpoint_path}: {e}"
                )

            initialization_stage = "moving model to device"
            # Try to move to device with OOM handling
            try:
                self.model = self.model.to(self.device)
                self.model.eval()
                logger.info(f"Loaded checkpoint from: {checkpoint_path}")
            except RuntimeError as e:
                if "out of memory" in str(e).lower():
                    logger.warning("GPU out of memory, falling back to CPU")
                    self.device = torch.device("cpu")
                    self.model = self.model.to(self.device)
                    self.model.eval()
                else:
                    raise

            initialization_stage = "creating dataset dictionaries"
            # Ensure all required dictionary files exist
            logger.info("Ensuring dataset dictionaries exist...")
            ensure_dataset_dictionaries(
                dataset_root=DEFAULT_CONFIG["dataset"]["root"],
                dataset_name=self.dataset.name,
                dataset=self.dataset,  # Pass the dataset to avoid reloading
            )

            initialization_stage = "loading ID dictionaries"
            # Load ID dictionaries from the created pickle files using helper function
            dataset_root = DEFAULT_CONFIG["dataset"]["root"]
            logger.info(f"Loading dictionaries from {dataset_root}")

            self.id2ent_dict = load_pickle(os.path.join(dataset_root, "id2ent_dict.pkl"))
            self.ent2id_dict = load_pickle(os.path.join(dataset_root, "ent2id_dict.pkl"))
            self.id2rel_dict = load_pickle(os.path.join(dataset_root, "id2rel_dict.pkl"))
            self.rel2id_dict = load_pickle(os.path.join(dataset_root, "rel2id_dict.pkl"))

            logger.info(
                f"Loaded dictionaries: {len(self.id2ent_dict)} entities, {len(self.id2rel_dict)} relations"
            )

            initialization_stage = "loading ent2name dictionary"
            # Load ent2name dictionary (located in data root, not dataset subdirectory)
            self.ent2name_dict = load_pickle(
                os.path.join(dataset_root, "ent2name_dict.pkl")
            )

            initialization_stage = "loading nodes dataframe"
            # Load nodes for validation
            dataset_dir = os.path.join(
                DEFAULT_CONFIG["dataset"]["root"], self.dataset.name
            )
            nodes_path = os.path.join(dataset_dir, "raw", "nodes.txt")
            # The actual format is: source_id, name, type, source, source_label
            self.nodes_df = pl.read_csv(
                nodes_path,
                separator="\t",
                has_header=True,
                schema={
                    "source_id": pl.Utf8,
                    "name": pl.Utf8,
                    "type": pl.Utf8,
                    "source": pl.Utf8,
                    "source_label": pl.Utf8,
                },
            )
            logger.info(f"Loaded nodes dataframe with {len(self.nodes_df)} nodes")

            initialization_stage = "building lookup caches"
            # Build lookup maps for O(1) access (performance optimization)
            logger.info("Building entity and relation lookup caches...")
            self.entity_type_map = dict(
                zip(self.nodes_df["source_label"], self.nodes_df["type"])
            )
            self.entity_id_list = self.nodes_df["source_label"].to_list()
            self.entity_name_list = self.nodes_df["name"].to_list()
            self.entity_name_to_id_map = dict(
                zip(self.nodes_df["name"], self.nodes_df["source_label"])
            )
            self.relation_list = list(self.rel2id_dict.keys())
            logger.info(
                f"Built lookup caches: {len(self.entity_type_map)} entities, "
                f"{len(self.relation_list)} relations"
            )

            initialization_stage = "creating filtered data"
            # Create filtered data for strict negative masking
            self.filtered_data = Data(
                edge_index=self.dataset._data.target_edge_index,
                edge_type=self.dataset._data.target_edge_type,
                num_nodes=self.dataset[0].num_nodes,
            ).to(self.device)

            logger.info("Model loading complete!")

        except Exception as e:
            logger.error(f"Failed to load model at stage '{initialization_stage}': {e}")
            import traceback

            traceback.print_exc()
            # Reset state on error so next attempt can retry
            self.model = None
            self.dataset = None
            self.nodes_df = None
            raise

    def get_model(self):
        """Get model, loading if necessary."""
        if self.model is None:
            self.load_model()
        return self.model


def validate_entity(entity: str, manager: ModelManager) -> tuple[str, str]:
    """
    Validate entity and return (entity_id, entity_name).

    Uses cached lookups for O(1) performance instead of O(n) DataFrame filters.

    Args:
        entity: Entity ID (e.g., "MONDO:5301") or name (e.g., "Crohn disease")
        manager: ModelManager instance with loaded nodes data

    Returns:
        Tuple of (entity_id, entity_name)

    Raises:
        ValueError: If entity not found, with suggestions
    """

    # Check if entity looks like an ID (contains colon)
    if ":" in entity:
        # Fast O(1) lookup using ent2name_dict
        entity_name = manager.ent2name_dict.get(entity)
        if entity_name:
            return entity, entity_name

        # Not found - suggest similar IDs using cached list
        suggestions = get_close_matches(entity, manager.entity_id_list, n=3, cutoff=0.6)

        if suggestions:
            suggestion_str = ", ".join(
                [f"{s} ({manager.ent2name_dict[s]})" for s in suggestions]
            )
            raise ValueError(
                f"Entity ID '{entity}' not found in PrimeKG. Did you mean: {suggestion_str}"
            )
        else:
            raise ValueError(
                f"Entity ID '{entity}' not found in PrimeKG and no similar IDs found."
            )
    else:
        # Fast O(1) lookup using entity_name_to_id_map
        entity_id = manager.entity_name_to_id_map.get(entity)
        if entity_id:
            return entity_id, entity

        # Case-insensitive search (still needs filter but only if exact match fails)
        result = manager.nodes_df.filter(
            pl.col("name").str.to_lowercase() == entity.lower()
        )
        if len(result) > 0:
            return result["source_label"][0], result["name"][0]

        # Not found - suggest similar names using cached list
        suggestions = get_close_matches(entity, manager.entity_name_list, n=3, cutoff=0.6)

        if suggestions:
            suggestion_str = ", ".join(
                [f"{manager.entity_name_to_id_map[s]} ({s})" for s in suggestions]
            )
            raise ValueError(
                f"Entity '{entity}' not found in PrimeKG. Did you mean: {suggestion_str}"
            )
        else:
            raise ValueError(
                f"Entity '{entity}' not found in PrimeKG and no similar entities found."
            )


def validate_relation(relation: str, manager: ModelManager) -> str:
    """
    Validate and normalize relation.

    Uses cached relation list for better performance.

    Args:
        relation: Relation label (e.g., "associated with", "associated_with")
        manager: ModelManager instance with loaded relation data

    Returns:
        Normalized relation string

    Raises:
        ValueError: If relation not found
    """
    # Use cached relation list for performance
    valid_relations = manager.relation_list

    # Try exact match
    if relation in valid_relations:
        return relation

    # Try with underscores replaced by spaces (e.g., "associated_with" -> "associated with")
    with_spaces = relation.replace("_", " ")
    if with_spaces in valid_relations:
        return with_spaces

    # Try lowercase version
    lower_relation = relation.lower()
    if lower_relation in valid_relations:
        return lower_relation

    # Try lowercase with underscores replaced by spaces
    normalized = relation.lower().replace("_", " ")
    if normalized in valid_relations:
        return normalized

    # Try with spaces replaced by underscores (for relations that use underscores)
    with_underscores = relation.replace(" ", "_")
    if with_underscores in valid_relations:
        return with_underscores

    # Not found - show available relations
    relation_list = ", ".join(valid_relations[:10])
    if len(valid_relations) > 10:
        relation_list += f", ... ({len(valid_relations)} total)"

    raise ValueError(
        f"Relation '{relation}' not found. Available relations include: {relation_list}"
    )


def load_graph_schema(manager: ModelManager) -> dict:
    """
    Load graph schema from static schema file.

    The schema file (primekg_schema.json) contains the complete graph schema
    built from all edges in PrimeKG. This provides fast validation of
    (head_type, relation) -> tail_type combinations.

    Returns:
        Dictionary mapping (head_type, relation) -> set of valid tail_types
    """
    # Return cached schema if available
    if manager.graph_schema is not None:
        return manager.graph_schema

    # Load schema from JSON file in data directory
    schema_path = os.path.join(MCP_ROOT, "data", "primekg_schema.json")

    if not os.path.exists(schema_path):
        logger.error(f"Schema file not found at {schema_path}")
        logger.error(
            "Please ensure primekg_schema.json exists in the data directory"
        )
        raise FileNotFoundError(
            f"Schema file not found at {schema_path}. "
            "Please ensure primekg_schema.json exists in the data directory."
        )

    logger.info(f"Loading graph schema from {schema_path}...")

    with open(schema_path, "r") as f:
        schema_data = json.load(f)

    # Convert schema from "head_type|relation" -> [tail_types] format
    # to (head_type, relation) -> {tail_types} format
    schema = {}
    for key, tail_types in schema_data["schema"].items():
        head_type, relation = key.split("|", 1)
        schema[(head_type, relation)] = set(tail_types)

    logger.info(f"Loaded schema with {len(schema)} (head_type, relation) pairs")

    # Cache the schema
    manager.graph_schema = schema

    return schema


def validate_query_schema(
    head_entity_id: str,
    relation: str,
    manager: ModelManager,
) -> set:
    """
    Validate that (head_type, relation) exists in PrimeKG schema.

    This function performs upfront validation to reject invalid queries
    before running inference, saving compute time and providing clear
    error messages to users.

    Uses cached lookups for O(1) performance.

    Args:
        head_entity_id: Entity ID (e.g., "MONDO:5301")
        relation: Relation label (e.g., "associated with")
        manager: ModelManager instance

    Returns:
        Set of valid tail entity types for this (head_type, relation) combination

    Raises:
        ValueError: If (head_type, relation) does not exist in PrimeKG schema
    """
    # Load schema
    schema = load_graph_schema(manager)

    # Get head entity type using O(1) cached lookup instead of O(n) filter
    head_type = manager.entity_type_map.get(head_entity_id)
    if not head_type:
        raise ValueError(f"Entity {head_entity_id} not found in nodes")

    # Get head entity name using O(1) cached lookup
    head_name = manager.ent2name_dict.get(head_entity_id, head_entity_id)

    # Check if (head_type, relation) exists in schema
    valid_tail_types = schema.get((head_type, relation))

    if not valid_tail_types:
        # This combination doesn't exist in PrimeKG
        # Provide a helpful error message
        raise ValueError(
            f"Invalid query: ({head_type}, {relation}) combination does not exist in PrimeKG. "
            f"Entity '{head_name}' ({head_entity_id}) has type '{head_type}', "
            f"but there are no edges with relation '{relation}' from this entity type in the graph."
        )

    logger.info(
        f"Schema validation passed: ({head_type}, {relation}) -> {valid_tail_types}"
    )

    return valid_tail_types


def create_single_query_data(dataset, head_entity_id: str, relation: str) -> Data:
    """
    Create a single-query Data object for zero-shot inference (not evaluation).

    Unlike util.inference_data_single() which extracts all existing edges for
    evaluation, this creates just one dummy query to predict new associations.

    Args:
        dataset: PyG dataset
        head_entity_id: Entity ID (e.g., "NCBI:7297")
        relation: Relation label (e.g., "associated with")

    Returns:
        PyG Data object with single target edge
    """
    # Translate entity and relation to integer indices
    trans_h_ent = util.translate_entity(dataset, head_entity_id)
    trans_t_ent_dummy = util.translate_entity(
        dataset, "NCBI:5340"
    )  # Plasminogen gene as dummy
    trans_rel = util.translate_relation(dataset, relation)

    # Create single dummy query
    target_edge_index = torch.tensor([[trans_h_ent], [trans_t_ent_dummy]])
    target_edge_type = torch.tensor([trans_rel])

    # Create Data object with full graph and single query
    return Data(
        edge_index=dataset[0].edge_index,
        edge_type=dataset[0].edge_type,
        target_edge_index=target_edge_index,
        target_edge_type=target_edge_type,
        num_relations=dataset[0].num_relations,
        num_nodes=dataset[0].num_nodes,
        relation_graph=dataset[0].relation_graph,
    )


@torch.no_grad()
def run_inference(
    head_entity_id: str,
    relation: str,
    manager: ModelManager,
) -> dict:
    """
    Run ULTRA model inference for a single (head, relation) query.

    Args:
        head_entity_id: Entity ID (e.g., "MONDO:5301")
        relation: Relation label (e.g., "associated with")
        manager: ModelManager instance

    Returns:
        Dictionary with raw predictions: h, r, t, scores, rankings, masks
    """
    logger.info(f"Running inference for ({head_entity_id}, {relation})")

    # Get model and data
    model = manager.get_model()
    dataset = manager.dataset
    device = manager.device

    # Create single-query inference data (not multi-edge evaluation data)
    infer_data = create_single_query_data(dataset, head_entity_id, relation)
    infer_data = infer_data.to(device)

    # Create batch (single query)
    test_triplets = torch.cat(
        [infer_data.target_edge_index, infer_data.target_edge_type.unsqueeze(0)]
    ).t()

    # Run inference
    t_batch, h_batch = tasks.all_negative(infer_data, test_triplets)
    t_pred = model(infer_data, t_batch)
    logger.info(
        f"t_pred shape: {t_pred.shape}, test_triplets shape: {test_triplets.shape}"
    )

    # Apply strict negative mask
    t_mask, h_mask = tasks.strict_negative_mask(manager.filtered_data, test_triplets)

    # Get rankings
    pos_h_index, pos_t_index, pos_r_index = test_triplets.t()
    t_ranking = tasks.compute_ranking(t_pred, pos_t_index, t_mask)
    t_unfilt_ranking = tasks.compute_ranking(t_pred, pos_t_index)

    # Get top-k predictions (filtered)
    t_predictions_filt = tasks.get_predictions(t_pred, t_mask)

    # Get unfiltered predictions sorted by score (all entities)
    # tasks.get_predictions returns entity IDs sorted by their prediction scores
    t_predictions_unfilt = tasks.get_predictions(pred=t_pred, top=None)

    # Return raw results
    return {
        "h": test_triplets[:, 0].cpu().tolist(),
        "r": test_triplets[:, 2].cpu().tolist(),
        "t": test_triplets[:, 1].cpu().tolist(),
        "t_filt_rank": t_ranking.cpu().tolist(),
        "t_unfilt_rank": t_unfilt_ranking.cpu().tolist(),
        "t_pred_filt": t_predictions_filt.cpu().tolist(),
        "t_pred_unfilt": t_predictions_unfilt.cpu().tolist(),
        "t_pred_score": t_pred.cpu().tolist(),
        "t_mask": t_mask.cpu().tolist(),
    }


def transform_predictions(
    raw_results: dict,
    manager: ModelManager,
    data_root: str,
    top_k: Optional[int] = None,
) -> "pl.DataFrame":
    """
    Transform raw inference results using data_util pipeline.

    Args:
        raw_results: Raw predictions from run_inference
        manager: ModelManager instance
        data_root: Path to dataset root directory (e.g., /path/to/data/primekg1)
        top_k: Number of top predictions to keep (None for all predictions)

    Returns:
        DataFrame with transformed predictions
    """
    logger.info("Transforming predictions...")

    # Create initial DataFrame
    df = pl.DataFrame(raw_results).unique()
    logger.info(f"Initial DataFrame shape after .unique(): {df.shape}")

    # Get parent directory where pickle files are located
    # If data_root is /path/to/data/primekg1, parent is /path/to/data
    parent_dir = os.path.dirname(data_root)

    # Step 1: Translate IDs to labels and names
    logger.info("Step 1: Translating IDs to labels and names")
    df = data_util.translate_hrt(df=df, data_path=parent_dir)
    logger.info(f"After translate_hrt: {df.shape}")

    # Step 2: Filter and process results (score sorting, novelty marking)
    logger.info("Step 2: Filtering and processing results")
    df = data_util.filter_process_results(df=df, results_path=parent_dir)
    logger.info(f"After filter_process_results: {df.shape}")

    # Step 3: Structure results (schema validation, type filtering)
    logger.info("Step 3: Structuring results with schema validation")
    node_path = os.path.join(data_root, "raw", "nodes.txt")
    graph_path = os.path.join(data_root, "raw")
    df = data_util.structure_results(df=df, node_path=node_path, graph_path=graph_path)
    logger.info(f"After structure_results: {len(df)} predictions")

    # Step 4: Add rank column and explode
    logger.info("Step 4: Adding rank column")
    df = df.with_columns(
        pl.int_ranges(1, pl.col("edge_in_primekg").list.len() + 1).alias("rank"),
        pl.col("r_label").str.replace(" ", "_"),
    ).explode(
        [
            "t_pred_label",
            "t_pred_name",
            "t_pred_score",
            "t_pred_type",
            "edge_in_primekg",
            "rank",
        ]
    )

    # Check if we have any predictions after filtering
    if len(df) == 0:
        logger.warning("No valid predictions after schema filtering")
        # Return empty dataframe with correct schema
        return pl.DataFrame(
            {
                "h_label": [],
                "h_name": [],
                "h_type": [],
                "r_label": [],
                "t_pred_label": [],
                "t_pred_name": [],
                "t_pred_score": [],
                "t_pred_type": [],
                "edge_in_primekg": [],
                "rank": [],
                "percentile_rank": [],
            }
        )

    # Step 5: Add percentile_rank column
    logger.info("Step 5: Adding percentile rank")
    total_predictions = len(df)
    df = df.with_columns(
        (1.0 - pl.col("rank") / total_predictions).alias("percentile_rank")
    )

    logger.info(f"Transformation complete: {total_predictions} predictions")
    # Step 6: Keep only top_k if specified
    if top_k is not None:
        df = df.filter(pl.col("rank") <= top_k)
        logger.info(f"Kept top {top_k} predictions")

    return df


def _predict_tail_entities_impl(
    head_entity: str, relation: str, top_k: Optional[int] = None
) -> dict:
    """
    Implementation of predict_tail_entities logic.

    Args:
        head_entity: Entity ID (e.g., "MONDO:5301") or name (e.g., "Crohn disease")
        relation: Relation label (e.g., "associated with", "associated_with")
        top_k: Number of top predictions to return (default: None for all predictions)

    Returns:
        Dictionary with:
        - success: bool
        - head_entity: str (normalized ID)
        - head_name: str
        - relation: str (normalized)
        - output_file: str (path to parquet file)
        - total_predictions: int
        - inference_time_seconds: float
        - preview: str (markdown table of top 10 results)
    """

    start_time = time.time()

    try:
        logger.info(
            f"Received request: head='{head_entity}', relation='{relation}', top_k={top_k}"
        )

        # Get model manager
        manager = ModelManager()
        manager.load_model()

        # Validate inputs
        logger.info("Validating inputs...")
        head_entity_id, head_name = validate_entity(head_entity, manager)
        relation_normalized = validate_relation(relation, manager)

        logger.info(
            f"Validated: {head_entity_id} ({head_name}) - {relation_normalized}"
        )

        # Validate query against schema BEFORE running inference
        logger.info("Validating query against graph schema...")
        expected_tail_types = validate_query_schema(
            head_entity_id=head_entity_id,
            relation=relation_normalized,
            manager=manager,
        )
        logger.info(
            f"Schema validation passed. Expected tail types: {expected_tail_types}"
        )

        # Run inference
        logger.info("Running inference...")
        raw_results = run_inference(
            head_entity_id=head_entity_id, relation=relation_normalized, manager=manager
        )

        # Transform predictions
        data_root = os.path.join(
            DEFAULT_CONFIG["dataset"]["root"], manager.dataset.name
        )
        df = transform_predictions(raw_results, manager, data_root, top_k=top_k)

        # Handle empty predictions
        if len(df) == 0:
            logger.warning(
                f"No valid predictions for ({head_entity_id}, {relation_normalized})"
            )
            return {
                "success": True,
                "head_entity": head_entity_id,
                "head_name": head_name,
                "relation": relation_normalized,
                "output_file": None,
                "total_predictions": 0,
                "inference_time_seconds": round(time.time() - start_time, 2),
                "preview": "No valid predictions found. This may indicate an unusual query or schema mismatch.",
                "warning": "No predictions after schema filtering",
            }

        # Create output directory
        output_dir = os.path.join(
            DEFAULT_CONFIG["output_dir"],
            head_entity_id.replace(":", "_"),
            relation_normalized,
        )
        os.makedirs(output_dir, exist_ok=True)

        # Write to parquet
        output_file = os.path.join(output_dir, "predictions.parquet")
        df.write_parquet(output_file)
        logger.info(f"Results written to: {output_file}")

        # Create preview (top 10 rows as markdown)
        preview_df = df.head(10)
        # Convert to string representation instead of markdown to avoid pyarrow dependency
        preview_markdown = str(preview_df)

        # Calculate inference time
        inference_time = time.time() - start_time

        # Return success response
        return {
            "success": True,
            "head_entity": head_entity_id,
            "head_name": head_name,
            "relation": relation_normalized,
            "output_file": output_file,
            "total_predictions": len(df),
            "inference_time_seconds": round(inference_time, 2),
            "preview": preview_markdown,
        }

    except ValueError as e:
        # Validation errors
        logger.error(f"Validation error: {e}")
        return {"success": False, "error": str(e), "error_type": "ValidationError"}
    except Exception as e:
        # Other errors
        logger.error(f"Inference error: {e}")
        import traceback

        traceback.print_exc()
        return {"success": False, "error": str(e), "error_type": type(e).__name__}


@mcp.tool
def predict_tail_entities(
    head_entity: str, relation: str, top_k: Optional[int] = None
) -> dict:
    """
    Predict tail entities for a given head entity and relation using ULTRA model.

    Args:
        head_entity: Entity ID (e.g., "MONDO:5301") or name (e.g., "Crohn disease")
        relation: Relation label (e.g., "associated with", "associated_with")
        top_k: Number of top predictions to return (default: None for all predictions)

    Returns:
        Dictionary with:
        - success: bool
        - head_entity: str (normalized ID)
        - head_name: str
        - relation: str (normalized)
        - output_file: str (path to parquet file)
        - total_predictions: int
        - inference_time_seconds: float
        - preview: str (markdown table of top 10 results)
    """
    return _predict_tail_entities_impl(head_entity, relation, top_k)


@mcp.tool
def setup_primekg_data(
    dataset_name: str = "PrimeKG1",
    force_redownload: bool = False,
    train_frac: float = 0.8,
    test_frac: float = 0.1,
    valid_frac: float = 0.1,
    seed: int = 42,
) -> dict:
    """
    Download and setup PrimeKG data with train/test/valid splits.

    This function will:
    1. Check if PrimeKG data already exists
    2. Download primekg.csv from Harvard Dataverse if needed
    3. Process the data into train/test/valid splits
    4. Create nodes file with entity metadata

    Args:
        dataset_name: Name of dataset (e.g., "PrimeKG1", "PrimeKG2")
        force_redownload: If True, re-download even if files exist
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
        - dataset_path: str (path to dataset directory)
    """
    try:
        # Construct dataset path
        data_root = DEFAULT_CONFIG["dataset"]["root"]
        dataset_path = os.path.join(data_root, dataset_name.lower())

        logger.info(f"Setting up PrimeKG data for {dataset_name} at {dataset_path}")

        # Run setup
        result = primekg_setup.setup_primekg(
            dataset_path=dataset_path,
            force_redownload=force_redownload,
            train_frac=train_frac,
            test_frac=test_frac,
            valid_frac=valid_frac,
            seed=seed,
        )

        return {"success": True, **result}

    except Exception as e:
        logger.error(f"Failed to setup PrimeKG: {e}")
        import traceback

        traceback.print_exc()
        return {"success": False, "error": str(e), "error_type": type(e).__name__}


if __name__ == "__main__":
    logger.info("Starting ULTRA Inference MCP Server...")
    mcp.run()
