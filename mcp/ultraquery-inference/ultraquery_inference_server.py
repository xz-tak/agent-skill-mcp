#!/usr/bin/env python3
"""
UltraQuery Inference MCP Server

On-demand inference for complex logical queries on PrimeKG using UltraQuery.
Supports multi-hop reasoning with intersection, union, and negation operations.
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
from typing import Optional, Union, List, Tuple, Any
import polars as pl
import time
from difflib import get_close_matches
from ultra.models import Ultra
from ultra.ultraquery import UltraQuery
from ultra.query_utils import Query
from ultra import util, tasks, datasets

# Import PrimeKG setup utilities
from src.primekg_setup import check_primekg_data, setup_primekg

# Initialize logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Initialize MCP server
mcp = FastMCP("ultraquery_inference", stateless_http=True)

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
            "class": "QueryNBFNet",
            "input_dim": 64,
            "hidden_dims": [64, 64, 64, 64, 64, 64],
            "message_func": "distmult",
            "aggregate_func": "sum",
            "short_cut": True,
            "layer_norm": True,
        },
    },
    "ultra_checkpoint": os.getenv(
        "ULTRA_CHECKPOINT_PATH",
        os.path.join(ULTRA_ROOT, "ckpts/ultra_primekg_50g_ft_epoch_1.pth"),
    ),
    "ultraquery_checkpoint": os.getenv(
        "ULTRAQUERY_CHECKPOINT_PATH",
        os.path.join(ULTRA_ROOT, "ckpts/ultraquery_primekg_ft_epoch_1.pth"),
    ),
    "logic": "product",  # product, godel, lukasiewicz
    "dropout_ratio": 0.25,  # Dropout ratio from config
    "threshold": 0.0,
    "more_dropout": 0.0,
    "output_dir": os.getenv("ULTRA_OUTPUT_DIR", os.path.join(MCP_ROOT, "output")),
    "device": "auto",
}

# Query type mapping (from datasets_query.py)
STRUCT2TYPE = {
    ("e", ("r",)): "1p",
    ("e", ("r", "r")): "2p",
    ("e", ("r", "r", "r")): "3p",
    (("e", ("r",)), ("e", ("r",))): "2i",
    (("e", ("r",)), ("e", ("r",)), ("e", ("r",))): "3i",
    ((("e", ("r",)), ("e", ("r",))), ("r",)): "ip",
    (("e", ("r", "r")), ("e", ("r",))): "pi",
    (("e", ("r",)), ("e", ("r", "n"))): "2in",
    (("e", ("r",)), ("e", ("r",)), ("e", ("r", "n"))): "3in",
    ((("e", ("r",)), ("e", ("r", "n"))), ("r",)): "inp",
    (("e", ("r", "r")), ("e", ("r", "n"))): "pin",
    (("e", ("r", "r", "n")), ("e", ("r",))): "pni",
    (("e", ("r",)), ("e", ("r",)), ("u",)): "2u-DNF",
    ((("e", ("r",)), ("e", ("r",)), ("u",)), ("r",)): "up-DNF",
    ((("e", ("r", "n")), ("e", ("r", "n"))), ("n",)): "2u-DM",
    ((("e", ("r", "n")), ("e", ("r", "n"))), ("n", "r")): "up-DM",
}


def ensure_dataset_dictionaries(
    dataset_root: str, dataset_name: str, dataset=None
) -> None:
    """
    Ensure all required pickle files exist for dataset translation.

    This function creates:
    - id2ent_dict.pkl, ent2id_dict.pkl, id2rel_dict.pkl, rel2id_dict.pkl in dataset_root/
    - ent2name_dict.pkl in dataset_root/

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
    """Singleton manager for lazy loading and caching UltraQuery model."""

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
            self.graph_schema = None  # Cache for graph schema
            # Additional caches for performance optimization
            self.entity_type_map = None  # source_label -> type
            self.entity_id_list = None  # List of all entity IDs for suggestions
            self.entity_name_list = None  # List of all entity names for suggestions
            self.entity_name_to_id_map = None  # name -> source_label for fast lookup
            self.relation_list = None  # List of all relations for validation
            ModelManager._initialized = True

    def load_model(self):
        """Load UltraQuery model, dataset, and dictionaries on first inference request."""
        if self.model is not None:
            logger.info("Model already loaded, using cached instance")
            return

        logger.info("Loading UltraQuery model and dataset for the first time...")

        # Track initialization stage for cleanup on error
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

            # Warn if using CPU (inference will be slow)
            if self.device.type == "cpu":
                logger.warning("=" * 80)
                logger.warning("WARNING: Using CPU for inference")
                logger.warning("Inference will take 10-15 minutes per query on CPU")
                logger.warning(
                    "For faster inference (1-3 minutes), use a GPU-enabled machine"
                )
                logger.warning("=" * 80)

            initialization_stage = "creating base Ultra model"
            # Create base Ultra model
            base_model = Ultra(
                rel_model_cfg=cfg.model.relation_model,
                entity_model_cfg=cfg.model.entity_model,
            )

            initialization_stage = "creating UltraQuery model"
            # Wrap in UltraQuery
            self.model = UltraQuery(
                model=base_model,
                logic=DEFAULT_CONFIG["logic"],
                dropout_ratio=DEFAULT_CONFIG["dropout_ratio"],
                threshold=DEFAULT_CONFIG["threshold"],
                more_dropout=DEFAULT_CONFIG["more_dropout"],
            )

            initialization_stage = "loading Ultra checkpoint"
            # Load base Ultra checkpoint if available
            ultra_checkpoint_path = DEFAULT_CONFIG["ultra_checkpoint"]
            if ultra_checkpoint_path and os.path.exists(ultra_checkpoint_path):
                try:
                    state = torch.load(
                        ultra_checkpoint_path, map_location="cpu", weights_only=False
                    )
                    if "model" not in state:
                        raise ValueError(
                            f"Invalid checkpoint format: 'model' key not found in {ultra_checkpoint_path}"
                        )
                    # Load into base model (model.model in UltraQuery)
                    self.model.model.model.load_state_dict(state["model"])
                    logger.info(
                        f"Loaded Ultra checkpoint from: {ultra_checkpoint_path}"
                    )
                except Exception as e:
                    logger.warning(
                        f"Could not load Ultra checkpoint from {ultra_checkpoint_path}: {e}. Continuing without it."
                    )

            initialization_stage = "loading UltraQuery checkpoint"
            # Load UltraQuery checkpoint if available
            ultraquery_checkpoint_path = DEFAULT_CONFIG["ultraquery_checkpoint"]
            if ultraquery_checkpoint_path and os.path.exists(
                ultraquery_checkpoint_path
            ):
                try:
                    state = torch.load(
                        ultraquery_checkpoint_path,
                        map_location="cpu",
                        weights_only=False,
                    )
                    if "model" not in state:
                        raise ValueError(
                            f"Invalid checkpoint format: 'model' key not found in {ultraquery_checkpoint_path}"
                        )
                    self.model.load_state_dict(state["model"])
                    logger.info(
                        f"Loaded UltraQuery checkpoint from: {ultraquery_checkpoint_path}"
                    )
                except Exception as e:
                    logger.warning(
                        f"Could not load UltraQuery checkpoint from {ultraquery_checkpoint_path}: {e}. Using base Ultra model only."
                    )

            initialization_stage = "moving model to device"
            # Try to move to device with OOM handling
            try:
                self.model = self.model.to(self.device)
                self.model.eval()
                logger.info(f"Model moved to device: {self.device}")
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
            # Load ID dictionaries from the created pickle files
            dataset_root = DEFAULT_CONFIG["dataset"]["root"]
            logger.info(f"Loading dictionaries from {dataset_root}")

            with open(os.path.join(dataset_root, "id2ent_dict.pkl"), "rb") as f:
                self.id2ent_dict = pickle.load(f)

            with open(os.path.join(dataset_root, "ent2id_dict.pkl"), "rb") as f:
                self.ent2id_dict = pickle.load(f)

            with open(os.path.join(dataset_root, "id2rel_dict.pkl"), "rb") as f:
                self.id2rel_dict = pickle.load(f)

            with open(os.path.join(dataset_root, "rel2id_dict.pkl"), "rb") as f:
                self.rel2id_dict = pickle.load(f)

            logger.info(
                f"Loaded dictionaries: {len(self.id2ent_dict)} entities, {len(self.id2rel_dict)} relations"
            )

            initialization_stage = "loading ent2name dictionary"
            # Load ent2name dictionary (located in data root, not dataset subdirectory)
            ent2name_path = os.path.join(
                DEFAULT_CONFIG["dataset"]["root"], "ent2name_dict.pkl"
            )
            with open(ent2name_path, "rb") as f:
                self.ent2name_dict = pickle.load(f)

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

            logger.info("UltraQuery model loading complete!")

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
        # Fast lookup using ent2name_dict (O(1) instead of O(n) filter)
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
        # Fast lookup using entity_name_to_id_map (O(1) instead of O(n) filter)
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

    # Try with underscores replaced by spaces
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

    # Try with spaces replaced by underscores
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


def parse_query_structure(
    query_structure: Any, manager: ModelManager
) -> Tuple[Any, dict, dict]:
    """
    Parse query structure and validate entities/relations.

    Converts entity names/IDs to integer indices and relation labels to integer indices.

    Args:
        query_structure: Nested tuple/list in BetaE format (e.g., ["e", ["r"]] for 1p)
                        Entities can be IDs or names, relations can be labels
                        Lists are automatically converted to tuples (for MCP JSON compatibility)
        manager: ModelManager instance

    Returns:
        Tuple of (parsed_structure_with_indices, entity_map, relation_map)
        - parsed_structure_with_indices: Nested tuple with integer indices
        - entity_map: Dict mapping index -> (entity_id, entity_name)
        - relation_map: Dict mapping index -> relation_label

    Raises:
        ValueError: If structure is invalid or entities/relations not found
    """
    entity_map = {}
    relation_map = {}
    entity_counter = 0

    def parse_recursive(struct):
        nonlocal entity_counter

        if isinstance(struct, str):
            # This is an entity
            entity_id, entity_name = validate_entity(struct, manager)
            # Get integer index from dataset
            entity_idx = manager.ent2id_dict[entity_id]
            entity_map[entity_idx] = (entity_id, entity_name)
            return entity_idx

        elif isinstance(struct, int):
            # Already an integer (relation index or entity index)
            # Check if it's a valid entity or relation index
            if struct in manager.id2ent_dict:
                entity_id = manager.id2ent_dict[struct]
                entity_name = manager.ent2name_dict.get(entity_id, entity_id)
                entity_map[struct] = (entity_id, entity_name)
                return struct
            elif struct in manager.id2rel_dict:
                relation_label = manager.id2rel_dict[struct]
                relation_map[struct] = relation_label
                return struct
            else:
                raise ValueError(
                    f"Integer {struct} is neither a valid entity nor relation index"
                )

        elif isinstance(struct, (tuple, list)):
            # Accept both tuples and lists (MCP sends JSON which converts tuples to lists)
            if len(struct) == 0:
                raise ValueError("Empty tuple/list in query structure")

            # Check if this is a relation tuple/list
            if all(isinstance(x, (str, int)) or x == "n" for x in struct):
                # This is a tuple/list of relations/operations
                parsed_ops = []
                for op in struct:
                    if op == "n":
                        # Negation operator
                        parsed_ops.append(-2)
                    elif op == "u":
                        # Union operator (kept as-is)
                        parsed_ops.append(op)
                    elif isinstance(op, str):
                        # Relation label
                        relation_label = validate_relation(op, manager)
                        relation_idx = manager.rel2id_dict[relation_label]
                        relation_map[relation_idx] = relation_label
                        parsed_ops.append(relation_idx)
                    elif isinstance(op, int):
                        # Already an integer relation index
                        if op == -2:
                            # Negation
                            parsed_ops.append(op)
                        elif op in manager.id2rel_dict:
                            relation_label = manager.id2rel_dict[op]
                            relation_map[op] = relation_label
                            parsed_ops.append(op)
                        else:
                            raise ValueError(f"Invalid relation index: {op}")
                    else:
                        raise ValueError(
                            f"Invalid operation in relation tuple/list: {op}"
                        )
                return tuple(parsed_ops)
            else:
                # This is a nested structure
                return tuple(parse_recursive(x) for x in struct)

        else:
            raise ValueError(f"Invalid type in query structure: {type(struct)}")

    parsed_structure = parse_recursive(query_structure)
    return parsed_structure, entity_map, relation_map


def identify_query_type(query_structure: tuple) -> str:
    """
    Identify the query type from the structure.

    Args:
        query_structure: Nested tuple structure

    Returns:
        Query type string (e.g., "1p", "2i", "ip")

    Raises:
        ValueError: If query type cannot be identified
    """

    def normalize_structure(struct):
        """Normalize structure to match STRUCT2TYPE keys."""
        if isinstance(struct, int):
            return "e"
        elif isinstance(struct, tuple):
            if len(struct) == 0:
                raise ValueError("Empty tuple in query structure")

            # Check if all elements are integers or "n" or "u"
            if all(
                isinstance(x, (int, str)) and (isinstance(x, int) or x in ["n", "u"])
                for x in struct
            ):
                # This is a relation tuple
                normalized = []
                for x in struct:
                    if x == -2:
                        normalized.append("n")
                    elif x == "u":
                        normalized.append("u")
                    elif isinstance(x, int):
                        normalized.append("r")
                    else:
                        normalized.append(x)
                return tuple(normalized)
            else:
                # Nested structure
                return tuple(normalize_structure(x) for x in struct)
        else:
            return struct

    normalized = normalize_structure(query_structure)

    # Look up in STRUCT2TYPE
    if normalized in STRUCT2TYPE:
        return STRUCT2TYPE[normalized]
    else:
        raise ValueError(
            f"Unknown query type. Normalized structure: {normalized}. "
            f"Supported types: {list(STRUCT2TYPE.values())}"
        )


def load_graph_schema(manager: ModelManager) -> dict:
    """
    Load graph schema from static schema file.

    The schema file (primekg_schema.json) contains the complete graph schema
    built from all edges in PrimeKG. This is much faster than computing from
    edges at inference time.

    Returns:
        Dictionary mapping [h_type, relation] -> set of valid t_types
    """
    # Return cached schema if available
    if manager.graph_schema is not None:
        return manager.graph_schema

    # Load schema from JSON file in data directory
    schema_path = os.path.join(MCP_ROOT, "data", "primekg_schema.json")

    if not os.path.exists(schema_path):
        logger.error(f"Schema file not found at {schema_path}")
        logger.error(
            "Please run 'pixi run python src/build_schema.py' to generate the schema file"
        )
        raise FileNotFoundError(
            f"Schema file not found at {schema_path}. "
            "Run 'pixi run python src/build_schema.py' to generate it."
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


class SchemaHelper:
    """
    Helper class for schema operations to avoid code duplication.
    Consolidates logic for entity type lookup, path following, and validation.
    """

    def __init__(self, manager: ModelManager):
        self.manager = manager
        self.schema = load_graph_schema(manager)
        self.nodes_df = manager.nodes_df

    def get_entity_type(self, entity_id: str) -> str:
        """Get entity type from entity ID."""
        row = self.nodes_df.filter(pl.col("source_label") == entity_id)
        return row["type"][0] if len(row) > 0 else None

    def get_tail_types(self, h_type: str, relation: str) -> set:
        """Get valid tail types for (h_type, relation) without validation."""
        return self.schema.get((h_type, relation), set())

    def validate_edge(
        self, h_type: str, relation: str, entity_name: str = ""
    ) -> set:
        """
        Validate that (h_type, relation) exists in schema.
        Returns valid tail types if valid, raises ValueError if invalid.
        """
        tail_types = self.get_tail_types(h_type, relation)
        if not tail_types:
            entity_info = f" (entity: {entity_name})" if entity_name else ""
            raise ValueError(
                f"Invalid query: No edges in graph schema for "
                f"({h_type}, {relation}){entity_info}. "
                f"This combination does not exist in PrimeKG."
            )
        return tail_types

    def follow_path(
        self, h_type: str, relations: list, validate: bool = False, entity_name: str = ""
    ) -> set:
        """
        Follow a path of relations from h_type.

        Args:
            h_type: Starting entity type
            relations: List of relations to follow
            validate: If True, raise errors for invalid paths
            entity_name: Entity name for error messages (used when validate=True)

        Returns:
            Set of reachable entity types at the end of the path
        """
        current_types = {h_type}
        for i, relation in enumerate(relations):
            next_types = set()
            for curr_type in current_types:
                if validate:
                    try:
                        tail_types = self.validate_edge(
                            curr_type, relation, entity_name if i == 0 else ""
                        )
                        next_types.update(tail_types)
                    except ValueError:
                        if i == 0:
                            raise  # First hop - use original error
                        else:
                            # Subsequent hop - no valid path from intermediate types
                            raise ValueError(
                                f"Invalid query: No valid path from {h_type} via {relations[:i+1]}. "
                                f"After {i} hop(s), reached types {current_types}, but none support "
                                f"relation '{relation}'."
                            )
                else:
                    # Non-validating mode: just accumulate tail types
                    next_types.update(self.get_tail_types(curr_type, relation))

            current_types = next_types
            if not current_types:
                if validate:
                    raise ValueError(
                        f"Invalid query: Path {h_type} → {' → '.join(relations[:i+1])} "
                        f"leads to no valid entity types."
                    )
                else:
                    break  # No valid path found
        return current_types

    def process_query(
        self,
        query_type: str,
        query_structure: tuple,
        entity_map: dict,
        relation_map: dict,
        validate: bool = True,
    ) -> set:
        """
        Process a query to compute expected tail types, optionally validating.

        Args:
            query_type: Query type (e.g., "1p", "2i", "ip")
            query_structure: Parsed query structure with indices
            entity_map: Dict mapping indices to (entity_id, entity_name)
            relation_map: Dict mapping indices to relation_label
            validate: If True, raise ValueError for invalid queries

        Returns:
            Set of valid tail types, or None if not applicable
        """
        try:
            if query_type == "1p":
                # Simple projection: (entity, (relation,))
                entity_idx, (rel_idx,) = query_structure
                entity_id, entity_name = entity_map[entity_idx]
                relation = relation_map[rel_idx]
                h_type = self.get_entity_type(entity_id)
                if h_type:
                    if validate:
                        return self.validate_edge(h_type, relation, entity_name)
                    else:
                        return self.get_tail_types(h_type, relation)

            elif query_type == "2p":
                # Two-hop projection: (entity, (r1, r2))
                entity_idx, (rel1_idx, rel2_idx) = query_structure
                entity_id, entity_name = entity_map[entity_idx]
                rel1 = relation_map[rel1_idx]
                rel2 = relation_map[rel2_idx]
                h_type = self.get_entity_type(entity_id)
                if h_type:
                    return self.follow_path(h_type, [rel1, rel2], validate, entity_name)

            elif query_type == "3p":
                # Three-hop projection: (entity, (r1, r2, r3))
                entity_idx, (rel1_idx, rel2_idx, rel3_idx) = query_structure
                entity_id, entity_name = entity_map[entity_idx]
                rel1 = relation_map[rel1_idx]
                rel2 = relation_map[rel2_idx]
                rel3 = relation_map[rel3_idx]
                h_type = self.get_entity_type(entity_id)
                if h_type:
                    return self.follow_path(h_type, [rel1, rel2, rel3], validate, entity_name)

            elif query_type == "2i":
                # Intersection: ((e1, (r1,)), (e2, (r2,)))
                (e1_idx, (r1_idx,)), (e2_idx, (r2_idx,)) = query_structure
                e1_id, e1_name = entity_map[e1_idx]
                e2_id, e2_name = entity_map[e2_idx]
                r1 = relation_map[r1_idx]
                r2 = relation_map[r2_idx]

                h1_type = self.get_entity_type(e1_id)
                h2_type = self.get_entity_type(e2_id)

                if h1_type and h2_type:
                    if validate:
                        types1 = self.validate_edge(h1_type, r1, e1_name)
                        types2 = self.validate_edge(h2_type, r2, e2_name)
                    else:
                        types1 = self.get_tail_types(h1_type, r1)
                        types2 = self.get_tail_types(h2_type, r2)

                    intersection = types1.intersection(types2)

                    if validate and not intersection:
                        raise ValueError(
                            f"Invalid query: No common entity types in intersection. "
                            f"Branch 1 ({h1_type}, {r1}) → {sorted(types1)}, "
                            f"Branch 2 ({h2_type}, {r2}) → {sorted(types2)}. "
                            f"These sets have no overlap."
                        )
                    elif validate:
                        logger.info(f"Intersection validation passed: common types = {intersection}")

                    return intersection

            elif query_type == "3i":
                # Three-way intersection: ((e1, (r1,)), (e2, (r2,)), (e3, (r3,)))
                (e1_idx, (r1_idx,)), (e2_idx, (r2_idx,)), (e3_idx, (r3_idx,)) = query_structure
                e1_id, e1_name = entity_map[e1_idx]
                e2_id, e2_name = entity_map[e2_idx]
                e3_id, e3_name = entity_map[e3_idx]
                r1 = relation_map[r1_idx]
                r2 = relation_map[r2_idx]
                r3 = relation_map[r3_idx]

                h1_type = self.get_entity_type(e1_id)
                h2_type = self.get_entity_type(e2_id)
                h3_type = self.get_entity_type(e3_id)

                if h1_type and h2_type and h3_type:
                    if validate:
                        types1 = self.validate_edge(h1_type, r1, e1_name)
                        types2 = self.validate_edge(h2_type, r2, e2_name)
                        types3 = self.validate_edge(h3_type, r3, e3_name)
                    else:
                        types1 = self.get_tail_types(h1_type, r1)
                        types2 = self.get_tail_types(h2_type, r2)
                        types3 = self.get_tail_types(h3_type, r3)

                    intersection = types1.intersection(types2).intersection(types3)

                    if validate and not intersection:
                        raise ValueError(
                            f"Invalid query: No common entity types in three-way intersection. "
                            f"Branch 1 ({h1_type}, {r1}) → {sorted(types1)}, "
                            f"Branch 2 ({h2_type}, {r2}) → {sorted(types2)}, "
                            f"Branch 3 ({h3_type}, {r3}) → {sorted(types3)}. "
                            f"These sets have no common overlap."
                        )
                    elif validate:
                        logger.info(f"Three-way intersection validation passed: common types = {intersection}")

                    return intersection

            elif query_type == "ip":
                # Intersection then projection: (((e1, (r1,)), (e2, (r2,))), (r3,))
                ((e1_idx, (r1_idx,)), (e2_idx, (r2_idx,))), (r3_idx,) = query_structure
                e1_id, e1_name = entity_map[e1_idx]
                e2_id, e2_name = entity_map[e2_idx]
                r1 = relation_map[r1_idx]
                r2 = relation_map[r2_idx]
                r3 = relation_map[r3_idx]

                h1_type = self.get_entity_type(e1_id)
                h2_type = self.get_entity_type(e2_id)

                if h1_type and h2_type:
                    # Validate/compute intersection branches
                    if validate:
                        types1 = self.validate_edge(h1_type, r1, e1_name)
                        types2 = self.validate_edge(h2_type, r2, e2_name)
                    else:
                        types1 = self.get_tail_types(h1_type, r1)
                        types2 = self.get_tail_types(h2_type, r2)

                    inter_types = types1.intersection(types2)

                    if validate and not inter_types:
                        raise ValueError(
                            f"Invalid query: No common entity types in intersection. "
                            f"Branch 1 ({h1_type}, {r1}) → {sorted(types1)}, "
                            f"Branch 2 ({h2_type}, {r2}) → {sorted(types2)}. "
                            f"These sets have no overlap."
                        )

                    if validate:
                        logger.info(f"Intersection produces types: {inter_types}")

                    # Project from intersection
                    final_types = set()
                    valid_projection = False
                    for inter_type in inter_types:
                        tail_types = self.get_tail_types(inter_type, r3)
                        if tail_types:
                            valid_projection = True
                            final_types.update(tail_types)

                    if validate and not valid_projection:
                        raise ValueError(
                            f"Invalid query: Intersection produces types {sorted(inter_types)}, "
                            f"but none of them support relation '{r3}'."
                        )

                    return final_types

            elif query_type == "pi":
                # Projection then intersection: ((e1, (r1, r2)), (e2, (r3,)))
                (e1_idx, (r1_idx, r2_idx)), (e2_idx, (r3_idx,)) = query_structure
                e1_id, e1_name = entity_map[e1_idx]
                e2_id, e2_name = entity_map[e2_idx]
                r1 = relation_map[r1_idx]
                r2 = relation_map[r2_idx]
                r3 = relation_map[r3_idx]

                h1_type = self.get_entity_type(e1_id)
                h2_type = self.get_entity_type(e2_id)

                if h1_type and h2_type:
                    # Validate/compute 2-hop projection
                    proj_types = self.follow_path(h1_type, [r1, r2], validate, e1_name)

                    # Validate/compute other branch
                    if validate:
                        other_types = self.validate_edge(h2_type, r3, e2_name)
                    else:
                        other_types = self.get_tail_types(h2_type, r3)

                    # Check intersection
                    intersection = proj_types.intersection(other_types)

                    if validate and not intersection:
                        raise ValueError(
                            f"Invalid query: No common entity types after projection and intersection. "
                            f"Path ({h1_type}, {r1}, {r2}) → {sorted(proj_types)}, "
                            f"Branch ({h2_type}, {r3}) → {sorted(other_types)}. "
                            f"These sets have no overlap."
                        )
                    elif validate:
                        logger.info(f"PI intersection validation passed: common types = {intersection}")

                    return intersection

            # For other query types, return None
            else:
                if validate:
                    logger.info(f"Schema validation not implemented for query type: {query_type}")
                else:
                    logger.info(f"Schema filtering not implemented for query type: {query_type}")
                return None

        except ValueError:
            raise  # Re-raise validation errors
        except Exception as e:
            if validate:
                logger.warning(f"Error during schema validation: {e}")
            else:
                logger.warning(f"Error computing expected tail types: {e}")
            return None

        return None


def validate_query_schema(
    query_structure: tuple,
    entity_map: dict,
    relation_map: dict,
    manager: ModelManager,
    query_type: str,
) -> None:
    """
    Validate that the query structure is valid according to the graph schema.
    Raises ValueError if the query contains invalid (h_type, relation) combinations.

    Args:
        query_structure: Parsed query structure with indices
        entity_map: Dict mapping indices to (entity_id, entity_name)
        relation_map: Dict mapping indices to relation_label
        manager: ModelManager instance
        query_type: Query type string (e.g., "2i", "ip")

    Raises:
        ValueError: If query contains invalid schema combinations
    """
    helper = SchemaHelper(manager)
    helper.process_query(query_type, query_structure, entity_map, relation_map, validate=True)


def compute_expected_tail_types(
    query_structure: tuple,
    entity_map: dict,
    relation_map: dict,
    manager: ModelManager,
    query_type: str,
) -> set:
    """
    Compute expected tail types based on query structure and graph schema.

    Args:
        query_structure: Parsed query structure with indices
        entity_map: Dict mapping indices to (entity_id, entity_name)
        relation_map: Dict mapping indices to relation_label
        manager: ModelManager instance
        query_type: Query type string (e.g., "2i", "ip")

    Returns:
        Set of valid tail types, or None if schema filtering is not applicable
    """
    helper = SchemaHelper(manager)
    return helper.process_query(query_type, query_structure, entity_map, relation_map, validate=False)


@torch.no_grad()
def run_ultraquery_inference(
    query_structure: Any,
    manager: ModelManager,
    top_k: Optional[int] = None,
) -> dict:
    """
    Run UltraQuery model inference for a complex logical query.

    IMPORTANT: Query is validated against graph schema BEFORE inference.
    Invalid queries (e.g., drug → "expression present") will raise ValueError.

    All predictions are automatically saved to parquet files in ./output/ directory:
    - predictions_all.parquet: All predictions with schema_match flag
    - predictions_filtered.parquet: Only predictions matching graph schema

    Schema validation checks:
    - 1p/2p/3p: Validates that each (h_type, relation) step exists in graph
    - 2i: Validates both branches exist and have overlapping tail types
    - 3i: Validates all three branches exist and have common tail types
    - ip: Validates intersection branches + projection from intermediate types
    - pi: Validates 2-hop path + intersection with other branch

    Schema filtering logic by query type:
    - 2i: t_pred_type must match schema for BOTH [h_type1, r1] AND [h_type2, r2]
    - 2p: t_pred_type must be reachable via h_type → r1 → intermediate → r2
    - 3p: t_pred_type must be reachable via h_type → r1 → i1 → r2 → i2 → r3
    - ip: First intersect to get intermediate types, then project via r3
    - pi: First project via 2-hop path, then intersect with other branch

    Args:
        query_structure: Nested tuple in BetaE format
        manager: ModelManager instance
        top_k: Number of top filtered predictions to return in response (None for all)

    Returns:
        Dictionary with:
        - query_type: str
        - query_structure: parsed structure with indices
        - entity_map: Dict mapping indices to (entity_id, entity_name)
        - relation_map: Dict mapping indices to relation_label
        - predictions: List of filtered predictions (limited by top_k)
        - total_predictions: Total number of predictions (before filtering)
        - filtered_predictions: Number of predictions after schema filtering
        - output_file_all: Path to parquet with all predictions
        - output_file_filtered: Path to parquet with filtered predictions
        - expected_tail_types: List of valid tail types from schema
    """
    logger.info(f"Running UltraQuery inference for structure: {query_structure}")

    # Parse and validate query
    parsed_structure, entity_map, relation_map = parse_query_structure(
        query_structure, manager
    )
    logger.info(f"Parsed structure: {parsed_structure}")
    logger.info(f"Entity map: {entity_map}")
    logger.info(f"Relation map: {relation_map}")

    # Identify query type
    query_type = identify_query_type(parsed_structure)
    logger.info(f"Query type: {query_type}")

    # Validate query against schema BEFORE running inference
    logger.info("Validating query against graph schema...")
    validate_query_schema(
        parsed_structure, entity_map, relation_map, manager, query_type
    )
    logger.info("Schema validation passed")

    # Convert to Query object
    query = Query.from_nested(parsed_structure, binary_op=True)
    logger.info(f"Query postfix notation: {query.tolist()}")
    logger.info(f"Query readable:\n{query.to_readable()}")

    # Get model and dataset
    model = manager.get_model()
    dataset = manager.dataset
    device = manager.device

    # Get graph data
    graph = dataset[0].to(device)

    # Create batch (single query)
    query_batch = query.unsqueeze(0).to(device)  # Shape: (1, query_length)

    # Run inference (no symbolic traversal during inference)
    pred = model(graph, query_batch, symbolic_traversal=False)
    logger.info(f"Prediction shape: {pred.shape}")  # Should be (1, num_entities)

    # Get scores and sort
    scores = pred[0].cpu()  # Shape: (num_entities,)
    sorted_indices = torch.argsort(scores, descending=True)
    total_entities = len(sorted_indices)
    logger.info(f"Processing {total_entities} entities...")

    # Compute expected tail types based on schema and query structure
    logger.info("Computing expected tail types from schema...")
    expected_tail_types = compute_expected_tail_types(
        parsed_structure, entity_map, relation_map, manager, query_type
    )
    logger.info(f"Expected tail types: {expected_tail_types}")

    # Use cached entity_type_map for O(1) lookups (already built in ModelManager)
    entity_type_map = manager.entity_type_map

    # Collect ALL predictions for parquet export
    all_predictions = []
    filtered_predictions = []
    for rank, entity_idx in enumerate(sorted_indices.tolist(), start=1):
        entity_id = manager.id2ent_dict[entity_idx]
        entity_name = manager.ent2name_dict.get(entity_id, entity_id)
        score = scores[entity_idx].item()

        # Get entity type from pre-built lookup map (O(1) instead of O(n))
        entity_type = entity_type_map.get(entity_id, "unknown")

        # Calculate percentile rank: 1 - (rank / total_entities)
        percentile_rank = 1.0 - (rank / total_entities)

        prediction = {
            "rank": rank,
            "entity_id": entity_id,
            "entity_name": entity_name,
            "entity_type": entity_type,
            "score": score,
            "percentile_rank": percentile_rank,
        }

        all_predictions.append(prediction)

        # Apply schema filtering
        if expected_tail_types is None or entity_type in expected_tail_types:
            prediction["schema_match"] = True
            filtered_predictions.append(prediction)
        else:
            prediction["schema_match"] = False

    # Re-rank filtered predictions
    for new_rank, pred in enumerate(filtered_predictions, start=1):
        pred["filtered_rank"] = new_rank
        # Recalculate percentile rank based on filtered results
        pred["filtered_percentile_rank"] = (
            1.0 - (new_rank / len(filtered_predictions))
            if filtered_predictions
            else 0.0
        )

    logger.info(
        f"Schema filtering: {len(all_predictions)} total → {len(filtered_predictions)} valid"
    )

    # Save all predictions to parquet (including both filtered and unfiltered)
    # Create output directory based on query type and entity names
    entity_names = "_".join(
        [v[1].replace(" ", "-").replace("/", "-") for v in entity_map.values()]
    )
    entity_names = entity_names[:100]  # Limit length

    output_dir = os.path.join(
        DEFAULT_CONFIG["output_dir"],
        query_type,
        entity_names,
    )
    os.makedirs(output_dir, exist_ok=True)

    # Save all predictions (unfiltered)
    output_file_all = os.path.join(output_dir, "predictions_all.parquet")
    logger.info(f"Saving {len(all_predictions)} predictions to {output_file_all}...")
    df_all = pl.DataFrame(all_predictions)
    df_all.write_parquet(output_file_all)
    logger.info(f"All results written to: {output_file_all}")

    # Save filtered predictions
    output_file_filtered = os.path.join(output_dir, "predictions_filtered.parquet")
    if filtered_predictions:
        logger.info(
            f"Saving {len(filtered_predictions)} filtered predictions to {output_file_filtered}..."
        )
        df_filtered = pl.DataFrame(filtered_predictions)
        df_filtered.write_parquet(output_file_filtered)
        logger.info(f"Filtered results written to: {output_file_filtered}")

    # Return only top_k filtered predictions in response to save context
    predictions_for_response = (
        filtered_predictions[:top_k] if top_k else filtered_predictions
    )

    return {
        "query_type": query_type,
        "query_structure": parsed_structure,
        "entity_map": {
            k: {"entity_id": v[0], "entity_name": v[1]} for k, v in entity_map.items()
        },
        "relation_map": relation_map,
        "predictions": predictions_for_response,
        "total_predictions": len(all_predictions),
        "filtered_predictions": len(filtered_predictions),
        "output_file_all": output_file_all,
        "output_file_filtered": output_file_filtered if filtered_predictions else None,
        "expected_tail_types": (
            list(expected_tail_types) if expected_tail_types else None
        ),
    }


@mcp.tool
def answer_complex_query(
    query_structure: Any,
    top_k: Optional[int] = 25,
) -> dict:
    """
    Answer a complex logical query using UltraQuery model.

    Supports multi-hop reasoning with intersection, union, and negation operations.

    IMPORTANT: All queries are validated against PrimeKG schema before inference.
    Invalid queries will be rejected with a clear error message explaining why.
    Example of invalid query: ["DRUGBANK:DB00001", ["expression present"]]
    This is invalid because drug entities don't have "expression present" relations in PrimeKG.

    Query Format (BetaE nested list format):
    - Entities: Use entity IDs (e.g., "MONDO:5301") or names (e.g., "Crohn disease")
    - Relations: Use relation labels (e.g., "associated with", "ppi")
    - Operations:
      - Projection: ["entity", ["relation"]] or multi-hop ["entity", ["r1", "r2"]]
      - Intersection: [["e1", ["r1"]], ["e2", ["r2"]]]
      - Union: [["e1", ["r1"]], ["e2", ["r2"]], ["u"]]
      - Negation: Use "n" in relation list, e.g., ["entity", ["relation", "n"]]

    Supported Query Types:
    - 1p: One projection - ["entity", ["relation"]]
    - 2p: Two projections - ["entity", ["relation1", "relation2"]]
    - 3p: Three projections - ["entity", ["relation1", "relation2", "relation3"]]
    - 2i: Two intersections - [["entity1", ["relation1"]], ["entity2", ["relation2"]]]
    - 3i: Three intersections - [["e1", ["r1"]], ["e2", ["r2"]], ["e3", ["r3"]]]
    - ip: Intersection then projection - [[["e1", ["r1"]], ["e2", ["r2"]]], ["r3"]]
    - pi: Projection then intersection - [["entity", ["r1", "r2"]], ["entity2", ["r3"]]]
    - 2in: Two intersections with negation - [["e1", ["r1"]], ["e2", ["r2", "n"]]]
    - 3in, inp, pin, pni: More complex patterns with negation
    - 2u-DNF, up-DNF: Union queries in Disjunctive Normal Form
    - 2u-DM, up-DM: Union queries in De Morgan form

    Examples:
    Note: For simple single-target analysis (1p queries), consider using ultra-inference-mcp
    (pipeline tool) instead, which is optimized for single-hop queries.

    1. Two-hop path (2p) - Target pathway analysis:
       ["GREM1", ["ppi", "associated with"]]
       → What diseases are associated with proteins that interact with GREM1?
       Use case: Find downstream disease associations through protein interactions

    2. Three-hop path (3p) - Extended target pathway analysis:
       ["ITGA4", ["ppi", "ppi", "associated with"]]
       → What diseases are associated with proteins that interact with ITGA4's interaction partners?
       Use case: Explore deeper pathway connections for combination therapy targets

    3. Intersection (2i) - Bispecific target analysis:
       [["TYK2", ["associated with"]], ["JAK1", ["associated with"]]]
       → What diseases are associated with both TYK2 AND JAK1?
       Use case: Identify shared disease indications for bispecific antibodies

    4. Intersection (3i) - Trispecific target analysis:
       [["ITGA4", ["associated with"]], ["ITGB7", ["associated with"]], ["TYK2", ["associated with"]]]
       → What diseases are associated with both ITGA4 AND ITGB7 AND TYK2?
       Use case: Identify shared disease indications for trispecific antibodies

    5. Intersection then projection (ip) - Pathway convergence:
       [[["GREM1", ["ppi"]], ["IL11", ["ppi"]]], ["associated with"]]
       → What diseases are associated with proteins that interact with both GREM1 and IL11?
       Use case: Find disease targets where two pathways converge

    Args:
        query_structure: Nested tuple in BetaE format (see examples above)
        top_k: Number of top filtered predictions to return (default: 25)

    Returns:
        Dictionary with:
        - success: bool
        - query_type: str (e.g., "2i", "ip")
        - query_readable: str (human-readable query)
        - entity_map: Dict of entities used in query
        - relation_map: Dict of relations used in query
        - predictions: List of filtered predictions (limited by top_k) with:
          - rank: Original rank among all predictions
          - filtered_rank: Rank among schema-filtered predictions
          - entity_id, entity_name, entity_type
          - score: Model prediction score
          - percentile_rank: 1 - (rank / total_predictions)
          - filtered_percentile_rank: 1 - (filtered_rank / filtered_predictions)
          - schema_match: bool
        - total_predictions: Total number of predictions (before filtering)
        - filtered_predictions: Number of predictions after schema filtering
        - output_file_all: Path to parquet file with all predictions
        - output_file_filtered: Path to parquet file with filtered predictions
        - expected_tail_types: List of valid tail types from graph schema
        - inference_time_seconds: float
    """
    start_time = time.time()

    try:
        logger.info(f"Received query: {query_structure}")

        # Get model manager
        manager = ModelManager()
        manager.load_model()

        # Run inference
        result = run_ultraquery_inference(
            query_structure=query_structure,
            manager=manager,
            top_k=top_k,
        )

        # Convert query to readable format
        query = Query.from_nested(result["query_structure"], binary_op=True)
        query_readable = query.to_readable()

        # Calculate inference time
        inference_time = time.time() - start_time

        return {
            "success": True,
            "query_type": result["query_type"],
            "query_readable": query_readable,
            "query_structure": str(result["query_structure"]),
            "entity_map": result["entity_map"],
            "relation_map": result["relation_map"],
            "predictions": result["predictions"],
            "total_predictions": result["total_predictions"],
            "filtered_predictions": result["filtered_predictions"],
            "output_file_all": result["output_file_all"],
            "output_file_filtered": result["output_file_filtered"],
            "expected_tail_types": result["expected_tail_types"],
            "inference_time_seconds": round(inference_time, 2),
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
def list_query_types() -> dict:
    """
    List all supported query types and their structure formats.

    Returns:
        Dictionary with:
        - query_types: List of query type codes
        - structures: Dict mapping query types to their structure templates
        - descriptions: Dict mapping query types to human-readable descriptions
    """
    descriptions = {
        "1p": "One-hop projection: Find entities connected by one relation",
        "2p": "Two-hop projection: Find entities connected by a path of two relations",
        "3p": "Three-hop projection: Find entities connected by a path of three relations",
        "2i": "Two-way intersection: Find entities satisfying both conditions",
        "3i": "Three-way intersection: Find entities satisfying all three conditions",
        "ip": "Intersection then projection: Apply intersection, then follow a relation",
        "pi": "Projection then intersection: Follow a path, then intersect with another condition",
        "2in": "Intersection with negation: Find entities satisfying first condition but NOT second",
        "3in": "Three-way intersection with negation",
        "inp": "Intersection with negation, then projection",
        "pin": "Projection, then intersection with negation",
        "pni": "Projection with negation, then intersection",
        "2u-DNF": "Union of two conditions (Disjunctive Normal Form)",
        "up-DNF": "Union then projection (DNF)",
        "2u-DM": "Union with negation (De Morgan form)",
        "up-DM": "Union with negation then projection (De Morgan)",
    }

    structure_templates = {
        "1p": '["entity", ["relation",]]',
        "2p": '["entity", ["relation1", "relation2"]]',
        "3p": '["entity", ["relation1", "relation2", "relation3"]]',
        "2i": '[["entity1", ["relation1"]], ["entity2", ["relation2"]]]',
        "3i": '[["e1", ["r1"]], ["e2", ["r2"]], ["e3", ["r3"]]]',
        "ip": '[[["e1", ["r1"]], ["e2", ["r2"]]], ["r3"]]',
        "pi": '[["entity", ["r1", "r2"]], ["entity2", ["r3"]]]',
        "2in": '[["e1", ["r1"]], ["e2", ["r2", "n"]]]',
        "3in": '[["e1", ["r1"]], ["e2", ["r2"]], ["e3", ["r3", "n"]]]',
        "inp": '[[["e1", ["r1"]], ["e2", ["r2", "n"]]], ["r3"]]',
        "pin": '[["entity", ["r1", "r2"]], ["entity2", ["r3", "n"]]]',
        "pni": '[["entity", ["r1", "r2", "n"]], ["entity2", ["r3"]]]',
        "2u-DNF": '[["e1", ["r1"]], ["e2", ["r2"]], ["u"]]',
        "up-DNF": '[[["e1", ["r1"]], ["e2", ["r2"]], ["u"]], ["r3"]]',
        "2u-DM": '[[["e1", ["r1", "n"]], ["e2", ["r2", "n"]]], ["n"]]',
        "up-DM": '[[["e1", ["r1", "n"]], ["e2", ["r2", "n"]]], ["n", "r3"]]',
    }

    return {
        "success": True,
        "query_types": list(STRUCT2TYPE.values()),
        "structures": structure_templates,
        "descriptions": descriptions,
    }


@mcp.tool
def setup_primekg_data(
    dataset_path: Optional[str] = None,
    force_redownload: bool = False,
    check_only: bool = False,
) -> dict:
    """
    Check if PrimeKG data is available and optionally install it if missing.

    This tool ensures that PrimeKG data required for UltraQuery inference is properly
    set up. It can check the status of existing data or download and process the
    complete PrimeKG dataset from Harvard Dataverse.

    Args:
        dataset_path: Path to dataset directory (e.g., /path/to/data/primekg1).
                     If None, uses the default path from server config.
        force_redownload: If True, re-download and reprocess even if data exists.
                         Use this to refresh the dataset. Default: False.
        check_only: If True, only check data status without installing.
                   Useful to verify setup before running queries. Default: False.

    Returns:
        Dictionary with status information:
        - success: bool - Whether the operation succeeded
        - status: str - Current status ("available", "missing", "installed", "error")
        - data_path: str - Path to the dataset directory
        - exists: dict - Status of each required file (if check_only=True)
        - nodes_count: int - Number of nodes (if installed)
        - train_edges: int - Number of training edges (if installed)
        - test_edges: int - Number of test edges (if installed)
        - valid_edges: int - Number of validation edges (if installed)
        - message: str - Human-readable status message
        - error: str - Error message (if status="error")

    Examples:
        # Check if PrimeKG data is available
        result = setup_primekg_data(check_only=True)
        if result["status"] == "missing":
            print("PrimeKG data needs to be installed")

        # Install PrimeKG data if missing
        result = setup_primekg_data()
        if result["status"] == "installed":
            print(f"Installed {result['nodes_count']} nodes")

        # Force re-download and reprocess
        result = setup_primekg_data(force_redownload=True)

    Notes:
        - Download size: ~500 MB (primekg.csv)
        - Processing time: ~2-5 minutes depending on system
        - Disk space required: ~2 GB for complete dataset
        - Data source: Harvard Dataverse (https://dataverse.harvard.edu/)
    """
    try:
        # Determine dataset path
        if dataset_path is None:
            dataset_path = DEFAULT_CONFIG["dataset"]["root"]
            # Append dataset name if not already in path
            if not dataset_path.endswith("primekg1"):
                dataset_path = os.path.join(dataset_path, "primekg1")

        logger.info(f"PrimeKG data path: {dataset_path}")

        # Check current status
        status_info = check_primekg_data(dataset_path)

        if check_only:
            # Only return status, don't install
            if status_info["all_present"]:
                return {
                    "success": True,
                    "status": "available",
                    "data_path": dataset_path,
                    "exists": status_info["exists"],
                    "message": "PrimeKG data is available and ready to use",
                }
            else:
                missing_files = [
                    name for name, exists in status_info["exists"].items() if not exists
                ]
                return {
                    "success": True,
                    "status": "missing",
                    "data_path": dataset_path,
                    "exists": status_info["exists"],
                    "missing_files": missing_files,
                    "message": f"PrimeKG data is missing. Missing files: {', '.join(missing_files)}",
                }

        # Install or update data
        if status_info["all_present"] and not force_redownload:
            return {
                "success": True,
                "status": "available",
                "data_path": dataset_path,
                "message": "PrimeKG data is already available. Use force_redownload=True to reinstall.",
            }

        # Run setup
        logger.info("Starting PrimeKG setup...")
        setup_result = setup_primekg(
            dataset_path=dataset_path,
            force_redownload=force_redownload,
            train_frac=0.8,
            test_frac=0.1,
            valid_frac=0.1,
            seed=42,
        )

        if setup_result["status"] == "completed":
            return {
                "success": True,
                "status": "installed",
                "data_path": dataset_path,
                "nodes_count": setup_result["nodes_count"],
                "train_edges": setup_result["train_edges"],
                "test_edges": setup_result["test_edges"],
                "valid_edges": setup_result["valid_edges"],
                "message": f"PrimeKG data successfully installed with {setup_result['nodes_count']:,} nodes and {setup_result['train_edges']:,} training edges",
            }
        else:
            return {
                "success": True,
                "status": "available",
                "data_path": dataset_path,
                "message": "PrimeKG data was already available",
            }

    except Exception as e:
        logger.error(f"Error in setup_primekg_data: {e}", exc_info=True)
        return {
            "success": False,
            "status": "error",
            "error": str(e),
            "message": f"Failed to setup PrimeKG data: {str(e)}",
        }


if __name__ == "__main__":
    logger.info("Starting UltraQuery Inference MCP Server...")
    mcp.run()
