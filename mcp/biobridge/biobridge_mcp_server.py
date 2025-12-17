#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ============================== Instruction ==============================
mcp_instruction ="""
BioBridge Multimodal KG — MCP Usage Instruction

You can call predict_associations to infer biomedical links. Use it whenever you need predict association between biomedical entities (e.g., “genes associated with disease X”, “what drugs relate to Y?”, “is A linked to B via relation R?”).

What the tool does (high level)
- Parses your context text, proposes a head entity, selects a relation family and tail type, and runs neural retrieval over BioBridge embeddings.
- Returns a compact list of high-confidence tail hits with node identifiers and scores.

When to call it
- Call predict_associations if the user asks for:
- associations/links between entities (gene↔disease, drug↔phenotype, etc.),
- candidates/shortlists for a head entity’s likely connections,
- validation of a specific pair (“Is GREM1 associated with Crohn disease?”).
- Do not call it for non-biomedical or non-KG tasks.

Parameters to pass
- context (required): A concise sentence capturing the task and key entities, including useful synonyms. Example:
  "Find diseases associated with GREM1 (aka gremlin 1). Focus on inflammatory bowel disease, Crohn disease variants."
- topk (optional): Desired number of returned results (the server may cap it globally). Reasonable default: 25.
- override_head_name / override_head_names / override_head_type (optional):
  * Use override_head_name for a single head entity (e.g., "GREM1")
  * Use override_head_names for multiple head entities as a list (e.g., ["IL11", "GREM1"]) - computes mean embedding
  * override_head_type is required when using override_head_names
  * If the head is implied, omit and rely on context.
- override_tail_name / override_tail_type (optional): Use only if the user specifies an exact tail target (e.g., "Crohn disease", "disease"). Otherwise omit and let retrieval surface candidates.
- relation_hint (optional): If the user names the relation family (e.g., "associated with", "treats", "interacts with"), pass it to narrow the relation.
- slidewindow: Omit unless you explicitly need to override the KG default.
- include_relation_catalog: False for normal user queries; True if the user asks to inspect relation IDs.
- include_debug: False normally; True only for troubleshooting (adds warnings/priors).

How to choose overrides
- Single head: If the user states a specific head (e.g., "Show diseases linked to GREM1"), set override_head_name="GREM1" and override_head_type="gene/protein".
- Multiple heads (gene signature): If the user mentions multiple genes/proteins (e.g., "IL11 and GREM1 signature"), set override_head_names=["IL11", "GREM1"] and override_head_type="gene/protein".
- If the user fixes the tail (e.g., "Is GREM1 linked to Crohn disease?"), also set override_tail_name="Crohn disease" and override_tail_type="disease".
- If the user names a relation (e.g., "associated with"), set relation_hint="associated with".
- Prefer not to invent types; use the canonical types used by the KG (e.g., "gene/protein", "disease", "drug", "phenotype").

How to present results
- If the tool returns an error, briefly explain and suggest a tighter query or adjusted overrides.
- Otherwise, read:
    - resolved.head_name, resolved.head_type
    - resolved.relation_family, resolved.tail_type
    - results (a list of hits with node_index, node_id/mondo_id, node_name/mondo_name, cos_sim, pct_rank)
- Show the top hits in a short list (3–10 items). Include the name and identifier plus a light explanation of the score (cosine similarity; higher is stronger) or percentile rank.
- Avoid clinical claims; say this is a model-based retrieval from a knowledge graph.

Do/Don’t
- Keep context compact and entity-rich (names + common synonyms).
- Use overrides only when the user is explicit.
- If users want a yes/no about a specific pair, set both head and tail overrides and topk=1.
- Don’t fabricate entities, relations, or IDs.
- Don’t rely on external web search to answer what the tool is designed to compute.
- Don’t expose raw embeddings or internal file paths.

Examples
- Head-only exploration
{
  "tool_name": "predict_associations",
  "arguments": {
    "context": "Find diseases associated with GREM1 (gremlin 1). Emphasize IBD and Crohn variants.",
    "override_head_name": "GREM1",
    "override_head_type": "gene/protein",
    "relation_hint": "associated with",
    "topk": 25,
    "include_relation_catalog": false,
    "include_debug": false
  }
}

- Pair validation (yes/no-ish)
{
  "tool_name": "predict_associations",
  "arguments": {
    "context": "Is GREM1 associated with Crohn disease (including esophageal involvement)?",
    "override_head_name": "GREM1",
    "override_head_type": "gene/protein",
    "override_tail_name": "Crohn disease",
    "override_tail_type": "disease",
    "relation_hint": "associated with",
    "topk": 1
  }
}

- Drug → phenotype
{
  "tool_name": "predict_associations",
  "arguments": {
    "context": "Which phenotypes are associated with infliximab (anti-TNF)? Focus on adverse effects.",
    "override_head_name": "Infliximab",
    "override_head_type": "drug",
    "relation_hint": "side effect",
    "topk": 15
  }
}

Interpreting low/empty results
- If results is empty or low-rank, suggest:
    - broadening the relation family (omit relation_hint),
    - relaxing specificity in context (fewer modifiers),
    - removing tail overrides to let the model retrieve candidates,
    - increasing topk.

Notes on identifiers
- Retrieval joins on node_index (authoritative key).
- Names (node_name/mondo_name) are display labels; identical names can differ by index.
- relation_family is mapped internally to one or more numeric relation IDs; you may show the family name and (optionally) the first relation ID from resolved.chosen_relations.
"""

# ============================== Standard Library ==============================
import json
import logging
import os
import pickle
import re
import sys
import difflib
from contextlib import asynccontextmanager
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Dict, List, Optional, Tuple

# ============================== Third-Party ==============================
import numpy as np
import pandas as pd
import torch
from mcp.server.fastmcp import Context, FastMCP, Icon
from dotenv import load_dotenv
import boto3
from botocore.exceptions import ClientError
from pathlib import Path

load_dotenv()

# ==============================================================================
# Logging
# ==============================================================================
LOG_PATH = "./biobridge_mcp.log"


def _setup_logging() -> None:
    """Configure file-based logging for the MCP server."""
    os.makedirs(os.path.dirname(LOG_PATH) or ".", exist_ok=True)
    fmt = "%(asctime)s [BioBridge MCP pid=%(process)d] %(message)s"
    logging.basicConfig(
        level=logging.INFO,
        format=fmt,
        handlers=[logging.FileHandler(LOG_PATH, mode="a", encoding="utf-8")],
        force=True,
    )


def _log(msg: str) -> None:
    """Write a message to the configured log file."""
    logging.info(msg)


_setup_logging()
_log("[main] starting MCP stdio server...")

# ==============================================================================
# S3 Utilities
# ==============================================================================
def _is_s3_uri(path: str) -> bool:
    """Check if a path is an S3 URI."""
    return isinstance(path, str) and path.startswith("s3://")


def _parse_s3_uri(s3_uri: str) -> Tuple[str, str]:
    """Parse S3 URI into bucket and key.

    Args:
        s3_uri: S3 URI in format s3://bucket/key/path

    Returns:
        Tuple of (bucket, key)
    """
    if not s3_uri.startswith("s3://"):
        raise ValueError(f"Invalid S3 URI: {s3_uri}")

    parts = s3_uri[5:].split("/", 1)
    bucket = parts[0]
    key = parts[1] if len(parts) > 1 else ""
    return bucket, key


def _get_s3_client():
    """Get or create boto3 S3 client."""
    return boto3.client("s3")


def _s3_file_exists(s3_uri: str) -> bool:
    """Check if a file exists in S3."""
    try:
        bucket, key = _parse_s3_uri(s3_uri)
        s3_client = _get_s3_client()
        s3_client.head_object(Bucket=bucket, Key=key)
        return True
    except ClientError:
        return False


def _read_s3_file_object(s3_uri: str):
    """Read file from S3 and return file-like object.

    Args:
        s3_uri: S3 URI of the file

    Returns:
        File-like object (BytesIO)
    """
    import io
    bucket, key = _parse_s3_uri(s3_uri)
    s3_client = _get_s3_client()

    _log(f"[s3] reading {s3_uri} on the fly")
    response = s3_client.get_object(Bucket=bucket, Key=key)
    return io.BytesIO(response['Body'].read())


def _read_s3_json(s3_uri: str) -> dict:
    """Read JSON file directly from S3.

    Args:
        s3_uri: S3 URI of the JSON file

    Returns:
        Parsed JSON as dict
    """
    file_obj = _read_s3_file_object(s3_uri)
    return json.load(file_obj)


def _read_s3_pickle(s3_uri: str) -> Any:
    """Read pickle file directly from S3.

    Args:
        s3_uri: S3 URI of the pickle file

    Returns:
        Unpickled object
    """
    file_obj = _read_s3_file_object(s3_uri)
    return pickle.load(file_obj)


def _read_s3_csv(s3_uri: str, **kwargs) -> pd.DataFrame:
    """Read CSV file directly from S3.

    Args:
        s3_uri: S3 URI of the CSV file
        **kwargs: Additional arguments to pass to pd.read_csv

    Returns:
        DataFrame
    """
    file_obj = _read_s3_file_object(s3_uri)
    return pd.read_csv(file_obj, **kwargs)


def _resolve_path(path: str, base_dir: str) -> str:
    """Resolve a path, returning S3 URI or local path.

    Args:
        path: Path to resolve (can be local or S3 URI)
        base_dir: Base directory (can be local or S3 URI)

    Returns:
        Full path (S3 URI if base is S3, local path otherwise)
    """
    # If path is already absolute and local, return as-is
    if os.path.isabs(path) and not _is_s3_uri(path):
        return path

    # Construct full path
    if _is_s3_uri(base_dir):
        # Base is S3, construct S3 URI and return it directly
        full_s3_uri = base_dir.rstrip("/") + "/" + path.lstrip("/")
        return full_s3_uri
    else:
        # Base is local, use normal path join
        return os.path.join(base_dir, path)


# ==============================================================================
# Paths & Constants
# ==============================================================================
# Support both local path and S3 URI
_default_src_dir = "s3://tec-rnd-sci-dev-gi2/gi2-xz/biobridge/" # os.path.dirname(os.path.abspath(__file__))
BIOBRIDGE_SRC_DIR = os.getenv("BIOBRIDGE_SRC_DIR", _default_src_dir)

from src.inference import BridgeInference
from src.model import BindingModel

KG = "primekg_bulk_excelraopiintegration-lfc_slidewindow_weightedclass-11112024"

# Helper function to build paths that work with both local and S3
def _build_path(*parts: str) -> str:
    """Build a path relative to BIOBRIDGE_SRC_DIR, handling both local and S3."""
    if _is_s3_uri(BIOBRIDGE_SRC_DIR):
        # For S3, join with forward slashes
        path = "/".join(parts)
        return _resolve_path(path, BIOBRIDGE_SRC_DIR)
    else:
        # For local, use os.path.join
        return os.path.join(BIOBRIDGE_SRC_DIR, *parts)

DATA_BASE = _build_path("data")
EMB_SUBDIR = _build_path("data", "embeddings", "esm2b_unimo_pubmedbert")

BIND_CONFIG_PATH = _build_path("data", "BindData", KG, "data_config.json")
MODEL_CONFIG_PATH = _build_path("ckpt", KG, "model_6layer_100epoch", "model_config.json")
CKPT_PATH = _build_path("ckpt", KG, "model_6layer_100epoch", "model.bin")
KG_CSV_PATH = _build_path("data", "PrimeKG", "kg.csv")
NODES_CSV_PATH = _build_path("data", "PrimeKG", "nodes.csv")

# Log source directory type
if _is_s3_uri(BIOBRIDGE_SRC_DIR):
    _log(f"[init] Using S3 source directory: {BIOBRIDGE_SRC_DIR}")
else:
    _log(f"[init] Using local source directory: {BIOBRIDGE_SRC_DIR}")


def _pick_device() -> str:
    """Pick CUDA device by env var (CUDA_DEVICE_INDEX) or CPU if unavailable."""
    if torch.cuda.is_available():
        try:
            n = torch.cuda.device_count()
            idx = int(os.getenv("CUDA_DEVICE_INDEX", "0"))
            if idx < 0 or idx >= n:
                idx = 0
            return f"cuda:{idx}"
        except Exception:
            return "cuda:0"
    return "cpu"


DEFAULT_DEVICE = _pick_device()
slidewindow_default = True if "_slidewindow" in KG else False

# Canonical entity types and mapping to short codes for file naming.
ENTITY_DICT: Dict[str, str] = {
    "biological_process": "bp",
    "molecular_function": "mf",
    "cellular_component": "cc",
    "gene/protein": "protein",
    "disease": "disease",
    "drug": "drug",
    "pathway": "pathway",
    "effect/phenotype": "phenotype",
    "anatomy": "anatomy",
    "exposure": "exposure",
    "biologics_drug": "biologics_drug",
}
CANON_ENTITY_TYPES = set(ENTITY_DICT.keys())

# Global cap (logic preserved — `predict_associations` assigns from this).
TOPK_HARD_LIMIT: Optional[int] = 100

# ==============================================================================
# MCP Server
# ==============================================================================
mcp = FastMCP(
    name="BioBridge MCP",
    instructions=mcp_instruction,
    icons=[
        Icon(
            src="https://raw.githubusercontent.com/ryanwangzf/biobridge/main/unimodal/esm2_logo.png",
            mimeType="image/png",
        )
    ],
)

# ==============================================================================
# String Normalization Utilities
# ==============================================================================
_WORD_RE = re.compile(r"[A-Za-z0-9]+")


def _tokens(s: str) -> List[str]:
    """Split a string into lowercase alphanumeric tokens."""
    return _WORD_RE.findall(str(s or "").lower())


def _norm(s: str) -> str:
    """Strong normalization (letters/digits only), e.g., Crohn’s -> crohns."""
    return re.sub(r"[^a-z0-9]+", "", str(s or "").lower())


def _normalize_type_alias(t: Optional[str]) -> Optional[str]:
    """Normalize synonymous type labels to canonical ones."""
    if not t:
        return t
    t = t.strip().lower()
    if t in {"gene", "protein", "gene/protein"}:
        return "gene/protein"
    if t in {"phenotype", "effect", "effect/phenotype"}:
        return "effect/phenotype"
    return t


def _emb_file_name(entity_type: str, slidewindow: bool) -> str:
    """Return base embedding filename for an entity type (slidewindow-aware)."""
    base = ENTITY_DICT[entity_type]
    if entity_type in ["gene/protein", "biologics_drug"] and slidewindow:
        return f"{base}_slidewindow"
    return base


def _load_pickle(path: str) -> Dict[str, Any]:
    """Load a pickle file and return its content (supports S3)."""
    if _is_s3_uri(path):
        return _read_s3_pickle(path)
    with open(path, "rb") as f:
        return pickle.load(f)


def _safe_int(v: Any) -> Optional[int]:
    """Safely cast a value to int, returning None on failure/NaN."""
    try:
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return None
        return int(v)
    except Exception:
        try:
            return int(float(v))
        except Exception:
            return None


def _avg_embeddings_for_nodes(
    raw_index: np.ndarray, raw_emb: np.ndarray, node_indices: np.ndarray
) -> np.ndarray:
    """Average embedding rows for each node index."""
    out: List[np.ndarray] = []
    for i in node_indices:
        sel = raw_emb[raw_index == i]
        out.append(sel[0] if sel.shape[0] == 1 else sel.mean(axis=0))
    return np.stack(out, axis=0)


def _df_best_match_by_name(df: pd.DataFrame, name: str) -> pd.DataFrame:
    """Best-effort match by exact, substring, or normalized containment."""
    if df.empty:
        return df
    s = str(name).strip().lower()
    node_col = "node_name" if "node_name" in df else "mondo_name"
    df_exact = df[df[node_col].astype(str).str.lower() == s]
    if not df_exact.empty:
        return df_exact
    df_cont = df[df[node_col].astype(str).str.lower().str.contains(re.escape(s))]
    if not df_cont.empty:
        return df_cont
    s2 = _norm(s)
    df_norm = df[
        df[node_col].apply(
            lambda x: _norm(str(x)).find(s2) >= 0 if pd.notna(x) else False
        )
    ]
    return df_norm


def _jsonify(obj: Any) -> Any:
    """Make an object JSON-serializable where possible (structure-preserving)."""
    try:
        json.dumps(obj)
        return obj
    except Exception:
        if isinstance(obj, dict):
            return {str(k): _jsonify(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_jsonify(x) for x in obj]
        if hasattr(obj, "tolist"):
            return obj.tolist()
        return str(obj)


# ==============================================================================
# nodes.csv Helpers
# ==============================================================================
@lru_cache(maxsize=4096)
def _context_tokens(context: str) -> set[str]:
    """Tokenize context text and return a cached set of tokens."""
    return set(_tokens(context or ""))


def _aliases_for_row(row: pd.Series, max_syn: int = 3) -> List[str]:
    """Extract up to `max_syn` unique synonyms from a nodes row."""
    if "synonyms" not in row.index:
        return []
    syn = row.get("synonyms", None)
    if pd.isna(syn):
        return []
    vals = [x.strip() for x in str(syn).split("|") if x.strip()]
    out: List[str] = []
    seen: set[str] = set()
    for n in vals:
        k = n.lower()
        if k not in seen:
            out.append(n)
            seen.add(k)
        if len(out) >= max_syn:
            break
    return out


def _build_node_index(
    nodes: Optional[pd.DataFrame],
) -> Tuple[Dict[str, List[int]], Dict[str, pd.DataFrame]]:
    """Build normalized name->indices and type->df indexes from nodes DataFrame.

    Optimization improvements:
    1. Pre-compute lowercase names for faster case-insensitive lookup
    2. Build token-based index for partial matching
    3. Index common abbreviations and acronyms
    4. Store node_index as the actual row identifier for O(1) lookups
    """
    if nodes is None or nodes.empty or not {"node_name", "node_type"}.issubset(nodes.columns):
        return {}, {}

    nodes["__norm_name__"] = nodes["node_name"].astype(str).map(_norm)
    nodes["__lower_name__"] = nodes["node_name"].astype(str).str.lower()
    nodes["__tokens__"] = nodes["node_name"].astype(str).map(_tokens)

    name_index: Dict[str, List[int]] = {}
    by_type: Dict[str, pd.DataFrame] = {}

    # Build type index
    for t in nodes["node_type"].dropna().unique():
        by_type[str(t)] = nodes[nodes["node_type"] == t].copy()

    # Build name index with multiple keys per entity
    for i, row in nodes.iterrows():
        nm = str(row.get("node_name", "")).strip()

        # 1. Normalized key (alphanumeric only)
        key = _norm(nm)
        if key:
            name_index.setdefault(key, []).append(int(i))

        # 2. Lowercase key (preserves structure)
        key_lower = nm.lower()
        if key_lower and key_lower != key:
            name_index.setdefault(key_lower, []).append(int(i))

        # 3. Acronym key (for multi-word entities like "Crohn Disease" -> "cd")
        tokens = _tokens(nm)
        if len(tokens) > 1:
            acronym = "".join([t[0] for t in tokens if t])
            if acronym and len(acronym) >= 2:
                name_index.setdefault(acronym.lower(), []).append(int(i))

        # 4. Individual significant tokens (for partial matching)
        for token in tokens:
            if len(token) >= 4:  # Only index meaningful tokens (4+ chars)
                token_key = _norm(token)
                if token_key:
                    # Use a special prefix to distinguish token-based keys
                    name_index.setdefault(f"_tok_{token_key}", []).append(int(i))

        # 5. Synonyms and aliases
        if "synonyms" in nodes.columns:
            for alias in _aliases_for_row(row, max_syn=5):
                ak = _norm(alias)
                if ak and ak not in name_index:
                    name_index.setdefault(ak, []).append(int(i))

                # Also index lowercase version of synonym
                ak_lower = alias.lower()
                if ak_lower != ak:
                    name_index.setdefault(ak_lower, []).append(int(i))

    _log(f"[index] built name_index with {len(name_index)} keys for {len(nodes)} nodes")
    return name_index, by_type


def _score_entity_for_context(
    name: str, ctx: str, ctx_toks: set[str]
) -> Tuple[bool, int, float]:
    """Score a candidate entity against context using exact/overlap heuristics."""
    ln = str(name or "").lower()
    exact = ln in (ctx or "").lower()
    ntoks = set(_tokens(ln))
    inter = len(ntoks & ctx_toks)
    uni = max(1, len(ntoks | ctx_toks))
    jacc = inter / uni
    score = (3.0 if exact else 0.0) + (1.0 if inter > 0 else 0.0) + jacc
    return exact, inter, score


def _candidate_rows_by_context(
    nodes: Optional[pd.DataFrame],
    context: str,
    type_hint: Optional[str] = None,
    limit: int = 40,
) -> pd.DataFrame:
    """Context-guided candidate selection (prefers exact/substr matches)."""
    if nodes is None or nodes.empty:
        return nodes.iloc[0:0] if nodes is not None else pd.DataFrame()

    df = nodes if not type_hint else nodes[nodes["node_type"] == type_hint]
    if df.empty:
        return df

    ctx = (context or "")
    ctx_toks = _context_tokens(ctx)

    exact_mask = df["node_name"].astype(str).str.lower().apply(lambda x: x in ctx.lower())
    if ctx_toks:
        patt = "|".join(map(re.escape, sorted(ctx_toks, key=len, reverse=True)))
        sub_mask = df["node_name"].astype(str).str.lower().str.contains(patt, na=False)
    else:
        sub_mask = False

    cand = df[exact_mask | sub_mask].copy()
    if cand.empty:
        cand = df.copy()

    cand[["_exact", "_overlap", "_score"]] = cand["node_name"].apply(
        lambda n: pd.Series(_score_entity_for_context(n, ctx, ctx_toks))
    )
    cand = cand.sort_values(
        ["_score", "_overlap", "node_name"], ascending=[False, False, True]
    ).head(limit)
    return cand.drop(columns=["_score"], errors="ignore")


@lru_cache(maxsize=1024)
def _fuzzy_match_cached(query: str, choices_tuple: tuple, cutoff: float = 0.6) -> Optional[str]:
    """Cached fuzzy matching to avoid redundant difflib calls."""
    best = difflib.get_close_matches(query, list(choices_tuple), n=1, cutoff=cutoff)
    return best[0] if best else None


def _lookup_entity(
    nodes: Optional[pd.DataFrame],
    name_index: Dict[str, List[int]],
    query_name: Optional[str],
    type_hint: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """Resolve an entity by normalized exact match, case-insensitive match, or fuzzy match (type-aware).

    Three-stage matching:
    1. Exact normalized match via index (fast)
    2. Case-insensitive exact match
    3. Partial substring or fuzzy match
    """
    if nodes is None or not query_name:
        return None

    query_name_clean = str(query_name).strip()
    df_scan = nodes if not type_hint else nodes[nodes["node_type"] == type_hint]

    # Stage 1: Exact normalized match via index
    key = _norm(query_name_clean)
    indices = name_index.get(key, [])

    if indices:
        rows = nodes.loc[indices]
        if type_hint and "node_type" in rows.columns:
            # Prioritize type-matching rows
            type_matches = rows[rows["node_type"] == type_hint]
            if not type_matches.empty:
                r = type_matches.iloc[0]
            else:
                r = rows.iloc[0]
        else:
            r = rows.iloc[0]
        return {
            "name": str(r.get("node_name")),
            "type": str(r.get("node_type")),
            "node_index": _safe_int(r.get("node_index")),
        }

    # Stage 2: Case-insensitive exact match
    query_lower = query_name_clean.lower()
    exact_match = df_scan[df_scan["node_name"].astype(str).str.lower() == query_lower]
    if not exact_match.empty:
        r = exact_match.iloc[0]
        return {
            "name": str(r.get("node_name")),
            "type": str(r.get("node_type")),
            "node_index": _safe_int(r.get("node_index")),
        }

    # Stage 3: Partial substring match, then fuzzy match as fallback
    # Try partial contains first (fast)
    contains_match = df_scan[df_scan["node_name"].astype(str).str.lower().str.contains(re.escape(query_lower), na=False)]
    if not contains_match.empty:
        r = contains_match.iloc[0]
        return {
            "name": str(r.get("node_name")),
            "type": str(r.get("node_type")),
            "node_index": _safe_int(r.get("node_index")),
        }

    # Fuzzy match as last resort
    choices = df_scan["node_name"].astype(str).tolist()
    if len(choices) > 1000:
        # For large datasets, sample for performance
        sample_df = df_scan.sample(n=min(1000, len(df_scan)), random_state=42)
        choices = sample_df["node_name"].astype(str).tolist()
        df_scan = sample_df

    if choices:
        choices_tuple = tuple(choices)
        best = _fuzzy_match_cached(query_name_clean, choices_tuple, cutoff=0.6)
        if best:
            hit = df_scan[df_scan["node_name"].astype(str) == best].iloc[0]
            return {
                "name": str(hit.get("node_name")),
                "type": str(hit.get("node_type")),
                "node_index": _safe_int(hit.get("node_index")),
            }

    return None


def _make_entity_payload(df: pd.DataFrame, *, context: str) -> List[dict]:
    """Build a lightweight payload for LLM selection with hints/aliases."""
    ctx_toks = _context_tokens(context)
    out: List[dict] = []
    need_cols = {"node_name", "node_type", "node_index"}
    if not need_cols.issubset(df.columns):
        return out
    for _, r in df.iterrows():
        name = str(r.get("node_name", ""))
        typ = str(r.get("node_type", ""))
        node_idx = _safe_int(r.get("node_index"))
        exact, overlap, _ = _score_entity_for_context(name, context, ctx_toks)
        aliases = _aliases_for_row(r, max_syn=3)
        out.append(
            {
                "name": name,
                "type": typ,
                "node_index": node_idx if node_idx is not None else -1,
                "aliases": aliases,
                "hint": {"exact": bool(exact), "overlap": int(overlap)},
            }
        )
    return out


# ==============================================================================
# Type-Aware Helpers
# ==============================================================================

def _allowed_head_types_from_examples(example_guide: List[Dict[str, str]]) -> List[str]:
    """Return allowed head types from the example guide (normalized & deduped)."""
    return sorted(
        {_normalize_type_alias(ex["head_type"]) for ex in example_guide if ex.get("head_type")}
    )


def _tail_types_for_head_from_examples(
    example_guide: List[Dict[str, str]], head_type: str
) -> List[str]:
    """List tail types that are valid for a given head type (from examples)."""
    ht = _normalize_type_alias(head_type)
    return sorted(
        {
            _normalize_type_alias(ex["tail_type"])
            for ex in example_guide
            if _normalize_type_alias(ex.get("head_type")) == ht
        }
    )


def _build_type_token_vocab(
    by_type: Dict[str, pd.DataFrame], max_rows: int = 2000
) -> Dict[str, set]:
    """Build a lightweight token vocabulary per type for type-prior estimation."""
    vocab: Dict[str, set] = {}
    for t, df in by_type.items():
        t_norm = _normalize_type_alias(t)
        sample = df.head(max_rows)
        toks: set = set()
        for n in sample["node_name"].astype(str):
            toks |= set(_tokens(n))
        vocab[t_norm] = toks
    return vocab


def _derive_type_priors(context: str, type_token_vocab: Dict[str, set], max_types: int = 4) -> List[str]:
    """Estimate top-k likely types from context using token overlap."""
    ctx = _context_tokens(context)
    scored: List[Tuple[str, int]] = []
    for t, vocab in type_token_vocab.items():
        inter = len(ctx & vocab)
        if inter > 0:
            scored.append((t, inter))
    scored.sort(key=lambda x: x[1], reverse=True)
    return [t for t, _ in scored][:max_types]


def _head_candidates_stratified(
    nodes: pd.DataFrame, context: str, head_types_scope: List[str], per_type: int = 12, total: int = 30
) -> pd.DataFrame:
    """Balanced head-candidate frame: top-k per type; fallback to all-types list."""
    parts: List[pd.DataFrame] = []
    for t in head_types_scope:
        df_t = _candidate_rows_by_context(nodes, context, type_hint=t, limit=per_type)
        if not df_t.empty:
            parts.append(df_t)
    if not parts:
        return _candidate_rows_by_context(nodes, context, type_hint=None, limit=total)
    df = pd.concat(parts, ignore_index=True)
    cols = [c for c in ["_exact", "_overlap", "node_name", "mondo_name"] if c in df.columns]
    df = df.drop_duplicates(subset=["node_index"]) if "node_index" in df.columns else df
    if cols:
        df = df.sort_values(cols, ascending=[False, False, True]).head(total)
    else:
        df = df.head(total)
    return df


# ==============================================================================
# LLM Selection (Head + Relation/Tail)
# ==============================================================================

def _safe_json_extract(txt: str) -> dict:
    """Extract a JSON object from a string; return empty dict on failure."""
    if not txt:
        return {}
    s, e = txt.find("{"), txt.rfind("}")
    if s == -1 or e == -1 or e <= s:
        return {}
    try:
        obj = json.loads(txt[s : e + 1])
        return obj if isinstance(obj, dict) else {}
    except Exception:
        return {}


def _chat_json(
    client: Any, *, model: str, system: str, user: str, fallback: dict
) -> dict:
    """Call a chat model and safely parse strict-JSON output with a fallback."""
    try:
        resp = client.chat.completions.create(
            model=model,
            temperature=0.0,
            max_tokens=192,
            messages=[{"role": "system", "content": system}, {"role": "user", "content": user}],
        )
        obj = _safe_json_extract(resp.choices[0].message.content)
        return obj if obj else fallback
    except Exception as e:
        _log(f"[llm] call failed; using fallback. {type(e).__name__}: {e}")
        return fallback


HEAD_PICK_SYS = (
    "You are selecting ONE head entity strictly from the provided list.\n"
    "- Output STRICT JSON only.\n"
    "- Do NOT invent values or types.\n"
    "- Prefer candidates whose TYPE is in allowed_head_types and that align with type_prior_hint.\n"
    "- If multiple look plausible, PICK index 0.\n"
    "- Prefer candidates with hint.exact=true; otherwise highest hint.overlap."
)


def _head_pick_user(
    context: str, candidates: List[dict], allowed_head_types: List[str], type_prior_hint: List[str]
) -> str:
    """Format the user prompt for the head-entity JSON picker."""
    return (
        f"Context:\n{context or ''}\n\n"
        f"Allowed head types: {allowed_head_types}\n"
        f"Type priors suggested by the context (descending priority): {type_prior_hint}\n\n"
        "Candidates (choose EXACTLY one by 0-based index):\n"
        + json.dumps(candidates, ensure_ascii=False)
        + '\nRespond ONLY with JSON object: {"head_idx": <int>}  (no prose).'
    )


COMBO_PICK_SYS = (
    "You will select ONE relation/tail_type combo and optionally one tail candidate.\n"
    "- Output STRICT JSON only.\n"
    "- Do NOT invent values.\n"
    "- Only pick from the menu (which is derived from example_guide).\n"
    "- Prefer relations and tail types whose words appear in the context or in type_prior_hint.\n"
    "- If unsure, PICK index 0 and tail_idx=null."
)


def _combo_pick_user(
    context: str,
    head_name: str,
    head_type: str,
    triplet_menu: List[dict],
    tails_by_type: Dict[str, List[dict]],
    allowed_tail_types: List[str],
    type_prior_hint: List[str],
) -> str:
    """Format the user prompt for the relation/tail JSON picker."""
    return (
        f"Context:\n{context or ''}\n\n"
        f"Head selected: {head_name} (type={head_type})\n"
        f"Allowed tail types for this head (from example_guide): {allowed_tail_types}\n"
        f"Tail type priors suggested by the context: {type_prior_hint}\n\n"
        "Pick ONE relation/tail_type combo (0-based index) from:\n"
        + json.dumps(triplet_menu, ensure_ascii=False)
        + "\n\nTail candidates by type (may be empty):\n"
        + json.dumps(tails_by_type, ensure_ascii=False)
        + '\nRespond ONLY with JSON object: {"combo_idx": <int>, "tail_idx": <int|null>}  (no prose).'
    )


def _enumerate_triplet_candidates(
    head_type: str,
    relation_hint: Optional[str],
    example_guide: List[Dict[str, str]],
    allowed_rel: List[str],
) -> List[Tuple[str, str, str]]:
    """Enumerate valid (head_type, relation, tail_type) combos filtered by hints."""
    combos = [
        (ex["head_type"], ex["rel"], ex["tail_type"])
        for ex in example_guide
        if _normalize_type_alias(ex.get("head_type")) == _normalize_type_alias(head_type)
        and ex.get("rel") in allowed_rel
    ]
    if relation_hint and relation_hint in allowed_rel:
        combos = [c for c in combos if c[1] == relation_hint] or combos
    seen: set = set()
    out: List[Tuple[str, str, str]] = []
    for c in combos:
        if c not in seen:
            out.append(c)
            seen.add(c)
    return out


def _llm_select_combo(
    *,
    context: str,
    nodes: Optional[pd.DataFrame],
    name_index: Dict[str, List[int]],
    relation_hint: Optional[str],
    example_guide: List[Dict[str, str]],
    allowed_rel: List[str],
    excluded_triples: set,
    override_head_name: Optional[str] = None,
    override_head_type: Optional[str] = None,
    type_priors: Optional[List[str]] = None,
) -> Optional[Dict[str, Any]]:
    """Select (head, relation family, tail_type) via LLM with strict-JSON protocol."""
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        return None
    try:
        from openai import OpenAI  # type: ignore
    except Exception:
        return None

    client = OpenAI(api_key=api_key, base_url=os.getenv("OPENAI_BASE_URL", None))
    model = os.getenv("MCP_LLM_MODEL", "gpt-4o-mini")

    # ---- Allowed types and priors ----
    allowed_head_types = _allowed_head_types_from_examples(example_guide)
    priors = [t for t in (type_priors or []) if t in allowed_head_types]
    head_types_scope = priors or allowed_head_types or list(CANON_ENTITY_TYPES)

    # -------------------- Step A: HEAD --------------------
    if override_head_name:
        hit = _lookup_entity(nodes, name_index, override_head_name, override_head_type)
        if not hit:
            return None
        head_name = hit["name"]
        head_type = _normalize_type_alias(override_head_type) or _normalize_type_alias(hit["type"])
        head_node_index = hit.get("node_index")
    else:
        if nodes is None:
            return None
        head_df = _head_candidates_stratified(nodes, context, head_types_scope, per_type=12, total=30)
        if head_df.empty:
            return None
        head_payload = _make_entity_payload(head_df, context=context)
        obj = _chat_json(
            client,
            model=model,
            system=HEAD_PICK_SYS,
            user=_head_pick_user(context, head_payload, allowed_head_types, priors),
            fallback={"head_idx": 0},
        )
        idx = obj.get("head_idx", 0)
        if not isinstance(idx, int) or not (0 <= idx < len(head_payload)):
            idx = 0
        head_sel = head_payload[idx]
        head_name = head_sel["name"]
        head_type = _normalize_type_alias(override_head_type) or _normalize_type_alias(head_sel["type"])
        head_node_index = _safe_int(head_sel["node_index"])

    if head_type not in CANON_ENTITY_TYPES:
        return None

    # -------------------- Step B: RELATION + TAIL --------------------
    triplets = _enumerate_triplet_candidates(head_type, relation_hint, example_guide, allowed_rel)
    triplets = [t for t in triplets if t not in excluded_triples]
    if not triplets:
        return None

    allowed_tail_types = _tail_types_for_head_from_examples(example_guide, head_type)
    tail_priors = [t for t in (type_priors or []) if t in allowed_tail_types]

    triplet_menu = [
        {
            "relation_family": rel,
            "tail_type": _normalize_type_alias(tt),
            "hint": {"prio": (_normalize_type_alias(tt) in tail_priors)},
        }
        for _, rel, tt in triplets
    ]

    tails_by_type: Dict[str, List[dict]] = {}
    if nodes is not None:
        ordered_tail_types = tail_priors + [t for t in allowed_tail_types if t not in tail_priors]
        limited_tail_types = ordered_tail_types[:6]
        for ttn in limited_tail_types:
            tdf = _candidate_rows_by_context(nodes, context, type_hint=ttn, limit=12)
            tails_by_type[ttn] = _make_entity_payload(tdf, context=context) if not tdf.empty else []

    obj2 = _chat_json(
        client,
        model=model,
        system=COMBO_PICK_SYS,
        user=_combo_pick_user(
            context, head_name, head_type, triplet_menu, tails_by_type, allowed_tail_types, tail_priors
        ),
        fallback={"combo_idx": 0, "tail_idx": None},
    )
    ci = obj2.get("combo_idx", 0)
    if not isinstance(ci, int) or not (0 <= ci < len(triplet_menu)):
        ci = 0
    sel_combo = triplet_menu[ci]
    rel_family = sel_combo["relation_family"]
    tail_type = sel_combo["tail_type"]

    tail_name: Optional[str] = None
    tail_node_index: Optional[int] = None
    ti = obj2.get("tail_idx", None)
    tails_list = tails_by_type.get(tail_type, [])
    if isinstance(ti, int) and 0 <= ti < len(tails_list):
        tail_name = tails_list[ti]["name"]
        tail_node_index = _safe_int(tails_list[ti]["node_index"])

    return {
        "head_name": head_name,
        "head_type": head_type,
        "head_node_index": head_node_index,
        "tail_name": tail_name,
        "tail_type": tail_type,
        "tail_node_index": tail_node_index,
        "relation_family": rel_family,
        "excluded_triples": list(excluded_triples),
    }


def _fallback_select(
    *,
    context: str,
    nodes: Optional[pd.DataFrame],
    relation_hint: Optional[str],
    example_guide: List[Dict[str, str]],
    allowed_rel: List[str],
) -> Optional[Dict[str, Any]]:
    """Deterministic fallback: best head candidate + first valid (rel, tail_type)."""
    if nodes is None or nodes.empty:
        return None
    head_df = _candidate_rows_by_context(nodes, context, type_hint=None, limit=1)
    if head_df.empty:
        return None
    r = head_df.iloc[0]
    head_name = str(r.get("node_name"))
    head_type = _normalize_type_alias(str(r.get("node_type", "")))
    head_node_index = _safe_int(r.get("node_index"))

    triplets = _enumerate_triplet_candidates(head_type, relation_hint, example_guide, allowed_rel)
    if not triplets:
        return None
    _, rel_family, tail_type = triplets[0]
    tail_type = _normalize_type_alias(tail_type)

    tail_df = _candidate_rows_by_context(nodes, context, type_hint=tail_type, limit=1)
    if not tail_df.empty:
        tr = tail_df.iloc[0]
        tail_name = str(tr.get("node_name"))
        tail_node_index = _safe_int(tr.get("node_index"))
    else:
        tail_name = None
        tail_node_index = None

    return {
        "head_name": head_name,
        "head_type": head_type,
        "head_node_index": head_node_index,
        "tail_name": tail_name,
        "tail_type": tail_type,
        "tail_node_index": tail_node_index,
        "relation_family": rel_family,
        "excluded_triples": [],
    }


# ==============================================================================
# Retrieval (NaN-safe, dedup-safe)
# ==============================================================================

def _retrieve_topk(
    tgt: torch.Tensor,
    tgt_emb: torch.Tensor,
    tgt_idx: torch.Tensor,
    df_tail: pd.DataFrame,
    topk: Optional[int],
) -> pd.DataFrame:
    """Cosine-similarity retrieval of top-k rows for the projected target vector."""
    # NOTE: logic preserved — 'all' or None => return all items
    topk = tgt_emb.shape[0] if topk in [None, "all"] else min(int(topk), tgt_emb.shape[0])
    cossim = torch.cosine_similarity(tgt, tgt_emb, dim=1)
    top_k = torch.topk(cossim, k=topk, dim=0)
    top_k_map = dict(zip(tgt_idx[top_k.indices].tolist(), top_k.values.tolist()))
    retrieved = df_tail.set_index("node_index").loc[tgt_idx[top_k.indices].cpu().numpy()]
    retrieved["cos_sim"] = retrieved.index.map(top_k_map)
    retrieved["pct_rank"] = retrieved["cos_sim"].rank(pct=True) * np.sign(retrieved["cos_sim"])
    return retrieved.reset_index().sort_values("cos_sim", ascending=False)


# ==============================================================================
# App State & Hydration
# ==============================================================================
@dataclass
class AppState:
    """Runtime container for configs, model, caches, and indices."""

    device: torch.device
    data_config: Optional[Dict[str, Any]]
    inference: Optional[BridgeInference]
    emb_cache: Dict[str, Dict[str, Any]]
    proj_tail_cache: Dict[str, torch.Tensor]
    nodes: Optional[pd.DataFrame]
    name_index: Dict[str, List[int]]
    by_type: Dict[str, pd.DataFrame]
    allowed_rel: List[str]
    example_guide: List[Dict[str, str]]
    init_errors: List[str]
    type_token_vocab: Dict[str, set]


APP_STATE: Optional[AppState] = None


def _safe_load_ckpt(device: torch.device) -> Tuple[Optional[BridgeInference], Optional[str]]:
    """Load model config + weights into a BridgeInference, returning (inference, error)."""
    # Check if checkpoint exists (S3 or local)
    if _is_s3_uri(CKPT_PATH):
        if not _s3_file_exists(CKPT_PATH):
            return None, f"Model checkpoint not found: {CKPT_PATH}"
    elif not os.path.exists(CKPT_PATH):
        return None, f"Model checkpoint not found: {CKPT_PATH}"

    try:
        # Load model config (supports S3)
        if _is_s3_uri(MODEL_CONFIG_PATH):
            model_config = _read_s3_json(MODEL_CONFIG_PATH)
        else:
            with open(MODEL_CONFIG_PATH, "r") as f:
                model_config = json.load(f)

        model = BindingModel(**model_config)

        # Load checkpoint (supports S3)
        if _is_s3_uri(CKPT_PATH):
            ckpt_file = _read_s3_file_object(CKPT_PATH)
            ckpt = torch.load(ckpt_file, map_location=device)
        else:
            ckpt = torch.load(CKPT_PATH, map_location=device)

        model.load_state_dict(ckpt)
        model = model.to(device).eval()
        return BridgeInference(model), None
    except Exception as e:
        return None, f"Checkpoint load failed: {e}"


def _load_embeddings_df(state: AppState, entity_type: str, slidewindow: bool) -> Dict[str, Any]:
    """Load (and cache) embeddings + metadata DataFrame for a given entity type."""
    key = f"{entity_type}|{slidewindow}"
    if key in state.emb_cache:
        return state.emb_cache[key]
    name = _emb_file_name(entity_type, slidewindow)

    # Build paths using S3-aware helper
    emb_pkl = _build_path("data", "embeddings", "esm2b_unimo_pubmedbert", f"{name}.pkl")

    # Check if file exists (S3 or local)
    if _is_s3_uri(emb_pkl):
        if not _s3_file_exists(emb_pkl):
            raise FileNotFoundError(f"Embeddings file not found: {emb_pkl}")
    elif not os.path.exists(emb_pkl):
        raise FileNotFoundError(f"Embeddings file not found: {emb_pkl}")

    raw = _load_pickle(emb_pkl)
    emb = torch.tensor(raw["embedding"], dtype=torch.float32, device=state.device)
    idx = torch.tensor(raw["node_index"], device=state.device)

    df_csv = _build_path("data", "Processed", f"{ENTITY_DICT[entity_type]}.csv")

    # Check if CSV exists (S3 or local)
    if _is_s3_uri(df_csv):
        if not _s3_file_exists(df_csv):
            raise FileNotFoundError(f"Processed CSV not found: {df_csv}")
        df = _read_s3_csv(df_csv, low_memory=False)
    elif not os.path.exists(df_csv):
        raise FileNotFoundError(f"Processed CSV not found: {df_csv}")
    else:
        df = pd.read_csv(df_csv, low_memory=False)

    obj = {"emb": emb, "idx": idx, "raw": raw, "df": df}
    state.emb_cache[key] = obj
    return obj


def _ensure_runtime(state: AppState) -> Optional[str]:
    """Ensure data_config and inference are initialized; return error message or None."""
    if state.data_config is None:
        # Check if config exists (S3 or local)
        if _is_s3_uri(BIND_CONFIG_PATH):
            if not _s3_file_exists(BIND_CONFIG_PATH):
                return f"BindData config not found: {BIND_CONFIG_PATH}"
        elif not os.path.exists(BIND_CONFIG_PATH):
            return f"BindData config not found: {BIND_CONFIG_PATH}"

        try:
            if _is_s3_uri(BIND_CONFIG_PATH):
                state.data_config = _read_s3_json(BIND_CONFIG_PATH)
            else:
                with open(BIND_CONFIG_PATH, "r") as f:
                    state.data_config = json.load(f)
        except Exception as e:
            return f"Failed to read data_config: {e}"
    if state.inference is None:
        inf, err = _safe_load_ckpt(state.device)
        if err:
            return err
        state.inference = inf
    return None


def _load_nodes_and_indices() -> Tuple[Optional[pd.DataFrame], Dict[str, List[int]], Dict[str, pd.DataFrame], List[str]]:
    """Load nodes.csv and build indices; return (nodes, name_index, by_type, init_errors)."""
    init_errors: List[str] = []
    nodes: Optional[pd.DataFrame] = None
    name_index: Dict[str, List[int]] = {}
    by_type: Dict[str, pd.DataFrame] = {}
    try:
        if _is_s3_uri(NODES_CSV_PATH):
            nodes = _read_s3_csv(NODES_CSV_PATH, low_memory=False)
        else:
            nodes = pd.read_csv(NODES_CSV_PATH, low_memory=False)
        _log(f"[hydrate] nodes columns={list(nodes.columns)}")
        name_index, by_type = _build_node_index(nodes)
    except Exception as e:
        init_errors.append(f"Failed to load nodes.csv: {e}")
    return nodes, name_index, by_type, init_errors


def _load_kg_pairs() -> Tuple[List[str], List[Dict[str, str]], Optional[str]]:
    """Load relation families and (head_type, rel, tail_type) examples from kg.csv."""
    allowed_rel: List[str] = []
    example_guide: List[Dict[str, str]] = []

    # Check if file exists (S3 or local)
    if _is_s3_uri(KG_CSV_PATH):
        if not _s3_file_exists(KG_CSV_PATH):
            return allowed_rel, example_guide, f"kg.csv not found at {KG_CSV_PATH}"
    elif not os.path.exists(KG_CSV_PATH):
        return allowed_rel, example_guide, f"kg.csv not found at {KG_CSV_PATH}"

    try:
        if _is_s3_uri(KG_CSV_PATH):
            kg_df = _read_s3_csv(KG_CSV_PATH, low_memory=False)
        else:
            kg_df = pd.read_csv(KG_CSV_PATH, low_memory=False)

        allowed_rel = kg_df.display_relation.unique().tolist()
        example_guide = [
            {"head_type": row["x_type"], "rel": row["display_relation"], "tail_type": row["y_type"]}
            for _, row in kg_df[["x_type", "y_type", "display_relation"]].drop_duplicates().iterrows()
        ]
        return allowed_rel, example_guide, None
    except Exception as e:
        return allowed_rel, example_guide, f"Failed to load kg.csv: {e}"


def _hydrate_app_if_needed(app: AppState) -> None:
    """Lazily populate app state with nodes, relations, and type vocab if missing."""
    if app.nodes is None or not isinstance(app.name_index, dict) or not app.name_index:
        _log("[hydrate] nodes/index missing — loading nodes.csv")
        nodes, name_index, by_type, errs = _load_nodes_and_indices()
        app.nodes, app.name_index, app.by_type = nodes, name_index, by_type
        app.init_errors.extend(errs)
        if nodes is not None:
            _log(f"[hydrate] nodes loaded shape={nodes.shape}; name_index={len(name_index)}")
    if not app.allowed_rel or not app.example_guide:
        _log("[hydrate] relations/example_guide missing — loading kg.csv")
        allowed_rel, example_guide, err = _load_kg_pairs()
        app.allowed_rel, app.example_guide = allowed_rel, example_guide
        if err:
            app.init_errors.append(err)
        _log(f"[hydrate] relations loaded allowed_rel={len(app.allowed_rel)} combos={len(app.example_guide)}")
    if app.nodes is not None and app.by_type and not getattr(app, "type_token_vocab", None):
        app.type_token_vocab = _build_type_token_vocab(app.by_type, max_rows=2000)
        _log(f"[hydrate] type-token vocab built for {len(app.type_token_vocab)} types")
    _log(f"[hydrate] AppState hydrated id={id(app)}")


@asynccontextmanager
async def lifespan(server: FastMCP):
    """Create and register a global AppState for the MCP server lifespan."""
    device = torch.device(DEFAULT_DEVICE)
    state = AppState(
        device=device,
        data_config=None,
        inference=None,
        emb_cache={},
        proj_tail_cache={},
        nodes=None,
        name_index={},
        by_type={},
        allowed_rel=[],
        example_guide=[],
        init_errors=[],
        type_token_vocab={},
    )
    global APP_STATE
    APP_STATE = state
    try:
        yield state
    finally:
        pass


def _build_state_now() -> AppState:
    """Construct a fresh (unhydrated) AppState."""
    device = torch.device(DEFAULT_DEVICE)
    return AppState(
        device=device,
        data_config=None,
        inference=None,
        emb_cache={},
        proj_tail_cache={},
        nodes=None,
        name_index={},
        by_type={},
        allowed_rel=[],
        example_guide=[],
        init_errors=[],
        type_token_vocab={},
    )


def _get_app() -> AppState:
    """Return a singleton AppState, hydrating it if necessary."""
    global APP_STATE
    if APP_STATE is None:
        APP_STATE = _build_state_now()
    _hydrate_app_if_needed(APP_STATE)
    return APP_STATE


mcp.settings.lifespan = lifespan

# ==============================================================================
# Tool: predict_associations
# ==============================================================================
@mcp.tool()
def predict_associations(
    ctx: Context,  # ctx unused but kept for MCP signature compatibility
    *,
    context: str,
    topk: int = 25,
    override_head_name: Optional[str] = None,
    override_head_names: Optional[List[str]] = None,
    override_head_type: Optional[str] = None,
    override_tail_name: Optional[str] = None,
    override_tail_type: Optional[str] = None,
    relation_hint: Optional[str] = None,
    slidewindow: Optional[bool] = None,
    include_relation_catalog: bool = True,
    include_debug: bool = False,
) -> Dict[str, Any]:
    """Main retrieval endpoint: LLM-assisted selection → projection → nearest retrieval.

    Parameters
    ----------
    context : str
        Natural-language context used to infer head/tail candidates and types.
    topk : int
        Desired number of results (overridden by TOPK_HARD_LIMIT if set).
    override_head_name : Optional[str]
        Single head entity name. Mutually exclusive with override_head_names.
    override_head_names : Optional[List[str]]
        Multiple head entity names for mean embedding calculation. Mutually exclusive with override_head_name.
    override_head_type : Optional[str]
        Head entity type (required when using override_head_names).
    override_tail_name : Optional[str]
        Optional explicit tail name override.
    override_tail_type : Optional[str]
        Optional explicit tail type override.
    relation_hint : Optional[str]
        Optional relation family hint.
    slidewindow : Optional[bool]
        Whether to use slidewindow variants for embeddings (default by KG name).
    include_relation_catalog : bool
        Include relation id catalog from BindData config in the response.
    include_debug : bool
        Include internal debug info (type priors, warnings, slidewindow flag).
    """
    try:
        app = _get_app()
        _log(
            f"[tool] using AppState id={id(app)}; nodes_ok={app.nodes is not None}; "
            f"rows={(app.nodes.shape[0] if isinstance(app.nodes, pd.DataFrame) else 0)}; "
            f"name_index={len(app.name_index)}"
        )

        maybe_err = _ensure_runtime(app)
        if maybe_err:
            return {"error": maybe_err, "startup_warnings": app.init_errors}

        data_config = app.data_config
        if data_config is None or app.inference is None:
            return {"error": "runtime not initialized", "startup_warnings": app.init_errors}

        # Validate head name parameters (mutually exclusive)
        if override_head_name and override_head_names:
            return {"error": "Cannot specify both override_head_name and override_head_names", "startup_warnings": app.init_errors}

        # If multiple heads provided, require head_type
        if override_head_names and not override_head_type:
            return {"error": "override_head_type is required when using override_head_names", "startup_warnings": app.init_errors}

        # Honor global cap if configured (keeps existing behavior)
        if TOPK_HARD_LIMIT is not None:
            topk = min(topk, TOPK_HARD_LIMIT)
        if slidewindow is None:
            slidewindow = slidewindow_default

        excluded_triples: set = set()
        device = app.device

        # ---------- Type priors from context ----------
        type_priors = _derive_type_priors(context, app.type_token_vocab, max_types=4)

        # ---------- OPTIMIZED: Fast-path selection when both head and tail are explicit ----------
        sel = None
        use_llm = True
        multi_head_mode = False

        # Fast path: Skip LLM if we have explicit head(s) and tail with types
        if (override_head_name or override_head_names) and override_head_type and override_tail_type:
            _log("[select] fast-path: explicit head and tail types provided, attempting direct resolution")

            # Normalize types
            norm_head_type = _normalize_type_alias(override_head_type)
            norm_tail_type = _normalize_type_alias(override_tail_type)

            # Handle single or multiple heads
            if override_head_names:
                # Multiple heads mode
                multi_head_mode = True
                head_hits = []
                for head_name in override_head_names:
                    hit = _lookup_entity(app.nodes, app.name_index, head_name, norm_head_type)
                    if hit:
                        head_hits.append(hit)
                    else:
                        _log(f"[select] Warning: head entity '{head_name}' not found, skipping")

                if not head_hits:
                    _log("[select] fast-path failed: no valid head entities found")
                    head_hit = None
                else:
                    # Create a combined representation for multiple heads
                    head_names = [h["name"] for h in head_hits]
                    head_node_indices = [h.get("node_index") for h in head_hits]
                    head_hit = {
                        "name": ", ".join(head_names),  # Combined name for display
                        "names": head_names,  # Individual names
                        "type": norm_head_type,
                        "node_index": None,  # Will compute mean embedding
                        "node_indices": head_node_indices,  # Multiple indices for mean
                    }
                    _log(f"[select] fast-path: found {len(head_hits)} head entities: {head_names}")
            else:
                # Single head mode
                head_hit = _lookup_entity(app.nodes, app.name_index, override_head_name, norm_head_type)

            if head_hit:
                # Find valid relation for this head-tail type combination
                valid_combos = [
                    ex for ex in app.example_guide
                    if _normalize_type_alias(ex.get("head_type")) == norm_head_type
                    and _normalize_type_alias(ex.get("tail_type")) == norm_tail_type
                ]

                # Filter by relation hint if provided
                if relation_hint and valid_combos:
                    rel_filtered = [c for c in valid_combos if relation_hint.lower() in c.get("rel", "").lower()]
                    if rel_filtered:
                        valid_combos = rel_filtered

                if valid_combos:
                    # Use first valid combo
                    chosen_combo = valid_combos[0]
                    rel_family = chosen_combo.get("rel")

                    # Optionally lookup tail if name provided
                    tail_name = None
                    tail_node_index = None
                    if override_tail_name:
                        tail_hit = _lookup_entity(app.nodes, app.name_index, override_tail_name, norm_tail_type)
                        if tail_hit:
                            tail_name = tail_hit["name"]
                            tail_node_index = tail_hit.get("node_index")

                    sel = {
                        "head_name": head_hit["name"],
                        "head_type": norm_head_type,
                        "head_node_index": head_hit.get("node_index"),
                        "head_node_indices": head_hit.get("node_indices"),  # For multiple heads
                        "multi_head_mode": multi_head_mode,
                        "tail_name": tail_name,
                        "tail_type": norm_tail_type,
                        "tail_node_index": tail_node_index,
                        "relation_family": rel_family,
                        "excluded_triples": [],
                    }
                    use_llm = False
                    _log(f"[select] fast-path success: head={sel['head_name']}, relation={rel_family}, tail_type={norm_tail_type}")

        # Standard LLM-based selection if fast-path didn't work
        if sel is None and use_llm:
            _log("[select] using LLM-based entity selection")
            sel = _llm_select_combo(
                context=context,
                nodes=app.nodes,
                name_index=app.name_index,
                relation_hint=relation_hint,
                example_guide=app.example_guide,
                allowed_rel=app.allowed_rel,
                excluded_triples=excluded_triples,
                override_head_name=override_head_name,
                override_head_type=override_head_type,
                type_priors=type_priors,
            )

        if sel is None:
            _log("[select] LLM selection failed; using last-ditch fallback.")
            sel = _fallback_select(
                context=context,
                nodes=app.nodes,
                relation_hint=relation_hint,
                example_guide=app.example_guide,
                allowed_rel=app.allowed_rel,
            )

        if sel is None:
            return {
                "error": "Selection failed to produce a usable combination.",
                "startup_warnings": app.init_errors,
            }

        # --------- Unpack selection ----------
        head_name = sel["head_name"]
        sel_head_type = _normalize_type_alias(sel["head_type"])
        sel_tail_type = _normalize_type_alias(sel["tail_type"])
        rel = sel["relation_family"]
        rel_low = (rel or "").lower()
        head_node_index = sel.get("head_node_index")
        head_node_indices = sel.get("head_node_indices")  # For multiple heads
        multi_head_mode = sel.get("multi_head_mode", False)
        tail_name = sel.get("tail_name")  # may be None
        tail_node_index = sel.get("tail_node_index")  # may be None

        if multi_head_mode:
            _log(
                f"[select] resolved {len(head_node_indices) if head_node_indices else 0} heads: {head_name}; "
                f"head_type={sel_head_type}; relation={rel}; tail_type={sel_tail_type}; tail_name={tail_name}"
            )
        else:
            _log(
                f"[select] resolved head={head_name} head_type={sel_head_type} node_index={head_node_index}; "
                f"relation={rel}; tail_type={sel_tail_type}; tail_name={tail_name}"
            )

        # ---------- Ensure chosen combo exists; if not, select next best ----------
        combo_valid = any(
            (
                _normalize_type_alias(ex.get("head_type")) == sel_head_type
                and ex.get("rel") == rel
                and _normalize_type_alias(ex.get("tail_type")) == sel_tail_type
            )
            for ex in app.example_guide
        )
        if not combo_valid:
            prefer = _enumerate_triplet_candidates(
                sel_head_type, relation_hint=None, example_guide=app.example_guide, allowed_rel=app.allowed_rel
            )
            # Exclude the current invalid combo and pick the next
            prefer = [
                c
                for c in prefer
                if not (
                    _normalize_type_alias(c[0]) == sel_head_type
                    and c[1] == rel
                    and _normalize_type_alias(c[2]) == sel_tail_type
                )
            ]
            if not prefer:
                return {
                    "error": "No valid (head_type, relation, tail_type) combos for selected head type.",
                    "startup_warnings": app.init_errors,
                }
            _, rel, sel_tail_type = prefer[0]
            sel_tail_type = _normalize_type_alias(sel_tail_type)
            rel_low = rel.lower()

            # Clear any stale tail selection from the previous (invalid) combo
            tail_name = None
            tail_node_index = None
            _log(
                f"[select-adjusted] head={head_name} head_type={sel_head_type} node_index={head_node_index}; "
                f"relation={rel}; tail_type={sel_tail_type}; tail_name={tail_name}"
            )

        # ---------- Enforce tail_type override if provided ----------
        if override_tail_type:
            otype = _normalize_type_alias(override_tail_type)
            if otype != sel_tail_type:
                combos = _enumerate_triplet_candidates(
                    sel_head_type, relation_hint=None, example_guide=app.example_guide, allowed_rel=app.allowed_rel
                )
                combos = [c for c in combos if _normalize_type_alias(c[2]) == otype]
                if not combos:
                    return {
                        "error": f"Override tail_type '{override_tail_type}' not valid for head_type '{sel_head_type}'.",
                        "startup_warnings": app.init_errors,
                    }
                _, rel, sel_tail_type = combos[0]
                sel_tail_type = _normalize_type_alias(sel_tail_type)
                rel_low = rel.lower()
                # Clear stale tail selection because tail_type changed
                tail_name = None
                tail_node_index = None
                _log(f"[select-adjusted-override] relation={rel}; tail_type={sel_tail_type}; tail_name={tail_name}")

        # ---------- Tail name override (resolve for selected tail type) ----------
        if override_tail_name and sel_tail_type:
            hit_tail = _lookup_entity(app.nodes, app.name_index, override_tail_name, sel_tail_type)
            tail_name = hit_tail["name"] if hit_tail is not None else None
            tail_node_index = hit_tail.get("node_index") if hit_tail is not None else None

        # ---------- Restrict relations to the chosen family ----------
        rel_pairs_exact = [(n, i) for n, i in data_config["relation_type"].items() if str(n).lower() == rel_low]
        rel_pairs = rel_pairs_exact or [
            (n, i) for n, i in data_config["relation_type"].items() if rel_low in str(n).lower()
        ]
        if not rel_pairs:
            return {"error": f"No relation ids matched for '{rel}'.", "startup_warnings": app.init_errors}

        # ---------- Load embeddings/frames ----------
        try:
            head_pkg = _load_embeddings_df(app, sel_head_type, slidewindow)
            tail_pkg = _load_embeddings_df(app, sel_tail_type, slidewindow)
        except Exception as e:
            _log(f"[emb] load failed: {e}")
            return {"error": f"Embeddings load failed: {e}", "startup_warnings": app.init_errors}

        # Head vector(s) - handle single or multiple heads
        head_node_indices_list = sel.get("head_node_indices")
        is_multi_head = sel.get("multi_head_mode", False)

        if is_multi_head and head_node_indices_list:
            # Multiple heads mode: get embeddings for all heads and calculate mean
            _log(f"[emb] multi-head mode: computing mean embedding for {len(head_node_indices_list)} entities")
            uniq_head_idx = np.array([int(idx) for idx in head_node_indices_list if idx is not None])
            if len(uniq_head_idx) == 0:
                return {
                    "error": "No valid head node indices found for multiple heads",
                    "startup_warnings": app.init_errors,
                }
        elif head_node_index is not None:
            # Single head with node_index
            uniq_head_idx = np.array([int(head_node_index)])
        else:
            # Single head without node_index - look up by name
            head_rows = _df_best_match_by_name(head_pkg["df"], head_name)
            if head_rows.empty:
                return {
                    "error": f"Head entity '{head_name}' not found in processed CSV.",
                    "startup_warnings": app.init_errors,
                }
            uniq_head_idx = head_rows["node_index"].unique()

        # Load embeddings and compute average
        pro_node_index = np.array(head_pkg["raw"]["node_index"])
        pro_node_emb = np.array(head_pkg["raw"]["embedding"])

        # Get embeddings for each head entity (handles multi-sequence proteins)
        head_embs = []
        for idx in uniq_head_idx:
            idx_emb = _avg_embeddings_for_nodes(pro_node_index, pro_node_emb, np.array([idx]))
            head_embs.append(idx_emb)

        # Calculate mean across all head entities
        head_avg = np.mean(head_embs, axis=0)
        head_avg_t = torch.tensor(head_avg, dtype=torch.float32, device=device)
        if head_avg_t.ndim == 1:
            head_avg_t = head_avg_t.unsqueeze(0)

        if is_multi_head:
            _log(f"[emb] multi-head mean embedding computed from {len(head_embs)} entities, shape={head_avg_t.shape}")

        # ---------- IDs & projection ----------
        try:
            src_type_id = data_config["node_type"][sel_head_type]
            tgt_type_id = data_config["node_type"][sel_tail_type]
        except KeyError as e:
            return {"error": f"Unknown node type in data_config: {e}", "startup_warnings": app.init_errors}

        proj_key = f"{sel_tail_type}|{slidewindow}|{tgt_type_id}"
        try:
            if proj_key not in app.proj_tail_cache:
                app.proj_tail_cache[proj_key] = app.inference.project(x=tail_pkg["emb"], src_type=tgt_type_id)
            tgt_emb_projected = app.proj_tail_cache[proj_key]
        except Exception as e:
            _log(f"[project] failed: {e}")
            return {"error": f"Projection failed: {e}", "startup_warnings": app.init_errors}

        # ---------- Transform & retrieve ----------
        results: List[Dict[str, Any]] = []
        try:
            # Determine if we should retrieve all results (for tail name filtering) or just topk
            retrieve_all = override_tail_name is not None
            retrieve_topk_param = None if retrieve_all else topk

            for _, rel_id in rel_pairs:
                hrt = app.inference.transform(
                    x=head_avg_t, src_type=src_type_id, tgt_type=tgt_type_id, rel_type=int(rel_id)
                )
                if hrt.ndim == 2 and hrt.shape[0] > 1:
                    hrt = hrt.mean(dim=0, keepdim=True)

                if hrt.device != tgt_emb_projected.device:
                    tgt_emb_projected = tgt_emb_projected.to(hrt.device)
                tail_idx_tensor = tail_pkg["idx"]
                if isinstance(tail_idx_tensor, torch.Tensor) and tail_idx_tensor.device != hrt.device:
                    tail_idx_tensor = tail_idx_tensor.to(hrt.device)

                retrieved = _retrieve_topk(
                    tgt=hrt,
                    tgt_emb=tgt_emb_projected,
                    tgt_idx=tail_idx_tensor,
                    df_tail=tail_pkg["df"],
                    topk=retrieve_topk_param,
                ).sort_values("cos_sim", ascending=False)
                retrieved["relation"] = rel

                # If override_tail_name is specified, filter to matching tail entities
                if override_tail_name:
                    # Use fuzzy/partial matching to find the specific tail
                    tail_name_lower = override_tail_name.lower()
                    name_col = "node_name" if "node_name" in retrieved.columns else "mondo_name"

                    # Try exact match first
                    matched = retrieved[retrieved[name_col].astype(str).str.lower() == tail_name_lower]

                    # If no exact match, try partial/contains match
                    if matched.empty:
                        matched = retrieved[retrieved[name_col].astype(str).str.lower().str.contains(re.escape(tail_name_lower), na=False)]

                    # If still no match, try normalized match
                    if matched.empty:
                        tail_norm = _norm(override_tail_name)
                        matched = retrieved[retrieved[name_col].apply(
                            lambda x: _norm(str(x)).find(tail_norm) >= 0 if pd.notna(x) else False
                        )]

                    if not matched.empty:
                        results.extend(matched.drop_duplicates().to_dict(orient="records"))
                        _log(f"[inference] matched tail '{override_tail_name}' -> '{matched.iloc[0][name_col]}' with pct_rank={matched.iloc[0].get('pct_rank')}")
                    else:
                        _log(f"[inference] no match found for tail name '{override_tail_name}'")
                else:
                    # No tail name specified, return topk results (no pct_rank calculated for this case)
                    results.extend(retrieved.drop_duplicates().to_dict(orient="records"))

        except Exception as e:
            _log(f"[inference] failure ({sel_head_type}, {rel}, {sel_tail_type}): {e}")
            return {"error": f"Inference failure: {e}", "startup_warnings": app.init_errors}

        # Sort by raw score desc; keep logic intact (topk was already applied inside retrieve)
        results = sorted(results, key=lambda x: float(x.get("cos_sim") or 0.0), reverse=True)

        if isinstance(topk, int) and len(results) > topk:
            results = results[:topk]

        if results:
            # target_hit returns selected tail if available else pass top results and drop seq,definition,feature to save token
            target_hit = ([r for r in results if r.get("node_index") == tail_node_index] if tail_node_index is not None
                          else results)
            target_hit = [{k: v for k, v in e_dict.items()
                           if k in ['node_index','node_id','mondo_id','node_name','mondo_name','cos_sim','pct_rank']} for e_dict in target_hit]
            if target_hit:
                _log(
                    f"[results] head={head_name} target={target_hit[0].get('node_name') or target_hit[0].get('mondo_name')} relation={rel} "
                    f"cos={target_hit[0].get('cos_sim')} pct_rank={target_hit[0].get('pct_rank')}"
                )
        else:
            target_hit = None

        resp: Dict[str, Any] = {
            "resolved": {
                "head_name": head_name,
                "head_type": sel_head_type,
                "tail_name": tail_name,
                "tail_type": sel_tail_type,
                "relation_family": rel,
                "chosen_relations": [{"name": rel, "id": rel_pairs[0][1]}] if rel_pairs else [],
                "excluded_triples": list(excluded_triples),
                "llm_used": bool(os.getenv("OPENAI_API_KEY")),
            },
            "results": target_hit,
        }
        if include_relation_catalog:
            resp["relation_catalog"] = _jsonify(data_config.get("relation_type", {}))
        if include_debug:
            resp["debug"] = {"slidewindow": slidewindow, "startup_warnings": app.init_errors, "type_priors": type_priors}
        return resp

    except Exception as e:
        _log(f"[tool] Unhandled exception: {type(e).__name__}: {e}")
        return {"error": f"Unhandled exception: {type(e).__name__}: {e}"}


# ==============================================================================
# Entrypoint
# ==============================================================================
if __name__ == "__main__":
    # IMPORTANT: never print to stdout; MCP uses stdout for protocol
    mcp.run(transport="stdio")
