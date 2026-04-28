# Consensus Annotation Workflow — Complete Code Patterns

## Table of Contents
1. [Setup & Environment](#setup--environment)
2. [CyteType Initialization](#cytetype-initialization)
3. [3 Independent Annotator Agents](#3-independent-annotator-agents)
4. [Extract Per-Cluster Annotations](#extract-per-cluster-annotations)
5. [Reviewer Harmonization](#reviewer-harmonization)
6. [Cache-Busted Re-runs](#cache-busted-re-runs)
7. [Save Results](#save-results)
8. [Config & Log Export](#config--log-export)

---

## Setup & Environment

```python
import csv
import gc
import json
import logging
import os
import sys
import time
from collections import Counter
from pathlib import Path
from uuid import uuid4

from dotenv import load_dotenv
load_dotenv(os.path.expanduser("~/.env"))

import h5py
import numpy as np
import scanpy as sc
from cytetype import CyteType
from openai import OpenAI

# Configuration
OUTPUT_DIR = Path("./cytetype")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

GROUP_KEY = "leiden"           # user-specified
RANK_KEY = "rank_genes_groups"  # user-specified
COORDINATES_KEY = "X_pca_harmony"  # user-specified
CONFIDENCE_THRESHOLD = 0.9
MAX_BATCHES = 10

STUDY_CONTEXT_BASE = "..."  # user-provided biological context

LLM_CONFIGS = [{
    "provider": "openai",
    "name": "gpt-5.2",
    "apiKey": os.environ["OPENAI_API_KEY"],
    "modelSettings": {"reasoning_effort": "high"},
}]

client = OpenAI(
    api_key=os.environ["OPENAI_API_KEY"],
    base_url=os.environ.get("OPENAI_API_BASE", "https://api.openai.com/v1"),
)
```

## CyteType Initialization

```python
adata = sc.read_h5ad(INPUT_PATH)

annotator = CyteType(
    adata,
    group_key=GROUP_KEY,
    rank_key=RANK_KEY,
    n_top_genes=200,
    coordinates_key=COORDINATES_KEY,
    pcent_batch_size=5000,
    max_cells_per_group=1000,
    max_metadata_categories=500,
    vars_h5_path=str(OUTPUT_DIR / "vars.h5"),
    obs_duckdb_path=str(OUTPUT_DIR / "obs.duckdb"),
)
```

## 3 Independent Annotator Agents

**Always cache-bust** with unique run_id per run:

```python
for run_idx in range(1, 4):
    run_id = uuid4().hex[:8]
    prefix = f"round1_run{run_idx}"
    study_context = f"{STUDY_CONTEXT_BASE} [run_id: {run_id}]"

    max_retries = 3
    for attempt in range(1, max_retries + 1):
        try:
            annotator.run(
                study_context=study_context,
                llm_configs=LLM_CONFIGS,
                results_prefix=prefix,
                n_parallel_clusters=5,
                save_query=True,
                query_filename=str(OUTPUT_DIR / f"query_{prefix}.json"),
                override_existing_results=True,
                show_progress=True,
                timeout_seconds=14400,
            )
            break
        except Exception as e:
            logging.warning(f"{prefix} attempt {attempt}/{max_retries} failed: {e}")
            if attempt == max_retries:
                logging.error(f"{prefix} FAILED after {max_retries} attempts")
            else:
                time.sleep(30)
```

## Extract Per-Cluster Annotations

```python
def extract_annotations_from_uns(adata, prefix):
    """Extract per-cluster annotations from CyteType results in adata.uns."""
    results_key = f"{prefix}_results"
    annotations = {}
    if results_key in adata.uns:
        result_data = adata.uns[results_key]
        if isinstance(result_data, dict) and "result" in result_data:
            result_str = result_data["result"]
            result = json.loads(result_str) if isinstance(result_str, str) else result_str
        else:
            result = {"annotations": []}

        for ann in result.get("annotations", []):
            cid = str(ann.get("clusterId", ""))
            annotations[cid] = {
                "run": prefix,
                "annotation": ann.get("annotation", ""),
                "ontologyTerm": ann.get("ontologyTerm", ""),
                "ontologyTermID": ann.get("ontologyTermID", ""),
                "granularAnnotation": ann.get("granularAnnotation", ""),
                "cellState": ann.get("cellState", ""),
                "justification": ann.get("justification", ""),
                "supportingMarkers": ann.get("supportingMarkers", []),
                "conflictingMarkers": ann.get("conflictingMarkers", []),
            }
    return annotations

# Collect all annotations per cluster
clusters = {}  # cluster_id → {"annotations": [...], "n_cells": N, "marker_genes": [...]}
for prefix in ["round1_run1", "round1_run2", "round1_run3"]:
    anns = extract_annotations_from_uns(adata, prefix)
    for cid, ann in anns.items():
        if cid not in clusters:
            clusters[cid] = {"annotations": []}
        clusters[cid]["annotations"].append(ann)
```

## Reviewer Harmonization

Read `reviewer_prompt.md` for the full system prompt and JSON schema.

```python
def call_reviewer(cluster_id, cluster_data):
    """Call GPT-5.2 reviewer for one cluster."""
    # Build prompt (see reviewer_prompt.md for template)
    # ...
    response = client.chat.completions.create(
        model="gpt-5.2",
        messages=[
            {"role": "system", "content": REVIEWER_SYSTEM_PROMPT},
            {"role": "user", "content": user_prompt},
        ],
        temperature=0.0,
        max_completion_tokens=1000,
        response_format={"type": "json_object"},
    )
    result = json.loads(response.choices[0].message.content)
    # Parse confidence (handle string fallback)
    conf = result.get("confidence", 0)
    if isinstance(conf, str):
        try: conf = float(conf)
        except ValueError: conf = {"high": 0.9, "medium": 0.6, "low": 0.3}.get(conf.lower(), 0.5)
    result["confidence"] = conf
    return result

# Harmonize all clusters
harmonized = {}
unresolved = []
for cid in sorted(clusters.keys(), key=lambda x: int(x)):
    result = call_reviewer(cid, clusters[cid])
    if result and result["confidence"] >= CONFIDENCE_THRESHOLD:
        harmonized[cid] = result
    else:
        harmonized[cid] = result  # keep best so far
        unresolved.append(cid)
```

## Cache-Busted Re-runs

For clusters with confidence < threshold:

```python
for batch_num in range(1, MAX_BATCHES + 1):
    if not unresolved:
        break

    # Run 3 CyteType with cache-busting
    for run_idx in range(1, 4):
        run_id = uuid4().hex[:8]
        prefix = f"batch{batch_num}_run{run_idx}"
        study_context = f"{STUDY_CONTEXT_BASE} [run_id: {run_id}]"
        annotator.run(
            study_context=study_context,
            llm_configs=LLM_CONFIGS,
            results_prefix=prefix,
            n_parallel_clusters=5,
            override_existing_results=True,
            timeout_seconds=14400,
        )
        # Append new annotations to cluster history
        new_anns = extract_annotations_from_uns(adata, prefix)
        for cid in unresolved:
            if cid in new_anns:
                clusters[cid]["annotations"].append(new_anns[cid])

    # Re-harmonize unresolved clusters with ALL accumulated annotations
    newly_resolved = []
    for cid in unresolved:
        result = call_reviewer(cid, clusters[cid])
        if result and result["confidence"] >= CONFIDENCE_THRESHOLD:
            harmonized[cid] = result
            newly_resolved.append(cid)
        elif result:
            prev = harmonized.get(cid, {}).get("confidence", 0)
            if result["confidence"] > prev:
                harmonized[cid] = result

    unresolved = [c for c in unresolved if c not in newly_resolved]

# Force-accept remaining
for cid in unresolved:
    harmonized[cid]["needs_review"] = True
```

## Save Results

```python
# Map to adata obs (5 columns only)
adata.obs["cytetype_cluster"] = adata.obs[GROUP_KEY].astype(str).map(
    {cid: h["annotation"] for cid, h in harmonized.items()})
adata.obs["cytetype_ontologyTerm"] = adata.obs[GROUP_KEY].astype(str).map(
    {cid: h["ontologyTerm"] for cid, h in harmonized.items()})
adata.obs["cytetype_ontologyTermID"] = adata.obs[GROUP_KEY].astype(str).map(
    {cid: h.get("ontologyTermID", "") for cid, h in harmonized.items()})
adata.obs["cytetype_cellState"] = adata.obs[GROUP_KEY].astype(str).map(
    {cid: h.get("cellState", "") for cid, h in harmonized.items()})
adata.obs["cytetype_confidence"] = adata.obs[GROUP_KEY].astype(str).map(
    {cid: h.get("confidence", 0) for cid, h in harmonized.items()})

# Save h5ad
adata.write(str(OUTPUT_DIR / "integration_cytetype.h5ad"), compression="gzip")

# Save harmonized annotations CSV
with open(OUTPUT_DIR / "harmonized_annotations.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["cluster_id", "n_cells", "annotation", "ontologyTerm",
                      "ontologyTermID", "cellState", "confidence", "agreement_level",
                      "reasoning", "needs_review", "total_annotations"])
    for cid in sorted(harmonized.keys(), key=lambda x: int(x)):
        h = harmonized[cid]
        writer.writerow([cid, clusters[cid].get("n_cells", 0), h["annotation"],
                          h["ontologyTerm"], h.get("ontologyTermID", ""),
                          h.get("cellState", ""), h["confidence"],
                          h.get("agreement_level", ""), h.get("reasoning", ""),
                          h.get("needs_review", False), len(clusters[cid]["annotations"])])

# Save full JSON
with open(OUTPUT_DIR / "harmonized_annotations.json", "w") as f:
    json.dump(harmonized, f, indent=2, default=str)

# Cleanup
annotator.cleanup()
```

## Config & Log Export

```python
import logging

# Logging setup (at start of script)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler(str(OUTPUT_DIR / "run.log"), mode="w"),
        logging.StreamHandler(sys.stdout),
    ],
)

# Config dump (at start)
config = {
    "input_file": str(INPUT_PATH),
    "output_dir": str(OUTPUT_DIR),
    "group_key": GROUP_KEY,
    "rank_key": RANK_KEY,
    "coordinates_key": COORDINATES_KEY,
    "llm_model": LLM_CONFIGS[0]["name"],
    "reasoning_effort": LLM_CONFIGS[0]["modelSettings"]["reasoning_effort"],
    "confidence_threshold": CONFIDENCE_THRESHOLD,
    "max_batches": MAX_BATCHES,
    "n_clusters": len(clusters),
    "timestamp_start": time.strftime("%Y-%m-%d %H:%M:%S"),
}
with open(OUTPUT_DIR / "config.json", "w") as f:
    json.dump(config, f, indent=2)

# Update config at end with results
config["timestamp_end"] = time.strftime("%Y-%m-%d %H:%M:%S")
config["clusters_resolved"] = sum(1 for h in harmonized.values() if h["confidence"] >= CONFIDENCE_THRESHOLD)
config["clusters_forced"] = sum(1 for h in harmonized.values() if h.get("needs_review"))
```
