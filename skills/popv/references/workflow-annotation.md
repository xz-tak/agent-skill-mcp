# Workflow 3: Query Annotation — Full Code Pattern

This is the primary popV workflow. It annotates a query h5ad dataset using pre-trained or custom reference models.

## Table of Contents
1. [Setup & Parameters](#setup--parameters)
2. [Load Query Data](#load-query-data)
3. [Gene ID Mapping](#gene-id-mapping)
4. [Cross-Species Homolog Mapping](#cross-species-homolog-mapping)
5. [Load Reference & Ontology](#load-reference--ontology)
6. [PopV Preprocessing](#popv-preprocessing)
7. [Annotation](#annotation)
8. [Post-Processing](#post-processing)
9. [Ontology ID Mapping](#ontology-id-mapping)
10. [Cluster-Level Majority Vote](#cluster-level-majority-vote)
11. [QC Visualizations](#qc-visualizations)
12. [Save Output](#save-output)

---

## Setup & Parameters

```python
import os
import boto3
import pandas as pd
import popv
import scanpy as sc
import numba
import sklearn
from pathlib import Path
from urllib.parse import urlparse
import matplotlib.pyplot as plt
import json
import mygene

mg = mygene.MyGeneInfo()

# === S3 helpers ===
S3_MODEL_BASE = "s3://tec-rnd-sci-dev-gi2/gi2-xz/models/popv/data/huggingface_data"
S3_ONTOLOGY = f"{S3_MODEL_BASE}/ontology"

def read_json_from_s3(s3_uri):
    """Read JSON file directly from S3 (no download)."""
    parsed = urlparse(s3_uri)
    obj = boto3.client("s3").get_object(Bucket=parsed.netloc, Key=parsed.path.lstrip("/"))
    return json.loads(obj["Body"].read().decode("utf-8"))

def download_from_s3(s3_dir, local_dir):
    """Download a directory from S3 to a local path."""
    parsed = urlparse(s3_dir)
    bucket, prefix = parsed.netloc, parsed.path.lstrip("/")
    if not prefix.endswith("/"): prefix += "/"
    s3 = boto3.client("s3")
    paginator = s3.get_paginator("list_objects_v2")
    file_count = 0
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            rel_path = obj["Key"][len(prefix):].lstrip("/")
            if not rel_path: continue
            local_file = os.path.join(local_dir, rel_path)
            os.makedirs(os.path.dirname(local_file), exist_ok=True)
            s3.download_file(bucket, obj["Key"], local_file)
            file_count += 1
    print(f"Downloaded {file_count} files from {s3_dir} -> {local_dir}")
    return local_dir

# === User-defined parameters (ask user for these; model is auto-selected + confirmed) ===
WORK_DIR = Path(".")                          # User's working directory
DATASET = ""                                  # Subdirectory (can be empty)
SHORTNAME = "my_dataset"                      # Short identifier for file naming
SOURCE = "huggingface_data"                   # Default for S3 models
COLLECTION = "tabula-sapiens"                 # Auto-selected from organism
MODEL = "popV_tabula_sapiens_Liver"           # Auto-selected from tissue, ALWAYS confirmed by user
PREDICTION_MODEL = "inference"                # "fast" or "inference"
N_SAMPLES_PER_LABEL = 300                     # Default, suggested by popV authors

# === Derived paths ===
base_data_dir = WORK_DIR / DATASET if DATASET else WORK_DIR

query_adata_path = base_data_dir / f"{SHORTNAME}.h5ad"
output_qc_dir = base_data_dir / "popV_qc" / SHORTNAME
output_qc_dir.mkdir(parents=True, exist_ok=True)
output_adata_path = base_data_dir / f"{SHORTNAME}_popv.h5ad"

# === Download model from S3 to /tmp ===
s3_model_path = f"{S3_MODEL_BASE}/{COLLECTION}/{MODEL}"
local_model_dir = f"/tmp/popv_models/{MODEL}"
if not os.path.exists(local_model_dir):
    download_from_s3(s3_model_path, local_model_dir)
else:
    print(f"Model already cached at {local_model_dir}")

model_dir = Path(local_model_dir)
ref_adata_path = model_dir / "minified_ref_adata.h5ad"
pretrained_scvi_path = str(model_dir / "scvi")

# === Download ontology to /tmp (popv Process_Query needs local cl_obo_folder) ===
local_ontology_dir = "/tmp/popv_ontology"
if not os.path.exists(local_ontology_dir):
    download_from_s3(S3_ONTOLOGY, local_ontology_dir)
else:
    print(f"Ontology already cached at {local_ontology_dir}")

ontology_dir = Path(local_ontology_dir)
```

## Load Query Data

```python
query_adata = sc.read_h5ad(query_adata_path)
print(query_adata)

# Drop any existing popv columns from prior runs
popv_cols = [col for col in query_adata.obs.columns if col.startswith("popv")]
if popv_cols:
    print("Dropping popv columns from previous annotation:", popv_cols)
    query_adata.obs.drop(columns=popv_cols, inplace=True)

# Set X to raw counts — CRITICAL: popV requires raw integer counts
# The layer name varies by dataset — common options: "raw_counts", "counts", "spliced"
query_adata.X = query_adata.layers["raw_counts"].copy()

# Build batch_key from user-specified columns
batch_list = ["batch"]  # Ask user which columns define batches
query_adata.obs['batch_key'] = query_adata.obs.apply(
    lambda row: "_".join([str(row[ele]) for ele in batch_list]), axis=1
)

# Optional: record which existing cluster column to compare against
cluster_group = "Cell_type"  # Ask user
```

## Gene ID Mapping

PopV reference models use Ensembl gene IDs. If query data uses gene symbols, map them:

```python
# Check if var_names are already Ensembl IDs
if query_adata.var_names[0].startswith("ENSG") or query_adata.var_names[0].startswith("ENSMUSG"):
    print("Gene IDs already in Ensembl format — skipping mapping.")
else:
    # Detect species for mygene query
    species = 'mouse' if COLLECTION in ["Mouse", "tabula-muris"] else 'human'

    geneid_query = mg.querymany(
        query_adata.var_names, scopes='symbol', fields='ensembl.gene', species=species
    )
    geneidmap_dict = {
        res["query"]: res["ensembl"]["gene"]
        if "ensembl" in res and "gene" in res["ensembl"]
        else "not_mapped"
        for res in geneid_query
    }
    query_adata.var["feature_name"] = query_adata.var.index.copy()
    query_adata.var["gene_ids"] = query_adata.var.index.map(geneidmap_dict)

    # Set var_names to Ensembl IDs
    query_adata.var_names = query_adata.var["gene_ids"]
    query_adata = query_adata[:, ~query_adata.var_names.duplicated()]
```

## Cross-Species Homolog Mapping

When annotating mouse data with a human reference model (or vice versa):

```python
if (COLLECTION in ["Human", "tabula-sapiens"]) and query_adata.var_names.str.startswith("ENSMUSG").any():
    homologs_df = pd.read_csv(
        f"{MODEL_DIR}/aux_files/mouse_human_homologs_one_to_one.txt", sep='\t'
    )
    homologs_df = homologs_df.drop_duplicates(subset='Human.gene.name', keep='first')

    genes_to_keep = homologs_df['Gene.stable.ID'].values
    mask = query_adata.var['gene_ids'].isin(genes_to_keep)
    query_adata = query_adata[:, mask].copy()

    homologs_df = homologs_df[homologs_df['Gene.stable.ID'].isin(query_adata.var['gene_ids'])]
    mouse_to_human_id = dict(zip(homologs_df['Gene.stable.ID'], homologs_df['Human.gene.stable.ID']))
    mouse_to_human_symbol = dict(zip(homologs_df['Gene.stable.ID'], homologs_df['Human.gene.name']))

    query_adata.var['human_gene_id'] = query_adata.var['gene_ids'].map(mouse_to_human_id)
    query_adata.var['human_gene_symbol'] = query_adata.var['gene_ids'].map(mouse_to_human_symbol)
    query_adata.var_names = query_adata.var['human_gene_id']
    print("Mouse genes successfully mapped to human homologs.")
else:
    print("No cross-species mapping needed.")
```

## Load Reference & Ontology

```python
# Remove precomputed PCs (if present) — downstream methods recompute these
if "PCs" in query_adata.varm.keys():
    del query_adata.varm["PCs"]

# Load reference
ref_adata = sc.read(ref_adata_path)

# For custom_scds2 models, raw counts are in ref_adata.raw.X
if SOURCE == "custom_scds2":
    ref_adata.X = ref_adata.raw.X
```

## PopV Preprocessing

```python
output_folder = model_dir
os.makedirs(output_folder, exist_ok=True)

ref_labels_key = "cell_type"
unknown_celltype_label = "unassigned"

adata = popv.preprocessing.Process_Query(
    query_adata,
    ref_adata,
    query_batch_key="batch_key",
    ref_labels_key=ref_labels_key,
    ref_batch_key="batch_key",
    unknown_celltype_label=unknown_celltype_label,
    save_path_trained_models=output_folder,
    cl_obo_folder=ontology_dir,
    prediction_mode=PREDICTION_MODEL,
    n_samples_per_label=N_SAMPLES_PER_LABEL,
    pretrained_scvi_path=pretrained_scvi_path,
    hvg=None  # Use all genes; leave as None for "fast" mode
).adata
```

## Annotation

```python
popv.annotation.annotate_data(adata)
```

## Post-Processing

```python
# In "inference" mode, cell count may change — filter to query cells only
if PREDICTION_MODEL != "fast":
    adata = adata[adata.obs_names.isin(query_adata.obs_names), :]

print(f"query adata cell number: {query_adata.n_obs}")
print(f"predicted adata cell number: {adata.n_obs}")

# Diagnostics: OnClass seen vs prediction
obs = adata.obs.copy()
if "popv_onclass_prediction" in obs.columns and "popv_onclass_seen" in obs.columns:
    onclass_same = (obs["popv_onclass_prediction"] == obs["popv_onclass_seen"])
    count_same = onclass_same.sum()
    total = len(obs)
    print(f"OnClass seen == prediction: {count_same}/{total} ({count_same/total:.2%})")

    # In fast mode these should be 100% identical — clean up if so
    if count_same == total:
        adata.obs.drop(columns=["popv_onclass_seen", "popv_onclass_seen_probabilities"],
                       inplace=True, errors="ignore")
        adata.uns["prediction_keys_seen"] = [
            "popv_onclass_prediction" if k == "popv_onclass_seen" else k
            for k in adata.uns["prediction_keys_seen"]
        ]

# Reload original query to copy annotations into a clean adata
query_adata2 = sc.read_h5ad(query_adata_path)
# Re-apply gene ID mapping for later use
query_adata2.var["gene_ids"] = query_adata2.var.index.map(geneidmap_dict)

popv_cols_existing = [col for col in query_adata2.obs.columns if col.startswith("popv")]
if popv_cols_existing:
    query_adata2.obs.drop(columns=popv_cols_existing, inplace=True)

# Copy popv columns
popv_cols = [col for col in adata.obs.columns if col.startswith("popv_")]
if query_adata2.n_obs == adata.n_obs:
    query_adata2.obs[popv_cols] = adata.obs[popv_cols].loc[query_adata2.obs.index]
else:
    for col in popv_cols:
        query_adata2.obs[col] = query_adata2.obs.index.map(adata.obs[col].to_dict())
        if col.endswith(('score', 'probabilities')):
            query_adata2.obs[col] = query_adata2.obs[col].fillna(0)
        elif not col.endswith('depth'):
            query_adata2.obs[col] = query_adata2.obs[col].fillna('unknown_celltype_label')

# Rename prediction columns to *_ontology_name
query_adata2.obs.rename(
    columns={
        c: f"{c}_ontology_name"
        for c in query_adata2.obs.columns
        if (c.startswith("popv") and c.endswith("prediction"))
           or c == "popv_parent"
           or c == "popv_onclass_seen"
    },
    inplace=True
)
```

## Ontology ID Mapping

```python
# Read cl.json directly from S3 (no download needed for this step)
cl_data = read_json_from_s3(f"{S3_ONTOLOGY}/cl.json")

name_to_id = {
    node["lbl"].lower(): node["id"].replace("http://purl.obolibrary.org/obo/", "").replace("_", ":")
    for node in cl_data["graphs"][0]["nodes"]
    if "lbl" in node and "id" in node
}

for col in query_adata2.obs.columns:
    if col.startswith("popv") and col.endswith("ontology_name"):
        new_col = col.replace("ontology_name", "ontology_id")
        mapped_series = query_adata2.obs[col].str.lower().map(name_to_id)
        cols = list(query_adata2.obs.columns)
        idx = cols.index(col)
        query_adata2.obs.insert(loc=idx + 1, column=new_col, value=mapped_series)

        unmatched = query_adata2.obs[col][mapped_series.isna()].unique()
        if len(unmatched) > 0:
            print(f"Unmatched in '{col}': {unmatched}")
```

## Cluster-Level Majority Vote

```python
counts = query_adata2.obs.groupby(cluster_group)['popv_prediction_ontology_name'].value_counts()
majority_vote = counts.groupby(level=0).idxmax().apply(lambda x: x[1])
majorityontology_vote = (
    query_adata2.obs.groupby(cluster_group)['popv_prediction_ontology_id']
    .value_counts().groupby(level=0).idxmax().apply(lambda x: x[1])
)

cluster_sizes = query_adata2.obs.groupby(cluster_group).size()
majority_score = counts.groupby(level=0).max() / cluster_sizes

query_adata2.obs['cluster_popv'] = query_adata2.obs[cluster_group].map(majority_vote)
query_adata2.obs['cluster_popv_id'] = query_adata2.obs[cluster_group].map(majorityontology_vote)
query_adata2.obs['cluster_popv_score'] = query_adata2.obs[cluster_group].map(majority_score)
```

## QC Visualizations

```python
# --- Individual classifier UMAPs ---
popv_prediction_cols = [
    'popv_celltypist_prediction_ontology_name',
    'popv_knn_on_scvi_prediction_ontology_name',
    'popv_onclass_prediction_ontology_name',
    'popv_scanvi_prediction_ontology_name',
    'popv_svm_prediction_ontology_name',
    'popv_xgboost_prediction_ontology_name',
]

for col in popv_prediction_cols:
    sc.pl.umap(query_adata2, color=col, show=False)
    plt.savefig(output_qc_dir / f"{col}_umap.png", bbox_inches="tight", dpi=300)
    plt.show()
    plt.close()

# --- Summary UMAPs ---
summary_cols = [
    "popv_prediction_ontology_id",
    "popv_prediction_ontology_name",
    "popv_prediction_score",
    cluster_group,
    "cluster_popv_id",
    "cluster_popv",
    "cluster_popv_score"
]
for col in summary_cols:
    sc.pl.umap(query_adata2, color=col, show=False)
    plt.savefig(output_qc_dir / f"{col}_umap.png", bbox_inches="tight", dpi=300)
    plt.show()
    plt.close()

# --- Agreement plots ---
popv.visualization.make_agreement_plots(
    adata,
    prediction_keys=adata.uns["prediction_keys"],
    save_folder=output_qc_dir
)

# --- Cluster vs popv agreement ---
popv.visualization.make_agreement_plots(
    query_adata2,
    prediction_keys=[cluster_group, "cluster_popv"],
    popv_prediction_key="popv_prediction_ontology_name",
    save_folder=output_qc_dir
)

# --- Score distribution ---
popv.visualization.prediction_score_bar_plot(
    adata,
    popv_prediction_score="popv_prediction_score",
    save_folder=output_qc_dir
)

# --- Marker gene UMAPs ---
MARKER_GENES = [
    "CD3E", "CD4", "CD8A", "MS4A1", "EPCAM", "PDGFRB", "ACTA2",
    "PTPRC", "PECAM1", "JCHAIN", "NCAM1", "KLRD1", "ZBTB16",
    "TPSAB1", "CD14", "LYZ", "FCGR3A", "CD68", "ITGAX", "FCGR3B",
    "SNAP25", "GFAP", "OLIG1", "OLIG2"
]

is_mouse = query_adata2.var.get("gene_ids", query_adata2.var.index).str.startswith("ENSMUSG").any()
marker_genes = [g.capitalize() for g in MARKER_GENES] if is_mouse else MARKER_GENES
marker_genes = [g for g in marker_genes if g in query_adata2.var_names]

if marker_genes:
    cols_per_row = 4
    n = len(marker_genes)
    rows = (n + cols_per_row - 1) // cols_per_row
    fig, axs = plt.subplots(rows, cols_per_row, figsize=(4 * cols_per_row, 4 * rows))
    axs = axs.flatten()
    for i, gene in enumerate(marker_genes):
        sc.pl.umap(query_adata2, color=[gene], ax=axs[i], show=False, frameon=False)
    for j in range(n, len(axs)):
        axs[j].axis("off")
    plt.tight_layout()
    plt.savefig(output_qc_dir / "marker_genes_umap.png", bbox_inches="tight", dpi=300)
    plt.show()
    plt.close()
```

## Config & Log Export

Always export configuration and execution logs for reproducibility. This section should be included in every PopV annotation run.

### Setup (at the start of every script)

```python
import json
import logging
import sys
import time

# --- Logging: file + stdout ---
log_file = output_qc_dir.parent / "run.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler(str(log_file), mode="w"),
        logging.StreamHandler(sys.stdout),
    ],
)
logging.info(f"PopV annotation starting for {SHORTNAME}")

# --- Config dict: build at start, update throughout ---
run_config = {
    "shortname": SHORTNAME,
    "model": MODEL,
    "collection": COLLECTION,
    "prediction_mode": PREDICTION_MODEL,
    "n_samples_per_label": N_SAMPLES_PER_LABEL,
    "batch_key_columns": batch_list,
    "cluster_group": cluster_group,
    "counts_layer": "raw_counts",   # or whichever layer is used
    "input_file": str(query_adata_path),
    "output_file": str(output_adata_path),
    "python_version": sys.version,
    "popv_version": popv.__version__ if hasattr(popv, "__version__") else "unknown",
    "scanpy_version": sc.__version__,
    "timestamp_start": time.strftime("%Y-%m-%d %H:%M:%S"),
}

# --- Step timing dict ---
step_times = {}
```

### Wrap each pipeline step with timing

```python
t0 = time.time()
# ... run step (e.g., Process_Query) ...
step_times["preprocessing"] = time.time() - t0
logging.info(f"Preprocessing done in {step_times['preprocessing']:.0f}s")
```

### After annotation: log cell counts

```python
logging.info(f"Input cells: {query_adata.n_obs:,}")
logging.info(f"Annotated cells: {adata.n_obs:,}")
logging.info(f"Output cells: {query_adata2.n_obs:,}")
run_config["n_cells_input"] = int(query_adata.n_obs)
run_config["n_cells_annotated"] = int(adata.n_obs)
run_config["n_cells_output"] = int(query_adata2.n_obs)
```

### Finalize config & save (at the end of every script)

```python
run_config["step_times"] = step_times
run_config["total_time_seconds"] = sum(step_times.values())
run_config["timestamp_end"] = time.strftime("%Y-%m-%d %H:%M:%S")
run_config["popv_columns"] = [
    c for c in query_adata2.obs.columns
    if c.startswith("popv") or c.startswith("cluster_popv")
]

config_path = output_qc_dir.parent / "config.json"
with open(config_path, "w") as f:
    json.dump(run_config, f, indent=2, default=str)
logging.info(f"Config saved to {config_path}")
logging.info(f"COMPLETE in {run_config['total_time_seconds']:.0f}s")
```

---

## Save Output

```python
# Store model metadata
query_adata2.uns["popv_model"] = MODEL
for item in cl_data["graphs"][0]["meta"].get("basicPropertyValues", []):
    if item.get("pred") == "http://www.w3.org/2002/07/owl#versionInfo":
        query_adata2.uns["popv_CL_ontology_version"] = item.get("val")
        break

# Save with gzip compression
query_adata2.write(output_adata_path, compression='gzip')
logging.info(f"popV annotated anndata saved to: {output_adata_path}")
```

## Output Columns Reference

The annotated adata will contain these `obs` columns:

| Column Pattern | Description |
|----------------|-------------|
| `popv_{method}_prediction_ontology_name` | Per-classifier cell type name |
| `popv_{method}_prediction_ontology_id` | Per-classifier CL ontology ID |
| `popv_{method}_prediction_probabilities` | Per-classifier confidence |
| `popv_majority_vote_prediction_ontology_name` | Raw majority vote |
| `popv_majority_vote_score` | Majority vote agreement fraction |
| `popv_prediction_ontology_name` | Final ontology-refined prediction |
| `popv_prediction_ontology_id` | Final prediction CL ID |
| `popv_prediction_score` | Final prediction confidence |
| `popv_prediction_depth` | Ontology tree depth |
| `popv_parent_ontology_name` | Parent cell type in ontology |
| `cluster_popv` | Cluster-level majority vote name |
| `cluster_popv_id` | Cluster-level majority vote CL ID |
| `cluster_popv_score` | Cluster-level vote confidence |
