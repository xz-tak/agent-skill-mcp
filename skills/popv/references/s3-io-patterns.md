# S3 I/O Patterns for PopV Pipeline

Use boto3 for all S3 access. Models are stored on S3 and downloaded to `/tmp` for popv to access. Ontology files can be read directly from S3 where possible.

## Default S3 Model Location

```python
S3_MODEL_BASE = "s3://tec-rnd-sci-dev-gi2/gi2-xz/models/popv/data/huggingface_data"
S3_ONTOLOGY = f"{S3_MODEL_BASE}/ontology"
```

## Model Selection

```python
def select_model(organism="human", tissue=None, compartment=None, assay="10x"):
    """Auto-select PopV model based on organism, tissue, and compartment.

    Returns (collection, model_name) tuple.
    IMPORTANT: Always confirm the selection with the user via AskUserQuestion before proceeding.
    """
    if organism == "human":
        collection = "tabula-sapiens"
        if tissue:
            TISSUE_MAP = {
                "large_intestine": "Large_Intestine", "colon": "Large_Intestine",
                "small_intestine": "Small_Intestine", "ileum": "Small_Intestine",
                "gut": "Large_Intestine",
                "liver": "Liver", "lung": "Lung", "blood": "Blood",
                "heart": "Heart", "kidney": "Kidney", "pancreas": "Pancreas",
                "skin": "Skin", "spleen": "Spleen", "stomach": "Stomach",
                "bladder": "Bladder", "bone_marrow": "Bone_Marrow",
                "lymph_node": "Lymph_Node", "mammary": "Mammary",
                "muscle": "Muscle", "thymus": "Thymus", "trachea": "Trachea",
                "eye": "Eye", "ear": "Ear", "fat": "Fat",
                "ovary": "Ovary", "prostate": "Prostate", "testis": "Testis",
                "tongue": "Tongue", "uterus": "Uterus",
                "salivary_gland": "Salivary_Gland", "vasculature": "Vasculature",
            }
            tissue_key = tissue.lower().replace(" ", "_")
            mapped = TISSUE_MAP.get(tissue_key, tissue.replace(" ", "_").title())
            model = f"popV_tabula_sapiens_{mapped}"
        elif compartment:
            model = f"popV_tabula_sapiens_{compartment}"
        else:
            model = "popV_tabula_sapiens_All_Cells"
    elif organism == "mouse":
        collection = "tabula-muris"
        suffix = f"_{assay}" if assay else ""
        if tissue:
            MOUSE_TISSUE_MAP = {
                "large_intestine": "Large_intestine", "colon": "Large_intestine",
                "liver": "Liver", "lung": "Lung", "heart": "Heart",
                "kidney": "Kidney", "pancreas": "Pancreas", "spleen": "Spleen",
                "skin": "Skin_of_body", "bone_marrow": "Bone_marrow",
                "mammary": "Mammary_gland", "muscle": "Limb_muscle",
                "thymus": "Thymus", "trachea": "Trachea", "tongue": "Tongue",
                "bladder": "Bladder_lumen",
            }
            tissue_key = tissue.lower().replace(" ", "_")
            mapped = MOUSE_TISSUE_MAP.get(tissue_key, tissue.replace(" ", "_"))
            model = f"popV_tabula_muris_{mapped}{suffix}"
        else:
            model = f"popV_tabula_muris_All{suffix}" if assay else "popV_tabula_muris_All"
    else:
        raise ValueError(f"Unsupported organism: {organism}. Use 'human' or 'mouse'.")

    return collection, model
```

## Helper Functions

```python
import boto3
import json
import os
import scanpy as sc
from urllib.parse import urlparse
from pathlib import Path

S3_MODEL_BASE = "s3://tec-rnd-sci-dev-gi2/gi2-xz/models/popv/data/huggingface_data"
S3_ONTOLOGY = f"{S3_MODEL_BASE}/ontology"


def is_s3_path(path: str) -> bool:
    """Check if a path is an S3 URI."""
    return str(path).startswith("s3://")


def read_json_from_s3(s3_uri: str) -> dict:
    """Read a JSON file directly from S3 without downloading.

    Use this for ontology files (cl.json, cl_popv.json) that only need
    to be parsed, not stored on disk.
    """
    parsed = urlparse(s3_uri)
    obj = boto3.client("s3").get_object(
        Bucket=parsed.netloc, Key=parsed.path.lstrip("/")
    )
    return json.loads(obj["Body"].read().decode("utf-8"))


def read_h5ad_from_s3(s3_uri: str) -> tuple:
    """Download a single h5ad from S3 to /tmp, return (AnnData, local_path)."""
    parsed = urlparse(s3_uri)
    bucket = parsed.netloc
    key = parsed.path.lstrip("/")
    s3 = boto3.client("s3")
    local_path = f"/tmp/{os.path.basename(key)}"
    print(f"Downloading {s3_uri} -> {local_path}")
    s3.download_file(bucket, key, local_path)
    return sc.read_h5ad(local_path), local_path


def download_from_s3(s3_dir: str, local_dir: str) -> str:
    """Download a directory from S3 to a local path.

    Used for:
    - Model directories → /tmp/popv_models/{MODEL}/
    - Ontology files → /tmp/popv_ontology/ (needed for popv Process_Query cl_obo_folder)
    """
    parsed = urlparse(s3_dir)
    bucket = parsed.netloc
    prefix = parsed.path.lstrip("/")
    if not prefix.endswith("/"):
        prefix += "/"

    s3 = boto3.client("s3")
    paginator = s3.get_paginator("list_objects_v2")

    file_count = 0
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            rel_path = obj["Key"][len(prefix):].lstrip("/")
            if not rel_path:
                continue
            local_file = os.path.join(local_dir, rel_path)
            os.makedirs(os.path.dirname(local_file), exist_ok=True)
            s3.download_file(bucket, obj["Key"], local_file)
            file_count += 1

    print(f"Downloaded {file_count} files from {s3_dir} -> {local_dir}")
    return local_dir


def upload_file_to_s3(local_path: str, s3_uri: str):
    """Upload a single local file to S3."""
    parsed = urlparse(s3_uri)
    bucket = parsed.netloc
    key = parsed.path.lstrip("/")
    print(f"Uploading {local_path} -> {s3_uri}")
    boto3.client("s3").upload_file(local_path, bucket, key)


def upload_directory_to_s3(local_dir: str, s3_uri: str):
    """Upload an entire local directory to S3, preserving structure."""
    parsed = urlparse(s3_uri)
    bucket = parsed.netloc
    prefix = parsed.path.lstrip("/")
    if not prefix.endswith("/"):
        prefix += "/"

    s3 = boto3.client("s3")
    file_count = 0
    for root, dirs, files in os.walk(local_dir):
        for fname in files:
            local_file = os.path.join(root, fname)
            rel_path = os.path.relpath(local_file, local_dir)
            s3_key = prefix + rel_path
            s3.upload_file(local_file, bucket, s3_key)
            file_count += 1

    print(f"Uploaded {file_count} files from {local_dir} -> {s3_uri}")
```

## Integration Pattern — S3 Model + Ontology Access

```python
# === Model: download from S3 to /tmp ===
COLLECTION, MODEL = select_model(organism="human", tissue="large_intestine")
s3_model_path = f"{S3_MODEL_BASE}/{COLLECTION}/{MODEL}"
local_model_dir = f"/tmp/popv_models/{MODEL}"

if not os.path.exists(local_model_dir):
    download_from_s3(s3_model_path, local_model_dir)
else:
    print(f"Model already cached at {local_model_dir}")

model_dir = Path(local_model_dir)
ref_adata_path = model_dir / "minified_ref_adata.h5ad"
pretrained_scvi_path = str(model_dir / "scvi")

# === Ontology for Process_Query: download to /tmp (popv needs local path) ===
local_ontology_dir = "/tmp/popv_ontology"
if not os.path.exists(local_ontology_dir):
    download_from_s3(S3_ONTOLOGY, local_ontology_dir)
else:
    print(f"Ontology already cached at {local_ontology_dir}")

ontology_dir = Path(local_ontology_dir)

# === Ontology cl.json for post-processing: read directly from S3 ===
cl_data = read_json_from_s3(f"{S3_ONTOLOGY}/cl.json")

# ... run standard popV pipeline with model_dir, ontology_dir, cl_data ...
```

## Input/Output with S3

```python
# --- Input h5ad ---
if is_s3_path(query_adata_path):
    query_adata, local_query_path = read_h5ad_from_s3(str(query_adata_path))
else:
    query_adata = sc.read_h5ad(query_adata_path)

# --- Output to S3 ---
if is_s3_path(output_adata_path):
    local_output = f"/tmp/{SHORTNAME}_popv.h5ad"
    query_adata2.write(local_output, compression='gzip')
    upload_file_to_s3(local_output, str(output_adata_path))

    if is_s3_path(str(output_qc_dir)):
        upload_directory_to_s3(str(local_qc_dir), str(output_qc_dir))
```

## S3 Path Examples

| Resource | S3 URI |
|----------|--------|
| Model base | `s3://tec-rnd-sci-dev-gi2/gi2-xz/models/popv/data/huggingface_data` |
| Specific model | `s3://tec-rnd-sci-dev-gi2/gi2-xz/models/popv/data/huggingface_data/tabula-sapiens/popV_tabula_sapiens_Large_Intestine/` |
| Ontology | `s3://tec-rnd-sci-dev-gi2/gi2-xz/models/popv/data/huggingface_data/ontology/` |
| cl.json (direct read) | `s3://tec-rnd-sci-dev-gi2/gi2-xz/models/popv/data/huggingface_data/ontology/cl.json` |

## Notes

- boto3 uses the default AWS credential chain (env vars, ~/.aws/credentials, IAM role)
- Models are downloaded to `/tmp/popv_models/{MODEL}/` — cached for reuse within session
- Ontology is downloaded to `/tmp/popv_ontology/` for `Process_Query(cl_obo_folder=...)`
- `cl.json` for post-processing ontology ID mapping is read directly from S3 (no download needed)
- For SageMaker environments, the IAM role typically provides S3 access automatically
