# Workflow 1: Download Reference Models from HuggingFace

Downloads pre-trained Tabula Muris (42 models), Tabula Sapiens (35 models), and Cell Ontology files from HuggingFace.

## Prerequisites

```bash
pip install huggingface_hub
```

## Full Code Pattern

```python
import os
from huggingface_hub import get_collection, snapshot_download
from pathlib import Path

# === User must specify where to save models ===
MODEL_DIR = Path("/path/to/models")  # ALWAYS ask user
BASE_DIR = MODEL_DIR / "huggingface_data"
BASE_DIR.mkdir(parents=True, exist_ok=True)

# --- Utility function ---
def download_collection(collection_url: str, local_subdir: str):
    """Download all models in a HuggingFace collection."""
    collection = get_collection(collection_url)
    target_dir = BASE_DIR / local_subdir
    target_dir.mkdir(parents=True, exist_ok=True)

    print(f"Downloading collection: {collection.title}")
    print(f"Saving to: {target_dir}\n")

    for item in collection.items:
        print(f"Found item: {item.item_id} (type: {item.item_type})")
        if item.item_type == "model":
            try:
                path = snapshot_download(
                    repo_id=item.item_id,
                    local_dir=target_dir / item.item_id.replace('/', '_')
                )
                print(f"Saved to {path}\n")
            except Exception as e:
                print(f"Error downloading {item.item_id}: {e}")


def download_dataset_with_snapshot(repo_id: str, local_subdir: str):
    """Download a HuggingFace dataset (used for ontology files)."""
    target_dir = BASE_DIR / local_subdir
    target_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nDownloading dataset: {repo_id} to {target_dir}")
    try:
        snapshot_download(repo_id=repo_id, repo_type="dataset", local_dir=target_dir)
        print(f"Dataset downloaded to {target_dir}")
    except Exception as e:
        print(f"Failed to download dataset {repo_id}: {e}")


# === Step 1: Download Tabula Muris (42 mouse tissue models) ===
download_collection(
    "popV/tabula-muris-6791cedaf0ecdeb1a8a4840e",
    "tabula-muris"
)

# === Step 2: Download Tabula Sapiens (35 human tissue models) ===
download_collection(
    "popV/tabula-sapiens-67627b2bd44ba09e9129589a",
    "tabula-sapiens"
)

# === Step 3: Download Cell Ontology files ===
download_dataset_with_snapshot("popV/ontology", "ontology")

# === Step 4: Verify structure ===
print("\nFinal structure:")
for root, dirs, files in os.walk(BASE_DIR):
    level = root.replace(str(BASE_DIR), "").count(os.sep)
    indent = "  " * level
    print(f"{indent}- {os.path.basename(root)}/")
    if level < 2:  # Only show first 2 levels to avoid clutter
        for f in files:
            print(f"{'  ' * (level + 1)}- {f}")
```

## Expected Directory Structure After Download

```
{MODEL_DIR}/huggingface_data/
├── tabula-muris/
│   ├── popV_tabula_muris_Adipose_tissue/
│   │   ├── minified_ref_adata.h5ad
│   │   ├── scvi/
│   │   ├── scanvi/
│   │   ├── celltypist.pkl
│   │   ├── xgboost_classifier.model
│   │   ├── svm_classifier.joblib
│   │   ├── OnClass.*
│   │   └── ...
│   ├── popV_tabula_muris_Aorta/
│   └── ... (42 models total)
├── tabula-sapiens/
│   ├── popV_tabula_sapiens_All_Cells/
│   ├── popV_tabula_sapiens_Bladder/
│   └── ... (35 models total)
└── ontology/
    ├── cl.json
    ├── cl.obo
    └── ...
```

## HuggingFace Collection URLs

| Collection | URL | Models |
|-----------|-----|--------|
| Tabula Muris | `popV/tabula-muris-6791cedaf0ecdeb1a8a4840e` | 42 |
| Tabula Sapiens | `popV/tabula-sapiens-67627b2bd44ba09e9129589a` | 35 |
| Ontology | `popV/ontology` (dataset) | — |

## Selective Download

To download a single model instead of the full collection:

```python
from huggingface_hub import snapshot_download

# Download just one model
model_name = "popV/tabula_sapiens_Liver"
target_dir = BASE_DIR / "tabula-sapiens" / model_name.replace('/', '_')
snapshot_download(repo_id=model_name, local_dir=target_dir)
```

## Disk Space

The full set of models is approximately 29 GB. Individual tissue models range from ~100 MB to ~2 GB depending on the number of reference cells. The `All_Cells` models are the largest.
