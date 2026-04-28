# Workflow 2: Create Custom Models from CellxGene Census

Train a custom popV reference model using a CellxGene Census dataset. The model trains all classifiers from scratch using `prediction_mode="retrain"`.

## Prerequisites

```bash
pip install cellxgene-census popv scanpy
```

## Full Code Pattern

```python
import pandas as pd
import popv
import scanpy as sc
import os
from pathlib import Path
import json
import cellxgene_census

# === User-defined parameters (ALWAYS ask for all of these) ===
MODEL_DIR = Path("/path/to/models")             # ALWAYS ask user
SOURCE = "custom_scds2"                          # Fixed — do not change
COLLECTION = "Human"                             # "Human", "Mouse", or "NHP"
TISSUE = "Heart"                                 # Tissue group name
MODEL = "popV_CellxGene_Human_Author_2022_Heart" # Descriptive model name
REF_DATASET_ID = "5500c673-1610-40a0-..."        # CellxGene dataset UUID
N_SAMPLES_PER_LABEL = 300                        # Default, suggested by popV authors
HVG = 4000                                       # Highly variable genes count

# === Derived paths ===
model_dir = MODEL_DIR / SOURCE / COLLECTION / TISSUE / MODEL
model_dir.mkdir(parents=True, exist_ok=True)

ontology_dir = MODEL_DIR / SOURCE / "ontology"
query_dummy_dir = MODEL_DIR / SOURCE / "query_dummy"

# === Step 1: Load query dummy ===
# A small dummy dataset is used as the "query" during model training
query_adata = sc.read(f"{query_dummy_dir}/{COLLECTION}_query_dummy.h5ad")
query_adata.X = query_adata.raw.X  # CellxGene format: counts in raw.X

# === Step 2: Download reference dataset from CellxGene Census ===
cellxgene_census.download_source_h5ad(
    dataset_id=REF_DATASET_ID,
    census_version="stable",
    to_path=f"{model_dir}/{MODEL}.h5ad"
)

# === Step 3: Load and prepare reference ===
ref_adata = sc.read(f"{model_dir}/{MODEL}.h5ad")
ref_adata.X = ref_adata.raw.X  # Raw counts from CellxGene format

# Verify raw integer counts
print(ref_adata.X)  # Should show integer values
print(ref_adata)
print(ref_adata.obs['tissue'].value_counts())
print(ref_adata.obs['cell_type'].value_counts())

# === Step 4: Optional tissue filtering ===
tissues_to_keep = []  # e.g., ["heart left ventricle", "heart right ventricle"]
tissue_filter_applied = False
if tissues_to_keep and any(t.strip() for t in tissues_to_keep):
    tissue_filter_applied = True
    ref_adata = ref_adata[ref_adata.obs['tissue'].isin(tissues_to_keep)]
ref_adata.uns["tissues_to_keep"] = tissues_to_keep

# === Step 5: Cell type filtering ===
cells_to_discard = ["unknown", "cell", "native cell", "animal cell"]
cell_filter_applied = ref_adata.obs['cell_type'].isin(cells_to_discard).any()
ref_adata = ref_adata[~ref_adata.obs['cell_type'].isin(cells_to_discard)]

# Save filtered version if either filter was applied
if cell_filter_applied or tissue_filter_applied:
    ref_adata.write(f"{model_dir}/{MODEL}.h5ad", compression="gzip")
    print(f"Filtered .h5ad saved to {model_dir}/{MODEL}.h5ad")

# === Step 6: Train popV model ===
popv.settings.n_jobs = 48  # Adjust to server resources

ref_labels_key = "cell_type"
unknown_celltype_label = "unassigned"

# Build batch keys
query_adata.obs['batch_key'] = query_adata.obs.apply(
    lambda row: row['donor_id'] + '_' + row['assay'] + '_' + row['tissue'], axis=1
)
ref_adata.obs['batch_key'] = ref_adata.obs.apply(
    lambda row: row['donor_id'] + '_' + row['assay'] + '_' + row['tissue'], axis=1
)

# Process and train — this is the most time-consuming step
adata = popv.preprocessing.Process_Query(
    query_adata,
    ref_adata,
    query_batch_key="batch_key",
    ref_labels_key=ref_labels_key,
    ref_batch_key="batch_key",
    unknown_celltype_label=unknown_celltype_label,
    save_path_trained_models=model_dir,
    cl_obo_folder=ontology_dir,
    prediction_mode="retrain",  # Train all classifiers from scratch
    n_samples_per_label=N_SAMPLES_PER_LABEL,
    hvg=HVG
).adata

# === Step 7: Annotate and save model ===
popv.annotation.annotate_data(adata, save_path=f"{model_dir}")

# Clean up predictions file (not needed for the model)
(model_dir / "predictions.csv").unlink(missing_ok=True)

# === Step 8: Generate model summary ===
import sys
sys.path.insert(0, "/home/sagemaker-user/popv")
from generate_model_summaries import summarize_models
summarize_models(model_dir)
```

## Model Naming Convention

Suggested format: `popV_CellxGene_{Organism}_{Author}_{Year}_{Tissue}`

Examples:
- `popV_CellxGene_Human_Knight-Schrijver_2022_Heart`
- `popV_CellxGene_Mouse_Muraro_2016_Pancreas`
- `popV_CellxGene_Human_Tabula_2022_Lung`

## Finding CellxGene Dataset IDs

Browse the [CellxGene Discover portal](https://cellxgene.cziscience.com/) to find datasets. The dataset ID is the UUID in the dataset URL.

## Custom Model Directory Structure

```
{MODEL_DIR}/custom_scds2/
├── {COLLECTION}/
│   └── {TISSUE}/
│       └── {MODEL}/
│           ├── {MODEL}.h5ad          # Reference anndata
│           ├── scvi/                 # scVI model
│           ├── scanvi/               # ScanVI model
│           ├── celltypist.pkl        # CellTypist model
│           ├── xgboost_classifier.model
│           ├── svm_classifier.joblib
│           ├── OnClass.*             # OnClass model files
│           ├── model.txt             # Generated summary
│           └── ...
├── ontology/                         # Must be present
│   ├── cl.json
│   └── cl.obo
└── query_dummy/
    ├── Human_query_dummy.h5ad
    ├── Mouse_query_dummy.h5ad
    └── NHP_query_dummy.h5ad
```

## Training Time Estimates

Training time depends on the reference dataset size and available hardware:
- Small dataset (<10k cells): ~30 min
- Medium dataset (10-50k cells): ~1-3 hours
- Large dataset (>50k cells): ~3-8+ hours

GPU availability significantly speeds up scVI/ScanVI training.
