# CELLxGENE Census Skill

Comprehensive toolkit for working with CELLxGENE Census single-cell RNA-seq data.

## Structure

```
cellxgene/
├── SKILL.md                              # Main router skill
├── README.md                             # This file
│
├── query/                                # Query subskill
│   ├── SKILL.md                         # Query subskill documentation
│   └── scripts/
│       ├── query_cellxgene.py          # Main query and download tool
│       ├── inspect_metadata_fields.py   # Inspect available metadata
│       └── example_query.py             # Example usage
│
├── specificity/                          # Specificity subskill
│   ├── SKILL.md                         # Specificity subskill documentation
│   └── scripts/
│       ├── extract_markers.py           # Extract marker genes
│       └── visualize_markers.py         # Create visualizations
│
└── references/                           # Reference documentation
    ├── CELLXGENE_QUERY_QUICKSTART.md   # Quick reference for queries
    ├── QUERY_README.md                  # Comprehensive query guide
    ├── METADATA_REFERENCE.md            # Metadata field descriptions
    └── METADATA_FIELDS_ACTUAL.md        # Complete field listing
```

## Usage

The main `SKILL.md` routes to two subskills:

### Query Subskill (`query/`)
**Purpose**: Download expression data from CELLxGENE Census

**Use when**:
- Querying single-cell data with filters
- Downloading expression matrices
- Building custom cohorts
- Extracting metadata

**Example**:
```bash
python query/scripts/query_cellxgene.py \
  --tissue lung \
  --cell-type "T cell" \
  --output lung_tcells
```

### Specificity Subskill (`specificity/`)
**Purpose**: Extract and visualize cell type marker genes

**Use when**:
- Finding marker genes for cell types
- Analyzing gene expression specificity
- Creating marker visualizations
- Comparing markers across contexts

**Example**:
```bash
python specificity/scripts/extract_markers.py \
  --cell-type "B cell" \
  --organism "Homo sapiens" \
  --output bcell_markers

python specificity/scripts/visualize_markers.py \
  --input bcell_markers_computational.csv \
  --output bcell_dotplot.png
```

## Quick Start

See `SKILL.md` for:
- Decision guide (which subskill to use)
- Quick start examples
- Common usage patterns
- Complete reference

See individual subskill SKILL.md files for detailed documentation:
- `query/SKILL.md` - Query and download workflows
- `specificity/SKILL.md` - Marker extraction and visualization

## Reference Documentation

The `references/` directory contains detailed documentation:
- **CELLXGENE_QUERY_QUICKSTART.md** - Copy-paste examples for queries
- **QUERY_README.md** - Comprehensive query usage guide
- **METADATA_REFERENCE.md** - Metadata field descriptions and usage
- **METADATA_FIELDS_ACTUAL.md** - Complete listing of all metadata fields

## Environment

Python 3.7+ with:
- **Query**: `cellxgene-census anndata pandas numpy tiledbsoma`
- **Specificity**: `pandas requests matplotlib seaborn numpy`

Install:
```bash
pip install cellxgene-census anndata pandas numpy requests matplotlib seaborn
```

## Data Source

- **CELLxGENE Census**: https://chanzuckerberg.github.io/cellxgene-census/
- **CellGuide API**: https://cellguide.cellxgene.cziscience.com
