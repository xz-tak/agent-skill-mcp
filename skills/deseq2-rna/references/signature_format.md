# Custom Signature File Format Guide

This guide describes the supported formats for custom gene signature files used with the DESeq2 RNA-seq analysis skill.

## Overview

Custom gene signatures allow you to test your own gene sets in addition to MSigDB collections. The analysis script supports three formats with automatic detection.

## Format Detection

The script automatically detects the format based on file content:

| Indicator | Detected Format |
|-----------|-----------------|
| Tab characters in first line | GMT format |
| Lines starting with `>` | Multi-set text format |
| Plain gene list (no tabs, no `>`) | Single signature format |

---

## Format A: GMT (Gene Matrix Transposed)

GMT is the standard format used by MSigDB and most gene set databases.

### Structure

```
SET_NAME<tab>DESCRIPTION<tab>GENE1<tab>GENE2<tab>GENE3<tab>...
```

Each line represents one gene set with:
1. **Gene set name** (required)
2. **Description** (can be empty, but tab is required)
3. **Gene symbols** (tab-separated list)

### Example GMT File

```
HALLMARK_APOPTOSIS	http://www.gsea-msigdb.org/	CASP1	CASP3	CASP4	CASP6	CASP7	CASP8	BCL2	BAX	BAD
HALLMARK_EMT	http://www.gsea-msigdb.org/	VIM	CDH2	SNAI1	SNAI2	TWIST1	ZEB1	ZEB2	FN1
CUSTOM_TGF_RESPONSE	Custom TGF-beta responsive genes	SERPINE1	CTGF	COL1A1	FN1	ACTA2	TAGLN
CUSTOM_FIBROSIS	Fibrosis signature	COL1A1	COL1A2	COL3A1	COL4A1	FN1	ACTA2	TIMP1
```

### GMT Parsing Rules

1. Lines are split by tab character (`\t`)
2. First field = gene set name
3. Second field = description (can be empty)
4. Remaining fields = gene symbols
5. Empty strings are ignored
6. Whitespace is trimmed from gene symbols

---

## Format B: Multi-Set Text Format

A human-readable format for multiple gene sets with one gene per line.

### Structure

```
>SET_NAME1
GENE1
GENE2
GENE3

>SET_NAME2
GENEA
GENEB
GENEC
```

### Components

| Component | Description |
|-----------|-------------|
| `>` prefix | Marks start of new gene set |
| Set name | Text after `>` (no tabs or newlines) |
| Genes | One gene symbol per line |
| Blank lines | Optional separators between sets |

### Example Multi-Set Text File

```
>TGF_BETA_RESPONSE
SERPINE1
CTGF
COL1A1
FN1
ACTA2
TAGLN
TGFBI

>EMT_MARKERS
VIM
CDH2
SNAI1
SNAI2
TWIST1
ZEB1
ZEB2

>INFLAMMATORY_RESPONSE
IL6
IL8
CXCL1
CXCL2
CCL2
TNF
IL1B
```

---

## Format C: Single Signature (Simple Gene List)

The simplest format - just a list of genes, one per line. Perfect for a single custom signature.

### Structure

```
GENE1
GENE2
GENE3
GENE4
GENE5
```

### Example Single Signature File

```
SERPINE1
CTGF
COL1A1
FN1
ACTA2
TAGLN
TGFBI
COL3A1
LOX
```

### Single Signature Parsing Rules

1. Each non-empty line is a gene symbol
2. Blank lines are ignored
3. Whitespace is trimmed
4. The signature is named "CUSTOM_SIGNATURE" by default
5. No header required

### When to Use Single Signature Format

- You have one specific gene list to test
- Quick ad-hoc analysis with a gene list
- Pasted genes from a publication or database
- Simple pathway or marker gene list

---

## Gene Symbol Requirements

### Supported Identifiers

The analysis uses gene symbols for matching. Supported formats:

| Species | Primary Symbol | Example |
|---------|---------------|---------|
| Human | HGNC symbol | TP53, BRCA1, MYC |
| Mouse | MGI symbol | Trp53, Brca1, Myc |
| Rat | RGD symbol | Tp53, Brca1, Myc |

### Case Sensitivity

- Gene symbols are **case-sensitive**
- Use official nomenclature for your species
- Human genes are typically ALL CAPS
- Mouse/rat genes are typically First-letter capitalized

### Common Issues

| Issue | Solution |
|-------|----------|
| ENSEMBL IDs | Convert to symbols using biomaRt or annotation packages |
| Entrez IDs | Convert to symbols using org.*.eg.db |
| Mixed case | Standardize to species-appropriate case |
| Aliases | Use official current symbols |
| Deprecated symbols | Update to current nomenclature |

---

## Usage in Analysis

### Command-Line

```bash
# With GMT file
conda run -n r_env Rscript scripts/deseq2_analysis.R \
  --custom_sigs /path/to/signatures.gmt \
  ... other arguments ...

# With multi-set text file
conda run -n r_env Rscript scripts/deseq2_analysis.R \
  --custom_sigs /path/to/signatures.txt \
  ... other arguments ...

# With single gene list
conda run -n r_env Rscript scripts/deseq2_analysis.R \
  --custom_sigs /path/to/my_genes.txt \
  ... other arguments ...
```

### JSON Config

```json
{
  "counts": "counts.txt",
  "metadata": "metadata.txt",
  "custom_sigs": "/path/to/signatures.txt",
  "msigdb": "H,C2"
}
```

### Combining with MSigDB

Custom signatures are combined with MSigDB collections:

```bash
--msigdb "H,C2" --custom_sigs my_signatures.gmt
```

This will include:
- MSigDB Hallmarks (H)
- MSigDB Curated gene sets (C2)
- Your custom signatures

---

## Creating Signature Files

### From Excel (GMT)

1. Create spreadsheet with columns: Name, Description, Gene1, Gene2, ...
2. Save as "Tab-delimited Text (.txt)"

### From R (Multi-Set Text)

```r
gene_sets <- list(
  "MY_SIGNATURE" = c("GENE1", "GENE2", "GENE3"),
  "ANOTHER_SET" = c("GENEA", "GENEB")
)

con <- file("signatures.txt", "w")
for (name in names(gene_sets)) {
  writeLines(paste0(">", name), con)
  writeLines(gene_sets[[name]], con)
  writeLines("", con)
}
close(con)
```

### From R (Single Signature)

```r
my_genes <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5")
writeLines(my_genes, "my_signature.txt")
```

### From Python (Single Signature)

```python
genes = ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"]
with open("my_signature.txt", "w") as f:
    f.write("\n".join(genes))
```

---

## Validation Checklist

Before using custom signatures, verify:

- [ ] File is plain text (not Word, not PDF)
- [ ] Correct format (GMT, multi-set text, or single signature)
- [ ] Gene symbols match your species
- [ ] No duplicate gene set names (for multi-set formats)
- [ ] At least 5 genes per set (recommended minimum)
- [ ] No more than 500 genes per set (recommended maximum)
- [ ] No special characters in gene set names

---

## Quick Examples

### Single Signature (Easiest)

File: `fibrosis_genes.txt`
```
COL1A1
COL3A1
ACTA2
FN1
CTGF
SERPINE1
TIMP1
LOX
```

### Multi-Set Text

File: `pathway_signatures.txt`
```
>WNT_TARGETS
AXIN2
MYC
CCND1
LEF1

>NOTCH_TARGETS
HES1
HEY1
NRARP

>HEDGEHOG_TARGETS
GLI1
PTCH1
HHIP
```

### GMT Format

File: `signatures.gmt`
```
WNT_TARGETS	Wnt pathway targets	AXIN2	MYC	CCND1	LEF1
NOTCH_TARGETS	Notch pathway targets	HES1	HEY1	NRARP
HEDGEHOG_TARGETS	Hedgehog targets	GLI1	PTCH1	HHIP
```

---

## Recommendations

| Use Case | Recommended Format |
|----------|-------------------|
| Single gene list from paper | Single signature |
| Multiple related signatures | Multi-set text |
| Integration with databases | GMT |
| Quick analysis | Single signature |
| Permanent signature collection | GMT |
