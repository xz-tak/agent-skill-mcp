# ARCHS4 Filter Reference

**Version:** 1.1
**Last Updated:** 2025-12-17

Complete reference for all available tissue, cell type, disease, and treatment filters in ARCHS4.

---

## Tissue & Cell Type Filters

### Overview

- **Total tissues available:** 72
- **Hierarchical format:** `System.Organ.Cell_Type`
- **Matching methods:** Substring (primary), fuzzy (fallback), smart suggestions (help)
- **Multiple filters:** Supported with OR logic (v1.1)

### Filtering Strategies

Use `tissue_filter` parameter in `gene_expression_analysis()`:

```python
# Single filter
tissue_filter='liver'

# Multiple filters (v1.1)
tissue_filter=['liver', 'hepatocyte', 'bile']

# Partial terms
tissue_filter='hepat'  # Matches HEPATOCYTE, HEPATIC STELLATE CELL

# System-level
tissue_filter='immune'  # Matches all 13 immune system tissues
```

---

## Complete Tissue Catalog (72 tissues)

### Nervous System (14 tissues)

```
System.Nervous System.CNS.CEREBRAL CORTEX
System.Nervous System.CNS.NEURON
System.Nervous System.CNS.ASTROCYTE
System.Nervous System.CNS.OLIGODENDROCYTE
System.Nervous System.CNS.MICROGLIA
System.Nervous System.CNS.CEREBELLUM
System.Nervous System.CNS.HYPOTHALAMUS
System.Nervous System.CNS.THALAMUS
System.Nervous System.CNS.MIDBRAIN
System.Nervous System.CNS.PONS
System.Nervous System.CNS.SPINAL CORD
System.Nervous System.CNS.RETINA
System.Nervous System.PNS.MOTOR NEURON
System.Nervous System.PNS.SENSORY NEURON
```

**Common filters:**
- `'brain'` - All brain tissues
- `'neuron'` - Neural cells
- `'astrocyte'` - Neural support cells
- `'oligodendrocyte'` - Myelin-producing cells
- `'microglia'` - Immune cells of CNS

---

### Immune System (13 tissues)

```
System.Immune System.Lymphoid.TLYMPHOCYTE
System.Immune System.Lymphoid.BLYMPHOCYTE
System.Immune System.Lymphoid.NK CELL
System.Immune System.Thymus.THYMUS
System.Immune System.Thymus.THYMOCYTE
System.Immune System.Spleen.SPLEEN
System.Immune System.Lymph Node.LYMPH NODE
System.Immune System.Myeloid.MONOCYTE
System.Immune System.Myeloid.MACROPHAGE
System.Immune System.Myeloid.DENDRITIC CELL
System.Immune System.Myeloid.PLASMACYTOID DENDRITIC CELL
System.Immune System.Myeloid.NEUTROPHIL
System.Immune System.Bone Marrow.PLASMA CELL
```

**Common filters:**
- `'immune'` - All immune tissues
- `'T cell'` or `'tlymphocyte'` - T cells
- `'B cell'` or `'blymphocyte'` - B cells
- `'NK cell'` - Natural killer cells
- `'monocyte'` - Monocytes
- `'macrophage'` - Macrophages
- `'dendritic'` - Dendritic cells
- `'neutrophil'` - Neutrophils

---

### Digestive System (12 tissues)

```
System.Digestive System.Liver.HEPATOCYTE
System.Digestive System.Liver.LIVER
System.Digestive System.Liver.HEPATIC STELLATE CELL
System.Digestive System.Liver.KUPFFER CELL
System.Digestive System.Pancreas.PANCREATIC ISLET
System.Digestive System.Pancreas.ALPHA CELL
System.Digestive System.Pancreas.BETA CELL
System.Digestive System.Stomach.GASTRIC EPITHELIAL CELL
System.Digestive System.Stomach.STOMACH
System.Digestive System.Intestine.INTESTINAL EPITHELIAL CELL
System.Digestive System.Intestine.SMALL INTESTINE
System.Digestive System.Colon.COLON
```

**Common filters:**
- `'liver'` - Liver tissues
- `'hepatocyte'` - Liver cells
- `'hepatic'` - Hepatic stellate cells
- `'kupffer'` - Liver macrophages
- `'pancreas'` - Pancreatic tissues
- `'beta cell'` - Insulin-producing cells
- `'alpha cell'` - Glucagon-producing cells
- `'intestine'` - Intestinal tissues
- `'colon'` - Colon tissues
- `'stomach'` - Gastric tissues

---

### Urogenital/Reproductive System (9 tissues)

```
System.Urogenital/Reproductive System.Kidney.PODOCYTE
System.Urogenital/Reproductive System.Kidney.RENAL CORTEX
System.Urogenital/Reproductive System.Kidney.KIDNEY
System.Urogenital/Reproductive System.Breast.BREAST
System.Urogenital/Reproductive System.Breast.MAMMARY GLAND
System.Urogenital/Reproductive System.Ovary.OVARY
System.Urogenital/Reproductive System.Ovary.GRANULOSA
System.Urogenital/Reproductive System.Testis.TESTIS
System.Urogenital/Reproductive System.Bladder.BLADDER
```

**Common filters:**
- `'kidney'` or `'renal'` - Kidney tissues
- `'podocyte'` - Kidney filtering cells
- `'breast'` or `'mammary'` - Breast tissues
- `'ovary'` - Ovarian tissues
- `'granulosa'` - Ovarian follicle cells
- `'testis'` - Testicular tissues
- `'bladder'` - Bladder tissues

---

### Connective Tissue (6 tissues)

```
System.Connective Tissue.FIBROBLAST
System.Connective Tissue.ADIPOCYTE
System.Connective Tissue.ADIPOSE
System.Connective Tissue.STROMAL CELL
System.Connective Tissue.CHONDROCYTE
System.Connective Tissue.BONE
```

**Common filters:**
- `'connective'` - All connective tissues
- `'fibroblast'` - Fibroblasts
- `'adipocyte'` or `'adipose'` - Fat cells/tissue
- `'stromal'` - Stromal cells
- `'chondrocyte'` - Cartilage cells
- `'bone'` - Bone tissue

---

### Integumentary System (6 tissues)

```
System.Integumentary System.Skin.KERATINOCYTE
System.Integumentary System.Skin.SKIN
System.Integumentary System.Skin.BASAL CELL
System.Integumentary System.Skin.MELANOCYTE
System.Integumentary System.Hair Follicle.HAIR FOLLICLE
System.Integumentary System.Nail.NAIL
```

**Common filters:**
- `'skin'` - Skin tissues
- `'keratinocyte'` - Skin epithelial cells
- `'basal cell'` - Basal cells
- `'melanocyte'` - Pigment cells
- `'hair'` - Hair follicle

---

### Muscular System (5 tissues)

```
System.Muscular System.SKELETAL MUSCLE
System.Muscular System.SMOOTH MUSCLE
System.Muscular System.MYOBLAST
System.Muscular System.MYOFIBROBLAST
System.Muscular System.AIRWAY SMOOTH MUSCLE
```

**Common filters:**
- `'muscle'` - All muscle tissues
- `'skeletal muscle'` - Skeletal muscle
- `'smooth muscle'` - Smooth muscle
- `'myoblast'` - Muscle precursors
- `'myofibroblast'` - Myofibroblasts

---

### Cardiovascular System (4 tissues)

```
System.Cardiovascular System.HEART
System.Cardiovascular System.ENDOTHELIAL CELL
System.Cardiovascular System.VASCULAR SMOOTH MUSCLE
System.Cardiovascular System.VALVE
```

**Common filters:**
- `'heart'` or `'cardiac'` - Heart tissues
- `'endothelial'` - Blood vessel lining
- `'vascular'` - Vascular tissues
- `'valve'` - Heart valve

---

### Respiratory System (3 tissues)

```
System.Respiratory System.Lung.LUNG EPITHELIAL CELL
System.Respiratory System.Lung.LUNG
System.Respiratory System.Trachea.TRACHEA
```

**Common filters:**
- `'lung'` - Lung tissues
- `'respiratory'` - All respiratory tissues
- `'trachea'` - Tracheal tissues

---

## Common Filter Combinations

Pre-built filter sets for comprehensive tissue coverage:

### Liver Comprehensive
```python
tissue_filter=['liver', 'hepatocyte', 'hepatic', 'kupffer']
```
Matches: LIVER, HEPATOCYTE, HEPATIC STELLATE CELL, KUPFFER CELL

### Immune Comprehensive
```python
tissue_filter=['T cell', 'B cell', 'NK cell', 'lymphocyte', 'monocyte', 'dendritic', 'neutrophil']
```
Matches: All major immune cell types

### Neural Comprehensive
```python
tissue_filter=['brain', 'neuron', 'astrocyte', 'oligodendrocyte', 'neural']
```
Matches: All brain and neural tissues

### Skin Comprehensive
```python
tissue_filter=['skin', 'keratinocyte', 'melanocyte', 'basal']
```
Matches: All skin-related tissues

### Kidney Comprehensive
```python
tissue_filter=['kidney', 'podocyte', 'renal']
```
Matches: All kidney tissues

### Heart Comprehensive
```python
tissue_filter=['heart', 'cardiac', 'endothelial', 'vascular']
```
Matches: All cardiac tissues

### Lung Comprehensive
```python
tissue_filter=['lung', 'respiratory', 'alveolar', 'bronchial']
```
Matches: All lung tissues

---

## Disease & Condition Filters

**Note:** Disease filters require two-step process using `quicksearch_metadata()`.

### Cancer-Related (~77,760 samples)

**Broad terms:**
- `'cancer'` - All cancer types
- `'tumor'` - Tumor samples
- `'carcinoma'` - Epithelial cancers
- `'adenocarcinoma'` - Glandular cancers

**Specific cancer types:**
- `'breast cancer'`
- `'lung cancer'`
- `'liver cancer'`
- `'colon cancer'`
- `'prostate cancer'`
- `'lymphoma'` - Lymphoid cancers
- `'leukemia'` - Blood cancers
- `'melanoma'` - Skin cancer

**Example:**
```python
# Step 1: Find cancer samples
samples = quicksearch_metadata('breast cancer', species='human')
sample_ids = samples.iloc[:, 0].tolist()

# Step 2: Use in analysis
results = gene_correlation('TP53', 'breast cancer', samples=sample_ids[:500])
```

### Neurological Diseases

- `'alzheimer'` - Alzheimer's disease
- `'parkinson'` - Parkinson's disease
- `'multiple sclerosis'`
- `'epilepsy'`
- `'stroke'`

### Metabolic Diseases

- `'diabetes'` - Diabetes-related
- `'obesity'`
- `'metabolic syndrome'`

### Infectious Diseases

- `'covid'` or `'COVID-19'` - COVID-19 related
- `'infection'` - Infectious diseases
- `'viral'` - Viral infections
- `'bacterial'` - Bacterial infections

### Inflammatory Conditions

- `'inflammation'` - Inflammatory conditions
- `'fibrosis'` - Fibrotic conditions
- `'arthritis'`
- `'colitis'`

### Controls

- `'normal'` - Normal/healthy controls
- `'healthy'` - Healthy controls
- `'control'` - Control samples

---

## Treatment Filters

**Note:** Treatment filters require two-step process using `quicksearch_metadata()`.

### General Treatment (~296,703 samples)

- `'treatment'` - Any treatment
- `'treated'` - Treated samples
- `'untreated'` - Untreated/control samples
- `'control'` - Control conditions

### Drug Treatments

- `'drug'` - Drug treatments
- `'inhibitor'` - Inhibitor treatments
- `'compound'` - Compound treatments
- `'chemotherapy'` - Chemotherapy

### Stimulation/Activation

- `'stimulation'` - Stimulated samples
- `'activation'` - Activated samples
- `'lipopolysaccharide'` or `'LPS'` - Immune stimulation

### Genetic Perturbations

- `'knockout'` - Gene knockout
- `'knockdown'` - Gene knockdown
- `'overexpression'` - Gene overexpression
- `'CRISPR'` - CRISPR modifications
- `'siRNA'` - siRNA knockdown
- `'shRNA'` - shRNA knockdown

### Specific Drug Classes

- `'interferon'`
- `'dexamethasone'`
- `'tamoxifen'`
- `'metformin'`

**Example:**
```python
# Find samples treated with a drug
samples = quicksearch_metadata('dexamethasone treatment', species='human')
sample_ids = samples.iloc[:, 0].tolist()

# Analyze gene correlation in treated samples
results = gene_correlation('NR3C1', 'dexamethasone', samples=sample_ids)
```

---

## GEO Series Filters

Query specific GEO series directly:

```python
# Find all samples from a GEO series
samples = quicksearch_metadata('GSE147507', species='human')
```

**Format:** `GSE` followed by series number (e.g., GSE147507, GSE12345)

---

## Smart Filter Suggestions (v1.1)

When filters don't match, system provides intelligent suggestions:

### Example 1: Typo
**Input:** `'livar'`
**Output:**
```
'livar': No matches found
  💡 Did you mean: liver, ovary, muscular?
  💡 Try one of these filters: liver, ovary, muscular, microglia, valve
```

### Example 2: Close Match
**Input:** `'nuerons'`
**Output:**
```
'nuerons': No matches found
  💡 Did you mean: neuron, nervous, pons?
  💡 Try one of these filters: neuron, nervous, pons, motor neuron, granulosa
```

### Example 3: Wrong Term
**Input:** `'pancake'`
**Output:**
```
'pancake': No matches found
  💡 Did you mean: pancreas, melanocyte, macrophage?
  💡 Try one of these filters: pancreas, melanocyte, macrophage, hepatocyte
```

**Suggestion sources:**
- All 72 available tissues
- Extracted terms (system, organ, cell type)
- Curated common filter terms

---

## Filtering Tips

### 1. Start Broad
Begin with organ or system-level filters to see what's available:
```python
tissue_filter='immune'  # See all immune tissues
```

### 2. Use Substrings
Short substrings catch more variations:
```python
tissue_filter='hepat'  # Matches HEPATOCYTE, HEPATIC STELLATE CELL
```

### 3. Combine Related Terms
Use multiple filters for comprehensive coverage:
```python
tissue_filter=['liver', 'hepatocyte', 'bile']
```

### 4. Check Results
Always verify what tissues matched:
```python
results = gene_expression_analysis(genes=['TP53'], tissue_filter='liver')
print(results['top_tissues'])
```

### 5. Iterate
Refine filters based on what you see in results.

---

## Common Mistakes

### 1. Exact Case Matching
❌ Wrong: Expecting exact case
```python
tissue_filter='HEPATOCYTE'  # Works but unnecessary
```

✅ Correct: Use lowercase (case-insensitive)
```python
tissue_filter='hepatocyte'  # Same result
```

### 2. Full Hierarchical Path
❌ Wrong: Using full path
```python
tissue_filter='System.Digestive System.Liver.HEPATOCYTE'
```

✅ Correct: Use simple term
```python
tissue_filter='hepatocyte'  # Substring matching works
```

### 3. Spaces in Cell Types
❌ Wrong: Missing space
```python
tissue_filter='Tcell'  # Won't match TLYMPHOCYTE
```

✅ Correct: Use substring
```python
tissue_filter='lymphocyte'  # Matches TLYMPHOCYTE, BLYMPHOCYTE
```

### 4. Not Using Multiple Filters
❌ Inefficient: Multiple separate queries
```python
# Multiple separate calls
result1 = gene_expression_analysis(genes=['ALB'], tissue_filter='liver')
result2 = gene_expression_analysis(genes=['ALB'], tissue_filter='hepatocyte')
```

✅ Efficient: Single query with multiple filters
```python
# Single call with OR logic
results = gene_expression_analysis(
    genes=['ALB'],
    tissue_filter=['liver', 'hepatocyte', 'bile']
)
```

---

## Performance Notes

- Filtering is fast (< 0.1 seconds)
- Substring matching is primary method
- Fuzzy matching adds negligible overhead
- Smart suggestions computed only when needed

---

## Version History

### v1.1 (2025-12-17)
- ✅ Multiple filter support with OR logic
- ✅ Smart filter suggestions
- ✅ Schema-aware vocabulary

### v1.0 (2025-12-16)
- ✅ 72 tissues catalogued
- ✅ Substring and fuzzy matching
- ✅ Disease and treatment filters documented
