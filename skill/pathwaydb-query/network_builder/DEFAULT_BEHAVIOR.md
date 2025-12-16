# Default Behavior - Pathway Network Builder

## Updated Defaults ✅

The `pathway_network_builder.py` script now has the following default behavior:

### 1. Databases (Default: ALL THREE)
```bash
--databases kegg msigdb reactome
```
**Included by default:**
- ✅ KEGG human pathways (~370 pathways)
- ✅ MSigDB all collections (~50,000 gene sets)
- ✅ Reactome human pathways (~2,500 pathways)

### 2. MSigDB Collections (Default: ALL)
```bash
--msigdb-collections None  # None means ALL
```
**Collections included by default:**
- ✅ H - Hallmark gene sets (~50)
- ✅ C1 - Positional gene sets (~300)
- ✅ C2 - Curated gene sets (~6,000)
- ✅ C3 - Regulatory target gene sets (~3,700)
- ✅ C4 - Computational gene sets (~850)
- ✅ C5 - Ontology gene sets (~15,000)
- ✅ C6 - Oncogenic signatures (~200)
- ✅ C7 - Immunologic signatures (~5,000)
- ✅ C8 - Cell type signatures (~700)

**Total: ~50,000 MSigDB gene sets**

### 3. Jaccard Index Filter (Default: 0.0)
```bash
--min-jaccard 0.0
```
- ✅ No filtering - includes all edges regardless of similarity

### 4. Organism/Species (Default: Human)
```bash
--kegg-org hsa
--species 9606
```
- ✅ Human for both KEGG and Reactome

### 5. Output File (Default)
```bash
--output pathway_network_edges.csv
```

## What This Means

### Running with minimal arguments:
```bash
python pathway_network_builder.py --output my_network.csv
```

### Is equivalent to:
```bash
python pathway_network_builder.py \
    --databases kegg msigdb reactome \
    --kegg-org hsa \
    --msigdb-collections H C1 C2 C3 C4 C5 C6 C7 C8 \
    --msigdb-version 2025.1.Hs \
    --species 9606 \
    --min-jaccard 0.0 \
    --output my_network.csv
```

## Expected Output Size

With default settings:

### Nodes (Pathways)
- KEGG: ~370 pathways
- MSigDB: ~50,000 gene sets
- Reactome: ~2,500 pathways
- **Total: ~53,000 nodes**

### Edges
- All pairwise comparisons: 53,000 × 52,999 / 2 = **~1.4 billion potential edges**
- With Jaccard > 0: **Millions of edges** (exact number depends on overlaps)

### Execution Time
- KEGG collection: ~40 seconds (API rate limiting)
- MSigDB collection: ~2 minutes (loading all GMT files)
- Reactome collection: ~5 minutes (API rate limiting)
- Network building: ~10-30 minutes (depends on number of edges)
- **Total: ~15-40 minutes for full network**

### File Size
- Edge list CSV: **100 MB - 10+ GB** (depends on min-jaccard threshold)
- With min-jaccard 0.0: Extremely large (may be impractical)
- With min-jaccard 0.05: More manageable (~500 MB - 2 GB)
- With min-jaccard 0.1: Smaller (~50-200 MB)

## Recommendations

### For Production Use:
```bash
# Use Jaccard filtering to reduce output size
python pathway_network_builder.py \
    --min-jaccard 0.05 \
    --output pathway_network.csv
```

### For Testing:
```bash
# Limit pathways per database
python pathway_network_builder.py \
    --max-pathways-per-db 50 \
    --min-jaccard 0.05 \
    --output test_network.csv
```

### For Specific Use Cases:

#### Only MSigDB Hallmark:
```bash
python pathway_network_builder.py \
    --databases msigdb \
    --msigdb-collections H \
    --output hallmark_network.csv
```

#### KEGG + MSigDB (no Reactome):
```bash
python pathway_network_builder.py \
    --databases kegg msigdb \
    --msigdb-collections H C2 \
    --min-jaccard 0.1 \
    --output kegg_msigdb_network.csv
```

#### Only curated pathways (KEGG + MSigDB C2):
```bash
python pathway_network_builder.py \
    --databases kegg msigdb \
    --msigdb-collections C2 \
    --min-jaccard 0.05 \
    --output curated_pathways.csv
```

## Override Defaults

You can override any default by explicitly specifying the parameter:

### Use only one database:
```bash
--databases msigdb
```

### Use specific MSigDB collections:
```bash
--msigdb-collections H C2
```

### Filter by Jaccard threshold:
```bash
--min-jaccard 0.1
```

### Use different organism:
```bash
--kegg-org mmu --species 10090  # Mouse
```

### Limit pathways for testing:
```bash
--max-pathways-per-db 100
```

## Viewing Current Defaults

To see all defaults, run:
```bash
python pathway_network_builder.py --help
```

Look for the "(default: ...)" text in the help output.

## Example Workflows

### Workflow 1: Quick test (1 minute)
```bash
python pathway_network_builder.py \
    --databases msigdb \
    --msigdb-collections H \
    --output quick_test.csv

python analyze_network_example.py quick_test.csv
```

### Workflow 2: Medium test with defaults (5 minutes)
```bash
python pathway_network_builder.py \
    --max-pathways-per-db 100 \
    --min-jaccard 0.05 \
    --output medium_test.csv

python analyze_network_example.py medium_test.csv
```

### Workflow 3: Full comprehensive network (30+ minutes)
```bash
python pathway_network_builder.py \
    --min-jaccard 0.05 \
    --output full_network.csv

# Analyze the large network
python analyze_network_example.py full_network.csv
```

## Performance Considerations

### Memory Usage
- Loading all MSigDB collections: ~500 MB RAM
- Building network structure: ~1-2 GB RAM
- Writing large CSV: Additional disk I/O

### Disk Space
- With min-jaccard 0.0: May need 10+ GB
- With min-jaccard 0.05: ~500 MB - 2 GB
- With min-jaccard 0.1: ~50-200 MB

### CPU Usage
- Single-threaded execution
- Compute-intensive during pairwise comparisons
- CPU usage varies with number of pathways

## Troubleshooting

### Issue: Script runs too long
**Solution:** Use `--max-pathways-per-db 50` for testing

### Issue: Output file too large
**Solution:** Increase `--min-jaccard` to 0.05 or 0.1

### Issue: Out of memory
**Solution:** Process databases separately or reduce pathways

### Issue: Want faster results
**Solution:** Use only one database or specific collections

## Summary

✅ **Default: Comprehensive network from ALL sources**
- All three databases
- All MSigDB collections
- Human organism
- No Jaccard filtering

⚠️ **Production: Add filtering for manageable output**
- Recommend: `--min-jaccard 0.05`
- Consider: `--msigdb-collections` to limit collections
- Optional: `--max-pathways-per-db` for testing

📊 **Expected scale:**
- Nodes: ~53,000 pathways
- Edges: Millions (with filtering)
- Time: 15-40 minutes
- Size: 100 MB - 10+ GB
