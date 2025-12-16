# Changelog - Interaction Database Query Skill

## Version 1.1.0 (November 2025)

### BioGRID Multi-Hop Enhancement

**New Features**:
- ✅ Multi-hop expansion (1→2→3 hops) now fully operational for BioGRID queries
- ✅ Automatic neighbor expansion to fill `top_n` requests for genes with sparse direct interactions
- ✅ Extended default timeout to 5 minutes (300 seconds) to accommodate BioGRID API latency
- ✅ Multi-hop works seamlessly in both standalone queries and unified interface

**Example Use Case**:
Gene GREM1 has only 4 direct neighbors. With `max_hops=3` and `top_n=30`:
- **Result**: 4 (1-hop) + 26 (2-hop) = 30 total neighbors found
- **Performance**: Completed in ~2 minutes with 5-minute timeout

**Technical Details**:
- Default timeout increased from 60s to 300s (`biogrid_api.py:60`)
- Multi-hop enabled by default in unified query interface (`unified_query.py:173`)
- BFS algorithm correctly expands through intermediate proteins

**Bug Fixes**:
- Fixed `TypeError` when BioGRID API returns integer `PUBMED_ID` values
  - **Issue**: BioGRID API sometimes returns `PUBMED_ID` as integer instead of string
  - **Impact**: Caused crashes during multi-hop expansion when accumulating publication IDs
  - **Fix**: Added explicit type conversion to string (`biogrid_api.py:295-296, 323-326`)

**Performance Impact**:
- **1-hop queries**: No change (~5-10 seconds)
- **Multi-hop queries**: Expect 1-5 minutes depending on gene connectivity
- **Sparse genes benefit most**: Genes with < 10 direct neighbors now automatically expand

**Validation**:
All changes validated with test suite:
- ✅ GREM1 multi-hop test (sparse gene, 4 direct neighbors)
- ✅ TP53 multi-hop test (highly connected gene)
- ✅ BioGRID shortest path algorithm unchanged
- ✅ No regression in STRING or IntAct functionality

---

## Version 1.0.0 (Initial Release)

### Core Features

**Three Database Support**:
- STRING: Fast, comprehensive, no API key required
- IntAct: Detailed experimental annotations, PSICQUIC integration
- BioGRID: Genetic interactions, manual curation focus

**Two Query Types**:
1. **Single-gene neighbor queries**: Find all interaction partners with comprehensive filtering
2. **Multi-gene shortest paths**: Dijkstra-based path finding between multiple genes

**Comprehensive Filtering**:
- STRING: 8 evidence channels (experimental, database, textmining, coexpression, etc.)
- IntAct: MI-score thresholds, organism filtering, experimental method types
- BioGRID: QUANTITATION scores, experimental system types, throughput tags

**Output Annotations**:
- Entity metadata: Gene names, IDs, organism information
- Edge annotations: Confidence scores, experimental methods, publications
- Path tracking: Full path reconstruction from seed to neighbors (1-hop, 2-hop, 3-hop)
- Algorithm metadata: Weight formulas, distance metrics

**Export Functionality**:
- CSV export for neighbor tables (all databases)
- JSON export for shortest path results
- Consistent format across databases for easy comparison

---

## Migration Guide

### Upgrading from v1.0 to v1.1

**No Breaking Changes** - All existing code continues to work.

**Optional: Leverage BioGRID Multi-Hop**

**Before (v1.0)** - Limited to direct neighbors:
```python
from biogrid_api import BioGRIDClient

client = BioGRIDClient(api_key)
neighbors = client.get_neighbors(
    seed_gene="GREM1",
    tax_id="9606",
    max_hops=1,  # Only direct
    max_neighbors=30
)
# Result: Only 4 neighbors found
```

**After (v1.1)** - Automatic multi-hop expansion:
```python
from biogrid_api import BioGRIDClient

client = BioGRIDClient(api_key)  # Now has 5-minute default timeout
neighbors = client.get_neighbors(
    seed_gene="GREM1",
    tax_id="9606",
    max_hops=3,  # Multi-hop enabled
    max_neighbors=30
)
# Result: 30 neighbors found (4 direct + 26 indirect)
```

**Unified Query Interface** - Now includes BioGRID multi-hop by default:
```python
from unified_query import query_single_gene_all_databases

results = query_single_gene_all_databases(
    gene="GREM1",
    species=9606,
    top_n=30
)

# All three databases now support multi-hop:
# - STRING: ✅ Multi-hop (fast)
# - IntAct: ✅ Multi-hop (moderate)
# - BioGRID: ✅ Multi-hop (5-minute timeout)
```

**Timeout Customization** (if needed):
```python
# For very large queries, increase timeout further
client = BioGRIDClient(api_key, timeout=600)  # 10 minutes

# For quick 1-hop only queries, reduce timeout
client = BioGRIDClient(api_key, timeout=30)   # 30 seconds
```

---

## Known Issues & Limitations

### BioGRID API Performance

**Issue**: BioGRID REST API has higher latency than STRING/IntAct
- **Impact**: Multi-hop queries may take 1-5 minutes
- **Mitigation**: Default 5-minute timeout accommodates most queries
- **Recommendation**: Use STRING for initial exploration, BioGRID for validation

**When BioGRID Multi-Hop Works Best**:
- ✅ Sparse genes (< 20 direct neighbors)
- ✅ max_hops ≤ 3
- ✅ max_neighbors ≤ 100
- ✅ When you can wait 1-5 minutes

**When to Use 1-Hop Only**:
- ⚠️ Highly connected genes (> 100 direct neighbors)
- ⚠️ Need results in < 30 seconds
- ⚠️ Running batch queries (many genes sequentially)

---

## Future Roadmap

**Planned Enhancements**:
- Asynchronous query support for parallel database queries
- Result caching layer for frequently queried genes
- Network visualization export helpers (Cytoscape, Gephi formats)
- Batch query optimization for multiple genes
- Additional evidence type filters across all databases

**Under Consideration**:
- Support for additional databases (MINT, HPRD, DIP)
- Graph database backend option (Neo4j) for complex path queries
- Integration with pathway databases (KEGG, Reactome)
- Machine learning confidence score harmonization across databases
