# DrugBank Query Examples

This directory contains comprehensive SQL query examples for the DrugBank database, organized by use case.

## Files Overview

### 01-basic-drug-search.sql
**Basic drug lookup queries**
- Search by name or DrugBank ID
- Include synonyms and molecular properties
- Filter by drug state and type
- Count drugs by category

**Use when:** You need to find or identify specific drugs.

### 02-drug-target-interactions.sql
**Drug-target relationship queries**
- Find targets for a drug
- Find drugs targeting a protein/gene
- Identify inhibitors, agonists, antagonists
- Explore enzyme interactions
- Count targets per drug

**Use when:** Analyzing mechanisms of action or finding drugs for specific targets.

### 03-clinical-trials.sql
**Clinical trial data queries**
- Get trials for a drug
- Filter by status (recruiting, completed, etc.)
- Search by trial phase
- Analyze trial timelines
- Generate trial statistics

**Use when:** Researching clinical development status or trial outcomes.

### 04-adverse-effects.sql
**Safety profile queries**
- Get adverse effects with frequency data
- Find common vs rare side effects
- Filter by severity
- Compare AE profiles across drugs
- Include medical coding (SNOMED, MedDRA)

**Use when:** Analyzing drug safety or comparing side effect profiles.

### 05-multi-target-discovery.sql
**Advanced polypharmacology queries**
- Find drugs targeting multiple proteins (AND logic)
- Find drugs targeting protein families (OR logic)
- Identify multi-kinase inhibitors
- Analyze pathway targeting (PI3K, MAPK, etc.)
- Compare target profiles
- Classify selectivity

**Use when:** Discovering multi-target drugs or analyzing drug selectivity.

### 06-indications-and-interactions.sql
**Indications and interaction queries**
- Get approved indications
- Find off-label uses
- Search by medical condition
- Analyze drug-drug interactions
- Find food interactions
- Include medical coding

**Use when:** Researching approved uses or checking for drug interactions.

## Query Pattern Guide

### Simple Lookups
For straightforward queries (single drug, single target), consider using the helper functions in `scripts/drugbank_helper.py`:

```python
from scripts.drugbank_helper import search_drugs, get_drug_targets

# Quick drug search
results = search_drugs("bevacizumab")

# Get targets
targets = get_drug_targets("DB00112")
```

### Complex Analysis
For multi-step queries or complex joins, use the SQL examples directly or adapt them:

```sql
-- Adapt from examples
WITH target_drugs AS (
    -- Use pattern from 02-drug-target-interactions.sql
    ...
)
SELECT ... FROM target_drugs WHERE ...
```

## Common Query Patterns

### 1. Flexible String Matching
Always use `ILIKE` with wildcards for robust searches:
```sql
WHERE d.name ILIKE '%search_term%'
```

### 2. Many-to-Many Relationships
Use `DISTINCT` to avoid duplicate rows:
```sql
SELECT DISTINCT d.drugbank_id, d.name, ...
```

### 3. Multi-Target Queries
Use CTEs (Common Table Expressions) for clarity:
```sql
WITH target1_drugs AS (...),
     target2_drugs AS (...)
SELECT ... WHERE drugbank_id IN (SELECT ... FROM target1_drugs)
  AND drugbank_id IN (SELECT ... FROM target2_drugs)
```

### 4. Aggregation
Use `STRING_AGG` to combine multiple values:
```sql
STRING_AGG(DISTINCT be.name, '; ' ORDER BY be.name) AS target_names
```

### 5. Null Handling
Use `LEFT JOIN` and handle nulls appropriately:
```sql
LEFT JOIN drug_calculated_properties dcp ON d.id = dcp.drug_id
...
ORDER BY aei.percent DESC NULLS LAST
```

## Example Use Cases

### Drug Discovery Workflow

**Step 1: Identify target drugs**
```sql
-- Use 02-drug-target-interactions.sql Example 2
-- Find drugs targeting your protein of interest
```

**Step 2: Analyze multi-target profile**
```sql
-- Use 05-multi-target-discovery.sql Example 6
-- Get complete target profile
```

**Step 3: Check clinical status**
```sql
-- Use 03-clinical-trials.sql Example 1
-- Review clinical trial data
```

**Step 4: Assess safety**
```sql
-- Use 04-adverse-effects.sql Example 1
-- Analyze adverse effect profile
```

### Competitive Analysis

**Compare similar drugs:**
```sql
-- Use 05-multi-target-discovery.sql Example 9
-- Compare target profiles of competing drugs
```

### Safety Assessment

**Comprehensive safety check:**
```sql
-- Use 04-adverse-effects.sql (multiple examples)
-- Then 06-indications-and-interactions.sql Example 13
-- Get complete safety profile including interactions
```

## Tips for Modifying Queries

### 1. Change Drug Selection
Replace the `WHERE` clause:
```sql
-- By name
WHERE d.name ILIKE '%your_drug%'

-- By ID
WHERE d.drugbank_id = 'DB00112'

-- By state
WHERE d.state = 'approved' AND d.name ILIKE '%pattern%'
```

### 2. Change Target Selection
Modify target matching:
```sql
-- By protein name
WHERE be.name ILIKE '%epidermal growth factor receptor%'

-- By gene name
WHERE p.gene_name ILIKE '%EGFR%'

-- By gene list
WHERE p.gene_name IN ('EGFR', 'ERBB2', 'ERBB3')
```

### 3. Add Filters
Restrict results:
```sql
-- Approved only
AND d.state = 'approved'

-- Small molecules only
AND d.type = 'small molecule'

-- Inhibitors only
AND b.inhibitor = true

-- Known pharmacological action
AND b.pharmacological_action IS NOT NULL
```

### 4. Limit Results
For exploratory queries:
```sql
-- Add at end of query
LIMIT 20
```

## Performance Considerations

### Efficient Queries
- Filter early (put restrictive WHERE conditions first)
- Use indexed columns (drugbank_id, gene_name)
- Add LIMIT for exploration

### Avoid
- `SELECT *` on large result sets
- Unnecessary `DISTINCT` (use only when needed)
- Missing WHERE clauses on large tables

## Getting Help

- **Schema details:** See `../references/schema.md`
- **Helper functions:** See `../../scripts/drugbank_helper.py`
- **Main documentation:** See `../SKILL.md`

## Query Execution

### Via MCP (Recommended)
Claude will execute these queries automatically via the PostgreSQL MCP server when analyzing your requests.

### Direct Execution
If running manually:
```bash
psql -h usvgarps11158-dev003.cm9aqaugy64i.us-east-1.rds.amazonaws.com \
     -p 5442 \
     -U your_username \
     -d drugbank \
     -f 01-basic-drug-search.sql
```

## Next Steps

1. **Browse examples** - Review queries in each file
2. **Adapt patterns** - Modify for your specific needs
3. **Combine queries** - Chain multiple queries for complex analysis
4. **Use helpers** - Leverage Python scripts for common tasks

For questions or to report issues, consult the main plugin documentation.
