---
name: database-gostar
description: Use to access and query the GOSTAR medicinal chemistry database for drug discovery. GOSTAR (Global Online Structure Activity Relationship) contains structure-activity relationship (SAR) data, bioassay/bioactivity information, and pharmacological data.
---

# GOSTAR Database

## Database Connection

- **Host**: usvgarps11158-dev003.cm9aqaugy64i.us-east-1.rds.amazonaws.com
- **Port**: 5442
- **Database**: gostar
- **User**: gostar_ro (read-only access)
- **Password**: Stored in AWS Secrets Manager as `gostar/ro-password` (or environment variable `GOSTAR_RO_PASSWORD` for GitHub Actions)

## Prerequisites

Install dependencies using pixi:
```bash
pixi install
```

Key dependencies:
- `psycopg2-binary` - PostgreSQL adapter for direct connections
- `polars` - High-performance data processing
- `connectorx` - Database connectivity for Polars
- `pyarrow` - Data operations for Polars
- `pandas` - Available for compatibility (if needed)
- `sqlalchemy` - SQL toolkit

## Usage Examples

### Basic Connection

```python
import polars as pl
from get_gostar_credentials import get_gostar_connection_string

# Get connection string (retrieves password from AWS Secrets Manager or environment)
conn_string = get_gostar_connection_string()

# Example query - note: no semicolon at the end
query = "SELECT * FROM information_schema.tables WHERE table_schema = 'v2026_01_client' LIMIT 10"

# Execute query and get results as Polars DataFrame
df = pl.read_database_uri(query, uri=conn_string)
print(df)
```

### Query with Polars

```python
import polars as pl
from get_gostar_credentials import get_gostar_connection_string

# Get connection string
conn_string = get_gostar_connection_string()

# Query directly into Polars DataFrame
# Note: Do not include semicolons at the end of SQL queries
df = pl.read_database_uri("SELECT * FROM v2026_01_client.all_activity_gostar LIMIT 100", uri=conn_string)
print(df.head())

# Polars provides efficient data manipulation for large datasets
# Example: filtering and aggregation
result = df.filter(pl.col("activity_value") > 0.5).group_by("target_name").agg([
    pl.len().alias("count"),
    pl.mean("activity_value").alias("mean_activity")
])
print(result)
```

### Connection Helper Functions

Use the provided credential helper module (`get_gostar_credentials.py`) which automatically retrieves the password from AWS Secrets Manager:

```python
from get_gostar_credentials import get_gostar_connection_params, get_gostar_connection_string

# For psycopg2 direct connections
import psycopg2
conn_params = get_gostar_connection_params()
conn = psycopg2.connect(**conn_params)
cursor = conn.cursor()
cursor.execute("SELECT * FROM v2026_01_client.structure_master LIMIT 5")
conn.close()

# For Polars (recommended for data analysis)
import polars as pl
conn_string = get_gostar_connection_string()
df = pl.read_database_uri("SELECT * FROM v2026_01_client.structure_master LIMIT 10", uri=conn_string)
```

**Credential Priority:**
1. Environment variable `GOSTAR_RO_PASSWORD` (GitHub Actions/CI)
2. AWS Secrets Manager `gostar/ro-password` (EC2/local development)
```

## Common Queries

### List All Tables

```python
import polars as pl
from get_gostar_credentials import get_gostar_connection_string

conn_string = get_gostar_connection_string()

query = """
SELECT table_schema, table_name
FROM information_schema.tables
WHERE table_schema NOT IN ('pg_catalog', 'information_schema')
ORDER BY table_schema, table_name
"""
df = pl.read_database_uri(query, uri=conn_string)
print(df)
```

### Get Table Schema

```python
import polars as pl
from get_gostar_credentials import get_gostar_connection_string

conn_string = get_gostar_connection_string()

query = """
SELECT column_name, data_type, character_maximum_length
FROM information_schema.columns
WHERE table_schema = 'v2026_01_client' AND table_name = 'structure_master'
ORDER BY ordinal_position
"""
df = pl.read_database_uri(query, uri=conn_string)
print(df)
```

### Search for Compounds by Structure

```python
import polars as pl
from get_gostar_credentials import get_gostar_connection_string

conn_string = get_gostar_connection_string()

# Search for compounds with activity against kinase targets
query = """
SELECT DISTINCT
    sm.smiles,
    sm.mol_formula,
    sm.mol_weight,
    aag.standard_value,
    aag.standard_type
FROM v2026_01_client.all_activity_gostar aag
JOIN v2026_01_client.structure_gvkid sg ON aag.gvk_id = sg.gvk_id
JOIN v2026_01_client.structure_master sm ON sg.str_id = sm.str_id
WHERE LOWER(aag.target_name) LIKE '%kinase%'
    AND sm.smiles IS NOT NULL
LIMIT 100
"""
df = pl.read_database_uri(query, uri=conn_string)

# Filter and sort the results using Polars
filtered = df.filter(pl.col("standard_value") > 1.0).sort("standard_value", descending=True)
print(filtered)
```

### Retrieve Compounds by Disease Indication

Find compounds for a specific disease (e.g., Acute Lymphoblastic Leukemia):

```python
import polars as pl
from get_gostar_credentials import get_gostar_connection_string

conn_string = get_gostar_connection_string()

query = """
SELECT
    dim.drug_name,
    dim.status,
    dim.highest_phase_all,
    sm.smiles,
    sm.mol_formula,
    sm.mol_weight
FROM v2026_01_client.disease_indication di
JOIN v2026_01_client.drug_information_master dim ON di.gvk_id = dim.gvk_id
LEFT JOIN v2026_01_client.structure_gvkid sg ON di.gvk_id = sg.gvk_id
LEFT JOIN v2026_01_client.structure_master sm ON sg.str_id = sm.str_id
WHERE LOWER(di.disease_indication) LIKE '%acute lymphoblastic leuk%'
ORDER BY dim.highest_phase_all, dim.drug_name
LIMIT 100
"""

df = pl.read_database_uri(query, uri=conn_string)
print(f"Found {len(df)} compounds")
print(df.head())
```

**Key joins:** `disease_indication` → `drug_information_master` (via `gvk_id`) → `structure_gvkid` → `structure_master` (via `str_id`)

## Best Practices

1. **Read-Only Access**: The provided credentials are read-only. No INSERT, UPDATE, or DELETE operations are permitted.

2. **Connection Management**: Always close connections when done:
   ```python
   from get_gostar_credentials import get_gostar_connection_params
   import psycopg2

   try:
       conn = psycopg2.connect(**get_gostar_connection_params())
       # Do work
   finally:
       conn.close()
   ```

3. **Query Optimization**: Use LIMIT clauses for exploratory queries to avoid loading large datasets.

4. **Pagination**: For large result sets, use OFFSET and LIMIT:
   ```python
   query = "SELECT * FROM table LIMIT 1000 OFFSET 0;"
   ```

5. **Connection Pooling**: For multiple queries, reuse connections or use SQLAlchemy's connection pooling.

## Security Notes

- Password is stored in AWS Secrets Manager as `gostar/ro-password`
- For GitHub Actions, use the `GOSTAR_RO_PASSWORD` environment variable
- This is a read-only user account
- Database is accessible only from authorized networks
- Always use the credential helper functions (`get_gostar_connection_string()` or `get_gostar_connection_params()`) to retrieve credentials

## Troubleshooting

### Connection Timeout

If you encounter connection timeouts:
```python
conn_params['connect_timeout'] = 30
```

### SSL Connection

If SSL is required:
```python
conn_params['sslmode'] = 'require'
```

### Network Access

Ensure your IP address has access to the RDS instance. Contact your database administrator if you cannot connect.

## Use Cases

- **SAR Analysis**: Query structure-activity relationships for drug discovery
- **Bioactivity Data**: Access pharmacological and biochemical assay data
- **Target Analysis**: Find compounds active against specific biological targets
- **Chemical Space Exploration**: Analyze compound libraries and chemical diversity
- **Lead Optimization**: Track compound modifications and potency improvements
- **Literature Data Mining**: Access curated data from medicinal chemistry publications

## Related Skills

- Skill(scientific-skills:chembl-database) - Alternative bioactivity database
- Skill(scientific-skills:pubchem-database) - Chemical structures and properties
- Skill(scientific-skills:biopython) - Sequence analysis and bioinformatics
- Skill(scientific-skills:rdkit) - Cheminformatics and molecular processing
- Skill(scientific-skills:polars) - High-performance data processing

## Documentation

For more information about GOSTAR database schema and features, consult your organization's internal GOSTAR documentation or contact the database administrator: Dan Myung (dan.myung@takeda.com)
