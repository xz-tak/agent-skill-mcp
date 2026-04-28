--
name: database-drugbank
description: Use when user asks to "query drugbank", "search drugbank database", "find drugs", "drug target", "clinical trial", "adverse effect", "drug indication", "pharmacology", or mentions drug names, protein targets, or pharmaceutical analysis. Provides comprehensive knowledge of the DrugBank medicinal database schema and query patterns.
version: 0.1.0
---

# DrugBank Database

## Overview

The DrugBank database is a comprehensive pharmaceutical database containing detailed information about drugs, drug targets, mechanisms of action, clinical trials, indications, adverse effects, and pharmacology. The database contains 91 tables with interconnected relationships covering all aspects of drug development and clinical use.

**Database connection details:**
- Host: usvgarps11158-dev003.cm9aqaugy64i.us-east-1.rds.amazonaws.com
- Port: 5442
- Database: drugbank
- Credentials: Retrieved from AWS Secrets Manager (profiles: `cmp-dev`, `sci-dev`)

**S3 fallback (when RDS is unreachable due to VPC/network issues):**
- CSV: `s3://tec-rnd-sci-dev-gi2/gi2-xz/insights/project/test/data/drugbank/drugbank.csv` (515 MB, 17,430 drugs)
- XML: `s3://tec-rnd-sci-dev-gi2/gi2-xz/insights/project/test/data/drugbank/drugbank_all_full_database_5.1.16.xml.zip` (186 MB → 2.5 GB unzipped)
- Access: SageMaker execution role (no profile needed) or `cmp-dev` profile
- CSV columns: `drugbank_id_primary, name, synonyms, target_genes, categories, mechanism_of_action, state, groups, ...`
- Use `aws s3 cp <s3_uri> /tmp/drugbank.csv` to download, then query locally with pandas

## When to Use This Skill

Use this skill when querying the DrugBank database for:
- Drug information (names, structures, properties)
- Drug-target interactions and mechanisms
- Multi-target drug discovery queries
- Clinical trial data
- Drug indications and approvals
- Adverse effects and safety profiles
- Drug-drug interactions
- Metabolic pathways and pharmacokinetics

## Core Database Structure

### Primary Tables

**Drugs Core:**
- `drugs` - Main drug entries (drugbank_id, name, type, state, description)
- `drug_calculated_properties` - Molecular properties (SMILES, InChI, molecular_formula, molecular_weight)
- `drug_synonyms` - Alternative drug names

**Drug-Target Interactions:**
- `bonds` - Drug-target relationships (drug_id, biodb_id, type, pharmacological_action, inhibitor, agonist, antagonist)
- `bio_entities` - Biological targets (biodb_id, name, type)
- `bio_entity_components` - Target components (uniprot_id links)
- `polypeptides` - Protein details (gene_name, general_function, specific_function)

**Clinical Data:**
- `clinical_trials` - Trial information (identifier, title, status, purpose, start_date, end_date)
- `clinical_trial_interventions` - Intervention details
- `clinical_trial_interventions_drugs` - Drug-intervention links

**Indications:**
- `structured_indications` - Drug indication records (kind, off_label, country)
- `indication_conditions` - Indication-condition relationships
- `conditions` - Medical conditions (title, snomed_id, meddra_id, icd10_id)

**Adverse Effects:**
- `structured_adverse_effects` - Adverse event records (severity)
- `adverse_effect_conditions` - AE-condition links
- `adverse_effect_incidences` - Frequency data (percent, kind)

**Pharmacology:**
- `pathways` - Metabolic pathways
- `metabolites` - Drug metabolites
- `structured_drug_interactions` - Drug-drug interactions

### Key Relationships

```
drugs (1) -----> (*) bonds -----> (1) bio_entities
                                        |
                                        v
                               bio_entity_components
                                        |
                                        v
                                   polypeptides

drugs (1) -----> (*) structured_indications -----> (*) indication_conditions -----> (1) conditions

drugs (1) -----> (*) structured_adverse_effects -----> (*) adverse_effect_conditions -----> (1) conditions

drugs (1) -----> (*) clinical_trial_interventions_drugs -----> (1) clinical_trial_interventions -----> (1) clinical_trials
```

## Query Patterns

### Drug Search

Search by name or DrugBank ID using ILIKE for flexible matching:

```sql
SELECT d.drugbank_id, d.name, d.type, d.state, d.description
FROM drugs d
WHERE d.name ILIKE '%search_term%' OR d.drugbank_id ILIKE '%search_term%';
```

Get molecular properties:

```sql
SELECT d.*, dcp.smiles, dcp.molecular_formula, dcp.molecular_weight
FROM drugs d
LEFT JOIN drug_calculated_properties dcp ON d.id = dcp.drug_id
WHERE d.drugbank_id = 'DB00112';
```

### Target-Based Queries

Find drugs targeting a specific protein:

```sql
SELECT DISTINCT
    d.drugbank_id, d.name, d.state,
    be.name AS target_name,
    p.gene_name,
    b.inhibitor, b.agonist, b.antagonist,
    b.pharmacological_action
FROM drugs d
JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
JOIN bio_entities be ON b.biodb_id = be.biodb_id
LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
WHERE be.name ILIKE '%target_name%' OR p.gene_name ILIKE '%gene_name%';
```

### Multi-Target Analysis

Find drugs affecting multiple targets (intersection):

```sql
-- Find drugs targeting both target1 AND target2
WITH target1_drugs AS (
    SELECT DISTINCT d.drugbank_id
    FROM drugs d
    JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
    JOIN bio_entities be ON b.biodb_id = be.biodb_id
    LEFT JOIN polypeptides p ON be.biodb_id = p.uniprot_id
    WHERE be.name ILIKE '%target1%' OR p.gene_name ILIKE '%target1%'
),
target2_drugs AS (
    SELECT DISTINCT d.drugbank_id
    FROM drugs d
    JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
    JOIN bio_entities be ON b.biodb_id = be.biodb_id
    LEFT JOIN polypeptides p ON be.biodb_id = p.uniprot_id
    WHERE be.name ILIKE '%target2%' OR p.gene_name ILIKE '%target2%'
)
SELECT d.drugbank_id, d.name, d.type, d.state, d.description
FROM drugs d
WHERE d.drugbank_id IN (SELECT drugbank_id FROM target1_drugs)
  AND d.drugbank_id IN (SELECT drugbank_id FROM target2_drugs);
```

### Clinical Trials

Query trials for a drug:

```sql
SELECT DISTINCT
    ct.identifier, ct.title, ct.status, ct.purpose,
    ct.start_date, ct.end_date,
    d.drugbank_id, d.name
FROM clinical_trials ct
JOIN clinical_trial_interventions cti ON ct.identifier = cti.trial_id
JOIN clinical_trial_interventions_drugs ctid ON cti.id = ctid.intervention_id
JOIN drugs d ON ctid.drug_id = d.id
WHERE d.name ILIKE '%drug_name%' OR d.drugbank_id ILIKE '%drug_id%'
ORDER BY ct.start_date DESC;
```

### Indications

Get approved indications:

```sql
SELECT DISTINCT
    d.drugbank_id, d.name,
    c.title AS indication,
    si.kind, si.off_label, si.country,
    c.snomed_id, c.meddra_id
FROM structured_indications si
JOIN drugs d ON si.drug_id = d.id
JOIN indication_conditions ic ON si.id = ic.indication_id AND ic.relationship = 'for_condition'
JOIN conditions c ON ic.condition_id = c.id
WHERE d.name ILIKE '%drug_name%' OR d.drugbank_id ILIKE '%drug_id%';
```

### Adverse Effects

Query adverse effects with frequency data:

```sql
SELECT
    d.drugbank_id, d.name,
    c.title AS adverse_effect,
    aei.percent AS frequency_percent,
    sae.severity,
    c.snomed_id, c.meddra_id
FROM structured_adverse_effects sae
JOIN drugs d ON sae.drug_id = d.id
JOIN adverse_effect_conditions aec ON sae.id = aec.adverse_effect_id AND aec.relationship = 'effect'
JOIN conditions c ON aec.condition_id = c.id
LEFT JOIN adverse_effect_incidences aei ON sae.id = aei.adverse_effect_id
WHERE d.name ILIKE '%drug_name%'
ORDER BY aei.percent DESC NULLS LAST;
```

## Using Helper Scripts

The plugin includes a Python helper module at `scripts/drugbank_helper.py` with pre-built query functions.

**Available functions:**
- `search_drugs(search_term)` - Search by name or ID
- `get_drug_by_id(drugbank_id)` - Get detailed drug info
- `get_drug_targets(drug_identifier)` - Get target interactions
- `get_clinical_trials(drug_identifier)` - Query trials
- `get_drug_indications(drug_identifier)` - Get indications
- `get_adverse_effects(drug_identifier, min_frequency)` - Get AEs
- `search_by_target(target_name)` - Find drugs by target
- `get_drug_interactions(drug_identifier)` - Get drug-drug interactions

**When to use helpers:**
- Standard single-drug queries
- Simple target searches
- Quick lookups

**When to use custom SQL:**
- Multi-target analysis
- Complex joins across multiple domains
- Custom filtering or aggregation
- Advanced analytical queries

## Query Strategy

### Approach Selection

1. **Simple lookups** → Use helper functions for token efficiency
2. **Complex analysis** → Write custom SQL with appropriate joins
3. **Multi-step queries** → Break into logical steps (find targets → find drugs → analyze results)

### Performance Considerations

- Always use ILIKE with wildcards for flexible string matching
- Filter early in joins to reduce result set size
- Use DISTINCT when joining many-to-many relationships
- Add LIMIT clauses for exploratory queries

### Data Interpretation

**Drug states:**
- `approved` - FDA/regulatory approved
- `investigational` - In clinical development
- `experimental` - Preclinical research
- `withdrawn` - Previously approved, now withdrawn

**Bond types:**
- `TargetBond` - Primary drug target
- `EnzymeBond` - Metabolic enzyme
- `CarrierBond` - Transport protein
- `TransporterBond` - Membrane transporter

**Pharmacological actions:**
- `inhibitor` - Inhibits target activity
- `agonist` - Activates target
- `antagonist` - Blocks target
- `substrate` - Metabolized by enzyme

## Additional Resources

### Reference Files

For detailed information, consult:
- **`references/schema.md`** - Complete table documentation (all 91 tables)
- **`references/query-patterns.md`** - Advanced SQL patterns and optimization techniques

### Example Queries

Working examples in `examples/`:
- **`common-queries.sql`** - Generic query templates for common use cases

### Helper Scripts

Utilities in `scripts/`:
- **`drugbank_helper.py`** - Python query utilities with pre-built functions

## Best Practices

1. **Flexible matching**: Always use ILIKE with wildcards for drug/target names
2. **ID vs name**: Use DrugBank IDs when available for exact matches
3. **Filter relationships**: Always filter bond.type = 'TargetBond' for drug targets
4. **Check nulls**: Use LEFT JOIN and handle NULL values in results
5. **Multi-target logic**: Use CTEs or subqueries for intersection/union queries
6. **Performance**: Add LIMIT for exploratory queries, remove for complete results
7. **Data quality**: Some fields may be NULL (e.g., pharmacological_action, frequency data)

## Error Handling

**Common issues:**

- **Permission denied**: Ensure database credentials are set in `.claude/drugbank.local.md`
- **No results**: Check spelling, try broader ILIKE patterns, verify target nomenclature
- **Ambiguous results**: Refine search terms or add filters (e.g., state = 'approved')
- **Missing data**: Some drugs lack complete target/clinical/AE data - this is expected

**Credential retrieval:**

Credentials are automatically retrieved from AWS Secrets Manager using profiles `cmp-dev` or `sci-dev`. If AWS Secrets Manager is not available, the system falls back to environment variables `DRUGBANK_USERNAME` and `DRUGBANK_PASSWORD`.

If credentials cannot be retrieved, the connection will fail.

## Workflow

For any DrugBank query:

1. **Understand the question** - Identify what data is needed (drugs, targets, clinical, etc.)
2. **Choose approach** - Helper function or custom SQL
3. **Construct query** - Use appropriate tables and joins
4. **Execute** - Via MCP PostgreSQL server or helper script
5. **Interpret results** - Present findings with context (state, pharmacological action, etc.)
6. **Follow-up** - Offer related analyses if helpful

Remember: This skill enables flexible, generic querying. Do not hardcode specific drug names or targets in queries - always use parameters from user requests.
