# Cortellis API Fields Reference

## Available Fields

The Cortellis gene query script supports the following data fields:

### annotation (always included)
Target information from the Cortellis Targets API:
- **Target Name**: Official name of the gene/protein target
- **Symbol**: Gene symbol
- **Gene ID**: NCBI Entrez Gene ID
- **UniProt ID**: UniProt accession
- **Target Type**: Classification (e.g., Protein, RNA)
- **Description**: Detailed description of the target's function
- **Organism**: Species (e.g., Homo sapiens)
- **Synonyms**: Alternative names
- **External Identifiers**: Cross-references to other databases (SwissProt, etc.)
- **Related Targets**: Related target IDs (isoforms, variants, complexes) - queried automatically for comprehensive drug coverage

### drug
Drug-target associations from both Targets API and Investigational Drugs API:

**From Targets API (basic records):**
- Associated drugs and their mechanisms
- Condition-drug associations
- Basic drug information

**From Investigational Drugs API (comprehensive records):**
- **Drug Name**: Official drug name
- **Drug Synonyms**: Alternative names, brand names, code names
- **Highest Phase**: Development phase (Preclinical, Phase 1/2/3, Approved, etc.) with phase ID
- **Originator**: Company that originated the drug
- **Primary Companies**: Companies involved in development
- **Primary Indications**: Therapeutic indications
- **Mechanisms of Action**: Drug's mechanism of action (ActionsPrimary)
- **Therapy Areas**: Therapeutic areas (e.g., Cancer, Immunology)
- **Technologies**: Drug modality (mAb, small molecule, biologic, etc.)
- **Regulatory Designations**: Special statuses (Orphan Drug, Fast Track, Breakthrough Therapy, etc.)
- **Structure (SMILES)**: Chemical structure representation
- **Clinical Trials**: Related clinical trial records with trial IDs, titles, indications

The script attempts to retrieve comprehensive drug records by:
1. Querying Targets API for basic drug associations
2. Searching Investigational Drugs API using drug names
3. Fetching full drug records with development information

Both basic (`Drug`) and comprehensive (`DrugRecord`) records are preserved in JSON output.

### biomarker
Biomarker information from the Biomarkers API:
- **Biomarker Uses**: Applications of the biomarker
- **Biomarker Application**: Type of application (diagnostic, prognostic, etc.)
- **Role**: Biomarker's role in disease/treatment
- **Indication**: Associated disease or condition
- **Indication Type**: Classification of the indication
- **Drugs Studied**: Drugs studied with this biomarker
- **Validity**: Validation status

### interaction
Protein-protein interactions from the Targets API:
- **Interaction Type**: Type of interaction (e.g., binding, phosphorylation)
- **Partner Name**: Name of the interacting protein
- **Partner Gene Symbol**: Gene symbol of interaction partner
- **Partner UniProt ID**: UniProt ID of partner
- **Evidence**: Supporting evidence for the interaction
- **PubMed ID**: Literature references

### association
Condition-gene associations (disease associations):
- **Condition-Gene Associations**: Links between genes and diseases
- **Condition Name**: Disease or condition name
- **Condition ID**: Identifier for the condition
- **Association Type**: Nature of the association (causal, risk factor, etc.)
- **Evidence**: Supporting evidence
- **Source**: Data source
- **PubMed ID**: Literature references
- **Condition-Gene-Variant Associations**: Genetic variant associations

## Default Behavior

- **Default fields**: If no fields specified, queries `drug` and `biomarker`
- **All fields**: Use `--all` flag to query all available fields
- **Annotation**: Always included regardless of field selection

## Output Formats

1. **JSON** (`{GENE}_cortellis_data.json`): Complete API response with nested structure
2. **Markdown** (`{GENE}_summary.md`): Human-readable summary of findings
3. **Excel** (`{GENE}_cortellis_data.xlsx`): Tabular data with separate sheets per field

## API Endpoints

The script uses multiple Cortellis API endpoints:

- **Targets API v2** (auth v3): Gene/target annotations, basic drug associations
- **Biomarkers API v3** (auth v4): Biomarker information
- **Investigational Drugs API v2.0** (auth v2): Comprehensive drug development data

## Data Completeness Notes

- **Comprehensive Drug Records**: Not all drugs from Targets API will have comprehensive records in Investigational Drugs API. Some research compounds may only have basic information.
- **Gene Matching**: Uses exact synonym matching to ensure correct gene identification
- **Search Results**: Drug name searches may not always return results if drug names differ between APIs
