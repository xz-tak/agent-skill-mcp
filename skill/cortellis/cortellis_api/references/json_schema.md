# Cortellis JSON Data Structure Reference

This document describes the JSON schema/structure of Cortellis data files from both API queries and Excel conversions, based on **Cortellis API v2.0**.

**Scope:** This documentation covers the most commonly used fields:
- **annotation** - Gene/target information with identifiers and related targets
- **drug** - Drug associations, development data, clinical trials
- **biomarker** - Biomarker data (when queried)

Other available fields (interaction, association) follow similar patterns. Refer to `api_fields.md` for complete field listing.

## File Naming Convention

- **API-generated:** `{GENE}_cortellis_data.json`
- **Excel-converted:** `{GENE}_cortellis_data.json` (compatible format)

---

## JSON Schema Overview

```json
{
  "annotation": {
    "@Id": "Target ID",
    "@namemain": "Target Name",
    "Symbol": "GENE",
    "GeneId": "Entrez Gene ID",
    "UniprotId": "UniProt ID",
    "TargetType": "Protein",
    "Description": "Protein description",
    "Organism": {"$": "Homo sapiens"},
    "RelatedTargets": {
      "Id": [{"@type": "Isoform", "$": "related_target_id"}]
    }
  },
  "drug": {
    "Target": {
      "@namemain": "Target name",
      "@id": "Target ID",
      "ConditionDrugAssociations": {...}
    },
    "Drug": [...],
    "Trial": [...],
    "DrugRecord": {
      "drug_id_1": {...},
      "drug_id_2": {...}
    }
  }
}
```

---

## 1. Annotation Section

Basic gene/target information including identifiers and related targets.

### Complete Structure

```json
"annotation": {
  "@Id": "156324211736763",
  "@namemain": "Integrin alpha-4",
  "Symbol": "ITGA4",
  "GeneId": "3676",
  "UniprotId": "P13612",
  "TargetType": "Protein",
  "Description": "Integrin alpha-4 is a receptor for...",
  "Organism": {
    "$": "Homo sapiens"
  },
  "EntrezgeneIdentifiers": {
    "Identifier": 3676
  },
  "ExternalIdentifiers": {
    "Identifier": [
      {
        "@type": "SwissProt",
        "@id": "P13612",
        "$": "P13612"
      }
    ]
  },
  "RelatedTargets": {
    "Id": [
      {
        "@type": "Isoform",
        "$": "related_target_id"
      }
    ]
  },
  "Synonyms": {
    "Synonym": ["CD49D", "VLA-4"]
  },
  "_source": "API",
  "_conversion_date": "2025-12-21T..."
}
```

### Key Fields

| Field | Type | Description |
|-------|------|-------------|
| `@Id` | String | Cortellis internal target ID |
| `@namemain` | String | Primary target name |
| `Symbol` | String | Gene symbol (e.g., ITGA4) |
| `GeneId` | String | Entrez Gene ID |
| `UniprotId` | String | UniProt accession |
| `TargetType` | String | Usually "Protein" |
| `Description` | String | Protein function description |
| `RelatedTargets` | Object | Related target IDs (isoforms, variants) |
| `Synonyms` | Object | Alternative gene/protein names |

### Related Targets

**Important for comprehensive drug queries:** The script automatically queries related targets (isoforms, variants) to capture all associated drugs.

```json
"RelatedTargets": {
  "Id": [
    {"@type": "Isoform", "$": "156324211736764"},
    {"@type": "Variant", "$": "156324211736765"}
  ]
}
```

---

## 2. Drug Section

The drug section contains multiple data sources with different levels of detail.

### Structure Overview

```json
"drug": {
  "Target": {...},           // Basic target info
  "Drug": [...],             // Basic drug list (Targets API)
  "Trial": [...],            // Clinical trials
  "DrugRecord": {...}        // Comprehensive drug records (Investigational Drugs API)
}
```

### 2.1 Target (Basic Info)

Condition-drug associations from Targets API:

```json
"Target": {
  "@namemain": "Integrin alpha-4",
  "@id": "156324211736763",
  "ConditionDrugAssociations": {
    "Condition": [
      {
        "@id": "291",
        "@name": "Multiple sclerosis",
        "DrugId": [
          {
            "@status": "Approved",
            "@highestphase": "Launched",
            "$": "4838"
          }
        ]
      }
    ]
  }
}
```

### 2.2 Drug (Basic Drug List)

Simple array of drug records from Targets API with basic information:

```json
"Drug": [
  {
    "@id": "4838",
    "@namemain": "natalizumab",
    "NamesChemicalAndDescriptions": {
      "Name": "Humanized anti-integrin alpha-4 monoclonal antibody"
    },
    "NamesCode": {
      "Name": "AN-100226"
    },
    "MechanismsMolecular": {
      "Mechanism": {
        "@id": "0",
        "$": "Integrin alpha-4 inhibitor"
      }
    }
  }
]
```

### 2.3 Trial (Clinical Trials)

Clinical trial records associated with the target:

```json
"Trial": [
  {
    "@id": "NCT01234567",
    "TitleDisplay": "Study of Drug X in Disease Y",
    "TitleOfficial": "A Phase 3, Randomized, Double-Blind...",
    "Indications": {
      "Indication": [
        {
          "@id": "213",
          "$": "Rheumatoid arthritis"
        }
      ]
    },
    "ProtocolAndOutcomes": {...}
  }
]
```

### 2.4 DrugRecord (Comprehensive Drug Data) ⭐

**Most Important Section for Analysis**

Comprehensive drug records from Investigational Drugs API v2.0. This is a **dictionary keyed by drug ID** (not drug name):

```json
"DrugRecord": {
  "106955": {
    "@id": "106955",
    "DrugName": "Liocyx-D",
    "DrugNamesKey": {
      "Name": {
        "@id": "106955",
        "$": "Liocyx-D"
      }
    },
    "DrugSynonyms": {
      "Name": [
        {"Value": "LIO-D"},
        {"Value": "Liocyx D"}
      ]
    },
    "PhaseHighest": {
      "@id": "C1",
      "$": "Phase 1 Clinical"
    },
    "CompanyOriginator": {
      "@id": "company_id",
      "$": "Liocell Biotherapeutics"
    },
    "CompaniesPrimary": {
      "Company": [
        {"@id": "1", "$": "Partner Company 1"},
        {"@id": "2", "$": "Partner Company 2"}
      ]
    },
    "IndicationsPrimary": {
      "Indication": {
        "@id": "1767",
        "$": "Hepatocellular carcinoma"
      }
    },
    "ActionsPrimary": {
      "Action": [
        {"@id": "0", "$": "Integrin alpha-4 antagonist"},
        {"@id": "1", "$": "Integrin beta-1 antagonist"}
      ]
    },
    "TherapyAreas": {
      "TherapyArea": "Cancer"
    },
    "Technologies": {
      "Technology": [
        {"@id": "0", "$": "Monoclonal antibody"}
      ]
    },
    "RegulatoryDesignations": {
      "RegulatoryDesignation": [
        {"@id": "0", "$": "Orphan Drug"}
      ]
    },
    "StructureSmiles": "...",
    "MolecularFormula": "C6H12O6",
    "MolecularWeight": "180.16",
    "CASNumber": "12345-67-8",
    "DevelopmentProfile": {
      "Summary": {
        "displayLabel": "Summary",
        "value": "<Summary><para>Drug description...</para></Summary>"
      }
    },
    "_source_drug_name": "Liocyx-D",
    "_source_display_name": "Liocyx-D",
    "_excel_row": 0
  },
  "4838": {
    "@id": "4838",
    "DrugName": "natalizumab",
    "PhaseHighest": {
      "@id": "LA",
      "$": "Launched"
    },
    "IndicationsPrimary": {
      "Indication": [
        {"@id": "213", "$": "Multiple sclerosis"},
        {"@id": "500", "$": "Crohn's disease"}
      ]
    },
    ...
  }
}
```

### DrugRecord Key Fields

| Field | Type | Description |
|-------|------|-------------|
| `@id` | String | Drug ID (also the dictionary key) |
| `DrugName` | String | Primary drug name |
| `DrugSynonyms` | Object | Alternative names/codes |
| `PhaseHighest` | Object | Highest development phase |
| `CompanyOriginator` | Object | Original developer |
| `CompaniesPrimary` | Object | Current developers/partners |
| `IndicationsPrimary` | Object | Primary therapeutic indications |
| `ActionsPrimary` | Object | Mechanisms of action |
| `TherapyAreas` | Object | Therapeutic areas |
| `Technologies` | Object | Drug modality (mAb, small molecule, etc.) |
| `RegulatoryDesignations` | Object | Special designations (Orphan, Fast Track, etc.) |
| `StructureSmiles` | String | Chemical structure (SMILES) |
| `_source_drug_name` | String | Original drug name from Targets API search |

---

## 3. Key Fields for Analysis

### PhaseHighest

Development phase with both ID code and human-readable value.

#### Format

```json
"PhaseHighest": {
  "@id": "C3",
  "$": "Phase 3 Clinical"
}
```

#### Phase Values and Scoring Weights

| @id | $ Value | Meaning | Disease-Specific Score | Non-Specific Score |
|-----|---------|---------|----------------------|-------------------|
| `LA` | `Launched` | FDA approved / On market | 7 | 4 |
| `APR` | `Approved` | Approved (alternate format) | 7 | 4 |
| `C3` | `Phase 3 Clinical` | Late-stage clinical | 3 | 2 |
| `C2` | `Phase 2 Clinical` | Mid-stage clinical | 2 | 1 |
| `C1` | `Phase 1 Clinical` | Early clinical | 1 | 0.5 |
| `PC` | `Preclinical` | Preclinical development | 0.1 | 0.1 |
| `DI` | `Discontinued` | Development stopped | 0 | 0 |
| `ND` | `No Development Reported` | No active development | 0 | 0 |

**Note:** Scoring weights are context-specific. See `scoring_framework.md` for detailed guidance.

### IndicationsPrimary

Therapeutic indications. **Can be a single object or array** - always check type!

#### Single Indication

```json
"IndicationsPrimary": {
  "Indication": {
    "@id": "213",
    "$": "Rheumatoid arthritis"
  }
}
```

#### Multiple Indications

```json
"IndicationsPrimary": {
  "Indication": [
    {
      "@id": "213",
      "$": "Rheumatoid arthritis"
    },
    {
      "@id": "500",
      "$": "Crohn's disease"
    },
    {
      "@id": "273",
      "$": "Ulcerative colitis"
    }
  ]
}
```

#### Empty Indications

```json
"IndicationsPrimary": {}
```

### ActionsPrimary

Mechanisms of action. Same single/array pattern as indications:

```json
"ActionsPrimary": {
  "Action": [
    {"@id": "0", "$": "Integrin alpha-4 antagonist"},
    {"@id": "1", "$": "Integrin beta-1 antagonist"}
  ]
}
```

### Company Information

#### Originator

```json
"CompanyOriginator": {
  "@id": "1011918",
  "$": "Pfizer Inc"
}
```

#### Primary Companies (Partners)

```json
"CompaniesPrimary": {
  "Company": [
    {"@id": "1", "$": "Pfizer Inc"},
    {"@id": "2", "$": "BioNTech SE"}
  ]
}
```

#### Single Company

```json
"CompaniesPrimary": {
  "Company": {
    "@id": "1",
    "$": "Single Company Name"
  }
}
```

### TherapyAreas

Can be a single string or array of strings:

#### Single Area

```json
"TherapyAreas": {
  "TherapyArea": "Cancer"
}
```

#### Multiple Areas

```json
"TherapyAreas": {
  "TherapyArea": ["Cancer", "Immunology", "Autoimmune diseases"]
}
```

### Technologies (Drug Modality)

Drug type/modality:

```json
"Technologies": {
  "Technology": [
    {"@id": "0", "$": "Monoclonal antibody"},
    {"@id": "1", "$": "Humanized antibody"}
  ]
}
```

### DrugSynonyms

Alternative drug names and codes:

```json
"DrugSynonyms": {
  "Name": [
    {"Value": "Tysabri"},
    {"Value": "AN-100226"},
    {"Value": "Antegren"}
  ]
}
```

---

## 4. Accessing Data in Python

### Loading JSON

```python
import json

with open('ITGA4_cortellis_data.json', 'r') as f:
    data = json.load(f)

# Access annotation
gene_name = data['annotation']['@namemain']
gene_symbol = data['annotation']['Symbol']
gene_id = data['annotation'].get('GeneId', '')

# Access drug records (IMPORTANT: nested under data['drug'])
drug_records = data['drug']['DrugRecord']  # Dictionary keyed by drug ID
```

### Iterating Through Drugs

```python
# DrugRecord is a dictionary keyed by drug ID, not a list!
for drug_id, drug_info in drug_records.items():
    # Get drug name
    name = drug_info.get('DrugName', '')

    # Get phase
    phase_data = drug_info.get('PhaseHighest', {})
    if isinstance(phase_data, dict):
        phase = phase_data.get('$', 'Unknown')
        phase_id = phase_data.get('@id', '')
    else:
        phase = str(phase_data)

    # Get indications (ALWAYS check if single or array!)
    ind_data = drug_info.get('IndicationsPrimary', {})
    indications = []
    if isinstance(ind_data, dict) and 'Indication' in ind_data:
        ind_list = ind_data['Indication']

        # Handle single indication vs array
        if not isinstance(ind_list, list):
            ind_list = [ind_list]

        # Extract indication names
        for ind in ind_list:
            if isinstance(ind, dict):
                ind_name = ind.get('$', '')
                if ind_name:
                    indications.append(ind_name)

    # Get mechanisms
    actions_data = drug_info.get('ActionsPrimary', {})
    mechanisms = []
    if isinstance(actions_data, dict) and 'Action' in actions_data:
        action_list = actions_data['Action']
        if not isinstance(action_list, list):
            action_list = [action_list]

        for action in action_list:
            if isinstance(action, dict):
                mechanism = action.get('$', '')
                if mechanism:
                    mechanisms.append(mechanism)

    # Get company
    originator = drug_info.get('CompanyOriginator', {})
    company = originator.get('$', 'Unknown') if isinstance(originator, dict) else ''

    # Get therapy areas
    therapy_data = drug_info.get('TherapyAreas', {})
    therapy_areas = []
    if isinstance(therapy_data, dict):
        areas = therapy_data.get('TherapyArea', [])
        if isinstance(areas, list):
            therapy_areas = areas
        elif areas:
            therapy_areas = [areas]

    print(f"{name} ({phase}) - {company}")
    print(f"  Indications: {', '.join(indications)}")
    print(f"  Mechanisms: {', '.join(mechanisms)}")
```

### Helper Function for Extracting Values

```python
def extract_value(obj, default=''):
    """Extract value from dict with '$' key or return the object itself."""
    if isinstance(obj, dict):
        return obj.get('$', default)
    return obj if obj is not None else default

def as_list(x):
    """Convert to list if not already."""
    if x is None:
        return []
    return x if isinstance(x, list) else [x]

# Usage
phase = extract_value(drug_info.get('PhaseHighest'))
indications = as_list(drug_info.get('IndicationsPrimary', {}).get('Indication'))
```

### Filtering Drugs by Disease

```python
def matches_disease(indications_list, target_diseases):
    """Check if any indication matches target diseases."""
    target_lower = [d.lower() for d in target_diseases]

    for indication in indications_list:
        ind_lower = indication.lower()
        if any(target in ind_lower for target in target_lower):
            return True
    return False

# Example: Find IBD drugs
target_diseases = ['crohn', 'colitis', 'inflammatory bowel']

ibd_drugs = []
for drug_id, drug_info in drug_records.items():
    ind_data = drug_info.get('IndicationsPrimary', {})
    ind_list = as_list(ind_data.get('Indication', []))

    indications = [extract_value(ind) for ind in ind_list]

    if matches_disease(indications, target_diseases):
        ibd_drugs.append({
            'name': drug_info.get('DrugName'),
            'phase': extract_value(drug_info.get('PhaseHighest')),
            'indications': indications
        })

print(f"Found {len(ibd_drugs)} IBD drugs")
```

---

## 5. Excel-Converted Data Format

When using `convert_excel_to_json.py`, the structure is **identical** to API format with these notes:

### Differences from API Data

1. **Annotation Section:**
   - Less comprehensive (some fields may be empty)
   - `_source` marked as `"Manual Excel Download"`
   - Missing RelatedTargets data

2. **Drug Section:**
   - All drugs placed in `DrugRecord` section
   - Basic `Drug` list also populated for compatibility
   - Phase strings normalized to match API format
   - `Trial` section empty (not available in Excel)

3. **Additional Fields:**
   - `RegionalDevelopment` - Development status by country/region (from Development Status sheet)
   - `Milestones` - Development milestones with dates (from Milestones sheet)
   - `_excel_row` - Original Excel row number for traceability

### Example Excel-Converted Record

```json
"DrugRecord": {
  "12345": {
    "@id": "12345",
    "DrugName": "Example Drug",
    "PhaseHighest": {
      "@id": "LA",
      "$": "Launched"
    },
    "RegionalDevelopment": [
      {
        "Country": "United States",
        "Phase": "Launched - 2022",
        "Organization": "Company Name",
        "Indication": "Disease X",
        "FormulationRoute": "Tablet (Oral)"
      }
    ],
    "Milestones": [
      {
        "Date": "2022-03-15",
        "Type": "FDA Approval",
        "Notes": "Approved for Disease X",
        "Organization": "FDA",
        "Country": "United States"
      }
    ],
    "_excel_row": 42
  }
}
```

---

## 6. Data Structure Variations

### Empty Sections

If no drugs are found:
```json
"drug": {
  "Target": {...},
  "Drug": [],
  "Trial": [],
  "DrugRecord": {}
}
```

### Single vs Multiple Items Pattern

**Critical:** Many fields can be either a single object or an array. Always check type:

```python
# WRONG - assumes array
for indication in drug_info['IndicationsPrimary']['Indication']:
    print(indication)  # ERROR if single object!

# CORRECT - handle both
ind_list = drug_info.get('IndicationsPrimary', {}).get('Indication', [])
if not isinstance(ind_list, list):
    ind_list = [ind_list]

for indication in ind_list:
    print(indication)
```

**Fields with this pattern:**
- `IndicationsPrimary.Indication`
- `ActionsPrimary.Action`
- `CompaniesPrimary.Company`
- `Technologies.Technology`
- `RegulatoryDesignations.RegulatoryDesignation`
- `TherapyAreas.TherapyArea`
- `ExternalIdentifiers.Identifier`
- `RelatedTargets.Id`

---

## 7. Common Pitfalls

### 1. DrugRecord is nested

❌ **WRONG:**
```python
drug_records = data['DrugRecord']  # KeyError!
```

✅ **CORRECT:**
```python
drug_records = data['drug']['DrugRecord']
```

### 2. DrugRecord is keyed by drug ID, not name

❌ **WRONG:**
```python
for drug in drug_records:  # Only gets keys (IDs)
    print(drug)  # Prints: "106955", "4838", ...
```

✅ **CORRECT:**
```python
for drug_id, drug_info in drug_records.items():
    print(drug_info['DrugName'])  # Prints actual names
```

### 3. Phase format varies

❌ **WRONG:**
```python
phase = drug_info['PhaseHighest']  # Gets dict, not string!
```

✅ **CORRECT:**
```python
phase_data = drug_info.get('PhaseHighest', {})
phase = phase_data.get('$', 'Unknown') if isinstance(phase_data, dict) else str(phase_data)
```

### 4. Indications can be single or array

❌ **WRONG:**
```python
for indication in drug_info['IndicationsPrimary']['Indication']:
    # Fails if single object!
```

✅ **CORRECT:**
```python
ind_data = drug_info.get('IndicationsPrimary', {})
ind_list = ind_data.get('Indication', [])
if not isinstance(ind_list, list):
    ind_list = [ind_list]

for indication in ind_list:
    print(extract_value(indication))
```

### 5. Dollar sign fields

❌ **WRONG:**
```python
phase = drug_info['PhaseHighest']  # Gets: {"@id": "C3", "$": "Phase 3 Clinical"}
```

✅ **CORRECT:**
```python
phase = drug_info['PhaseHighest']['$']  # Gets: "Phase 3 Clinical"
# Or use helper function:
phase = extract_value(drug_info['PhaseHighest'])
```

### 6. Empty vs missing fields

❌ **WRONG:**
```python
if drug_info['IndicationsPrimary']:  # KeyError if missing!
```

✅ **CORRECT:**
```python
ind_data = drug_info.get('IndicationsPrimary', {})
if ind_data and 'Indication' in ind_data:
    # Process indications
```

### 7. Related targets

❌ **MISSING DATA:**
```python
# Only querying main target misses drugs from related targets (isoforms, variants)
result = target_drug(target_id)
```

✅ **COMPREHENSIVE:**
```python
# Script automatically includes related targets for complete drug coverage
result = target_drug(target_id, annotation)  # Pass annotation to get related targets
```

---

## 8. Validation Patterns

### Check Data Completeness

```python
def validate_drug_record(drug_info):
    """Validate essential fields are present."""
    issues = []

    if not drug_info.get('DrugName'):
        issues.append("Missing DrugName")

    if not drug_info.get('PhaseHighest'):
        issues.append("Missing PhaseHighest")

    if not drug_info.get('IndicationsPrimary'):
        issues.append("Missing IndicationsPrimary")

    return issues

# Usage
for drug_id, drug_info in drug_records.items():
    issues = validate_drug_record(drug_info)
    if issues:
        print(f"Drug {drug_id}: {', '.join(issues)}")
```

### Check for Data Loss

```python
# Check if API data vs Excel conversion
if data['annotation'].get('_source') == 'Manual Excel Download':
    print("Data from Excel conversion")
    print("Note: Trial data not available from Excel")
    print(f"Drug records: {len(data['drug']['DrugRecord'])}")
else:
    print("Data from API query")
    print(f"Drug records: {len(data['drug']['DrugRecord'])}")
    print(f"Trials: {len(data['drug'].get('Trial', []))}")
```

---

## 9. Example Analysis Workflows

See comprehensive reference implementations in [`examples/`](../examples/):

- **Excel Conversion:** `convert_excel_to_json.py` - Convert Cortellis Excel exports to API-compatible JSON
- **Target Scoring:** `analyze_targets_example.py` - Score targets based on drug development data
- **Report Generation:** `generate_report_example.py` - Generate clinical intelligence reports

---

## 10. Related Documentation

- **API Fields Reference:** [`api_fields.md`](./api_fields.md) - Complete field descriptions
- **API Endpoints:** [`api_reference.md`](./api_reference.md) - API endpoints and authentication
- **Scoring Framework:** [`scoring_framework.md`](./scoring_framework.md) - Clinical validation scoring methodology
- **Excel Conversion Guide:** [`../examples/EXCEL_TO_JSON_CONVERSION_GUIDE.md`](../examples/EXCEL_TO_JSON_CONVERSION_GUIDE.md) - Comprehensive Excel to JSON conversion guide

---

**Document Version:** 2.0 (Updated for Cortellis API v2.0)
**Last Updated:** 2025-12-21
