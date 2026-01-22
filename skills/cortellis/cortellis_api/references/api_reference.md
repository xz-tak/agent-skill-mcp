# Cortellis API Reference

Comprehensive reference for the Cortellis API endpoints used by this skill.

## Authentication

All Cortellis APIs use digest authentication with API key and secret.

### Authentication Endpoints

| API | Auth Version | Auth URL |
|-----|--------------|----------|
| Targets API | v3 | `https://api.cortellis.com/api-ws/ws/rs/auth-v3/token` |
| Biomarkers API | v4 | `https://api.cortellis.com/api-ws/ws/rs/auth-v4/token` |
| Investigational Drugs API | v2 | `https://api.cortellis.com/api-ws/ws/rs/auth-v2/token` |

**Authentication Process:**
1. Make GET request to auth endpoint with HTTPDigestAuth(API_KEY, API_SECRET)
2. Parse response (may be XML or plain text)
3. Extract token from XML `<token>` element or use plain text body
4. Include token in subsequent requests via `API-Token` header

**Example:**
```python
from requests.auth import HTTPDigestAuth
import xml.etree.ElementTree as ET

r = requests.get(AUTH_URL, auth=HTTPDigestAuth(API_KEY, API_SECRET), timeout=20)
root = ET.fromstring(r.text.strip())
token = root.findtext('token')
headers = {'API-Token': token}
```

## Targets API v2

Base URL: `https://api.cortellis.com/api-ws/ws/rs/targets-v2`

### Search for Targets

**Endpoint:** `/target/search`

**Parameters:**
- `query` (required): Query string using field:value syntax
- `hits`: Number of results (default: 30)
- `offset`: Result offset for pagination (default: 0)
- `sortBy`: Sort field (default: "targetNameMain")
- `sortDirection`: "ascending" or "descending"
- `fmt`: Output format ("json" or "xml")

**Query Syntax:**
- `targetSynonyms:BRCA1` - Search by gene symbol or synonym
- `targetId:12345` - Search by target ID

**Example Request:**
```
GET /target/search?query=targetSynonyms:BRCA1&hits=50&sortBy=targetNameMain&fmt=json
```

**Response Structure:**
```json
{
  "TargetResultsOutput": {
    "SearchResults": {
      "TargetResult": [
        {
          "@Id": "4897967856723",
          "targetNameMain": "Breast cancer type 1 susceptibility protein",
          "Symbol": "BRCA1",
          "targetSynonyms": {...}
        }
      ]
    }
  }
}
```

### Get Target Records

**Endpoint:** `/targets`

**Parameters:**
- `idList` (required): Comma-separated target IDs
- `fmt`: Output format ("json" or "xml")

**Example:**
```
GET /targets?idList=4897967856723&fmt=json
```

**Response:**
```json
{
  "TargetRecordsOutput": {
    "Targets": {
      "Target": {
        "@Id": "4897967856723",
        "NameMain": "Breast cancer type 1 susceptibility protein",
        "EntrezgeneIdentifiers": {...},
        "ExternalIdentifiers": {...},
        "Descriptions": {...}
      }
    }
  }
}
```

### Get Target Interactions

**Endpoint:** `/target/interactions`

**Parameters:**
- `idList` (required): Target ID
- `fmt`: Output format

**Response:**
```json
{
  "TargetRecordsOutput": {
    "Targets": {
      "Target": {
        "Interactions": {
          "Interaction": [
            {
              "Type": {"$": "Binding"},
              "PartnerName": {"$": "BARD1"},
              "PartnerGeneSymbol": {"$": "BARD1"},
              "Evidence": {"$": "Experimental"}
            }
          ]
        }
      }
    }
  }
}
```

### Get Condition-Drug Associations

**Endpoint:** `/target/conditionDrugAssociations`

**Parameters:**
- `idList` (required): Target ID
- `fmt`: Output format

**Response:**
```json
{
  "TargetRecordsOutput": {
    "Targets": {
      "Target": {
        "ConditionDrugAssociations": {
          "Condition": [
            {
              "ConditionName": {"$": "Breast cancer"},
              "DrugId": [{"$": "12345"}]
            }
          ]
        }
      }
    }
  }
}
```

### Get Drugs (Basic Records)

**Endpoint:** `/drugs`

**Parameters:**
- `idList` (required): Comma-separated drug IDs
- `fmt`: Output format

**Response:**
```json
{
  "drugRecordsOutput": {
    "Drug": [
      {
        "@id": "12345",
        "@namemain": "olaparib",
        "NamesChemicalAndDescriptions": {
          "Name": "4-[(3-{[4-(Cyclopropylcarbonyl)piperazin-1-yl]carbonyl}-4-fluorophenyl)methyl]phthalazin-1(2H)-one"
        },
        "MechanismsMolecular": {...}
      }
    ]
  }
}
```

### Get Trials

**Endpoint:** `/target/trials`

**Parameters:**
- `idList` (required): Comma-separated trial IDs
- `fmt`: Output format

### Get Condition-Gene Associations

**Endpoint:** `/target/conditionGeneAssociations`

**Parameters:**
- `idList` (required): Target ID
- `fmt`: Output format

### Get Condition-Gene-Variant Associations

**Endpoint:** `/target/conditionGeneVariantAssociations`

**Parameters:**
- `idList` (required): Target ID
- `fmt`: Output format

## Biomarkers API v3

Base URL: `https://api.cortellis.com/api-ws/ws/rs/biomarkers-v3`

### Search Biomarker Uses

**Endpoint:** `/biomarkerUse/search`

**Parameters:**
- `query` (required): Query string
- `hits`: Number of results
- `offset`: Result offset
- `sortDirection`: Sort direction
- `fmt`: Output format

**Query Syntax:**
- `biomarkerSynonyms:BRCA1` - Search by biomarker name or synonym

**Example:**
```
GET /biomarkerUse/search?query=biomarkerSynonyms:BRCA1&hits=1000&fmt=json
```

**Response:**
```json
{
  "biomarkerUseResultsOutput": {
    "SearchResults": {
      "BiomarkerUseResult": [
        {
          "@id": "123456",
          "Biomarker": {
            "@id": "789",
            "@mainName": "BRCA1"
          },
          "BiomarkerApplication": {"$": "Diagnostic"}
        }
      ]
    }
  }
}
```

### Get Biomarker Use Records

**Endpoint:** `/biomarkerUses`

**Parameters:**
- `idList` (required): Comma-separated biomarker use IDs
- `fmt`: Output format

### Get Biomarker Records

**Endpoint:** `/biomarkers`

**Parameters:**
- `idList` (required): Comma-separated biomarker IDs
- `fmt`: Output format

## Investigational Drugs API v2.0

Base URL: `https://api.cortellis.com/api-ws/ws/rs/drugs-v2`

### Search Drugs

**Endpoint:** `/drug/search`

**Parameters:**
- `query` (required): Query string using field:value syntax
- `hits`: Number of results
- `fmt`: Output format

**Query Syntax:**
- `drugNamesAll:olaparib` - Search across all drug name fields
- `drugId:12345` - Search by drug ID

**Important:** Use `drugNamesAll` field for name searches. The `drugName` field is not valid.

**Example:**
```
GET /drug/search?query=drugNamesAll:olaparib&hits=5&fmt=json
```

**Response:**
```json
{
  "drugResultsOutput": {
    "SearchResults": {
      "Drug": {
        "@id": "12004",
        "DrugName": "olaparib",
        "PhaseHighest": "Approved"
      }
    }
  }
}
```

**Note:** Response may contain `SearchResults.Drug` (single result as dict) or `SearchResults.DrugResult` (multiple results as list). The script handles both cases.

### Get Drug Record (Single)

**Endpoint:** `/drug/{drugId}`

**Parameters:**
- `drugId` (in path): Drug ID
- `fmt`: Output format (query parameter)

**Example:**
```
GET /drug/12004?fmt=json
```

**Response:**
```json
{
  "drugRecordOutput": {
    "@id": "12004",
    "DrugName": "olaparib",
    "PhaseHighest": {
      "@id": "LA",
      "$": "Launched"
    },
    "CompanyOriginator": {
      "@id": "company_id",
      "$": "AstraZeneca"
    },
    "CompaniesPrimary": {
      "Company": [
        {"@id": "1", "$": "AstraZeneca"},
        {"@id": "2", "$": "Merck & Co"}
      ]
    },
    "IndicationsPrimary": {
      "Indication": [
        {"@id": "1767", "$": "Ovarian cancer"},
        {"@id": "213", "$": "Breast cancer"}
      ]
    },
    "ActionsPrimary": {
      "Action": [
        {"@id": "0", "$": "PARP-1 inhibitor"},
        {"@id": "1", "$": "PARP-2 inhibitor"}
      ]
    },
    "TherapyAreas": {
      "TherapyArea": "Cancer"
    },
    "Technologies": {
      "Technology": [
        {"@id": "0", "$": "Small molecule"}
      ]
    },
    "DrugSynonyms": {
      "Name": [
        {"Value": "Lynparza"},
        {"Value": "AZD-2281"}
      ]
    }
  }
}
```

**Note:** Values are directly in the `$` field, not nested under field-specific names like `IndicationName`.

### Get Drug Records (Bulk)

**Endpoint:** `/drugs`

**Parameters:**
- `idList` (required): Comma-separated drug IDs (supports batching up to 30+ IDs)
- `fmt`: Output format

**Example:**
```
GET /drugs?idList=12004,4838,106955&fmt=json
```

**Response:**
```json
{
  "drugRecordsOutput": {
    "Drug": [
      {
        "@id": "12004",
        "DrugName": "olaparib",
        "PhaseHighest": {"@id": "LA", "$": "Launched"},
        ...
      },
      {
        "@id": "4838",
        "DrugName": "natalizumab",
        "PhaseHighest": {"@id": "LA", "$": "Launched"},
        ...
      }
    ]
  }
}
```

**Best Practice:** For querying multiple drugs (common when processing target results), use bulk endpoint with chunking (30 IDs per request) instead of individual `/drug/{drugId}` calls.

## Response Data Structures

### Value Extraction Pattern

Many fields use a `{"$": "value"}` pattern. The script includes an `extract_value()` helper:

```python
def extract_value(obj, default=''):
    if isinstance(obj, dict):
        return obj.get('$', default)
    return obj if obj is not None else default
```

### List Handling

Single vs. multiple results are inconsistent. The script uses `as_list()` helper:

```python
def as_list(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]
```

## Error Handling

### HTTP Status Codes

- **200**: Success
- **400**: Bad request (invalid query syntax)
- **403**: Invalid or expired token
- **404**: Resource not found
- **500**: Internal server error
- **Timeout**: Connection or read timeout (30 seconds)

### Common Issues

1. **Invalid Token (403)**
   - Re-authenticate to get fresh token
   - Token may have expired (typically 2-3 hours)

2. **Drug Not Found (404/500)**
   - Drug ID from Targets API may not exist in Investigational Drugs API
   - Research compounds may not be in commercial drug database
   - Fallback to basic records from Targets API

3. **Search Returns No Results**
   - Drug name from Targets API may not match Investigational Drugs API naming
   - Try alternative drug names or synonyms
   - Some drugs only available via basic Targets API records

4. **Timeout**
   - API may be slow or unavailable
   - Retry with exponential backoff
   - Consider reducing `hits` parameter for searches

## Rate Limiting

- No explicit rate limits documented
- Use reasonable delays between requests for large batch operations
- Connection pooling and session reuse recommended for multiple requests

## Best Practices

1. **Authentication**: Cache tokens (valid for ~2 hours), reuse across requests
2. **Chunking**: Batch drug/trial ID requests in chunks of 30 (optimal for Cortellis API)
3. **Error Handling**: Continue processing remaining items if one fails
4. **Pagination**: Use `offset` parameter for large result sets
5. **Search Optimization**: Use specific field searches (e.g., `targetSynonyms:`) rather than free text
6. **Related Targets**: Query related targets (isoforms, variants) for comprehensive drug coverage
7. **Async Requests**: Use concurrent requests with adaptive concurrency for optimal performance

## Performance Optimization

### Async Bulk Fetching

The `cortellis_gene_query.py` script implements async bulk fetching with adaptive concurrency:

**Key Features:**
- **Adaptive Concurrency**: Starts at 25 concurrent requests, backs off to 20→15→10→5 if failures detected
- **Exponential Backoff**: 2s → 5s → 10s → 15s delays between retry stages
- **Aggressive Retry**: 1% failure threshold triggers retry cascade
- **Chunk Size**: 30 items per batch for optimal throughput

**Performance Gains:**
- **11.8x faster** than sequential requests
- **99.5%+ data completeness** with retry logic
- Processes 1700+ drugs in ~18 seconds (vs ~3.5 minutes sequential)

**Implementation:**
```python
# Async fetch with adaptive concurrency
results = fetch_with_adaptive_concurrency(
    urls=chunk_urls,
    headers=DRUGS_HEADER,
    initial_concurrent=25,
    fallback_concurrent=20,
    timeout_seconds=120
)
```

See [`PERFORMANCE_ANALYSIS_REPORT.md`](../PERFORMANCE_ANALYSIS_REPORT.md) for detailed benchmarks.

### Related Targets Querying

**Important:** To capture all drugs associated with a target, the script automatically queries related targets:

**Related Target Types:**
- **Isoforms**: Alternative splice variants of the same gene
- **Variants**: Genetic variants or mutant forms
- **Complexes**: Multi-protein complexes containing the target

**Implementation:**
```python
# Extract related target IDs from annotation
related_targets = annotation.get('RelatedTargets', {})
related_ids = as_list(related_targets.get('Id', []))

# Query each related target for comprehensive drug coverage
for related_id in related_ids:
    # Fetch drugs for related target
    ...
```

**Example:** Querying ITGA4 also retrieves drugs for:
- ITGA4 isoforms
- ITGA4/ITGB1 heterodimer complex
- Related integrin family members

This ensures comprehensive drug coverage without missing drugs that target specific protein forms.

## API Limitations

1. **Drug ID Mismatch**: Drug IDs from Targets API are not always valid in Investigational Drugs API
2. **Name Variations**: Drug names may differ between APIs (use search to match)
3. **Incomplete Data**: Not all drugs have comprehensive development information
4. **Response Inconsistency**: Single results may be dict, multiple results may be list
5. **Field Naming**: Some fields use different names across APIs (e.g., `@Id` vs `@id`)
