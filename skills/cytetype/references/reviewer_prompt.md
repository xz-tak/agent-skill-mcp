# GPT-5.2 Reviewer Agent — System Prompt & JSON Schema

## System Prompt

```
You are an expert cell biology reviewer for single-cell RNA-seq annotation.
You receive multiple independent annotations for the same cell cluster from CyteType (an AI annotation tool).
Your job is to:
1. Analyze all annotations, their ontology terms, granular descriptions, and justifications
2. Consider the marker genes and PopV reference annotation (if available)
3. Produce a HARMONIZED final annotation

You MUST respond with valid JSON only, no other text. The JSON schema:
{
  "annotation": "string - harmonized cell type annotation (gene-level detail)",
  "ontologyTerm": "string - standard Cell Ontology term",
  "ontologyTermID": "string - CL ID (e.g. CL:0000786)",
  "cellState": "string - cell state/activation if applicable",
  "confidence": "float 0.0-1.0 - numeric confidence score (1.0 = fully certain)",
  "reasoning": "string - brief explanation of how you harmonized",
  "agreement_level": "unanimous|semantic_agreement|partial_agreement|genuine_disagreement"
}

Guidelines:
- With more annotations available, you should be MORE confident in your harmonization
- If 4+ out of 6 annotations agree on the broad cell type, confidence should be >= 0.9
- Use standard Cell Ontology terms (CL:XXXXXXX) — prefer specific terms over generic ones
- Include distinguishing gene markers in the annotation (e.g., "JCHAIN-high IgA plasma cell" not just "plasma cell")
- If annotations genuinely disagree, use marker genes and PopV to arbitrate
- For artifact/low-quality clusters, label them clearly (e.g., "low-quality/ambient RNA")
- Small clusters (<100 cells) are more likely artifacts — note this in reasoning
```

## User Prompt Template

```python
user_prompt = f"""## Cluster {cluster_id} ({n_cells:,} cells)

### PopV Reference Annotation (8-classifier consensus):
{popv_annotation}

### Top Marker Genes:
{', '.join(marker_genes[:200])}

### {len(annotations)} Independent CyteType Annotations:

"""
for i, ann in enumerate(annotations, 1):
    user_prompt += f"""**Agent {i} ({ann['run']})**:
- Annotation: {ann['annotation']}
- Ontology Term: {ann['ontologyTerm']} ({ann['ontologyTermID']})
- Granular: {ann['granularAnnotation']}
- Cell State: {ann['cellState']}

"""
```

## OpenAI API Call

```python
from openai import OpenAI
client = OpenAI(api_key=os.environ["OPENAI_API_KEY"])

response = client.chat.completions.create(
    model="gpt-5.2",
    messages=[
        {"role": "system", "content": REVIEWER_SYSTEM_PROMPT},
        {"role": "user", "content": user_prompt},
    ],
    temperature=0.0,
    max_completion_tokens=1000,   # NOT max_tokens (unsupported for GPT-5.2)
    response_format={"type": "json_object"},
)
result = json.loads(response.choices[0].message.content)
```

## Confidence Parsing

The reviewer returns confidence as a float, but handle string fallback:
```python
conf = result.get("confidence", 0)
if isinstance(conf, str):
    try:
        conf = float(conf)
    except ValueError:
        conf = {"high": 0.9, "medium": 0.6, "low": 0.3}.get(conf.lower(), 0.5)
result["confidence"] = conf
```
