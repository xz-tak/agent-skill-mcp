-- =====================================================
-- Basic Drug Search Examples
-- =====================================================
-- These examples show how to search for drugs by name,
-- ID, or other properties in the DrugBank database.
-- =====================================================

-- Example 1: Search for a drug by name (flexible matching)
-- Use ILIKE for case-insensitive substring matching
SELECT
    d.drugbank_id,
    d.name,
    d.type,
    d.state,
    d.description
FROM drugs d
WHERE d.name ILIKE '%bevacizumab%';

-- Example 2: Search by DrugBank ID
SELECT
    d.drugbank_id,
    d.name,
    d.type,
    d.state,
    d.description
FROM drugs d
WHERE d.drugbank_id = 'DB00112';

-- Example 3: Search with synonyms (alternative drug names)
SELECT DISTINCT
    d.drugbank_id,
    d.name AS official_name,
    ds.synonym,
    d.type,
    d.state
FROM drugs d
LEFT JOIN drug_synonyms ds ON d.id = ds.drug_id
WHERE d.name ILIKE '%aspirin%'
   OR ds.synonym ILIKE '%aspirin%';

-- Example 4: Get drug with molecular properties
SELECT
    d.drugbank_id,
    d.name,
    d.type,
    d.state,
    dcp.smiles,
    dcp.inchi,
    dcp.molecular_formula,
    dcp.molecular_weight
FROM drugs d
LEFT JOIN drug_calculated_properties dcp ON d.id = dcp.drug_id
WHERE d.name ILIKE '%imatinib%';

-- Example 5: Find approved drugs only
SELECT
    d.drugbank_id,
    d.name,
    d.type,
    d.description
FROM drugs d
WHERE d.state = 'approved'
  AND d.name ILIKE '%mab%'  -- Search for monoclonal antibodies
ORDER BY d.name
LIMIT 20;

-- Example 6: Find investigational drugs
SELECT
    d.drugbank_id,
    d.name,
    d.type,
    d.description
FROM drugs d
WHERE d.state = 'investigational'
  AND d.type = 'small molecule'
LIMIT 20;

-- Example 7: Search by drug type
SELECT
    d.drugbank_id,
    d.name,
    d.state,
    d.description
FROM drugs d
WHERE d.type = 'biotech'
  AND d.state = 'approved'
ORDER BY d.name
LIMIT 20;

-- Example 8: Broad search across names and IDs
SELECT
    d.drugbank_id,
    d.name,
    d.type,
    d.state
FROM drugs d
WHERE d.name ILIKE '%erlotinib%'
   OR d.drugbank_id ILIKE '%erlotinib%';

-- Example 9: Get multiple drugs by ID list
SELECT
    d.drugbank_id,
    d.name,
    d.type,
    d.state
FROM drugs d
WHERE d.drugbank_id IN ('DB00112', 'DB00619', 'DB01268')
ORDER BY d.name;

-- Example 10: Count drugs by state
SELECT
    d.state,
    COUNT(*) as drug_count
FROM drugs d
GROUP BY d.state
ORDER BY drug_count DESC;
