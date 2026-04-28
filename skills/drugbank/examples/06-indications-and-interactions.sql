-- =====================================================
-- Drug Indications and Drug Interactions Examples
-- =====================================================
-- These examples show how to query drug indications,
-- approved uses, drug-drug interactions, and food interactions.
-- =====================================================

-- ========== INDICATIONS ==========

-- Example 1: Get all indications for a drug
SELECT DISTINCT
    d.drugbank_id,
    d.name AS drug_name,
    c.title AS indication,
    si.kind,
    si.off_label,
    si.country,
    c.snomed_id,
    c.meddra_id,
    c.icd10_id
FROM structured_indications si
JOIN drugs d ON si.drug_id = d.id
JOIN indication_conditions ic ON si.id = ic.indication_id
    AND ic.relationship = 'for_condition'
JOIN conditions c ON ic.condition_id = c.id
WHERE d.name ILIKE '%bevacizumab%'
ORDER BY si.off_label, c.title;

-- Example 2: Find approved (non-off-label) indications
SELECT DISTINCT
    d.drugbank_id,
    d.name AS drug_name,
    c.title AS indication,
    si.country
FROM structured_indications si
JOIN drugs d ON si.drug_id = d.id
JOIN indication_conditions ic ON si.id = ic.indication_id
    AND ic.relationship = 'for_condition'
JOIN conditions c ON ic.condition_id = c.id
WHERE d.name ILIKE '%imatinib%'
  AND si.off_label = false
ORDER BY c.title;

-- Example 3: Find drugs indicated for a specific condition
SELECT DISTINCT
    d.drugbank_id,
    d.name AS drug_name,
    d.state,
    d.type,
    si.off_label,
    si.country
FROM structured_indications si
JOIN drugs d ON si.drug_id = d.id
JOIN indication_conditions ic ON si.id = ic.indication_id
    AND ic.relationship = 'for_condition'
JOIN conditions c ON ic.condition_id = c.id
WHERE c.title ILIKE '%lung cancer%'
  AND d.state = 'approved'
ORDER BY si.off_label, d.name;

-- Example 4: Find off-label uses
SELECT DISTINCT
    d.drugbank_id,
    d.name AS drug_name,
    c.title AS off_label_indication,
    si.country
FROM structured_indications si
JOIN drugs d ON si.drug_id = d.id
JOIN indication_conditions ic ON si.id = ic.indication_id
    AND ic.relationship = 'for_condition'
JOIN conditions c ON ic.condition_id = c.id
WHERE d.drugbank_id = 'DB00619'
  AND si.off_label = true
ORDER BY c.title;

-- Example 5: Count indications per drug
SELECT
    d.drugbank_id,
    d.name AS drug_name,
    d.state,
    COUNT(DISTINCT CASE WHEN si.off_label = false THEN c.id END) AS approved_indications,
    COUNT(DISTINCT CASE WHEN si.off_label = true THEN c.id END) AS off_label_indications,
    COUNT(DISTINCT c.id) AS total_indications
FROM drugs d
LEFT JOIN structured_indications si ON d.id = si.drug_id
LEFT JOIN indication_conditions ic ON si.id = ic.indication_id
    AND ic.relationship = 'for_condition'
LEFT JOIN conditions c ON ic.condition_id = c.id
WHERE d.name ILIKE '%pembrolizumab%'
GROUP BY d.drugbank_id, d.name, d.state;

-- Example 6: Find indications with medical coding (SNOMED, MedDRA, ICD-10)
SELECT DISTINCT
    d.drugbank_id,
    d.name AS drug_name,
    c.title AS indication,
    c.snomed_id,
    c.meddra_id,
    c.icd10_id
FROM structured_indications si
JOIN drugs d ON si.drug_id = d.id
JOIN indication_conditions ic ON si.id = ic.indication_id
    AND ic.relationship = 'for_condition'
JOIN conditions c ON ic.condition_id = c.id
WHERE d.name ILIKE '%nivolumab%'
  AND (c.snomed_id IS NOT NULL OR c.meddra_id IS NOT NULL OR c.icd10_id IS NOT NULL)
ORDER BY c.title;

-- ========== DRUG-DRUG INTERACTIONS ==========

-- Example 7: Get all drug-drug interactions
SELECT DISTINCT
    d1.drugbank_id AS drug1_id,
    d1.name AS drug1_name,
    d2.drugbank_id AS drug2_id,
    d2.name AS drug2_name,
    sdi.severity,
    sdi.description
FROM structured_drug_interactions sdi
JOIN drugs d1 ON sdi.drug_id = d1.id
JOIN drugs d2 ON sdi.affected_drug_id = d2.id
WHERE d1.name ILIKE '%warfarin%'
ORDER BY sdi.severity DESC, d2.name;

-- Example 8: Find severe drug interactions
SELECT DISTINCT
    d1.drugbank_id AS drug1_id,
    d1.name AS drug1_name,
    d2.drugbank_id AS drug2_id,
    d2.name AS drug2_name,
    sdi.severity,
    sdi.description
FROM structured_drug_interactions sdi
JOIN drugs d1 ON sdi.drug_id = d1.id
JOIN drugs d2 ON sdi.affected_drug_id = d2.id
WHERE d1.drugbank_id = 'DB00619'
  AND sdi.severity IN ('major', 'severe')
ORDER BY d2.name;

-- Example 9: Check for interactions between two specific drugs
SELECT
    d1.name AS drug1_name,
    d2.name AS drug2_name,
    sdi.severity,
    sdi.description
FROM structured_drug_interactions sdi
JOIN drugs d1 ON sdi.drug_id = d1.id
JOIN drugs d2 ON sdi.affected_drug_id = d2.id
WHERE (d1.name ILIKE '%aspirin%' AND d2.name ILIKE '%warfarin%')
   OR (d1.name ILIKE '%warfarin%' AND d2.name ILIKE '%aspirin%');

-- Example 10: Count drug interactions by severity
SELECT
    d.drugbank_id,
    d.name AS drug_name,
    sdi.severity,
    COUNT(*) AS interaction_count
FROM structured_drug_interactions sdi
JOIN drugs d ON sdi.drug_id = d.id
WHERE d.name ILIKE '%methotrexate%'
GROUP BY d.drugbank_id, d.name, sdi.severity
ORDER BY sdi.severity DESC, interaction_count DESC;

-- ========== FOOD INTERACTIONS ==========

-- Example 11: Get food interactions for a drug
SELECT
    d.drugbank_id,
    d.name AS drug_name,
    fi.description AS food_interaction
FROM food_interactions fi
JOIN drugs d ON fi.drug_id = d.id
WHERE d.name ILIKE '%imatinib%'
ORDER BY fi.description;

-- Example 12: Find drugs with grapefruit interactions
SELECT DISTINCT
    d.drugbank_id,
    d.name AS drug_name,
    d.state,
    fi.description AS food_interaction
FROM food_interactions fi
JOIN drugs d ON fi.drug_id = d.id
WHERE fi.description ILIKE '%grapefruit%'
  AND d.state = 'approved'
ORDER BY d.name;

-- Example 13: Get all interactions (drug-drug and food) for a drug
SELECT
    'Drug-Drug' AS interaction_type,
    d2.name AS interacting_substance,
    sdi.severity,
    sdi.description
FROM structured_drug_interactions sdi
JOIN drugs d1 ON sdi.drug_id = d1.id
JOIN drugs d2 ON sdi.affected_drug_id = d2.id
WHERE d1.name ILIKE '%erlotinib%'

UNION ALL

SELECT
    'Drug-Food' AS interaction_type,
    'Food' AS interacting_substance,
    NULL AS severity,
    fi.description
FROM food_interactions fi
JOIN drugs d ON fi.drug_id = d.id
WHERE d.name ILIKE '%erlotinib%'
ORDER BY interaction_type, severity DESC NULLS LAST;

-- Example 14: Find drugs with many interactions
SELECT
    d.drugbank_id,
    d.name AS drug_name,
    d.state,
    COUNT(DISTINCT sdi.affected_drug_id) AS drug_interaction_count,
    COUNT(DISTINCT fi.id) AS food_interaction_count
FROM drugs d
LEFT JOIN structured_drug_interactions sdi ON d.id = sdi.drug_id
LEFT JOIN food_interactions fi ON d.id = fi.drug_id
WHERE d.state = 'approved'
GROUP BY d.drugbank_id, d.name, d.state
HAVING COUNT(DISTINCT sdi.affected_drug_id) >= 50
ORDER BY drug_interaction_count DESC
LIMIT 20;

-- Example 15: Get comprehensive drug profile (indications + interactions)
WITH drug_info AS (
    SELECT id, drugbank_id, name, state
    FROM drugs
    WHERE name ILIKE '%cetuximab%'
),
indications AS (
    SELECT
        'Indication' AS data_type,
        c.title AS detail,
        si.off_label::text AS metadata
    FROM drug_info di
    JOIN structured_indications si ON di.id = si.drug_id
    JOIN indication_conditions ic ON si.id = ic.indication_id
        AND ic.relationship = 'for_condition'
    JOIN conditions c ON ic.condition_id = c.id
),
drug_interactions AS (
    SELECT
        'Drug Interaction' AS data_type,
        d2.name AS detail,
        sdi.severity AS metadata
    FROM drug_info di
    JOIN structured_drug_interactions sdi ON di.id = sdi.drug_id
    JOIN drugs d2 ON sdi.affected_drug_id = d2.id
),
food_interactions AS (
    SELECT
        'Food Interaction' AS data_type,
        fi.description AS detail,
        NULL AS metadata
    FROM drug_info di
    JOIN food_interactions fi ON di.id = fi.drug_id
)
SELECT * FROM indications
UNION ALL
SELECT * FROM drug_interactions
UNION ALL
SELECT * FROM food_interactions
ORDER BY data_type, detail;
