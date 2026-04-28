-- =====================================================
-- Adverse Effects Query Examples
-- =====================================================
-- These examples show how to query adverse effects data
-- including frequency, severity, and condition mappings.
-- =====================================================

-- Example 1: Get all adverse effects for a drug with frequency
SELECT
    d.drugbank_id,
    d.name AS drug_name,
    c.title AS adverse_effect,
    aei.percent AS frequency_percent,
    aei.kind AS frequency_kind,
    sae.severity,
    c.snomed_id,
    c.meddra_id
FROM structured_adverse_effects sae
JOIN drugs d ON sae.drug_id = d.id
JOIN adverse_effect_conditions aec ON sae.id = aec.adverse_effect_id
    AND aec.relationship = 'effect'
JOIN conditions c ON aec.condition_id = c.id
LEFT JOIN adverse_effect_incidences aei ON sae.id = aei.adverse_effect_id
WHERE d.name ILIKE '%erlotinib%'
ORDER BY aei.percent DESC NULLS LAST;

-- Example 2: Find common adverse effects (frequency >= 10%)
SELECT
    d.drugbank_id,
    d.name AS drug_name,
    c.title AS adverse_effect,
    aei.percent AS frequency_percent,
    sae.severity
FROM structured_adverse_effects sae
JOIN drugs d ON sae.drug_id = d.id
JOIN adverse_effect_conditions aec ON sae.id = aec.adverse_effect_id
    AND aec.relationship = 'effect'
JOIN conditions c ON aec.condition_id = c.id
LEFT JOIN adverse_effect_incidences aei ON sae.id = aei.adverse_effect_id
WHERE d.name ILIKE '%imatinib%'
  AND aei.percent >= 10.0
ORDER BY aei.percent DESC;

-- Example 3: Find severe adverse effects
SELECT
    d.drugbank_id,
    d.name AS drug_name,
    c.title AS adverse_effect,
    sae.severity,
    aei.percent AS frequency_percent
FROM structured_adverse_effects sae
JOIN drugs d ON sae.drug_id = d.id
JOIN adverse_effect_conditions aec ON sae.id = aec.adverse_effect_id
    AND aec.relationship = 'effect'
JOIN conditions c ON aec.condition_id = c.id
LEFT JOIN adverse_effect_incidences aei ON sae.id = aei.adverse_effect_id
WHERE d.name ILIKE '%bevacizumab%'
  AND sae.severity IN ('severe', 'major')
ORDER BY aei.percent DESC NULLS LAST;

-- Example 4: Compare adverse effects across drugs
SELECT
    d.name AS drug_name,
    c.title AS adverse_effect,
    aei.percent AS frequency_percent,
    sae.severity
FROM structured_adverse_effects sae
JOIN drugs d ON sae.drug_id = d.id
JOIN adverse_effect_conditions aec ON sae.id = aec.adverse_effect_id
    AND aec.relationship = 'effect'
JOIN conditions c ON aec.condition_id = c.id
LEFT JOIN adverse_effect_incidences aei ON sae.id = aei.adverse_effect_id
WHERE d.name IN ('Cetuximab', 'Panitumumab')
  AND c.title ILIKE '%rash%'
ORDER BY d.name, aei.percent DESC NULLS LAST;

-- Example 5: Get adverse effects with MedDRA codes
SELECT
    d.drugbank_id,
    d.name AS drug_name,
    c.title AS adverse_effect,
    c.meddra_id,
    aei.percent AS frequency_percent,
    sae.severity
FROM structured_adverse_effects sae
JOIN drugs d ON sae.drug_id = d.id
JOIN adverse_effect_conditions aec ON sae.id = aec.adverse_effect_id
    AND aec.relationship = 'effect'
JOIN conditions c ON aec.condition_id = c.id
LEFT JOIN adverse_effect_incidences aei ON sae.id = aei.adverse_effect_id
WHERE d.drugbank_id = 'DB00619'
  AND c.meddra_id IS NOT NULL
ORDER BY aei.percent DESC NULLS LAST;

-- Example 6: Find drugs with specific adverse effect
SELECT DISTINCT
    d.drugbank_id,
    d.name AS drug_name,
    d.state,
    c.title AS adverse_effect,
    aei.percent AS frequency_percent
FROM structured_adverse_effects sae
JOIN drugs d ON sae.drug_id = d.id
JOIN adverse_effect_conditions aec ON sae.id = aec.adverse_effect_id
    AND aec.relationship = 'effect'
JOIN conditions c ON aec.condition_id = c.id
LEFT JOIN adverse_effect_incidences aei ON sae.id = aei.adverse_effect_id
WHERE c.title ILIKE '%neutropenia%'
  AND d.state = 'approved'
ORDER BY aei.percent DESC NULLS LAST
LIMIT 20;

-- Example 7: Count adverse effects per drug
SELECT
    d.drugbank_id,
    d.name AS drug_name,
    d.state,
    COUNT(DISTINCT c.id) AS adverse_effect_count,
    COUNT(DISTINCT CASE WHEN sae.severity IN ('severe', 'major') THEN c.id END) AS severe_ae_count
FROM drugs d
LEFT JOIN structured_adverse_effects sae ON d.id = sae.drug_id
LEFT JOIN adverse_effect_conditions aec ON sae.id = aec.adverse_effect_id
    AND aec.relationship = 'effect'
LEFT JOIN conditions c ON aec.condition_id = c.id
WHERE d.name ILIKE '%paclitaxel%'
GROUP BY d.drugbank_id, d.name, d.state;

-- Example 8: Get adverse effects by frequency range
SELECT
    d.drugbank_id,
    d.name AS drug_name,
    c.title AS adverse_effect,
    aei.percent AS frequency_percent,
    sae.severity,
    CASE
        WHEN aei.percent >= 25 THEN 'Very Common (≥25%)'
        WHEN aei.percent >= 10 THEN 'Common (10-25%)'
        WHEN aei.percent >= 1 THEN 'Uncommon (1-10%)'
        WHEN aei.percent < 1 THEN 'Rare (<1%)'
        ELSE 'Unknown'
    END AS frequency_category
FROM structured_adverse_effects sae
JOIN drugs d ON sae.drug_id = d.id
JOIN adverse_effect_conditions aec ON sae.id = aec.adverse_effect_id
    AND aec.relationship = 'effect'
JOIN conditions c ON aec.condition_id = c.id
LEFT JOIN adverse_effect_incidences aei ON sae.id = aei.adverse_effect_id
WHERE d.name ILIKE '%pembrolizumab%'
ORDER BY aei.percent DESC NULLS LAST;

-- Example 9: Find adverse effects with SNOMED codes
SELECT
    d.drugbank_id,
    d.name AS drug_name,
    c.title AS adverse_effect,
    c.snomed_id,
    aei.percent AS frequency_percent
FROM structured_adverse_effects sae
JOIN drugs d ON sae.drug_id = d.id
JOIN adverse_effect_conditions aec ON sae.id = aec.adverse_effect_id
    AND aec.relationship = 'effect'
JOIN conditions c ON aec.condition_id = c.id
LEFT JOIN adverse_effect_incidences aei ON sae.id = aei.adverse_effect_id
WHERE d.drugbank_id = 'DB00112'
  AND c.snomed_id IS NOT NULL
ORDER BY aei.percent DESC NULLS LAST;

-- Example 10: Get adverse effect summary statistics
SELECT
    d.name AS drug_name,
    COUNT(DISTINCT c.id) AS total_adverse_effects,
    COUNT(DISTINCT CASE WHEN aei.percent IS NOT NULL THEN c.id END) AS with_frequency_data,
    COUNT(DISTINCT CASE WHEN sae.severity = 'severe' THEN c.id END) AS severe_effects,
    AVG(aei.percent) AS avg_frequency_percent,
    MAX(aei.percent) AS max_frequency_percent
FROM drugs d
LEFT JOIN structured_adverse_effects sae ON d.id = sae.drug_id
LEFT JOIN adverse_effect_conditions aec ON sae.id = aec.adverse_effect_id
    AND aec.relationship = 'effect'
LEFT JOIN conditions c ON aec.condition_id = c.id
LEFT JOIN adverse_effect_incidences aei ON sae.id = aei.adverse_effect_id
WHERE d.name ILIKE '%nivolumab%'
GROUP BY d.name;
