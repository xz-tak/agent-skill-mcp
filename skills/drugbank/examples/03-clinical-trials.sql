-- =====================================================
-- Clinical Trial Query Examples
-- =====================================================
-- These examples show how to query clinical trial data
-- including trial status, interventions, and timelines.
-- =====================================================

-- Example 1: Get all clinical trials for a drug
SELECT DISTINCT
    ct.identifier,
    ct.title,
    ct.status,
    ct.purpose,
    ct.start_date,
    ct.end_date,
    d.drugbank_id,
    d.name AS drug_name
FROM clinical_trials ct
JOIN clinical_trial_interventions cti ON ct.identifier = cti.trial_id
JOIN clinical_trial_interventions_drugs ctid ON cti.id = ctid.intervention_id
JOIN drugs d ON ctid.drug_id = d.id
WHERE d.name ILIKE '%bevacizumab%'
ORDER BY ct.start_date DESC;

-- Example 2: Find active/recruiting trials for a drug
SELECT DISTINCT
    ct.identifier,
    ct.title,
    ct.status,
    ct.purpose,
    ct.start_date,
    d.name AS drug_name
FROM clinical_trials ct
JOIN clinical_trial_interventions cti ON ct.identifier = cti.trial_id
JOIN clinical_trial_interventions_drugs ctid ON cti.id = ctid.intervention_id
JOIN drugs d ON ctid.drug_id = d.id
WHERE d.name ILIKE '%imatinib%'
  AND ct.status IN ('recruiting', 'active, not recruiting', 'enrolling by invitation')
ORDER BY ct.start_date DESC;

-- Example 3: Get trials by DrugBank ID
SELECT DISTINCT
    ct.identifier,
    ct.title,
    ct.status,
    ct.purpose,
    ct.start_date,
    ct.end_date,
    d.drugbank_id,
    d.name AS drug_name
FROM clinical_trials ct
JOIN clinical_trial_interventions cti ON ct.identifier = cti.trial_id
JOIN clinical_trial_interventions_drugs ctid ON cti.id = ctid.intervention_id
JOIN drugs d ON ctid.drug_id = d.id
WHERE d.drugbank_id = 'DB00619'
ORDER BY ct.start_date DESC;

-- Example 4: Find trials by purpose (phase)
SELECT DISTINCT
    ct.identifier,
    ct.title,
    ct.status,
    ct.purpose,
    ct.start_date,
    d.name AS drug_name
FROM clinical_trials ct
JOIN clinical_trial_interventions cti ON ct.identifier = cti.trial_id
JOIN clinical_trial_interventions_drugs ctid ON cti.id = ctid.intervention_id
JOIN drugs d ON ctid.drug_id = d.id
WHERE d.name ILIKE '%cetuximab%'
  AND ct.purpose ILIKE '%phase 3%'
ORDER BY ct.start_date DESC;

-- Example 5: Get completed trials with dates
SELECT DISTINCT
    ct.identifier,
    ct.title,
    ct.status,
    ct.purpose,
    ct.start_date,
    ct.end_date,
    d.name AS drug_name,
    (ct.end_date - ct.start_date) AS trial_duration_days
FROM clinical_trials ct
JOIN clinical_trial_interventions cti ON ct.identifier = cti.trial_id
JOIN clinical_trial_interventions_drugs ctid ON cti.id = ctid.intervention_id
JOIN drugs d ON ctid.drug_id = d.id
WHERE d.name ILIKE '%erlotinib%'
  AND ct.status = 'completed'
  AND ct.end_date IS NOT NULL
ORDER BY ct.end_date DESC;

-- Example 6: Find trials with intervention details
SELECT DISTINCT
    ct.identifier,
    ct.title,
    ct.status,
    cti.intervention_type,
    cti.intervention_name,
    d.name AS drug_name
FROM clinical_trials ct
JOIN clinical_trial_interventions cti ON ct.identifier = cti.trial_id
JOIN clinical_trial_interventions_drugs ctid ON cti.id = ctid.intervention_id
JOIN drugs d ON ctid.drug_id = d.id
WHERE d.drugbank_id = 'DB00112'
ORDER BY ct.start_date DESC
LIMIT 20;

-- Example 7: Count trials by status for a drug
SELECT
    ct.status,
    COUNT(DISTINCT ct.identifier) AS trial_count
FROM clinical_trials ct
JOIN clinical_trial_interventions cti ON ct.identifier = cti.trial_id
JOIN clinical_trial_interventions_drugs ctid ON cti.id = ctid.intervention_id
JOIN drugs d ON ctid.drug_id = d.id
WHERE d.name ILIKE '%paclitaxel%'
GROUP BY ct.status
ORDER BY trial_count DESC;

-- Example 8: Find recent trials (last 5 years)
SELECT DISTINCT
    ct.identifier,
    ct.title,
    ct.status,
    ct.purpose,
    ct.start_date,
    d.name AS drug_name
FROM clinical_trials ct
JOIN clinical_trial_interventions cti ON ct.identifier = cti.trial_id
JOIN clinical_trial_interventions_drugs ctid ON cti.id = ctid.intervention_id
JOIN drugs d ON ctid.drug_id = d.id
WHERE d.name ILIKE '%pembrolizumab%'
  AND ct.start_date >= CURRENT_DATE - INTERVAL '5 years'
ORDER BY ct.start_date DESC;

-- Example 9: Find trials by purpose keyword
SELECT DISTINCT
    ct.identifier,
    ct.title,
    ct.status,
    ct.purpose,
    ct.start_date,
    d.name AS drug_name
FROM clinical_trials ct
JOIN clinical_trial_interventions cti ON ct.identifier = cti.trial_id
JOIN clinical_trial_interventions_drugs ctid ON cti.id = ctid.intervention_id
JOIN drugs d ON ctid.drug_id = d.id
WHERE ct.purpose ILIKE '%lung cancer%'
  AND d.state = 'approved'
ORDER BY ct.start_date DESC
LIMIT 30;

-- Example 10: Get trial summary statistics
SELECT
    d.name AS drug_name,
    d.state,
    COUNT(DISTINCT ct.identifier) AS total_trials,
    COUNT(DISTINCT CASE WHEN ct.status = 'completed' THEN ct.identifier END) AS completed_trials,
    COUNT(DISTINCT CASE WHEN ct.status LIKE '%recruiting%' THEN ct.identifier END) AS recruiting_trials,
    MIN(ct.start_date) AS first_trial_date,
    MAX(ct.start_date) AS latest_trial_date
FROM drugs d
LEFT JOIN clinical_trial_interventions_drugs ctid ON d.id = ctid.drug_id
LEFT JOIN clinical_trial_interventions cti ON ctid.intervention_id = cti.id
LEFT JOIN clinical_trials ct ON cti.trial_id = ct.identifier
WHERE d.name ILIKE '%nivolumab%'
GROUP BY d.name, d.state;
