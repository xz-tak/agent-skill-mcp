-- =====================================================
-- Multi-Target Drug Discovery Examples
-- =====================================================
-- These examples show how to find drugs that target
-- multiple proteins simultaneously - crucial for
-- combination therapy and polypharmacology analysis.
-- =====================================================

-- Example 1: Find drugs targeting BOTH EGFR AND HER2
WITH egfr_drugs AS (
    SELECT DISTINCT d.drugbank_id
    FROM drugs d
    JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
    JOIN bio_entities be ON b.biodb_id = be.biodb_id
    LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
    LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
    WHERE be.name ILIKE '%epidermal growth factor receptor%'
       OR be.name ILIKE '%EGFR%'
       OR p.gene_name ILIKE '%EGFR%'
),
her2_drugs AS (
    SELECT DISTINCT d.drugbank_id
    FROM drugs d
    JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
    JOIN bio_entities be ON b.biodb_id = be.biodb_id
    LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
    LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
    WHERE be.name ILIKE '%receptor tyrosine-protein kinase erbB-2%'
       OR be.name ILIKE '%HER2%'
       OR p.gene_name ILIKE '%ERBB2%'
)
SELECT
    d.drugbank_id,
    d.name,
    d.type,
    d.state,
    d.description
FROM drugs d
WHERE d.drugbank_id IN (SELECT drugbank_id FROM egfr_drugs)
  AND d.drugbank_id IN (SELECT drugbank_id FROM her2_drugs)
ORDER BY d.state DESC, d.name;

-- Example 2: Find drugs targeting PD-1 OR PD-L1 (union)
WITH pd1_drugs AS (
    SELECT DISTINCT d.drugbank_id
    FROM drugs d
    JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
    JOIN bio_entities be ON b.biodb_id = be.biodb_id
    LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
    LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
    WHERE be.name ILIKE '%programmed cell death protein 1%'
       OR be.name ILIKE '%PD-1%'
       OR p.gene_name ILIKE '%PDCD1%'
),
pdl1_drugs AS (
    SELECT DISTINCT d.drugbank_id
    FROM drugs d
    JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
    JOIN bio_entities be ON b.biodb_id = be.biodb_id
    LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
    LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
    WHERE be.name ILIKE '%programmed cell death 1 ligand 1%'
       OR be.name ILIKE '%PD-L1%'
       OR p.gene_name ILIKE '%CD274%'
)
SELECT
    d.drugbank_id,
    d.name,
    d.type,
    d.state,
    d.description
FROM drugs d
WHERE d.drugbank_id IN (
    SELECT drugbank_id FROM pd1_drugs
    UNION
    SELECT drugbank_id FROM pdl1_drugs
)
ORDER BY d.state DESC, d.name;

-- Example 3: Find approved drugs targeting 3+ kinases
WITH kinase_targets AS (
    SELECT DISTINCT
        d.drugbank_id,
        d.name,
        d.state,
        COUNT(DISTINCT be.biodb_id) AS kinase_target_count
    FROM drugs d
    JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
    JOIN bio_entities be ON b.biodb_id = be.biodb_id
    LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
    LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
    WHERE (be.name ILIKE '%kinase%' OR p.gene_name ILIKE '%kinase%')
      AND d.state = 'approved'
    GROUP BY d.drugbank_id, d.name, d.state
    HAVING COUNT(DISTINCT be.biodb_id) >= 3
)
SELECT
    kt.drugbank_id,
    kt.name,
    kt.state,
    kt.kinase_target_count,
    STRING_AGG(DISTINCT be.name, '; ' ORDER BY be.name) AS target_names
FROM kinase_targets kt
JOIN drugs d ON kt.drugbank_id = d.drugbank_id
JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
JOIN bio_entities be ON b.biodb_id = be.biodb_id
LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
WHERE be.name ILIKE '%kinase%' OR p.gene_name ILIKE '%kinase%'
GROUP BY kt.drugbank_id, kt.name, kt.state, kt.kinase_target_count
ORDER BY kt.kinase_target_count DESC, kt.name;

-- Example 4: Find drugs targeting VEGFR family (VEGFR1, VEGFR2, VEGFR3)
WITH vegfr_drugs AS (
    SELECT DISTINCT
        d.drugbank_id,
        d.name,
        d.state,
        be.name AS target_name,
        p.gene_name
    FROM drugs d
    JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
    JOIN bio_entities be ON b.biodb_id = be.biodb_id
    LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
    LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
    WHERE be.name ILIKE '%vascular endothelial growth factor receptor%'
       OR p.gene_name IN ('FLT1', 'KDR', 'FLT4')
       OR p.gene_name ILIKE '%VEGFR%'
)
SELECT
    drugbank_id,
    name,
    state,
    COUNT(DISTINCT target_name) AS vegfr_count,
    STRING_AGG(DISTINCT gene_name, ', ' ORDER BY gene_name) AS gene_names
FROM vegfr_drugs
GROUP BY drugbank_id, name, state
ORDER BY vegfr_count DESC, state DESC, name;

-- Example 5: Find drugs targeting PI3K pathway components
WITH pi3k_targets AS (
    'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG',  -- PI3K isoforms
    'AKT1', 'AKT2', 'AKT3',                   -- AKT isoforms
    'MTOR'                                     -- mTOR
)
SELECT DISTINCT
    d.drugbank_id,
    d.name,
    d.state,
    p.gene_name,
    be.name AS target_name,
    b.inhibitor
FROM drugs d
JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
JOIN bio_entities be ON b.biodb_id = be.biodb_id
LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
WHERE p.gene_name IN ('PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'AKT1', 'AKT2', 'AKT3', 'MTOR')
  AND d.state IN ('approved', 'investigational')
ORDER BY d.state DESC, d.name;

-- Example 6: Find drugs with detailed multi-target profile
WITH drug_targets AS (
    SELECT
        d.drugbank_id,
        d.name,
        d.state,
        d.type,
        be.name AS target_name,
        p.gene_name,
        b.inhibitor,
        b.agonist,
        b.antagonist
    FROM drugs d
    JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
    JOIN bio_entities be ON b.biodb_id = be.biodb_id
    LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
    LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
    WHERE d.name ILIKE '%imatinib%'
)
SELECT
    drugbank_id,
    name,
    state,
    type,
    target_name,
    gene_name,
    inhibitor,
    agonist,
    antagonist
FROM drug_targets
ORDER BY target_name;

-- Example 7: Find drugs targeting both a kinase AND a checkpoint protein
WITH kinase_drugs AS (
    SELECT DISTINCT d.drugbank_id
    FROM drugs d
    JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
    JOIN bio_entities be ON b.biodb_id = be.biodb_id
    LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
    LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
    WHERE (be.name ILIKE '%kinase%' OR p.gene_name ILIKE '%kinase%')
      AND b.inhibitor = true
),
checkpoint_drugs AS (
    SELECT DISTINCT d.drugbank_id
    FROM drugs d
    JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
    JOIN bio_entities be ON b.biodb_id = be.biodb_id
    WHERE be.name ILIKE '%PD-1%'
       OR be.name ILIKE '%PD-L1%'
       OR be.name ILIKE '%CTLA%'
)
SELECT
    d.drugbank_id,
    d.name,
    d.state,
    d.description
FROM drugs d
WHERE d.drugbank_id IN (SELECT drugbank_id FROM kinase_drugs)
  AND d.drugbank_id IN (SELECT drugbank_id FROM checkpoint_drugs);

-- Example 8: Find drugs targeting specific combination (e.g., BRAF + MEK)
WITH braf_drugs AS (
    SELECT DISTINCT d.drugbank_id
    FROM drugs d
    JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
    JOIN bio_entities be ON b.biodb_id = be.biodb_id
    LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
    LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
    WHERE p.gene_name ILIKE '%BRAF%'
       OR be.name ILIKE '%serine/threonine-protein kinase B-raf%'
),
mek_drugs AS (
    SELECT DISTINCT d.drugbank_id
    FROM drugs d
    JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
    JOIN bio_entities be ON b.biodb_id = be.biodb_id
    LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
    LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
    WHERE p.gene_name IN ('MAP2K1', 'MAP2K2')
       OR be.name ILIKE '%MAP kinase kinase%'
)
SELECT
    d.drugbank_id,
    d.name,
    d.state,
    d.type,
    d.description
FROM drugs d
WHERE d.drugbank_id IN (SELECT drugbank_id FROM braf_drugs)
  AND d.drugbank_id IN (SELECT drugbank_id FROM mek_drugs)
ORDER BY d.state DESC, d.name;

-- Example 9: Compare target profiles of similar drugs
WITH drug_list AS (
    SELECT id, drugbank_id, name
    FROM drugs
    WHERE name IN ('Imatinib', 'Dasatinib', 'Nilotinib')
)
SELECT
    dl.name AS drug_name,
    dl.drugbank_id,
    p.gene_name,
    be.name AS target_name,
    b.inhibitor,
    b.pharmacological_action
FROM drug_list dl
JOIN bonds b ON dl.id = b.drug_id AND b.type = 'TargetBond'
JOIN bio_entities be ON b.biodb_id = be.biodb_id
LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
ORDER BY p.gene_name, dl.name;

-- Example 10: Find selective vs non-selective inhibitors
WITH target_counts AS (
    SELECT
        d.drugbank_id,
        d.name,
        d.state,
        COUNT(DISTINCT be.biodb_id) AS target_count
    FROM drugs d
    JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
    JOIN bio_entities be ON b.biodb_id = be.biodb_id
    WHERE b.inhibitor = true
      AND d.state = 'approved'
    GROUP BY d.drugbank_id, d.name, d.state
)
SELECT
    tc.drugbank_id,
    tc.name,
    tc.state,
    tc.target_count,
    CASE
        WHEN tc.target_count = 1 THEN 'Highly Selective'
        WHEN tc.target_count BETWEEN 2 AND 3 THEN 'Selective'
        WHEN tc.target_count BETWEEN 4 AND 10 THEN 'Multi-target'
        ELSE 'Promiscuous'
    END AS selectivity_class,
    STRING_AGG(DISTINCT be.name, '; ' ORDER BY be.name) AS targets
FROM target_counts tc
JOIN drugs d ON tc.drugbank_id = d.drugbank_id
JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
JOIN bio_entities be ON b.biodb_id = be.biodb_id
WHERE b.inhibitor = true
GROUP BY tc.drugbank_id, tc.name, tc.state, tc.target_count
ORDER BY tc.target_count, tc.name
LIMIT 50;
