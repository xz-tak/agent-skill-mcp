-- =====================================================
-- Drug-Target Interaction Examples
-- =====================================================
-- These examples show how to query drug-target relationships,
-- including mechanisms of action, target proteins, and genes.
-- =====================================================

-- Example 1: Get all targets for a specific drug
SELECT DISTINCT
    d.drugbank_id,
    d.name AS drug_name,
    be.name AS target_name,
    be.type AS target_type,
    p.gene_name,
    b.inhibitor,
    b.agonist,
    b.antagonist,
    b.pharmacological_action
FROM drugs d
JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
JOIN bio_entities be ON b.biodb_id = be.biodb_id
LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
WHERE d.name ILIKE '%imatinib%';

-- Example 2: Find drugs that target a specific protein
SELECT DISTINCT
    d.drugbank_id,
    d.name AS drug_name,
    d.type AS drug_type,
    d.state AS drug_state,
    be.name AS target_name,
    p.gene_name,
    b.inhibitor,
    b.pharmacological_action
FROM drugs d
JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
JOIN bio_entities be ON b.biodb_id = be.biodb_id
LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
WHERE be.name ILIKE '%epidermal growth factor receptor%'
   OR be.name ILIKE '%EGFR%'
ORDER BY d.state, d.name;

-- Example 3: Find drugs targeting a specific gene
SELECT DISTINCT
    d.drugbank_id,
    d.name AS drug_name,
    d.state,
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
WHERE p.gene_name ILIKE '%KRAS%'
ORDER BY d.state DESC, d.name;

-- Example 4: Find approved inhibitors for a target
SELECT DISTINCT
    d.drugbank_id,
    d.name AS drug_name,
    be.name AS target_name,
    p.gene_name,
    b.pharmacological_action
FROM drugs d
JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
JOIN bio_entities be ON b.biodb_id = be.biodb_id
LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
WHERE d.state = 'approved'
  AND b.inhibitor = true
  AND (be.name ILIKE '%vascular endothelial growth factor%' OR be.name ILIKE '%VEGF%')
ORDER BY d.name;

-- Example 5: Get target with detailed protein information
SELECT DISTINCT
    d.drugbank_id,
    d.name AS drug_name,
    be.name AS target_name,
    p.gene_name,
    p.general_function,
    p.specific_function,
    bec.uniprot_id,
    b.inhibitor,
    b.agonist,
    b.antagonist
FROM drugs d
JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
JOIN bio_entities be ON b.biodb_id = be.biodb_id
LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
WHERE d.name ILIKE '%cetuximab%';

-- Example 6: Find all agonists for a target class
SELECT DISTINCT
    d.drugbank_id,
    d.name AS drug_name,
    d.state,
    be.name AS target_name,
    p.gene_name
FROM drugs d
JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
JOIN bio_entities be ON b.biodb_id = be.biodb_id
LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
WHERE b.agonist = true
  AND be.name ILIKE '%dopamine receptor%'
ORDER BY d.state DESC, d.name;

-- Example 7: Find drugs with known pharmacological action
SELECT DISTINCT
    d.drugbank_id,
    d.name AS drug_name,
    d.state,
    be.name AS target_name,
    p.gene_name,
    b.pharmacological_action,
    b.inhibitor,
    b.agonist,
    b.antagonist
FROM drugs d
JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
JOIN bio_entities be ON b.biodb_id = be.biodb_id
LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
WHERE b.pharmacological_action IS NOT NULL
  AND d.name ILIKE '%erlotinib%';

-- Example 8: Count targets per drug
SELECT
    d.drugbank_id,
    d.name AS drug_name,
    d.state,
    COUNT(DISTINCT be.biodb_id) AS target_count
FROM drugs d
JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
JOIN bio_entities be ON b.biodb_id = be.biodb_id
WHERE d.state = 'approved'
GROUP BY d.drugbank_id, d.name, d.state
HAVING COUNT(DISTINCT be.biodb_id) >= 5
ORDER BY target_count DESC
LIMIT 20;

-- Example 9: Find drugs with enzyme interactions (metabolism)
SELECT DISTINCT
    d.drugbank_id,
    d.name AS drug_name,
    be.name AS enzyme_name,
    p.gene_name,
    b.pharmacological_action
FROM drugs d
JOIN bonds b ON d.id = b.drug_id AND b.type = 'EnzymeBond'
JOIN bio_entities be ON b.biodb_id = be.biodb_id
LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
WHERE d.name ILIKE '%warfarin%';

-- Example 10: Find all kinase inhibitors
SELECT DISTINCT
    d.drugbank_id,
    d.name AS drug_name,
    d.state,
    be.name AS target_name,
    p.gene_name
FROM drugs d
JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
JOIN bio_entities be ON b.biodb_id = be.biodb_id
LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
WHERE b.inhibitor = true
  AND (be.name ILIKE '%kinase%' OR p.gene_name ILIKE '%kinase%')
  AND d.state = 'approved'
ORDER BY d.name
LIMIT 30;
