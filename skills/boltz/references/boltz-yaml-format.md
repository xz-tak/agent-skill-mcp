# Boltz YAML Input Format Reference

## Basic Structure

Every Boltz input YAML must have `version: 1` and a `sequences` list. Optional sections: `constraints`, `templates`, `properties`, `modifications`.

## Entity Types

### Protein
```yaml
- protein:
    id: A
    sequence: MVTPEGNVSLVDES...
    msa: ./path/to/msa.a3m    # optional, omit with --use_msa_server
    cyclic: true               # optional, for cyclic proteins
```

### DNA
```yaml
- dna:
    id: B
    sequence: ATCGATCG...
```

### RNA
```yaml
- rna:
    id: C
    sequence: AUCGAUCG...
```

### Ligand (CCD code)
```yaml
- ligand:
    id: D
    ccd: SAH
```

### Ligand (SMILES)
```yaml
- ligand:
    id: E
    smiles: 'N[C@@H](Cc1ccc(O)cc1)C(=O)O'
```

## Constraints (optional)

### Bond constraint
```yaml
constraints:
  - bond:
      atom1: [A, 1, CA]     # [chain, residue, atom]
      atom2: [B, 1, CA]
```

### Pocket constraint
```yaml
constraints:
  - pocket:
      binder: D              # ligand chain ID
      contacts: [100, 200]   # residue indices
```

### Contact constraint
```yaml
constraints:
  - contact:
      atom1: [A, 50, CA]
      atom2: [B, 100, CA]
      max_distance: 8.0
```

## Templates (optional)

```yaml
templates:
  - cif: ./template.cif
```

## Properties (optional)

### Affinity prediction (ligand must be ≤128 atoms)
```yaml
properties:
  affinity:
    binder: D    # ligand chain ID
```

Returns:
- `affinity_probability_binary` (0-1): binder vs non-binder probability
- `affinity_pred_value` (log10 scale): binding strength as log10(IC50) in μM

## Modifications (optional)

```yaml
sequences:
  - protein:
      id: A
      sequence: MVTPEG...
      modifications:
        - position: 5
          ccd: SEP    # phosphoserine
```

## Complete Example: Protein-Ligand Complex with Affinity

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: MVTPEGNVSLVDESLLVGVTDEDRAVRSAHQFYERLIGLWAPAVMEAAHELGR
  - ligand:
      id: B
      smiles: 'c1ccc(cc1)C(=O)O'
constraints:
  - pocket:
      binder: B
      contacts: [10, 20, 30]
properties:
  affinity:
    binder: B
```

## Complete Example: Protein-Protein Multimer

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: MVTPEGNVSLVDES...
  - protein:
      id: B
      sequence: QLEDSEVEAVAKGL...
```

## Complete Example: Protein-DNA Complex

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: MVTPEGNVSLVDES...
  - dna:
      id: B
      sequence: ATCGATCGATCG
  - dna:
      id: C
      sequence: CGATCGATCGAT
```
