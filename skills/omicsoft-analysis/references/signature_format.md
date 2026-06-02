# Gene Signature Format Guide

This document explains how to provide gene signatures for omicsoft DEG analysis.

## Important: Gene Symbol Conventions

Gene symbols MUST match the organism in your h5ad file:

- **Human genes**: ALL UPPERCASE
  - Examples: `CGAS`, `TMEM173`, `TBK1`, `IRF3`, `IFNB1`

- **Mouse genes**: First letter uppercase, rest lowercase
  - Examples: `Cgas`, `Tmem173`, `Tbk1`, `Irf3`, `Ifnb1`

Using incorrect capitalization will result in genes not being found in the analysis.

## Signature Input Format

Users provide signatures as a comma-separated list of "SignatureName:Gene1,Gene2,Gene3" format.

### Single Signature Example

```
CGAS_STING:CGAS,MB21D1,TMEM173,STING1,TBK1,IKBKE,IRF3,IRF7
```

This will be converted to:
```python
{
    'CGAS_STING': ['CGAS', 'MB21D1', 'TMEM173', 'STING1', 'TBK1', 'IKBKE', 'IRF3', 'IRF7']
}
```

### Multiple Signatures Example

Separate multiple signatures with semicolons:

```
CGAS_STING:CGAS,TMEM173,TBK1,IRF3;TGFB:TGFB1,TGFB2,TGFB3,TGFBR1,TGFBR2;JAK_STAT:JAK1,JAK2,JAK3,STAT1,STAT3
```

This will be converted to:
```python
{
    'CGAS_STING': ['CGAS', 'TMEM173', 'TBK1', 'IRF3'],
    'TGFB': ['TGFB1', 'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2'],
    'JAK_STAT': ['JAK1', 'JAK2', 'JAK3', 'STAT1', 'STAT3']
}
```

## Naming Conventions

- Signature names should be descriptive and use underscores for spaces
- Avoid special characters in signature names (except underscores and forward slashes)
- Examples of good names:
  - `CGAS_STING`
  - `Type_I_IFN`
  - `Bcell_SLE/SSc`
  - `MASH_Fibrosis`

## Finding Gene Symbols

1. **For human genes**: Use HGNC (HUGO Gene Nomenclature Committee) official gene symbols
   - Website: https://www.genenames.org/
   - All uppercase format

2. **For mouse genes**: Use MGI (Mouse Genome Informatics) official gene symbols
   - Website: http://www.informatics.jax.org/
   - First letter uppercase, rest lowercase

3. **Cross-species conversion**:
   - Many human-mouse orthologs have the same name but different capitalization
   - Example: Human `CGAS` → Mouse `Cgas`
   - Always verify ortholog relationships for critical genes

## Verifying Gene Symbols

Before running analysis, verify your gene symbols are in the dataset by checking the h5ad file's gene list (adata.var_names). The analysis script will report which genes from your signatures are found in the dataset.

## Common Mistakes to Avoid

1. ❌ Using human capitalization for mouse data: `CGAS` (should be `Cgas`)
2. ❌ Using mouse capitalization for human data: `Cgas` (should be `CGAS`)
3. ❌ Mixing spaces and commas: `CGAS, TMEM173 ,TBK1` (use consistent comma separation)
4. ❌ Using gene aliases instead of official symbols: `STING` (should be `TMEM173` or use both)
5. ❌ Including whitespace around gene names: ` CGAS ` (should be `CGAS`)

## Alternative Gene Names

Some genes have multiple aliases. Include all common aliases in your signature if needed:

```
CGAS_STING:CGAS,MB21D1,TMEM173,STING1,STING,cGAS
```

Both `TMEM173` and `STING1` refer to STING protein - including both ensures the gene is found regardless of annotation.
