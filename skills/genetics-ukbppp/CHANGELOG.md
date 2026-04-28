# genetics-ukbppp Changelog

## 1.2.0 — 2026-04-21

- TSV filenames now include ICD-10 codes: `{GENE}_{ICD10_codes}_ukbppp.tsv` (e.g., `IL12B_K50_K51_ukbppp.tsv`)
- No TSV export for genes not found on the Olink panel (missing genes score 0, reported in summary only)
- Added Status column to Individual Gene Rankings table: `Significant`, `Not significant`, `Not on panel`
- Log entries for missing genes now say "no TSV exported" instead of referencing a file

## 1.1.0 — Initial

- Core query script with download, resolve-disease, and query actions
- Scoring system: FDR < 0.05 direction-based scoring across logit, ARD, CoxPH models
- Combo scoring with mean-of-components approach
- Summary markdown with individual gene rankings, detail tables, and combo sections
- Per-gene TSV export with enriched columns (disease name, protein name, Ensembl ID)
- Query log with column descriptions and model notes
