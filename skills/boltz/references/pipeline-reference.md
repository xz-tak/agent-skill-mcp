# Protein Structure Prediction Pipeline — Reference

## Prerequisites

- **Conda environment**: `boltz` (Python 3.12, Boltz 2.2.1)
- **GPU**: NVIDIA GPUs (A10G or similar). Pipeline auto-detects idle GPUs.
- **Network**: Access to UniProt, RCSB PDB, and Boltz MSA servers

Activate the environment:
```bash
conda activate boltz
```

---

## Quick Start

### Single gene (default: human)
```bash
cd protein_structure
python predict.py BRCA1
```

### Raw amino acid sequence
```bash
python predict.py --raw-sequence MVTPEGNVSLVDESLL...
```

### Batch mode
```bash
# genes.txt — one gene per line
python predict.py --input-file genes.txt
```

### Download PDB only (skip Boltz if structure exists)
```bash
python predict.py BRCA1 --skip-boltz-if-pdb
```

### High-accuracy prediction
```bash
python predict.py EGFR --recycling-steps 5 --diffusion-samples 3 --use-potentials
```

### Non-human organism
```bash
python predict.py Brca1 --organism "Mus musculus" --organism-id 10090
```

---

## CLI Reference

### Input (mutually exclusive, one required)

| Argument | Description |
|----------|-------------|
| `INPUT` | Positional: gene name or amino acid sequence |
| `--input-file FILE` | Batch file with one gene/sequence per line |
| `--raw-sequence SEQ` | Explicit raw amino acid sequence |

### Pipeline Control

| Flag | Default | Description |
|------|---------|-------------|
| `--skip-boltz-if-pdb` | false | Skip Boltz when PDB structure found |
| `--skip-pdb` | false | Skip PDB lookup entirely |

### Organism

| Flag | Default | Description |
|------|---------|-------------|
| `--organism NAME` | from config | Organism name (config default: Homo sapiens) |
| `--organism-id ID` | from config | NCBI taxonomy ID (config default: 9606) |

### Boltz Parameters

All Boltz parameters are **required** by the module functions and must come from config.yaml or CLI flags. There are no hidden function-level defaults.

| Flag | Config Key | Config Default | Description |
|------|-----------|----------------|-------------|
| `--recycling-steps N` | `boltz.recycling_steps` | 3 | Prediction refinement iterations |
| `--diffusion-samples N` | `boltz.diffusion_samples` | 1 | Structure samples to generate |
| `--sampling-steps N` | `boltz.sampling_steps` | 200 | Diffusion sampling steps |
| `--output-format FMT` | `boltz.output_format` | mmcif | Output: `pdb` or `mmcif` |
| `--use-potentials` | `boltz.use_potentials` | false | Physical plausibility potentials |
| `--no-msa-server` | `boltz.use_msa_server` | true | Disable MSA server |
| `--model MODEL` | `boltz.model` | boltz2 | Model: `boltz1` or `boltz2` |
| `--timeout SECS` | `boltz.timeout_seconds` | 3600 | Per-job timeout in seconds |
| `--rerun-threshold N` | `boltz.rerun_confidence_threshold` | 0.5 | Auto re-run predictions below this confidence with boosted params (recycling=10, samples=5, potentials) |

### Parallelism & Config

| Flag | Config Key | Default | Description |
|------|-----------|---------|-------------|
| `--max-parallel N` | `gpu.max_parallel` | auto | Max parallel jobs (auto = #available GPUs) |
| `--config PATH` | — | ./config.yaml | Configuration file path |
| `--output-dir PATH` | `output.base_dir` | . | Base output directory |

---

## Configuration File (`config.yaml`)

The config file is the **single source of truth** for default parameter values. CLI flags override config. Config overrides built-in Python fallbacks.

**Priority**: CLI flags > config.yaml > DEFAULT_CONFIG (Python fallback)

```yaml
organism: "Homo sapiens"
organism_id: 9606

boltz:
  use_msa_server: true
  recycling_steps: 3
  diffusion_samples: 1
  sampling_steps: 200
  output_format: mmcif
  use_potentials: false
  model: boltz2
  timeout_seconds: 3600
  rerun_confidence_threshold: 0.5   # Auto re-run below this confidence

output:
  pdb_dir_suffix: "_pdb"
  boltz_dir: "."
  summary_file: "results_summary.json"

gpu:
  utilization_threshold: 10.0   # GPUs below this % util are considered available
  max_parallel: null            # null = auto-detect from available GPUs
```

---

## Output Structure

```
protein_structure/
├── {GENE}_pdb/                          # PDB downloads
│   └── {PDB_ID}.pdb                    # Best-resolution experimental structure
├── {GENE}.yaml                          # Generated Boltz input YAML
├── boltz_results_{GENE}/                # Boltz predictions (flat, same level as _pdb dirs)
│   └── predictions/{GENE}/
│       ├── {GENE}_model_0.cif          # Predicted structure (mmcif)
│       ├── {GENE}_model_0.pdb          # Predicted structure (if converted)
│       └── confidence_{GENE}_model_0.json
└── results_summary.json                 # Run log with all metadata
```

---

## results_summary.json Fields

Each entry contains:

| Field | Type | Description |
|-------|------|-------------|
| `job_name` | string | Gene name or sequence identifier |
| `input_type` | string | `"gene"` or `"sequence"` |
| `sequence_length` | int | Amino acid count (from UniProt canonical) |
| `pdb_id` | string? | Best PDB ID found (null if none) |
| `pdb_resolution` | float? | Resolution in angstroms |
| `pdb_path` | string? | Downloaded PDB file path |
| `pdb_title` | string? | PDB entry title |
| `pdb_method` | string? | Experimental method (e.g., X-RAY DIFFRACTION) |
| `pdb_deposition_date` | string? | PDB deposition date (ISO 8601) |
| `pdb_entity_sequence_length` | int? | Chain length in the crystal structure |
| `pdb_entity_description` | string? | Protein description from PDB entity |
| `pdb_entity_gene_name` | string? | Gene name as recorded in PDB entity |
| `pdb_entity_organism` | string? | Source organism from PDB entity |
| `pdb_entity_weight_kda` | float? | Molecular weight in kDa |
| `uniprot_accession` | string? | UniProt accession (Swiss-Prot preferred) |
| `boltz_ran` | bool | Whether Boltz was executed |
| `boltz_rerun` | bool | Whether result came from an auto re-run (low confidence) |
| `boltz_skipped` | bool | Whether Boltz was skipped due to PDB |
| `boltz_status` | string? | `"completed"`, `"failed"`, `"timeout"` |
| `boltz_output_dir` | string? | Boltz prediction directory |
| `boltz_confidence` | float? | Overall confidence score (0-1) |
| `errors` | list | Any errors encountered |
| `timestamp` | string | ISO 8601 timestamp |

---

## Module API Reference

### `modules.fetch_sequence`

**`fetch_uniprot_sequence(gene_name, organism_id=9606) -> dict`**

Queries UniProt REST API for canonical protein sequence. Tries Swiss-Prot (`reviewed:true`) first, falls back to TrEMBL if no reviewed entry exists.

- Returns: `{accession, gene_name, sequence, organism, length}`
- Raises: `ValueError` if gene not found in either Swiss-Prot or TrEMBL
- Retries: 3 attempts with exponential backoff
- Known issue (fixed): Previously picked TrEMBL fragments (e.g., 34 aa for ITGA4) instead of canonical Swiss-Prot entries (1032 aa). Fixed by adding `reviewed:true` filter with fallback.

**`validate_sequence(sequence) -> bool`**

Checks if string contains only standard amino acid characters (A-Y + X).

### `modules.check_pdb`

**`search_pdb_by_gene(gene_name, organism_id=9606) -> list[dict]`**

Searches RCSB PDB for experimental structures, sorted by resolution. For the best hit, fetches full entry and entity metadata.

- Returns: `[{pdb_id, score, resolution, title, method, deposition_date, entity_sequence_length, entity_description, entity_gene_name, entity_organism, entity_weight_kda}]`
- Entity matching: iterates polymer entities 1-9 to find the one matching the query gene name (case-insensitive substring match on `pdbx_gene_src_gene`)
- Note: `entity_sequence_length` is the chain length in the crystal, which may be a domain fragment shorter than the full-length UniProt sequence

**`download_pdb(pdb_id, output_dir) -> Path`**

Downloads PDB file from RCSB.

**`check_and_download_best_pdb(gene_name, base_dir, organism_id=9606, pdb_dir_suffix="_pdb") -> dict | None`**

Orchestrator: search + download best resolution structure with full metadata.

- Returns: `{pdb_id, resolution, file_path, title, method, deposition_date, entity_sequence_length, entity_description, entity_gene_name, entity_organism, entity_weight_kda}` or None

### `modules.run_boltz`

**`run_boltz_prediction(sequence, job_name, output_dir, gpu_id, use_msa_server, recycling_steps, diffusion_samples, sampling_steps, output_format, use_potentials, model, timeout_seconds) -> dict`**

Runs Boltz prediction via subprocess with GPU isolation. All parameters are required — values come from config.yaml via predict.py.

- Returns: `{success, status, job_name, output_path, confidence, error_message, command, stderr}`
- CUDA setup: auto-detects `cu13` lib path from CONDA_PREFIX for LD_LIBRARY_PATH
- Handles timeouts: catches `TimeoutExpired`, returns status="timeout"
- Handles silent failures: checks for empty prediction dir even on exit code 0

**`parse_boltz_confidence(prediction_dir, job_name) -> dict | None`**

Reads Boltz confidence JSON from `confidence_{job_name}_model_0.json`.

### `modules.utils`

**`detect_available_gpus(threshold=10.0) -> list[int]`** — Returns GPU indices with utilization below threshold via `nvidia-smi`.

**`generate_boltz_yaml(sequence, chain_id, output_path) -> Path`** — Writes Boltz-compatible YAML input (`version: 1`).

**`load_config(config_path) -> dict`** — Loads config.yaml, deep-merges with DEFAULT_CONFIG fallback.

**`merge_config_with_args(config, cli_args) -> dict`** — Merges CLI overrides into config. Maps CLI arg names to nested config paths.

**`retry_request(func, max_retries=3, base_delay=1.0)`** — Exponential backoff wrapper (delays: 1s, 2s, 4s).

**`write_results_summary(summary_path, results)`** — Appends timestamped results to JSON summary file.

---

## Troubleshooting

### Gene not found in UniProt
```
ValueError: Gene 'XYZ' not found in UniProt for organism_id=9606
```
Check spelling. Use official HGNC gene symbols. Use `--organism-id` for non-human.

### UniProt returns fragment/short sequence
The pipeline queries Swiss-Prot (`reviewed:true`) first. If it still returns a short sequence, confirm the gene symbol is the primary name (not an alias). Compare `sequence_length` in results_summary.json against expected full-length from UniProt web.

### PDB entity sequence length differs from UniProt length
Normal. PDB structures are often domain fragments (e.g., TYK2 kinase domain = 318 aa vs full-length 1187 aa). The `pdb_entity_sequence_length` field reports the chain in the crystal, not the full protein. Compare against `sequence_length` (from UniProt) to identify partial coverage.

### No PDB structures found
Normal for many proteins. Boltz prediction will still run.

### Boltz timeout
Default timeout is 3600s (1 hour) per job, set in config.yaml. Increase with `--timeout 7200` or edit `boltz.timeout_seconds` in config.yaml. Long sequences (>1000 aa) may need more time.

### GPU not detected
Run `nvidia-smi` to verify. The pipeline looks for GPUs below the `gpu.utilization_threshold` (default 10%). If all GPUs are busy, it falls back to GPU 0.

### CUDA library errors / empty predictions
Boltz requires the `cu13` CUDA library. The pipeline auto-detects and sets `LD_LIBRARY_PATH` from `$CONDA_PREFIX/lib/python3.12/site-packages/nvidia/cu13/lib`. If predictions are empty despite exit code 0, this path may be wrong. Verify:
```bash
ls $CONDA_PREFIX/lib/python3.12/site-packages/nvidia/cu13/lib/
```

### MSA server issues
If the Boltz MSA server is unreachable, use `--no-msa-server` (slower, uses local MSA generation). Set `boltz.use_msa_server: false` in config.yaml for persistent change.

### Network errors
All API calls (UniProt, RCSB PDB) retry 3 times with exponential backoff (1s, 2s, 4s delays).

---

## Known Issues and Past Fixes

| Issue | Root Cause | Fix |
|-------|-----------|-----|
| UniProt returning 34 aa for ITGA4 (should be 1032) | Query lacked `reviewed:true`, picked TrEMBL fragment | Added Swiss-Prot filter with TrEMBL fallback in `fetch_sequence.py` |
| PDB metadata missing (only pdb_id + resolution) | `check_pdb.py` only called entry API, not entity API | Added `_get_entry_metadata()` and `_get_entity_metadata()` |
| Hardcoded defaults in function signatures masking config | Defaults in 3 places: functions, DEFAULT_CONFIG, config.yaml | Removed defaults from `run_boltz_prediction()` signature; kept only in config.yaml and DEFAULT_CONFIG |

---

## Boltz Background

**Boltz-2** is an open-source biomolecular foundation model for structure prediction. It predicts protein structures approaching AlphaFold2 accuracy, with additional capabilities for binding affinity estimation. Released under MIT license.

- Paper: Wohlwend et al. (2024) — Boltz-1 and Boltz-2
- GitHub: github.com/jwohlwend/boltz
- Model: Trained on PDB structures with diffusion-based architecture
