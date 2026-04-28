#!/usr/bin/env python3
"""
ReCoN Multinetwork Configuration

Shared configuration dataclass for all pipeline modules.
Consolidates ~40 parameters into a single serializable config.
"""

import argparse
import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional


@dataclass
class ReconConfig:
    """Configuration for the ReCoN multicellular network pipeline."""

    # === PATHS (user provides) ===
    h5ad_path: str = ""
    output_dir: str = "results/"
    data_dir: str = "data/"
    scatac_path: Optional[str] = None
    scatac_metadata_path: Optional[str] = None
    cellchat_data_path: Optional[str] = None

    # === Data columns ===
    condition_col: str = "condition"
    celltype_col: str = "cluster"
    disease_conditions: List[str] = field(default_factory=list)
    normal_condition: str = "normal"

    # === GRN parameters (M2) ===
    ref_genome: str = "hg38"
    circe_window: int = 500_000
    min_cells_grn: int = 20
    n_cpus: int = 8
    motif_fpr: float = 0.01
    tss_distance: int = 10_000
    min_features_scatac: int = 300
    min_cells_scatac: int = 50
    nb_features_selected: int = 100_000
    scatac_celltype_mapping: Dict[str, List[str]] = field(default_factory=dict)

    # === CCC compute parameters (M3) ===
    ccc_compute_methods: List[str] = field(
        default_factory=lambda: ["cellphonedb", "cellchat"]
    )
    resource_name: str = "consensus"
    expr_prop: float = 0.1
    min_lr_means: float = 0.5

    # === CCC source for downstream (M4-M7) ===
    ccc_source: str = "merged"

    # === ReCoN parameters (M4) ===
    restart_proba: float = 0.6
    alpha: float = 0.8
    min_grn_weight: float = 0.5
    n_jobs: int = 16
    extend_seeds: bool = True

    # === Seeds for ReCoN (M4) ===
    seeds: List[str] = field(default_factory=list)
    seeds_file: Optional[str] = None

    # === Seeds for Sankey visualization (M7) ===
    seed_categories: Dict[str, List[str]] = field(default_factory=dict)

    # === Focal cell types for Sankey (M7) ===
    focal_celltypes: List[str] = field(default_factory=list)
    ligand_source_cells: List[str] = field(default_factory=list)

    # === Target Prediction Parameters (M8) ===
    target_genes: List[str] = field(default_factory=list)
    focus_cell_types: List[str] = field(default_factory=list)
    seed_type: str = "receptor_activation"
    prediction_output_dir: Optional[str] = None
    top_tfs_sankey: int = 20
    top_grn_genes_sankey: int = 20
    min_rtf_weight: float = 0.01
    ppi_min_score: int = 400
    min_sankey_grn_weight: float = 0.5

    # === GSEA enrichment (integrated in M8) ===
    gsea_gene_sets: List[str] = field(default_factory=lambda: ["MSigDB_Hallmark_2020"])
    gsea_min_size: int = 3
    gsea_max_size: int = 5000
    gsea_permutations: int = 1000
    gsea_fdr_threshold: float = 0.05

    # === Cascade parameters (M6) ===
    n_permutations: int = 1000
    edge_weight_threshold: float = 0.5
    min_cascade_grn_weight: float = 0.01
    max_cascades_per_cellpair: int = 50_000
    fdr_threshold: float = 0.05

    # === Multinetwork parameters (M5) ===
    grn_score_threshold: float = 0.001
    module_file: Optional[str] = None

    def to_json(self, path: Optional[Path] = None) -> str:
        """Serialize config to JSON string. Optionally write to file."""
        data = asdict(self)
        json_str = json.dumps(data, indent=2)
        if path is not None:
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            Path(path).write_text(json_str)
        return json_str

    @classmethod
    def from_json(cls, path_or_str: str) -> "ReconConfig":
        """Load config from JSON file path or JSON string."""
        p = Path(path_or_str)
        if p.exists():
            data = json.loads(p.read_text())
        else:
            data = json.loads(path_or_str)
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})

    @classmethod
    def from_cli_args(cls, args: Optional[argparse.Namespace] = None) -> "ReconConfig":
        """Build config from CLI arguments, optionally layered on a JSON config."""
        parser = _build_parser()
        if args is None:
            args = parser.parse_args()

        # Start from JSON config if provided
        if hasattr(args, "config") and args.config:
            config = cls.from_json(args.config)
        else:
            config = cls()

        # Override with any explicitly provided CLI args
        cli_dict = vars(args)
        for key, value in cli_dict.items():
            if key == "config":
                continue
            if value is not None and hasattr(config, key):
                setattr(config, key, value)

        return config

    def validate(self) -> List[str]:
        """Validate config, return list of error messages (empty = valid)."""
        errors = []
        if self.h5ad_path and not Path(self.h5ad_path).exists():
            errors.append(f"h5ad_path not found: {self.h5ad_path}")
        if self.scatac_path and not Path(self.scatac_path).exists():
            errors.append(f"scatac_path not found: {self.scatac_path}")
        if not self.disease_conditions:
            errors.append("disease_conditions is empty")
        if not self.normal_condition:
            errors.append("normal_condition is empty")
        if self.ccc_source not in ("merged", "cellphonedb", "cellchat"):
            errors.append(f"Invalid ccc_source: {self.ccc_source}")
        if self.restart_proba <= 0 or self.restart_proba >= 1:
            errors.append(f"restart_proba must be in (0,1): {self.restart_proba}")
        if self.alpha < 0 or self.alpha > 1:
            errors.append(f"alpha must be in [0,1]: {self.alpha}")
        if self.ppi_min_score < 0 or self.ppi_min_score > 1000:
            errors.append(f"ppi_min_score must be in [0,1000]: {self.ppi_min_score}")
        if self.top_tfs_sankey < 1:
            errors.append(f"top_tfs_sankey must be >= 1: {self.top_tfs_sankey}")
        if self.min_rtf_weight < 0:
            errors.append(f"min_rtf_weight must be >= 0: {self.min_rtf_weight}")
        return errors

    @property
    def conditions(self) -> List[str]:
        """All conditions (disease + normal)."""
        return list(self.disease_conditions) + [self.normal_condition]

    def get_output_path(self, *parts: str) -> Path:
        """Get a path under output_dir, creating parent dirs."""
        p = Path(self.output_dir).joinpath(*parts)
        p.parent.mkdir(parents=True, exist_ok=True)
        return p

    def get_data_path(self, *parts: str) -> Path:
        """Get a path under data_dir."""
        return Path(self.data_dir).joinpath(*parts)

    def get_ccc_dir(self) -> Path:
        """Get CCC results directory."""
        return self.get_output_path("ccc")

    def get_grn_dir(self) -> Path:
        """Get GRN results directory."""
        return self.get_output_path("grn")

    def get_recon_dir(self) -> Path:
        """Get ReCoN results directory (ccc_source-specific)."""
        return self.get_output_path("recon", f"{self.ccc_source}_ccc")

    def get_multinetwork_dir(self) -> Path:
        """Get multinetwork output directory (ccc_source-specific)."""
        return self.get_output_path("multinetwork", f"{self.ccc_source}_ccc")

    def get_figures_dir(self) -> Path:
        """Get figures output directory (ccc_source-specific)."""
        return self.get_output_path("figures", f"{self.ccc_source}_ccc")

    def get_cascade_dir(self) -> Path:
        """Get differential cascades directory (ccc_source-specific)."""
        return self.get_output_path("differential_cascades", f"{self.ccc_source}_ccc")

    def get_prediction_dir(self) -> Path:
        """Get target prediction output directory."""
        if self.prediction_output_dir:
            p = Path(self.prediction_output_dir)
            p.mkdir(parents=True, exist_ok=True)
            return p
        return self.get_output_path("target_prediction", f"{self.ccc_source}_ccc")

    def load_seeds(self) -> List[str]:
        """Load seed genes from seeds list or seeds_file."""
        if self.seeds:
            return list(self.seeds)
        if self.seeds_file and Path(self.seeds_file).exists():
            with open(self.seeds_file) as f:
                return [line.strip() for line in f if line.strip()]
        return []


def _build_parser() -> argparse.ArgumentParser:
    """Build CLI argument parser with all config parameters."""
    parser = argparse.ArgumentParser(
        description="ReCoN Multicellular Network Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", type=str, help="Path to JSON config file")
    parser.add_argument("--h5ad", dest="h5ad_path", type=str)
    parser.add_argument("--output-dir", type=str)
    parser.add_argument("--data-dir", type=str)
    parser.add_argument("--scatac-path", type=str)
    parser.add_argument("--scatac-metadata-path", type=str)
    parser.add_argument("--cellchat-data-path", type=str)
    parser.add_argument("--condition-col", type=str)
    parser.add_argument("--celltype-col", type=str)
    parser.add_argument("--disease-conditions", nargs="+", type=str)
    parser.add_argument("--normal-condition", type=str)
    parser.add_argument("--ccc-source", choices=["merged", "cellphonedb", "cellchat"])
    parser.add_argument("--ccc-compute-methods", nargs="+", type=str)
    parser.add_argument("--seeds", nargs="+", type=str)
    parser.add_argument("--seeds-file", type=str)
    parser.add_argument("--restart-proba", type=float)
    parser.add_argument("--alpha", type=float)
    parser.add_argument("--min-grn-weight", type=float)
    parser.add_argument("--n-jobs", type=int)
    parser.add_argument("--n-cpus", type=int)
    parser.add_argument("--n-permutations", type=int)
    parser.add_argument("--fdr-threshold", type=float)
    parser.add_argument("--max-cascades-per-cellpair", type=int)
    parser.add_argument("--grn-score-threshold", type=float)
    parser.add_argument("--module-file", type=str)
    parser.add_argument("--ref-genome", type=str)
    parser.add_argument("--resource-name", type=str)
    parser.add_argument("--expr-prop", type=float)
    # M8 target prediction
    parser.add_argument("--target-genes", nargs="+", type=str)
    parser.add_argument("--focus-cell-types", nargs="+", type=str)
    parser.add_argument("--seed-type", type=str)
    parser.add_argument("--prediction-output-dir", type=str)
    parser.add_argument("--top-tfs-sankey", type=int)
    parser.add_argument("--top-grn-genes-sankey", type=int)
    parser.add_argument("--min-rtf-weight", type=float)
    parser.add_argument("--ppi-min-score", type=int)
    parser.add_argument("--min-sankey-grn-weight", type=float)
    # GSEA enrichment
    parser.add_argument("--gsea-gene-sets", nargs="+", type=str)
    parser.add_argument("--gsea-min-size", type=int)
    parser.add_argument("--gsea-max-size", type=int)
    parser.add_argument("--gsea-permutations", type=int)
    parser.add_argument("--gsea-fdr-threshold", type=float)
    return parser


def get_config(config: Optional["ReconConfig"] = None) -> "ReconConfig":
    """Get config: use provided config, or build from CLI args."""
    if config is not None:
        return config
    return ReconConfig.from_cli_args()


if __name__ == "__main__":
    config = ReconConfig.from_cli_args()
    errors = config.validate()
    if errors:
        print("Validation errors:")
        for e in errors:
            print(f"  - {e}")
    else:
        print("Config valid.")
    print(config.to_json())
