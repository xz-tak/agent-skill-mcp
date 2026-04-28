import copy
import json
import logging
import re
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Optional

import yaml

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWXY")
GENE_NAME_PATTERN = re.compile(r"^[A-Za-z0-9_.\-]{1,50}$")
MIN_SEQUENCE_LENGTH = 30

DEFAULT_CONFIG = {
    "organism": "Homo sapiens",
    "organism_id": 9606,
    "boltz": {
        "use_msa_server": True,
        "recycling_steps": 3,
        "diffusion_samples": 1,
        "sampling_steps": 200,
        "output_format": "mmcif",
        "use_potentials": False,
        "model": "boltz2",
        "timeout_seconds": 3600,
        "rerun_confidence_threshold": 0.5,
    },
    "output": {
        "pdb_dir_suffix": "_pdb",
        "boltz_dir": ".",
        "summary_file": "results_summary.json",
    },
    "gpu": {
        "utilization_threshold": 10.0,
        "max_parallel": None,
    },
}


def setup_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s [%(name)s] %(levelname)s: %(message)s")
        )
        logger.addHandler(handler)
    logger.setLevel(level)
    return logger


def load_config(config_path: Optional[Path] = None) -> dict:
    if config_path and config_path.exists():
        with open(config_path) as f:
            loaded = yaml.safe_load(f) or {}
        return _deep_merge(DEFAULT_CONFIG, loaded)
    return copy.deepcopy(DEFAULT_CONFIG)


def _deep_merge(base: dict, override: dict) -> dict:
    result = {}
    for key in set(base) | set(override):
        if key in override and key in base:
            if isinstance(base[key], dict) and isinstance(override[key], dict):
                result[key] = _deep_merge(base[key], override[key])
            else:
                result[key] = override[key]
        elif key in override:
            result[key] = override[key]
        else:
            result[key] = base[key]
    return result


def merge_config_with_args(config: dict, cli_args: dict) -> dict:
    merged = copy.deepcopy(config)
    arg_to_config = {
        "organism": ("organism",),
        "organism_id": ("organism_id",),
        "recycling_steps": ("boltz", "recycling_steps"),
        "diffusion_samples": ("boltz", "diffusion_samples"),
        "sampling_steps": ("boltz", "sampling_steps"),
        "output_format": ("boltz", "output_format"),
        "use_potentials": ("boltz", "use_potentials"),
        "model": ("boltz", "model"),
        "timeout": ("boltz", "timeout_seconds"),
        "rerun_threshold": ("boltz", "rerun_confidence_threshold"),
        "max_parallel": ("gpu", "max_parallel"),
        "output_dir": ("output", "base_dir"),
    }
    for arg_name, config_path in arg_to_config.items():
        value = cli_args.get(arg_name)
        if value is None:
            continue
        target = merged
        for key in config_path[:-1]:
            target = target.setdefault(key, {})
        target[config_path[-1]] = value

    if cli_args.get("no_msa_server"):
        merged.setdefault("boltz", {})["use_msa_server"] = False

    return merged


def generate_boltz_yaml(
    sequence: str,
    chain_id: str,
    output_path: Path,
) -> Path:
    data = {
        "version": 1,
        "sequences": [
            {"protein": {"id": chain_id, "sequence": sequence}}
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)
    return output_path


def detect_available_gpus(utilization_threshold: float = 10.0) -> list[int]:
    try:
        result = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=index,utilization.gpu",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode != 0:
            return []
        available = []
        for line in result.stdout.strip().split("\n"):
            if not line.strip():
                continue
            parts = line.split(",")
            if len(parts) < 2:
                continue
            try:
                idx = int(parts[0].strip())
                util = float(parts[1].strip())
            except (ValueError, IndexError):
                continue
            if util < utilization_threshold:
                available.append(idx)
        return available
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return []


def retry_request(
    func: Callable,
    max_retries: int = 3,
    base_delay: float = 1.0,
) -> Any:
    last_error = None
    for attempt in range(max_retries):
        try:
            return func()
        except Exception as e:
            last_error = e
            if attempt < max_retries - 1:
                delay = base_delay * (2 ** attempt)
                time.sleep(delay)
    raise last_error


def write_results_summary(
    summary_path: Path,
    results: list[dict],
) -> None:
    existing = []
    if summary_path.exists():
        with open(summary_path) as f:
            existing = json.load(f)

    timestamped = [
        {**r, "timestamp": r.get("timestamp", datetime.now(timezone.utc).isoformat())}
        for r in results
    ]

    combined = existing + timestamped
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_path, "w") as f:
        json.dump(combined, f, indent=2, default=str)


def is_amino_acid_sequence(s: str) -> bool:
    if len(s) <= MIN_SEQUENCE_LENGTH:
        return False
    return all(c in STANDARD_AA for c in s.upper())


def sanitize_gene_name(name: str) -> str:
    if not GENE_NAME_PATTERN.match(name):
        raise ValueError(
            f"Invalid gene name '{name}': must contain only alphanumeric "
            "characters, underscores, hyphens, and dots (max 50 chars)"
        )
    return name
