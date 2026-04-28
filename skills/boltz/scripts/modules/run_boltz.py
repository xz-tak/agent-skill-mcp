import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional

from modules.utils import generate_boltz_yaml, setup_logger

logger = setup_logger("run_boltz")


def run_boltz_prediction(
    sequence: str,
    job_name: str,
    output_dir: Path,
    gpu_id: int,
    use_msa_server: bool,
    recycling_steps: int,
    diffusion_samples: int,
    sampling_steps: int,
    output_format: str,
    use_potentials: bool,
    model: str,
    timeout_seconds: int,
) -> dict:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    yaml_path = output_dir / f"{job_name}.yaml"
    generate_boltz_yaml(sequence, "A", yaml_path)

    cmd = [
        "boltz", "predict",
        str(yaml_path),
        "--out_dir", str(output_dir),
        "--devices", "1",
        "--recycling_steps", str(recycling_steps),
        "--diffusion_samples", str(diffusion_samples),
        "--sampling_steps", str(sampling_steps),
        "--output_format", output_format,
        "--model", model,
    ]

    if use_msa_server:
        cmd.append("--use_msa_server")
    if use_potentials:
        cmd.append("--use_potentials")

    env = {**os.environ, "CUDA_VISIBLE_DEVICES": str(gpu_id)}
    pyver = f"python{sys.version_info.major}.{sys.version_info.minor}"
    cu13_lib = Path(os.environ.get("CONDA_PREFIX", "")) / f"lib/{pyver}/site-packages/nvidia/cu13/lib"
    if cu13_lib.exists():
        existing = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = f"{cu13_lib}:{existing}" if existing else str(cu13_lib)

    logger.info(
        "Running Boltz for %s on GPU %d (timeout=%ds)",
        job_name, gpu_id, timeout_seconds,
    )
    logger.info("Command: %s", " ".join(cmd))

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            env=env,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as e:
        logger.error("Boltz timed out after %ds for %s", timeout_seconds, job_name)
        return {
            "success": False,
            "status": "timeout",
            "job_name": job_name,
            "output_path": None,
            "confidence": None,
            "error_message": f"Boltz timed out after {timeout_seconds}s for {job_name}",
            "command": " ".join(cmd),
            "stderr": (e.stderr or "")[:2000] if e.stderr else "",
        }

    if result.returncode != 0:
        logger.error("Boltz failed for %s: %s", job_name, result.stderr[:500])
        return {
            "success": False,
            "status": "failed",
            "job_name": job_name,
            "output_path": None,
            "confidence": None,
            "error_message": result.stderr[:2000],
            "command": " ".join(cmd),
            "stderr": result.stderr[:2000],
        }

    prediction_dir = output_dir / f"boltz_results_{job_name}" / "predictions" / job_name
    has_predictions = prediction_dir.exists() and any(prediction_dir.iterdir())

    if not has_predictions:
        logger.error(
            "Boltz exited 0 but no predictions found for %s (likely CUDA/GPU error)",
            job_name,
        )
        return {
            "success": False,
            "status": "failed",
            "job_name": job_name,
            "output_path": None,
            "confidence": None,
            "error_message": "Boltz completed but produced no prediction files (check GPU/CUDA setup)",
            "command": " ".join(cmd),
            "stderr": result.stderr[:2000] if result.stderr else "",
        }

    confidence = parse_all_model_confidences(prediction_dir, job_name)
    output_path = str(prediction_dir)

    best_idx = confidence.get("best_model_index", 0) if confidence else 0
    best_score = confidence.get("confidence_score") if confidence else None
    num_models = confidence.get("num_models", 1) if confidence else 1

    logger.info(
        "Boltz completed for %s — best model: %d/%d, confidence: %s, output: %s",
        job_name, best_idx, num_models,
        f"{best_score:.3f}" if best_score is not None else "N/A",
        output_path,
    )

    return {
        "success": True,
        "status": "completed",
        "job_name": job_name,
        "output_path": output_path,
        "confidence": confidence,
        "best_model_index": best_idx,
        "error_message": None,
        "command": " ".join(cmd),
        "stderr": "",
    }


def parse_all_model_confidences(
    prediction_dir: Path,
    job_name: str,
) -> Optional[dict]:
    pattern = f"confidence_{job_name}_model_*.json"
    confidence_files = sorted(prediction_dir.glob(pattern))

    if not confidence_files:
        logger.warning("No confidence files matching %s in %s", pattern, prediction_dir)
        return None

    all_models = []
    for cf in confidence_files:
        try:
            stem = cf.stem
            model_idx = int(stem.rsplit("_", 1)[-1])
            with open(cf) as f:
                data = json.load(f)
            score = data.get("confidence_score")
            all_models.append({
                "model_index": model_idx,
                "confidence_score": score,
                "confidence_data": data,
                "file": str(cf),
            })
        except (json.JSONDecodeError, OSError, ValueError) as e:
            logger.warning("Failed to parse confidence file %s: %s", cf, e)

    if not all_models:
        return None

    all_models.sort(key=lambda m: m["confidence_score"] or -1, reverse=True)
    best = all_models[0]

    return {
        "best_model_index": best["model_index"],
        "confidence_score": best["confidence_score"],
        "best_confidence": best["confidence_data"],
        "all_models": [
            {"model_index": m["model_index"], "confidence_score": m["confidence_score"]}
            for m in all_models
        ],
        "num_models": len(all_models),
    }
