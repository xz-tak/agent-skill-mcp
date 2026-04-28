#!/usr/bin/env python3
"""Protein structure prediction pipeline: PDB lookup + Boltz AI prediction."""

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from modules.check_pdb import check_and_download_best_pdb
from modules.fetch_sequence import fetch_uniprot_sequence, validate_sequence
from modules.run_boltz import run_boltz_prediction
from modules.utils import (
    detect_available_gpus,
    is_amino_acid_sequence,
    load_config,
    merge_config_with_args,
    sanitize_gene_name,
    setup_logger,
    write_results_summary,
)

logger = setup_logger("predict")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Protein structure prediction: PDB lookup + Boltz prediction",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python predict.py BRCA1                              # predict by gene name
  python predict.py --raw-sequence MVTPEGNVSLVDES...   # predict from sequence
  python predict.py --input-file genes.txt             # batch mode
  python predict.py BRCA1 --skip-boltz-if-pdb          # PDB only if available
  python predict.py EGFR --recycling-steps 5           # custom Boltz params
        """,
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "input",
        nargs="?",
        help="Gene name (e.g. BRCA1) or raw amino acid sequence",
    )
    input_group.add_argument(
        "--input-file",
        type=str,
        help="Path to file with one gene/sequence per line (batch mode)",
    )
    input_group.add_argument(
        "--raw-sequence",
        type=str,
        help="Raw amino acid sequence (skips PDB and UniProt lookup)",
    )

    parser.add_argument(
        "--skip-boltz-if-pdb",
        action="store_true",
        default=False,
        help="Skip Boltz prediction when PDB structure is found",
    )
    parser.add_argument(
        "--skip-pdb",
        action="store_true",
        default=False,
        help="Skip PDB lookup entirely",
    )

    parser.add_argument("--organism", type=str, default=None)
    parser.add_argument("--organism-id", type=int, default=None)

    parser.add_argument("--recycling-steps", type=int, default=None)
    parser.add_argument("--diffusion-samples", type=int, default=None)
    parser.add_argument("--sampling-steps", type=int, default=None)
    parser.add_argument("--output-format", choices=["pdb", "mmcif"], default=None)
    parser.add_argument("--use-potentials", action="store_true", default=None)
    parser.add_argument("--no-msa-server", action="store_true", default=False)
    parser.add_argument("--model", choices=["boltz1", "boltz2"], default=None)
    parser.add_argument("--timeout", type=int, default=None, help="Boltz timeout in seconds")

    parser.add_argument("--rerun-threshold", type=float, default=None,
        help="Re-run predictions below this confidence with boosted params (default: 0.5)")
    parser.add_argument("--max-parallel", type=int, default=None)
    parser.add_argument("--config", type=str, default=None)
    parser.add_argument("--output-dir", type=str, default=None)

    return parser


def classify_input(value: str, force_sequence: bool = False) -> str:
    if force_sequence:
        return "sequence"
    if is_amino_acid_sequence(value):
        return "sequence"
    return "gene"


def build_jobs(args, config: dict) -> list[dict]:
    jobs = []

    if args.raw_sequence:
        if not validate_sequence(args.raw_sequence):
            logger.error("Invalid amino acid sequence: contains non-standard characters")
            sys.exit(1)
        jobs.append({
            "input_value": args.raw_sequence,
            "input_type": "sequence",
            "job_name": "raw_sequence",
        })
    elif args.input_file:
        input_path = Path(args.input_file)
        if not input_path.exists():
            logger.error("Input file not found: %s", args.input_file)
            sys.exit(1)
        for line in input_path.read_text().strip().split("\n"):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            input_type = classify_input(line)
            if input_type == "gene":
                sanitize_gene_name(line)
            job_name = line if input_type == "gene" else f"seq_{len(jobs)}"
            jobs.append({
                "input_value": line,
                "input_type": input_type,
                "job_name": job_name,
            })
    else:
        input_type = classify_input(args.input)
        if input_type == "gene":
            sanitize_gene_name(args.input)
        job_name = args.input if input_type == "gene" else "input_sequence"
        jobs.append({
            "input_value": args.input,
            "input_type": input_type,
            "job_name": job_name,
        })

    return jobs


def process_job_io(job: dict, config: dict, base_dir: Path, skip_pdb: bool) -> dict:
    """Phase A: sequential I/O — PDB check + UniProt fetch for gene inputs."""
    result = {
        **job,
        "pdb_result": None,
        "sequence": job["input_value"] if job["input_type"] == "sequence" else None,
        "uniprot_info": None,
        "errors": [],
    }

    if job["input_type"] != "gene":
        return result

    gene_name = job["input_value"]
    organism_id = config.get("organism_id", 9606)
    pdb_dir_suffix = config.get("output", {}).get("pdb_dir_suffix", "_pdb")

    if not skip_pdb:
        try:
            pdb_result = check_and_download_best_pdb(
                gene_name, base_dir, organism_id, pdb_dir_suffix,
            )
            result["pdb_result"] = pdb_result
        except Exception as e:
            logger.error("PDB check failed for %s: %s", gene_name, e)
            result["errors"].append(f"PDB check failed: {e}")

    try:
        uniprot_info = fetch_uniprot_sequence(gene_name, organism_id)
        result["sequence"] = uniprot_info["sequence"]
        result["uniprot_info"] = uniprot_info
    except ValueError as e:
        logger.error("UniProt lookup failed for %s: %s", gene_name, e)
        result["errors"].append(str(e))
    except Exception as e:
        logger.error("UniProt fetch error for %s: %s", gene_name, e)
        result["errors"].append(f"UniProt error: {e}")

    return result


def _run_single_boltz(kwargs: dict) -> dict:
    return run_boltz_prediction(**kwargs)


def run_boltz_parallel(
    boltz_jobs: list[dict],
    available_gpus: list[int],
    max_parallel: int,
) -> list[dict]:
    if not boltz_jobs:
        return []

    results = []
    workers = min(len(boltz_jobs), max_parallel) if max_parallel > 0 else 1

    for i, job in enumerate(boltz_jobs):
        job["gpu_id"] = available_gpus[i % len(available_gpus)] if available_gpus else 0

    if workers <= 1:
        for job in boltz_jobs:
            results.append(_run_single_boltz(job))
        return results

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(_run_single_boltz, job): job["job_name"]
            for job in boltz_jobs
        }
        for future in as_completed(futures):
            job_name = futures[future]
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                logger.error("Boltz job %s raised exception: %s", job_name, e)
                results.append({
                    "success": False,
                    "status": "exception",
                    "job_name": job_name,
                    "output_path": None,
                    "confidence": None,
                    "error_message": str(e),
                    "command": "",
                    "stderr": "",
                })

    return results


def main():
    parser = build_parser()
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    config_path = Path(args.config) if args.config else script_dir / "config.yaml"
    config = load_config(config_path)

    cli_args = {
        "organism": args.organism,
        "organism_id": args.organism_id,
        "recycling_steps": args.recycling_steps,
        "diffusion_samples": args.diffusion_samples,
        "sampling_steps": args.sampling_steps,
        "output_format": args.output_format,
        "use_potentials": args.use_potentials,
        "no_msa_server": args.no_msa_server,
        "model": args.model,
        "timeout": args.timeout,
        "rerun_threshold": args.rerun_threshold,
        "max_parallel": args.max_parallel,
        "output_dir": args.output_dir,
    }
    config = merge_config_with_args(config, cli_args)

    base_dir = Path(config.get("output", {}).get("base_dir", str(script_dir)))
    boltz_dir = base_dir / config.get("output", {}).get("boltz_dir", ".")
    summary_path = base_dir / config.get("output", {}).get("summary_file", "results_summary.json")

    jobs = build_jobs(args, config)
    logger.info("Processing %d job(s)", len(jobs))

    processed_jobs = []
    for job in jobs:
        processed = process_job_io(job, config, base_dir, args.skip_pdb)
        processed_jobs.append(processed)

    boltz_config = config.get("boltz", {})
    boltz_job_args = []
    for job in processed_jobs:
        if not job["sequence"]:
            logger.warning("Skipping %s — no sequence available", job["job_name"])
            continue

        if args.skip_boltz_if_pdb and job["pdb_result"] is not None:
            logger.info(
                "Skipping Boltz for %s — PDB structure found (%s)",
                job["job_name"], job["pdb_result"]["pdb_id"],
            )
            job["boltz_skipped"] = True
            continue

        boltz_job_args.append({
            "sequence": job["sequence"],
            "job_name": job["job_name"],
            "output_dir": boltz_dir,
            "use_msa_server": boltz_config.get("use_msa_server", True),
            "recycling_steps": boltz_config.get("recycling_steps", 3),
            "diffusion_samples": boltz_config.get("diffusion_samples", 1),
            "sampling_steps": boltz_config.get("sampling_steps", 200),
            "output_format": boltz_config.get("output_format", "mmcif"),
            "use_potentials": boltz_config.get("use_potentials", False),
            "model": boltz_config.get("model", "boltz2"),
            "timeout_seconds": boltz_config.get("timeout_seconds", 3600),
        })

    available_gpus = detect_available_gpus(
        config.get("gpu", {}).get("utilization_threshold", 10.0)
    )
    max_parallel = config.get("gpu", {}).get("max_parallel") or len(available_gpus) or 1
    logger.info("Available GPUs: %s (max_parallel=%d)", available_gpus, max_parallel)

    boltz_results = run_boltz_parallel(boltz_job_args, available_gpus, max_parallel)

    boltz_by_name = {r["job_name"]: r for r in boltz_results}

    rerun_threshold = boltz_config.get("rerun_confidence_threshold", 0.5)
    rerun_jobs = []
    for result in boltz_results:
        conf = (
            result.get("confidence", {}).get("confidence_score")
            if result.get("confidence") else None
        )
        if result.get("success") and conf is not None and conf < rerun_threshold:
            original = next(
                (j for j in boltz_job_args if j["job_name"] == result["job_name"]),
                None,
            )
            if original:
                rerun_jobs.append({
                    **original,
                    "recycling_steps": 10,
                    "diffusion_samples": 5,
                    "use_potentials": True,
                    "job_name": f"{original['job_name']}_rerun",
                })
                logger.info(
                    "Re-running %s (confidence=%.3f < %.1f) with boosted params",
                    result["job_name"], conf, rerun_threshold,
                )

    if rerun_jobs:
        logger.info("Auto re-running %d low-confidence predictions", len(rerun_jobs))
        rerun_results = run_boltz_parallel(rerun_jobs, available_gpus, max_parallel)
        for rerun_result in rerun_results:
            original_name = rerun_result["job_name"].removesuffix("_rerun")
            rerun_conf = (
                rerun_result.get("confidence", {}).get("confidence_score")
                if rerun_result.get("confidence") else None
            )
            original_result = boltz_by_name.get(original_name)
            original_conf = (
                original_result.get("confidence", {}).get("confidence_score")
                if original_result and original_result.get("confidence") else None
            )
            if rerun_conf is not None and (original_conf is None or rerun_conf > original_conf):
                rerun_result["job_name"] = original_name
                rerun_result["rerun"] = True
                boltz_by_name[original_name] = rerun_result
                logger.info(
                    "Replaced %s: %.3f -> %.3f (rerun improved)",
                    original_name, original_conf or 0, rerun_conf,
                )
            else:
                logger.info(
                    "Kept original %s: %.3f (rerun %.3f did not improve)",
                    original_name, original_conf or 0, rerun_conf or 0,
                )

    summary_results = []
    for job in processed_jobs:
        boltz_result = boltz_by_name.get(job["job_name"])
        pdb = job["pdb_result"] or {}
        summary_entry = {
            "job_name": job["job_name"],
            "input_type": job["input_type"],
            "input_value": job["input_value"][:50],
            "sequence_length": len(job["sequence"]) if job["sequence"] else 0,
            "pdb_id": pdb.get("pdb_id"),
            "pdb_resolution": pdb.get("resolution"),
            "pdb_path": pdb.get("file_path"),
            "pdb_title": pdb.get("title"),
            "pdb_method": pdb.get("method"),
            "pdb_deposition_date": pdb.get("deposition_date"),
            "pdb_entity_sequence_length": pdb.get("entity_sequence_length"),
            "pdb_entity_description": pdb.get("entity_description"),
            "pdb_entity_gene_name": pdb.get("entity_gene_name"),
            "pdb_entity_organism": pdb.get("entity_organism"),
            "pdb_entity_weight_kda": pdb.get("entity_weight_kda"),
            "uniprot_accession": (
                job["uniprot_info"]["accession"]
                if job.get("uniprot_info") else None
            ),
            "boltz_ran": boltz_result is not None,
            "boltz_rerun": boltz_result.get("rerun", False) if boltz_result else False,
            "boltz_skipped": job.get("boltz_skipped", False),
            "boltz_status": boltz_result["status"] if boltz_result else None,
            "boltz_output_dir": boltz_result["output_path"] if boltz_result else None,
            "boltz_confidence": (
                boltz_result["confidence"].get("confidence_score")
                if boltz_result and boltz_result.get("confidence")
                else None
            ),
            "boltz_best_model": (
                boltz_result.get("best_model_index", 0)
                if boltz_result else None
            ),
            "boltz_num_models": (
                boltz_result["confidence"].get("num_models", 1)
                if boltz_result and boltz_result.get("confidence")
                else None
            ),
            "boltz_all_confidences": (
                boltz_result["confidence"].get("all_models")
                if boltz_result and boltz_result.get("confidence")
                else None
            ),
            "errors": job["errors"] + (
                [boltz_result["error_message"]]
                if boltz_result and boltz_result.get("error_message")
                else []
            ),
        }
        summary_results.append(summary_entry)

    write_results_summary(summary_path, summary_results)
    logger.info("Results written to %s", summary_path)

    _print_summary(summary_results)


def _print_summary(results: list[dict]) -> None:
    total = len(results)
    succeeded = sum(1 for r in results if r.get("boltz_status") == "completed")
    failed = sum(1 for r in results if r.get("boltz_status") in ("failed", "timeout", "exception"))
    skipped = sum(1 for r in results if r.get("boltz_skipped"))
    no_seq = sum(1 for r in results if not r.get("boltz_ran") and not r.get("boltz_skipped"))

    print("\n" + "=" * 60)
    print("PREDICTION SUMMARY")
    print("=" * 60)

    for r in results:
        status_icon = "OK" if r.get("boltz_status") == "completed" else (
            "SKIP" if r.get("boltz_skipped") else
            "FAIL" if r.get("boltz_status") in ("failed", "timeout") else
            "ERR"
        )
        pdb_info = f"PDB:{r['pdb_id']}" if r.get("pdb_id") else "no PDB"
        pdb_len = f"pdb_len={r['pdb_entity_sequence_length']}" if r.get("pdb_entity_sequence_length") else ""
        method = r.get("pdb_method", "") or ""

        model_info = ""
        all_conf = r.get("boltz_all_confidences") or []
        if r.get("boltz_confidence") is not None:
            best_idx = r.get("boltz_best_model", 0)
            n_models = r.get("boltz_num_models", 1)
            scores = [m["confidence_score"] for m in all_conf if m.get("confidence_score") is not None]
            if len(scores) > 1:
                model_info = f"best=model_{best_idx} conf={r['boltz_confidence']:.3f} ({min(scores):.3f}-{max(scores):.3f}, n={n_models})"
            else:
                model_info = f"conf={r['boltz_confidence']:.3f}"

        print(f"  [{status_icon}] {r['job_name']:20s} | {pdb_info:15s} | {pdb_len:14s} | {method:20s} | {model_info}")

    print("-" * 60)
    print(f"Total: {total} | Completed: {succeeded} | Failed: {failed} | Skipped: {skipped}")
    if no_seq > 0:
        print(f"No sequence available: {no_seq}")
    if failed > 0:
        print("\nFailed jobs:")
        for r in results:
            if r.get("boltz_status") in ("failed", "timeout", "exception"):
                errors = r.get("errors", [])
                msg = errors[-1] if errors else "unknown error"
                print(f"  - {r['job_name']}: {msg[:100]}")
    print("=" * 60)


if __name__ == "__main__":
    main()
