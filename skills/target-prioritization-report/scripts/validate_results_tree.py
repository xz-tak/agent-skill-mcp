#!/usr/bin/env python3
"""Validate that required artifacts exist for an offline prioritization report run.

This is an existence check only (no parsing, no scoring). It is safe to run early.

It validates one primary artifact root:
- `<WORKDIR>/results/` (Cortellis/DEG/OFF-X/CI + optional analysis tabs + kgpred_<indication>/)

Legacy layouts are supported as fallbacks:
- `<WORKDIR>/kgpred_<indication>/` (BioBridge/ULTRA/PrimeKG Step1/Step2 reports)
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Optional


def _section_candidates(base: str, indication: Optional[str]) -> list[str]:
    ind = (indication or "").strip().lower()
    if ind:
        return [f"results/{base}_{ind}", f"results/{base}"]
    return [f"results/{base}"]


def _resolve_first_existing(workdir: Path, candidates: list[str]) -> str:
    """
    Return the first candidate that exists; otherwise return the first candidate
    (so missing-path output stays explicit and actionable).
    """
    for rel in candidates:
        if (workdir / rel).exists():
            return rel
    return candidates[0]


def _resolve_results_dir_specs(indication: Optional[str]) -> tuple[list[list[str]], list[list[str]]]:
    """
    Return (required_specs, recommended_specs) where each spec is a list of
    candidate relative paths, ordered by preference.
    """
    ind = (indication or "").strip()
    ind_lower = ind.lower() if ind else ""

    def cand(base: str) -> list[str]:
        # Preferred contract: results/<base>_<indication>/
        # Legacy contract:    results/<base>/
        if ind_lower:
            return [f"results/{base}_{ind_lower}", f"results/{base}"]
        return [f"results/{base}"]

    required = [cand("cortellis"), cand("deg_results"), cand("ci")]
    recommended = [
        cand("offx"),
        cand("pathwaydb"),
        cand("interactdb"),
        cand("bulk_coexpression"),
        cand("sc_coexp"),
    ]
    return required, recommended


def _dir_has_any_glob(root: Path, rel_dir: str, globs: list[str]) -> bool:
    p = root / rel_dir
    if not p.exists() or not p.is_dir():
        return False
    for pat in globs:
        if any(p.glob(pat)):
            return True
    return False


def _dir_has_named_md_anywhere(root: Path, rel_dir: str, names: list[str]) -> bool:
    p = root / rel_dir
    if not p.exists() or not p.is_dir():
        return False
    wanted = {n.lower() for n in names}
    for md in p.rglob("*.md"):
        if md.name.lower() in wanted:
            return True
    return False


def _dir_has_named_file_anywhere(root: Path, rel_dir: str, names: list[str]) -> bool:
    p = root / rel_dir
    if not p.exists() or not p.is_dir():
        return False
    wanted = {n.lower() for n in names}
    for f in p.rglob("*"):
        if f.is_file() and f.name.lower() in wanted:
            return True
    return False


def _iter_missing_any(root: Path, candidates: list[str]) -> list[str]:
    """
    For a single expected artifact, accept any of the candidate paths as satisfying it.
    If none exist, return the first candidate (for explicit missing-path output).
    """
    for rel in candidates:
        if (root / rel).exists():
            return []
    return [candidates[0]]


def _resolve_pathway_report_candidates(indication: Optional[str]) -> list[str]:
    ind = (indication or "").strip().lower()
    if ind:
        base = f"results/pathwaydb_{ind}"
    else:
        base = "results/pathwaydb"
    # Prefer the explicit contract filename, but accept common capitalization variants.
    return [
        f"{base}/COMPREHENSIVE_PATHWAY_ANALYSIS_REPORT.md",
        f"{base}/Comprehensive_Pathway_Analysis_Report.md",
    ]


def _iter_missing(root: Path, rel_paths: Iterable[str]) -> list[str]:
    missing: list[str] = []
    for rel in rel_paths:
        if not (root / rel).exists():
            missing.append(rel)
    return missing


def _resolve_kgpred_root(workdir: Path, indication: str, override: Optional[str]) -> Path:
    if override:
        p = Path(override)
        if p.is_absolute():
            return p
        # Support both `<WORKDIR>/<override>` and `<WORKDIR>/results/<override>`.
        if (workdir / p).exists():
            return workdir / p
        if (workdir / "results" / p).exists():
            return workdir / "results" / p
        return workdir / p

    ind = (indication or "").strip()
    if not ind:
        raise ValueError("indication required when --kgpred is not provided")

    candidates = []
    # Preferred contract: `<WORKDIR>/results/kgpred_<indication>/`
    candidates += [
        workdir / "results" / f"kgpred_{ind}",
        workdir / "results" / f"kgpred_{ind.lower()}",
        workdir / "results" / f"kgpred_{ind.upper()}",
    ]
    # Legacy contract: `<WORKDIR>/kgpred_<indication>/`
    candidates += [
        workdir / f"kgpred_{ind}",
        workdir / f"kgpred_{ind.lower()}",
        workdir / f"kgpred_{ind.upper()}",
    ]
    for c in candidates:
        if c.exists():
            return c

    # Default expected location (even if missing, so output is explicit).
    return workdir / "results" / f"kgpred_{ind.lower()}"


def _autodetect_candidates(workdir: Path) -> dict[str, list[str]]:
    """
    Best-effort discovery of non-standard artifact roots in the working directory.
    This does not change validation behavior; it only prints actionable suggestions.
    """

    def dirs(globs: list[str]) -> list[str]:
        out: list[str] = []
        for pat in globs:
            for p in sorted(workdir.glob(pat)):
                if p.is_dir():
                    try:
                        out.append(str(p.relative_to(workdir)))
                    except Exception:
                        out.append(str(p))
        # de-dupe while preserving order
        seen = set()
        uniq = []
        for x in out:
            if x in seen:
                continue
            seen.add(x)
            uniq.append(x)
        return uniq

    def files(globs: list[str]) -> list[str]:
        out: list[str] = []
        for pat in globs:
            for p in sorted(workdir.glob(pat)):
                if p.is_file():
                    try:
                        out.append(str(p.relative_to(workdir)))
                    except Exception:
                        out.append(str(p))
        seen = set()
        uniq = []
        for x in out:
            if x in seen:
                continue
            seen.add(x)
            uniq.append(x)
        return uniq

    return {
        "kgpred_roots": dirs(["kgpred_*", "results/kgpred_*"]),
        "cortellis_like_dirs": dirs(["*cortellis*", "results/*cortellis*"]),
        "cortellis_xlsx": files(["**/gene_cortellis_data.xlsx"]),
        "deg_like_dirs": dirs(["*deg*", "results/*deg*"]),
        "ci_like_dirs": dirs(["*ci*", "results/*ci*"]),
        "offx_like_dirs": dirs(["*offx*", "results/*offx*"]),
        "pathway_like_dirs": dirs(["*pathway*", "results/*pathway*"]),
        "ppi_like_dirs": dirs(["*interact*", "results/*interact*"]),
        "bulk_coexp_like_dirs": dirs(["*bulk*coexp*", "*bulk_coexp*", "results/*bulk*coexp*", "results/*bulk_coexp*"]),
        "sc_coexp_like_dirs": dirs(["*sc*coexp*", "*sc_coexp*", "results/*sc*coexp*", "results/*sc_coexp*"]),
        "ci_dashboard_html": files(["**/*_CI_Dashboard.html"]),
        "kgpred_part1_md": files(
            [
                "kgpred_*/*/Part1_Individual_Analysis.md",
                "results/kgpred_*/*/Part1_Individual_Analysis.md",
            ]
        ),
        "kgpred_part2_md": files(
            [
                "kgpred_*/*/Part2_Combo_Analysis.md",
                "results/kgpred_*/*/Part2_Combo_Analysis.md",
            ]
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True, help="Abstract project root (<WORKDIR>)")
    ap.add_argument(
        "--indication",
        default=None,
        help=(
            "Indication token used to resolve results/*_<indication>/ and results/kgpred_<indication>/ "
            "(e.g., hs, ibd). Required unless --kgpred is provided."
        ),
    )
    ap.add_argument(
        "--kgpred",
        default=None,
        help="Override kgpred root (absolute path or relative to <WORKDIR>), e.g. kgpred_ibd",
    )
    ap.add_argument(
        "--mode",
        choices=["individual", "combo", "both"],
        default="both",
        help="Which kgpred reports must exist (Part1=individual, Part2=combo).",
    )
    ap.add_argument(
        "--require-ultra",
        action="store_true",
        help="Fail if ULTRA kgpred reports are missing (otherwise warn).",
    )
    args = ap.parse_args()

    workdir = Path(args.workdir).resolve()

    # Resolve section dirs (prefer results/<section>_<indication>/; fallback to results/<section>/).
    cortellis_dir = _resolve_first_existing(workdir, _section_candidates("cortellis", args.indication))
    deg_dir = _resolve_first_existing(workdir, _section_candidates("deg_results", args.indication))
    ci_dir = _resolve_first_existing(workdir, _section_candidates("ci", args.indication))

    offx_dir = _resolve_first_existing(workdir, _section_candidates("offx", args.indication))
    pathway_dir = _resolve_first_existing(workdir, _section_candidates("pathwaydb", args.indication))
    interact_dir = _resolve_first_existing(workdir, _section_candidates("interactdb", args.indication))
    bulk_dir = _resolve_first_existing(
        workdir,
        _section_candidates("bulk_coexpression", args.indication) + _section_candidates("bulk_coexp", args.indication),
    )
    sc_dir = _resolve_first_existing(workdir, _section_candidates("sc_coexp", args.indication))

    missing_required = _iter_missing(workdir, [cortellis_dir, deg_dir, ci_dir])
    missing_recommended = _iter_missing(workdir, [offx_dir, pathway_dir, interact_dir, bulk_dir, sc_dir])

    # File-level checks (recursive within the corresponding section dir).
    if (workdir / cortellis_dir).exists() and not _dir_has_any_glob(
        workdir, cortellis_dir, ["**/*.md", "**/*.json"]
    ):
        missing_required.append(f"{cortellis_dir} (no .md/.json found under this directory)")

    # Cortellis XLSX sidecar (used for gene summary stacked bars; optional but recommended).
    # Accept either a pooled file (`gene_cortellis_data.xlsx`) or per-gene files (`*_cortellis_data.xlsx`).
    if (workdir / cortellis_dir).exists():
        has_pooled = _dir_has_named_file_anywhere(workdir, cortellis_dir, ["gene_cortellis_data.xlsx"])
        has_per_gene = _dir_has_any_glob(workdir, cortellis_dir, ["**/*_cortellis_data.xlsx"])
        if not (has_pooled or has_per_gene):
            missing_recommended.append(f"{cortellis_dir}/**/gene_cortellis_data.xlsx (or {cortellis_dir}/**/*_cortellis_data.xlsx)")

    if (workdir / deg_dir).exists() and not _dir_has_any_glob(workdir, deg_dir, ["**/*.md", "**/*.csv"]):
        missing_required.append(f"{deg_dir} (no .md/.csv found under this directory)")

    if (workdir / ci_dir).exists() and not _dir_has_any_glob(workdir, ci_dir, ["**/*.html"]):
        missing_required.append(f"{ci_dir} (no .html found under this directory)")

    if (workdir / offx_dir).exists() and not _dir_has_any_glob(workdir, offx_dir, ["**/*.md", "**/*.json"]):
        missing_recommended.append(f"{offx_dir} (no .md/.json found under this directory)")

    # Pathways: prefer the comprehensive report, but search recursively under the section dir.
    if (workdir / pathway_dir).exists() and not _dir_has_named_md_anywhere(
        workdir,
        pathway_dir,
        ["COMPREHENSIVE_PATHWAY_ANALYSIS_REPORT.md", "Comprehensive_Pathway_Analysis_Report.md"],
    ):
        # Keep a concrete “missing path” string for readability.
        missing_recommended.append(f"{pathway_dir}/**/COMPREHENSIVE_PATHWAY_ANALYSIS_REPORT.md")

    # PPI: accept SUMMARY_REPORT.md anywhere under the interactdb section dir.
    if (workdir / interact_dir).exists() and not _dir_has_named_md_anywhere(workdir, interact_dir, ["SUMMARY_REPORT.md"]):
        missing_recommended.append(f"{interact_dir}/**/SUMMARY_REPORT.md")

    # kgpred checks
    kgpred_missing_required: list[str] = []
    kgpred_missing_recommended: list[str] = []

    if args.kgpred or args.indication:
        kgpred_root = _resolve_kgpred_root(workdir, args.indication or "", args.kgpred)

        if not kgpred_root.exists():
            kgpred_missing_required.append(str(kgpred_root.relative_to(workdir)))
        else:
            required_methods = ["biobridge", "primekg"]
            recommended_methods = ["ultra"]

            def need_part1() -> bool:
                return args.mode in ("individual", "both")

            def need_part2() -> bool:
                return args.mode in ("combo", "both")

            def has_named_md_anywhere(base: Path, expected: str) -> bool:
                exp = expected.lower()
                for md in base.rglob("*.md"):
                    if md.name.lower() == exp:
                        return True
                return False

            for method in required_methods:
                method_dir = kgpred_root / method
                if need_part1():
                    if not has_named_md_anywhere(method_dir, "Part1_Individual_Analysis.md"):
                        kgpred_missing_required.append(
                            str((method_dir).relative_to(workdir)) + "/**/Part1_Individual_Analysis.md"
                        )
                if need_part2():
                    if not has_named_md_anywhere(method_dir, "Part2_Combo_Analysis.md"):
                        kgpred_missing_required.append(
                            str((method_dir).relative_to(workdir)) + "/**/Part2_Combo_Analysis.md"
                        )

            for method in recommended_methods:
                method_dir = kgpred_root / method
                missing_here: list[str] = []
                if need_part1() and not has_named_md_anywhere(method_dir, "Part1_Individual_Analysis.md"):
                    missing_here.append(str((method_dir).relative_to(workdir)) + "/**/Part1_Individual_Analysis.md")
                if need_part2() and not has_named_md_anywhere(method_dir, "Part2_Combo_Analysis.md"):
                    missing_here.append(str((method_dir).relative_to(workdir)) + "/**/Part2_Combo_Analysis.md")
                if args.require_ultra:
                    kgpred_missing_required += missing_here
                else:
                    kgpred_missing_recommended += missing_here

    if missing_required or kgpred_missing_required:
        print("Missing required paths:")
        for rel in missing_required + kgpred_missing_required:
            print(f"- {rel}")
        print("\nAutodetect (best-effort) to help fill PLAN.md path mappings:")
        cand = _autodetect_candidates(workdir)
        for k, vals in cand.items():
            if not vals:
                continue
            print(f"- {k}:")
            for v in vals[:8]:
                print(f"  - {v}")
            if len(vals) > 8:
                print(f"  - ... (+{len(vals) - 8} more)")
        if missing_recommended or kgpred_missing_recommended:
            print("\nMissing recommended paths (will degrade sections):")
            for rel in missing_recommended + kgpred_missing_recommended:
                print(f"- {rel}")
        return 2

    if missing_recommended or kgpred_missing_recommended:
        print("OK: required artifacts exist.")
        print("Missing recommended paths (will degrade sections):")
        for rel in missing_recommended + kgpred_missing_recommended:
            print(f"- {rel}")
        return 0

    print("OK: required artifacts exist.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
