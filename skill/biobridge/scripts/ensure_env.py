#!/usr/bin/env python3
"""
Environment Check and Setup Wrapper

This script checks if the biobridge conda environment exists and is properly configured.
If not, it triggers the setup process. Use this before running any BioBridge tasks.

Usage:
    python ensure_env.py
    python ensure_env.py --check-only
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path


ENV_NAME = "biobridge"
SCRIPT_DIR = Path(__file__).parent
SKILL_DIR = SCRIPT_DIR.parent
SETUP_SCRIPT = SCRIPT_DIR / "setup_env.sh"


def run_command(cmd, capture_output=True):
    """Run a shell command and return the result."""
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=capture_output,
            text=True,
            timeout=300
        )
        return result.returncode == 0, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return False, "", "Command timed out"
    except Exception as e:
        return False, "", str(e)


def check_conda_available():
    """Check if conda is available."""
    success, _, _ = run_command("command -v conda")
    return success


def check_env_exists():
    """Check if the conda environment exists."""
    success, stdout, _ = run_command("conda env list")
    if success:
        for line in stdout.split('\n'):
            if line.strip().startswith(ENV_NAME + ' ') or line.strip().startswith(ENV_NAME + '\t'):
                return True
    return False


def check_packages_installed():
    """Check if required packages are installed in the environment."""
    required_packages = [
        "torch", "numpy", "pandas", "openai", "transformers",
        "scipy", "openpyxl", "scikit-learn"
    ]

    success, stdout, _ = run_command(f"conda run -n {ENV_NAME} pip list")
    if not success:
        return False

    installed = stdout.lower()
    for package in required_packages:
        if package.lower() not in installed:
            print(f"✗ Missing package: {package}", file=sys.stderr)
            return False

    return True


def setup_environment():
    """Run the environment setup script."""
    print("Setting up BioBridge environment...", file=sys.stderr)
    print("This may take several minutes...", file=sys.stderr)
    print("", file=sys.stderr)

    # Run setup script
    success, stdout, stderr = run_command(f"bash {SETUP_SCRIPT}", capture_output=False)

    if success:
        print("", file=sys.stderr)
        print("✓ Environment setup completed successfully", file=sys.stderr)
        return True
    else:
        print("", file=sys.stderr)
        print("✗ Environment setup failed", file=sys.stderr)
        if stderr:
            print(stderr, file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Ensure BioBridge environment is set up correctly"
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Only check environment status, don't set up if missing"
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress output (only return exit code)"
    )

    args = parser.parse_args()

    # Check conda
    if not check_conda_available():
        if not args.quiet:
            print("✗ Conda is not available", file=sys.stderr)
            print("Please install Miniconda or Anaconda", file=sys.stderr)
        sys.exit(1)

    if not args.quiet:
        print("✓ Conda is available", file=sys.stderr)

    # Check environment exists
    if not check_env_exists():
        if not args.quiet:
            print(f"✗ Conda environment '{ENV_NAME}' does not exist", file=sys.stderr)

        if args.check_only:
            sys.exit(1)

        # Setup environment
        if not args.quiet:
            print("", file=sys.stderr)

        if setup_environment():
            sys.exit(0)
        else:
            sys.exit(1)

    if not args.quiet:
        print(f"✓ Conda environment '{ENV_NAME}' exists", file=sys.stderr)

    # Check packages
    if not check_packages_installed():
        if not args.quiet:
            print(f"✗ Required packages are missing in '{ENV_NAME}'", file=sys.stderr)

        if args.check_only:
            sys.exit(1)

        if not args.quiet:
            print("Re-running setup to install missing packages...", file=sys.stderr)

        if setup_environment():
            sys.exit(0)
        else:
            sys.exit(1)

    if not args.quiet:
        print(f"✓ All required packages are installed", file=sys.stderr)
        print("", file=sys.stderr)
        print("Environment is ready!", file=sys.stderr)

    sys.exit(0)


if __name__ == "__main__":
    main()
