#!/bin/bash
# BioBridge Environment Setup Script
# This script checks if the biobridge conda environment exists and sets it up if needed

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILL_DIR="$(dirname "$SCRIPT_DIR")"
ENV_NAME="biobridge"
PYTHON_VERSION="3.11"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "================================================"
echo "BioBridge Environment Setup"
echo "================================================"
echo ""

# Function to check if conda is available
check_conda() {
    if command -v conda &> /dev/null; then
        return 0
    else
        return 1
    fi
}

# Function to check if conda environment exists
check_conda_env() {
    if check_conda; then
        if conda env list | grep -q "^${ENV_NAME} "; then
            return 0
        fi
    fi
    return 1
}

# Function to install uv
install_uv() {
    echo -e "${YELLOW}Installing uv (fast Python package installer)...${NC}"
    if command -v uv &> /dev/null; then
        echo -e "${GREEN}✓ uv is already installed${NC}"
        return 0
    fi

    # Install uv using the official installer
    curl -LsSf https://astral.sh/uv/install.sh | sh

    # Add uv to PATH for current session
    export PATH="$HOME/.cargo/bin:$PATH"

    if command -v uv &> /dev/null; then
        echo -e "${GREEN}✓ uv installed successfully${NC}"
        return 0
    else
        echo -e "${RED}✗ Failed to install uv${NC}"
        return 1
    fi
}

# Function to create conda environment with uv
create_conda_env_with_uv() {
    echo -e "${YELLOW}Creating conda environment: ${ENV_NAME}${NC}"

    # Create conda environment with Python only
    conda create -n "${ENV_NAME}" python="${PYTHON_VERSION}" -y

    echo -e "${GREEN}✓ Conda environment created${NC}"

    # Activate environment and install packages with uv
    echo -e "${YELLOW}Installing packages with uv...${NC}"

    # Get the environment's Python path
    CONDA_PREFIX=$(conda run -n "${ENV_NAME}" python -c "import sys; print(sys.prefix)")

    # Use uv to install packages from pyproject.toml
    cd "${SKILL_DIR}"
    conda run -n "${ENV_NAME}" bash -c "
        # Ensure uv is available
        export PATH=\"\$HOME/.cargo/bin:\$PATH\"

        # Install packages using uv
        uv pip install -e .
    "

    echo -e "${GREEN}✓ Packages installed successfully${NC}"
}

# Function to create environment with pip (fallback)
create_conda_env_fallback() {
    echo -e "${YELLOW}Creating conda environment with pip (fallback)...${NC}"

    # Create conda environment with Python
    conda create -n "${ENV_NAME}" python="${PYTHON_VERSION}" -y

    # Install packages with pip
    cd "${SKILL_DIR}"
    conda run -n "${ENV_NAME}" pip install -e .

    echo -e "${GREEN}✓ Environment created with pip${NC}"
}

# Main setup logic
main() {
    # Check if conda is available
    if ! check_conda; then
        echo -e "${RED}✗ Conda is not available${NC}"
        echo "Please install Miniconda or Anaconda first."
        echo "Visit: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi

    echo -e "${GREEN}✓ Conda is available${NC}"

    # Check if environment already exists
    if check_conda_env; then
        echo -e "${GREEN}✓ Conda environment '${ENV_NAME}' already exists${NC}"
        echo ""
        echo "To use the environment, run:"
        echo "  conda activate ${ENV_NAME}"
        echo ""
        echo "To recreate the environment, first remove it with:"
        echo "  conda env remove -n ${ENV_NAME}"
        echo "Then run this script again."
        exit 0
    fi

    echo -e "${YELLOW}Environment '${ENV_NAME}' not found. Setting up...${NC}"
    echo ""

    # Install uv
    if install_uv; then
        # Try to create environment with uv
        if create_conda_env_with_uv; then
            echo ""
            echo -e "${GREEN}================================================${NC}"
            echo -e "${GREEN}✓ BioBridge environment setup complete!${NC}"
            echo -e "${GREEN}================================================${NC}"
            echo ""
            echo "To activate the environment, run:"
            echo "  conda activate ${ENV_NAME}"
            echo ""
            return 0
        else
            echo -e "${YELLOW}uv installation failed, falling back to pip...${NC}"
            create_conda_env_fallback
        fi
    else
        # Fallback to pip if uv installation fails
        echo -e "${YELLOW}Could not install uv, using pip instead...${NC}"
        create_conda_env_fallback
    fi

    echo ""
    echo -e "${GREEN}================================================${NC}"
    echo -e "${GREEN}✓ BioBridge environment setup complete!${NC}"
    echo -e "${GREEN}================================================${NC}"
    echo ""
    echo "To activate the environment, run:"
    echo "  conda activate ${ENV_NAME}"
    echo ""
}

# Run main function
main
