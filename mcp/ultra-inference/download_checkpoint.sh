#!/bin/bash
# Download ULTRA checkpoint from GitHub if it doesn't already exist

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CKPTS_DIR="${SCRIPT_DIR}/ckpts"
CHECKPOINT_FILE="ultra_primekg_50g_ft_epoch_1.pth"
CHECKPOINT_PATH="${CKPTS_DIR}/${CHECKPOINT_FILE}"

# Direct GitHub raw URL for the checkpoint
GITHUB_URL="https://raw.githubusercontent.com/roger-tu/ULTRA/b149f9d42921047b58475d2d18929d864a2321a7/ckpts/ultra_primekg_50g_ft_epoch_1.pth"

echo "============================================================"
echo "ULTRA MCP Server - Checkpoint Setup"
echo "============================================================"
echo ""

# Create ckpts directory if it doesn't exist
mkdir -p "${CKPTS_DIR}"

# Check if checkpoint already exists
if [ -f "${CHECKPOINT_PATH}" ]; then
    echo "✓ Checkpoint already exists:"
    echo "  ${CHECKPOINT_PATH}"
    ls -lh "${CHECKPOINT_PATH}"
    echo ""
    echo "No download needed."
    echo "============================================================"
    exit 0
fi

# Checkpoint doesn't exist, download it
echo "Checkpoint not found. Downloading from GitHub..."
echo "Source: ${GITHUB_URL}"
echo ""

# Download with wget (with progress bar)
if command -v wget &> /dev/null; then
    wget --progress=bar:force -O "${CHECKPOINT_PATH}" "${GITHUB_URL}"
# Fallback to curl if wget is not available
elif command -v curl &> /dev/null; then
    curl -L --progress-bar -o "${CHECKPOINT_PATH}" "${GITHUB_URL}"
else
    echo "✗ Error: Neither wget nor curl is installed."
    echo "  Please install wget or curl and try again."
    exit 1
fi

# Verify download
if [ -f "${CHECKPOINT_PATH}" ]; then
    echo ""
    echo "✓ Checkpoint successfully downloaded!"
    echo "  Location: ${CHECKPOINT_PATH}"
    ls -lh "${CHECKPOINT_PATH}"
    echo ""

    # Check file size (checkpoints are typically 100MB+)
    FILE_SIZE=$(stat -f%z "${CHECKPOINT_PATH}" 2>/dev/null || stat -c%s "${CHECKPOINT_PATH}" 2>/dev/null)
    FILE_SIZE_MB=$((FILE_SIZE / 1024 / 1024))

    echo "  File size: ${FILE_SIZE_MB} MB"

    if [ ${FILE_SIZE} -lt 10000000 ]; then
        echo ""
        echo "⚠ Warning: File size seems small (${FILE_SIZE_MB} MB)."
        echo "   This might not be a valid checkpoint file."
        echo "   Expected size: 100-500 MB"
        echo ""
        echo "   The download may have failed. Please check:"
        echo "   - Network connection"
        echo "   - GitHub repository access"
        echo "   - File URL: ${GITHUB_URL}"
        exit 1
    else
        echo ""
        echo "✓ MCP server is ready to use!"
    fi
else
    echo ""
    echo "✗ Download failed."
    echo "   Please check network connection and try again."
    echo ""
    echo "   If the problem persists, you can:"
    echo "   1. Download manually from:"
    echo "      https://github.com/roger-tu/ULTRA/blob/b149f9d42921047b58475d2d18929d864a2321a7/ckpts/ultra_primekg_50g_ft_epoch_1.pth"
    echo ""
    echo "   2. Copy from parent ULTRA directory:"
    echo "      cp ../ckpts/ultra_primekg_50g_ft_epoch_1.pth ./ckpts/"
    exit 1
fi

echo "============================================================"
