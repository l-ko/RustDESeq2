#!/bin/bash

# Install script for RustDESeq2
# This script installs the package locally using remotes

echo "=== RustDESeq2 Install Script ==="
echo "Directory: $(pwd)"
echo ""

# Ensure we're in the package root directory
cd "$(dirname "${BASH_SOURCE[0]}")"

echo "Changed to package directory: $(pwd)"
echo ""

# Run the R install script
echo "Running R install script..."
Rscript tests/install.R

echo ""
echo "=== Install Complete ==="
