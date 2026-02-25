#!/bin/bash

# Test script for RustDESeq2 interfaces
# This script runs the R test script to verify all functions

echo "=== Running RustDESeq2 Test ==="
echo "Directory: $(pwd)"
echo ""

# Ensure we're in the package root directory
cd "$(dirname "${BASH_SOURCE[0]}")"

echo "Changed to package directory: $(pwd)"
echo ""

# Run the R test script
Rscript tests/test_interfaces.R

echo ""
echo "=== Test Finished ==="
