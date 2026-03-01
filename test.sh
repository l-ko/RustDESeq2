#!/bin/bash

# Test script for RustDESeq2 interfaces
# This script tests all functions with presence, type, and functionality checks

echo "=== RustDESeq2 Test Script ==="
echo "Directory: $(pwd)"
echo ""

# Ensure we're in the package root directory
cd "$(dirname "${BASH_SOURCE[0]}")"

echo "Changed to package directory: $(pwd)"
echo ""

# Run the R test script
echo "Running R test script..."
Rscript tests/test.R

echo ""
echo "=== Test Complete ==="
