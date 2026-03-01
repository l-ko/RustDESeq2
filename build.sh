#!/bin/bash

# Build script for RustDESeq2
# This script regenerates the R wrappers from Rust code

echo "=== RustDESeq2 Build Script ==="
echo "Directory: $(pwd)"
echo ""

# Ensure we're in the package root directory
cd "$(dirname "${BASH_SOURCE[0]}")"

echo "Changed to package directory: $(pwd)"
echo ""

# Run the R build script and check for errors
echo "Running R build script..."
if Rscript tests/build.R; then
    echo ""
    echo "✅ Build completed successfully!"
    echo "📦 R wrappers are ready for installation"
else
    echo ""
    echo "❌ Build failed!"
    echo "🔥 Please check the error messages above"
    exit 1
fi

echo ""
echo "=== Build Complete ==="
