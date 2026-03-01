#!/usr/bin/env Rscript

# Build script for RustDESeq2
# This script regenerates the R wrappers from Rust code

cat("=== RustDESeq2 Build Script ===\n")
cat("Current working directory:", getwd(), "\n")

# Regenerate wrappers using the exact command from original working script
cat("Running: Rscript -e \"rextendr::document()\"\n")
result <- system("Rscript -e \"rextendr::document()\"", intern = TRUE)

# Check if compilation failed
if (any(grepl("ERROR: compilation failed", result))) {
  cat("❌ Compilation failed!\n")
  cat("🔥 Build errors detected:\n")
  # Show the error lines
  error_lines <- result[grepl("error\\[", result) | grepl("ERROR:", result)]
  if (length(error_lines) > 0) {
    for (line in error_lines) {
      cat("   ", line, "\n")
    }
  }
  cat("\n❌ Build failed!\n")
  quit(status = 1)
} else {
  cat("✓ R wrappers regenerated\n")
}

# Check if wrapper files were created
wrapper_files <- c('R/extendr-wrappers.R', 'man/RustDESeq2.Rd')
missing_files <- wrapper_files[!file.exists(wrapper_files)]

if (length(missing_files) > 0) {
  cat("❌ Missing wrapper files:\n")
  for (file in missing_files) {
    cat(sprintf("   - %s\n", file))
  }
  cat("\n❌ Build failed!\n")
  quit(status = 1)
} else {
  for (file in wrapper_files) {
    cat(sprintf("✓ %s created\n", file))
  }
}

cat("\n🎉 Build completed successfully!\n")
cat("📦 R wrappers are ready for installation\n")
cat("=== Build Complete ===\n")
