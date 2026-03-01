#!/usr/bin/env Rscript

# Install script for RustDESeq2
# This script installs the package locally using remotes

cat("=== RustDESeq2 Install Script ===\n")
cat("Current working directory:", getwd(), "\n")

# Step 1: Install package locally (exact command from original working script)
cat("Step 1: Installing package locally...\n")

if (!require('remotes', quietly = TRUE)) {
  install.packages('remotes', repos='https://cran.rstudio.com/')
}

# Install from current directory (exact command from original working script)
remotes::install_local('.', force = TRUE)
cat("✓ Package installed locally\n")

# Step 2: Verify installation
cat("Step 2: Verifying installation...\n")

tryCatch({
  # Try to load the package
  library(RustDESeq2)
  cat("✓ RustDESeq2 loaded successfully\n")
  
  # Check if main functions are available
  main_functions <- c('deseq_dataset_from_matrix', 'deseq', 'results', 'counts', 'vst', 'assay')
  available_functions <- main_functions[sapply(main_functions, function(f) exists(f, mode = "function"))]
  
  cat(sprintf("✓ %d/%d main functions available\n", length(available_functions), length(main_functions)))
  
  if (length(available_functions) == length(main_functions)) {
    cat("✅ All main functions are available\n")
  } else {
    missing_functions <- setdiff(main_functions, available_functions)
    cat("⚠️  Missing functions:", paste(missing_functions, collapse=", "), "\n")
  }
  
}, error = function(e) {
  cat('❌ Error loading package:', e$message, '\n')
  quit(status = 1)
})

cat("\n🎉 Installation completed successfully!\n")
cat("📦 RustDESeq2 is ready for use\n")
cat("=== Install Complete ===\n")
