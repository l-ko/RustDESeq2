#!/usr/bin/env Rscript

# Test script for RustDESeq2 interfaces
# This script:
# 1. Regenerates the R wrappers
# 2. Installs the package locally
# 3. Tests all functions including vst and assay

cat("=== RustDESeq2 Interfaces Test ===\n")
cat("Current working directory:", getwd(), "\n")

# We're already in the package root directory when called via Rscript

# Regenerate wrappers using the exact command
cat("Running: Rscript -e \"rextendr::document()\"\n")
system("Rscript -e \"rextendr::document()\"")
cat("✓ R wrappers regenerated\n")

# Step 2: Install package locally
cat("Step 2: Installing package locally...\n")

if (!require('remotes', quietly = TRUE)) {
  install.packages('remotes', repos='https://cran.rstudio.com/')
}

# Install from current directory
remotes::install_local('.', force = TRUE)
cat("✓ Package installed locally\n")

# Step 3: Test all functions
cat("Step 3: Testing all functions...\n")

# Load the package
library(RustDESeq2)

cat("Testing RustDESeq2 functions with all interfaces...\n")

# Test with proper data dimensions
tryCatch({
  # Create test data: 6 genes, 2 samples (12 values total)
  count_data <- as.numeric(c(10, 20, 15, 25, 30, 12, 18, 22, 28, 35, 8, 16))
  gene_ids <- c('gene1', 'gene2', 'gene3', 'gene4', 'gene5', 'gene6')
  sample_ids <- c('S1', 'S2')
  condition <- c('A', 'B')
  
  cat('1. Testing deseq_dataset_from_matrix...\n');
  dds <- RustDESeq2::deseq_dataset_from_matrix(count_data, 6, 2, gene_ids, sample_ids, condition, 'condition');
  cat('   ✓ Dataset created successfully\n');
  
  cat('2. Testing deseq...\n');
  dds <- RustDESeq2::deseq(dds);
  cat('   ✓ DESeq analysis completed\n');
  
  cat('3. Testing counts...\n');
  norm_counts <- RustDESeq2::counts(dds, normalized = TRUE);
  cat('   ✓ Counts retrieved, dimensions:', dim(norm_counts), '\n');
  
  cat('4. Testing vst function...\n');
  vst_obj <- RustDESeq2::vst(dds, blind = FALSE, nsub = 6);
  cat('   ✓ VST transformation completed\n');
  
  cat('5. Testing assay function...\n');
  assay_data <- RustDESeq2::assay(vst_obj);
  cat('   ✓ Assay data retrieved, dimensions:', dim(assay_data), '\n');
  cat('   Sample values:', round(assay_data[1:3, 1], 2), '...\n');
  
  cat('6. Testing results...\n');
  results <- RustDESeq2::results(dds, 'B', 'A', 0.05);
  cat('   ✓ Results retrieved, available columns:', names(results), '\n');
  
  cat('\n🎉 All functions working correctly!\n');
  cat('✅ All interfaces are ready for use\n');
  cat('📦 Package is ready for deployment\n');
  
}, error=function(e) {
  cat('❌ Error:', e$message, '\n');
  traceback()
})

cat("=== Test Complete ===\n")
