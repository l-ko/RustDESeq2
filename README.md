# RustDESeq2

This package is a wrapper around the rust_deseq2 crate. I recommend to check out the description of the crate author @necoli1822's repository for more information: https://github.com/necoli1822/rust_deseq2.

The benefits of using rust_deseq2 against using R implememntation of DESeq2 are enormous. Interface has minor differences to classic R DESeq2 but implements the same model based on the same study. I decided to create a quick wrapper to be loaded as an R module.

## Installation

### Prerequisites
```r
# Install required packages with CRAN mirror
install.packages(c("remotes", "rextendr"), repos="https://cran.rstudio.com/")
```

### Install Package
```r
# Install from GitHub
remotes::install_github("l-ko/RustDESeq2")
```

## Quick Start

```r
library(RustDESeq2)

# Create DESeq dataset from count matrix
dds <- RustDESeq2::deseq_dataset_from_matrix(
  count_data = your_count_vector,
  nrows = n_genes,
  ncols = n_samples,
  gene_ids = gene_names,
  sample_ids = sample_names,
  condition = sample_conditions,
  design_variable = "condition"
)

# Run DESeq2 analysis
dds <- RustDESeq2::deseq(dds)

# Get results
results <- RustDESeq2::results(dds, "treated", "control", alpha = 0.05)

# Get normalized counts
norm_counts <- RustDESeq2::counts(dds, normalized = TRUE)
```

## Features

Currently supports interfaces:

- `deseq_dataset_from_matrix`
- `deseq`
- `results`
- `counts`
- `vst` - Variance stabilizing transformation
- `assay` - Extract transformed data from vst objects

to be used in classic R DEseq2 style (eg. you want to switch your R DEseq2 implementation to Rust for performance reasons).

## System Requirements

- R (>= 4.2)
- Rust (>= 1.65.0) - automatically installed if needed
- Cargo (Rust package manager) - automatically installed if needed

## License

MIT License
