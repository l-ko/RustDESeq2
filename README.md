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
  count_data = your_count_matrix,
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

# Get VST transformed data
vst_obj <- RustDESeq2::vst_transform(dds, blind = FALSE)
vst_data <- RustDESeq2::assay(vst_obj)
```

## Features

Currently supports functions:

- `deseq_dataset_from_matrix` - Create DESeqDataSet from count matrix
- `deseq` - Run DESeq2 differential expression analysis
- `results` - Extract differential expression results with contrast specification

- `counts` - Get raw or normalized counts from DESeqDataSet
- `assay` - Extract data from DESeqDataSet or VST objects (supports "counts", "normalized", etc.)
- `vst_transform` - Variance stabilizing transformation using proper DESeq2 algorithm

- `dispersions` - Get dispersion estimates for each gene
- `dispersion_function` - Get dispersion function coefficients
- `base_mean` - Get mean of normalized counts for each gene
- `size_factors` - Get size factors for each sample

- `design_matrix` - Get the design matrix used in the analysis
- `coefficients` - Get model coefficients (log2 fold changes)
- `coefficient_se` - Get standard errors for model coefficients

## System Requirements

- R (>= 4.2)
- Rust (>= 1.65.0)
- Cargo (Rust package manager)

## License

MIT License
