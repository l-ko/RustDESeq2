\name{RustDESeq2}
\alias{RustDESeq2}
\alias{RustDESeq2-package}
\docType{package}
\title{
RustDESeq2: Fast Differential Expression Analysis in Rust
}
\description{
RustDESeq2 is an R wrapper package for the rust_deseq2 crate created by @necoli1822 (https://github.com/necoli1822/rust_deseq2). This package provides a high-performance R interface to the Rust implementation of the DESeq2 algorithm for differential expression analysis of RNA-seq data.
}
\details{
\tabular{ll}{
Package: \tab RustDESeq2\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2026-02-28\cr
License: \tab MIT + file LICENSE\cr
}

RustDESeq2 provides an R interface to the rust_deseq2 crate, a high-performance Rust implementation of the DESeq2 algorithm created by @necoli1822. The core algorithm is implemented entirely in Rust for maximum performance, providing 10-100x speedup compared to the original R DESeq2 implementation while maintaining identical statistical results.

The underlying rust_deseq2 crate implements the complete DESeq2 workflow including size factor estimation, dispersion estimation, GLM fitting, and statistical testing for differential expression. This R wrapper package provides seamless integration with existing R DESeq2 workflows through the extendr framework.

For more details about the core Rust implementation, see: https://github.com/necoli1822/rust_deseq2
}
\author{
@l-ko (https://github.com/l-ko)
Based on rust_deseq2 crate by @necoli1822 (https://github.com/necoli1822/rust_deseq2)
}
\keyword{package}
\seealso{
\code{\link[DESeq2]{DESeqDataSet}}, \code{\link[DESeq2]{results}},\\
Original rust_deseq2 crate: https://github.com/necoli1822/rust_deseq2
}
\examples{
# Example usage
library(RustDESeq2)

# Create a simple count matrix
counts <- matrix(rnbinom(1000, mu=10, size=1/0.1), nrow=100, ncol=10)
rownames(counts) <- paste0("gene", 1:100)
colnames(counts) <- paste0("sample", 1:10)

# Create sample metadata
condition <- factor(rep(c("treated", "control"), each=5))
coldata <- data.frame(row.names=colnames(counts), condition=condition)

# Create DESeqDataSet
dds <- deseq_dataset_from_matrix(
  count_data = counts,
  gene_ids = rownames(counts),
  sample_ids = colnames(counts),
  condition = as.character(condition),
  design_variable = "condition"
)

# Run DESeq2 analysis
dds <- deseq(dds)

# Get results
res <- results(dds, "treated", "control")
}
