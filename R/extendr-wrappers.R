# nolint start

#' @docType package
#' @usage NULL
#' @useDynLib RustDESeq2, .registration = TRUE
NULL

# Check for required dependencies on package load
.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("rextendr", quietly = TRUE)) {
    stop("Package 'rextendr' is required but not installed. Please install with: install.packages('rextendr')")
  }
}

#' Create a DESeq dataset (Rust-backed)
#' @export
deseq_dataset_from_matrix <- function(count_data, nrows, ncols, gene_ids, sample_ids, condition, design_variable) .Call(wrap__deseq_dataset_from_matrix, count_data, nrows, ncols, gene_ids, sample_ids, condition, design_variable)

#' Run DESeq pipeline (Rust-backed)
#' @export
deseq <- function(dds) .Call(wrap__deseq, dds)

#' Get results for numerator/denominator contrast (Rust-backed)
#' @export
results <- function(dds, numerator, denominator, alpha) .Call(wrap__results, dds, numerator, denominator, alpha)

#' Get counts (Rust-backed)
#' @export
counts <- function(dds, normalized = FALSE) .Call(wrap__counts, dds, normalized)

# nolint end
