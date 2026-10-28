suppressPackageStartupMessages({
  library("BPCells")
  library("Matrix")
  library("data.table")
  library("hdf5r")
})


#' Convert a joint 10x HDF5 count matrix to native BPCells storage.
#'
#' @param counts_h5_f Path to joint counts HDF5 (10x-format CSC matrix)
#' @param out_dir     Output BPCells matrix directory
write_join_bpcells_matrix <- function(counts_h5_f, out_dir) {
  message("Loading count matrix from HDF5")
  h5_readonly <- .open_join_counts_readonly(counts_h5_f)
  on.exit(h5_readonly$close_all(), add = TRUE)
  mat <- open_matrix_10x_hdf5(counts_h5_f)
  message(sprintf("  Matrix dimensions: %d genes x %d cells", nrow(mat), ncol(mat)))

  message("Writing uint32 native BPCells matrix")
  dir.create(dirname(out_dir), recursive = TRUE, showWarnings = FALSE)
  write_matrix_dir(
    convert_matrix_type(mat, "uint32_t"),
    dir = out_dir,
    compress = TRUE
  )
  message(sprintf("Done. Wrote native BPCells matrix to %s", out_dir))
}


#' Compute PCA on a native BPCells joint count matrix.
#'
#' @param counts_bp_dir Path to native BPCells count matrix directory
#' @param n_dims       Number of PCA dimensions to compute
#' @param n_threads    Number of threads for BPCells matrix-vector products
#' @param out_pca_f    Output path for PCA embeddings (CSV.gz)
run_join_pca <- function(counts_bp_dir, n_dims, n_threads, out_pca_f,
                         tricycle_f = NULL, harmonics = 2L, ridge_lambda = 0.1,
                         correction_centering = "within_sample",
                         out_regression_f = NULL, out_regression_genes_f = NULL) {
  correction_centering <- match.arg(correction_centering, c("within_sample", "global"))
  if (length(n_threads) != 1L || !is.finite(n_threads) || n_threads < 1 ||
      n_threads != floor(n_threads)) stop("n_threads must be a positive integer")
  n_dims <- as.integer(n_dims)
  n_threads <- as.integer(n_threads)

  message("Loading native BPCells count matrix (disk-backed)")
  mat <- open_matrix_dir(counts_bp_dir)
  message(sprintf("  Matrix dimensions: %d genes x %d cells", nrow(mat), ncol(mat)))
  if (min(dim(mat)) < 3L || length(n_dims) != 1L || is.na(n_dims) ||
      n_dims < 1L || n_dims >= min(dim(mat))) {
    stop("n_dims must be positive and smaller than both matrix dimensions")
  }

  message("Normalizing: library size -> scale 1e4 -> log1p")
  lib_sizes <- BPCells::colSums(mat)
  if (any(!is.finite(lib_sizes) | lib_sizes <= 0)) stop("Invalid or zero count library size")
  mat_norm <- log1p(multiply_cols(mat, 1e4 / lib_sizes))

  message("Computing per-gene mean and variance (streaming)")

  # BPCells 0.3.x accepts one row statistic; requesting variance also returns mean.
  stats <- matrix_stats(mat_norm, row_stats = "variance")
  gene_means <- stats$row_stats["mean", ]
  gene_vars  <- stats$row_stats["variance", ]
  gene_sds   <- sqrt(gene_vars)

  # Replace zero sds with 1 to avoid division by zero (constant genes)
  gene_sds[gene_sds == 0 | is.na(gene_sds)] <- 1

  message("Queueing delayed per-gene centering and scaling")
  mat_scaled <- (mat_norm - gene_means) / gene_sds

  # Join PCA retains its uncapped scaling; cell-cycle correction is optional.
  pca_matrix <- mat_scaled
  ridge <- NULL
  design_info <- NULL
  if (!is.null(tricycle_f) && nzchar(tricycle_f)) {
    if (!exists(".write_cycle_regression_provenance", mode = "function")) {
      stop("Source scripts/bpcells_pca.R before using join cell-cycle regression")
    }
    message("Regressing source-project tricycle scores; correction centring: ", correction_centering)
    design_info <- .cyclic_design(tricycle_f, colnames(mat), harmonics, correction_centering)
    ridge <- .make_ridge_operator(t(mat_scaled), design_info$fit,
                                  design_info$correction, ridge_lambda)
    pca_matrix <- ridge$matrix
  } else if (!is.null(out_regression_f) || !is.null(out_regression_genes_f)) {
    stop("Regression diagnostics require tricycle_f")
  }

  message(sprintf(
    "Running threaded BPCells truncated SVD (n_dims = %d, threads = %d)",
    n_dims, n_threads))
  svd_result <- BPCells::svds(
    pca_matrix,
    k = n_dims,
    opts = list(tol = 1e-5),
    threads = n_threads
  )
  if (length(svd_result$nconv) != 1L || !is.finite(svd_result$nconv) ||
      svd_result$nconv != n_dims || length(svd_result$d) != n_dims ||
      any(!is.finite(svd_result$d))) {
    stop("BPCells::svds did not converge to all requested PCs")
  }

  message("Computing PC scores (V * D)")

  # Input is genes x cells, so V contains cell embeddings.
  # The regression operator is cells x genes, so its cell embeddings are U.
  cell_vectors <- if (is.null(ridge)) svd_result$v else svd_result$u
  pc_scores <- sweep(cell_vectors, 2, svd_result$d, "*")

  message("Writing PCA embeddings to CSV.gz")
  pca_colnames <- sprintf("pca_%d", seq_len(n_dims))
  pca_dt <- data.table(cell_id = colnames(mat))
  pca_dt[, (pca_colnames) := as.data.table(pc_scores)]

  dir.create(dirname(out_pca_f), recursive = TRUE, showWarnings = FALSE)
  fwrite(pca_dt, out_pca_f)
  if (!is.null(ridge)) {
    .write_cycle_regression_provenance(mat, ridge, design_info, harmonics,
      ridge_lambda, out_regression_f, out_regression_genes_f)
  }

  message(sprintf("Done. Wrote %d cells x %d PCs to %s", nrow(pca_dt), n_dims, out_pca_f))
  invisible(pca_dt)
}


# Validate without modifying an upstream output. Return the read-only handle
# so the caller can retain it throughout BPCells' lazy reads.
.open_join_counts_readonly <- function(h5_path) {
  h5 <- H5File$new(h5_path, mode = "r")
  tryCatch({
    if (!h5$exists("matrix/features/id")) {
      stop("Missing matrix/features/id in ", h5_path,
           "; regenerate this matrix with join_build_matrix. Input left unchanged.")
    }
    h5
  }, error = function(e) {
    h5$close_all()
    stop(e)
  })
}
