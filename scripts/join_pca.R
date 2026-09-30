suppressPackageStartupMessages({
  library("BPCells")
  library("Matrix")
  library("data.table")
  library("hdf5r")
})


# Validate without modifying an upstream output. Return the read-only handle
# so the caller can retain it throughout BPCells' lazy reads.
.ensure_10x_hdf5_compat <- function(h5_path) {
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


#' Convert a joint 10x HDF5 count matrix to native BPCells storage.
#'
#' @param counts_h5_f Path to joint counts HDF5 (10x-format CSC matrix)
#' @param out_dir     Output BPCells matrix directory
write_join_bpcells_matrix <- function(counts_h5_f, out_dir) {
  message("Loading count matrix from HDF5")
  h5_readonly <- .ensure_10x_hdf5_compat(counts_h5_f)
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
run_join_pca <- function(counts_bp_dir, n_dims, n_threads, out_pca_f) {
  n_dims <- as.integer(n_dims)
  n_threads <- as.integer(n_threads)

  message("Loading native BPCells count matrix (disk-backed)")
  mat <- open_matrix_dir(counts_bp_dir)
  message(sprintf("  Matrix dimensions: %d genes x %d cells", nrow(mat), ncol(mat)))

  message("Normalizing: library size -> scale 1e4 -> log1p")
  lib_sizes <- BPCells::colSums(mat)
  mat_norm <- multiply_cols(mat, 1e4 / lib_sizes) %>% log1p()

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

  message(sprintf(
    "Running threaded BPCells truncated SVD (n_dims = %d, threads = %d)",
    n_dims, n_threads))
  svd_result <- BPCells::svds(
    mat_scaled,
    k = n_dims,
    opts = list(tol = 1e-5),
    threads = n_threads
  )

  message("Computing PC scores (V * D)")
  # Input is genes x cells, so V contains cell embeddings.
  pc_scores <- sweep(svd_result$v, 2, svd_result$d, "*")

  message("Writing PCA embeddings to CSV.gz")
  pca_colnames <- sprintf("pca_%d", seq_len(n_dims))
  pca_dt <- data.table(cell_id = colnames(mat))
  pca_dt[, (pca_colnames) := as.data.table(pc_scores)]

  dir.create(dirname(out_pca_f), recursive = TRUE, showWarnings = FALSE)
  fwrite(pca_dt, out_pca_f)

  message(sprintf("Done. Wrote %d cells x %d PCs to %s", nrow(pca_dt), n_dims, out_pca_f))
}
