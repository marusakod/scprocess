suppressPackageStartupMessages({
  library("BPCells")
  library("Matrix")
  library("irlba")
  library("data.table")
  library("hdf5r")
})


# Ensure HDF5 has all datasets required by BPCells open_matrix_10x_hdf5()
.ensure_10x_hdf5_compat <- function(h5_path) {
  h5 <- H5File$new(h5_path, mode = "r+")
  on.exit(h5$close_all())

  grp <- h5[["matrix/features"]]
  if (!"id" %in% names(grp)) {
    message("  Adding missing 'matrix/features/id' dataset to HDF5")
    grp[["id"]] <- grp[["name"]]$read()
  }
  invisible(NULL)
}


#' Compute PCA on a joint count matrix using BPCells disk-backed streaming.
#'
#' @param counts_h5_f  Path to joint counts HDF5 (10x-format CSC matrix)
#' @param n_dims       Number of PCA dimensions to compute
#' @param out_pca_f    Output path for PCA embeddings (CSV.gz)
run_join_pca <- function(counts_h5_f, n_dims, out_pca_f) {
  n_dims <- as.integer(n_dims)

  message("Loading count matrix from HDF5 (disk-backed)")
  .ensure_10x_hdf5_compat(counts_h5_f)
  mat <- open_matrix_10x_hdf5(counts_h5_f)
  message(sprintf("  Matrix dimensions: %d genes x %d cells", nrow(mat), ncol(mat)))

  message("Normalizing: library size -> scale 1e4 -> log1p")
  lib_sizes <- BPCells::colSums(mat)
  mat_norm <- multiply_cols(mat, 1e4 / lib_sizes) %>% log1p()

  # Transpose to cells x genes for PCA (irlba center/scale apply to columns)
  mat_t <- t(mat_norm)

  message("Computing per-gene mean and variance (streaming)")
  stats <- matrix_stats(mat_norm, row_stats = "variance")
  gene_means <- stats$row_stats["mean", ]
  gene_vars  <- stats$row_stats["variance", ]
  gene_sds   <- sqrt(gene_vars)

  # Replace zero sds with 1 to avoid division by zero (constant genes)
  gene_sds[gene_sds == 0] <- 1

  message(sprintf("Running truncated SVD (n_dims = %d)", n_dims))
  svd_result <- irlba(mat_t, nv = n_dims, center = gene_means, scale = gene_sds)

  message("Computing PC scores (U * D)")
  # svd_result$u is (n_cells x n_dims), d is singular values
  pc_scores <- sweep(svd_result$u, 2, svd_result$d, "*")

  message("Writing PCA embeddings to CSV.gz")
  pca_colnames <- sprintf("pca_%d", seq_len(n_dims))
  pca_dt <- data.table(cell_id = colnames(mat))
  pca_dt[, (pca_colnames) := as.data.table(pc_scores)]

  dir.create(dirname(out_pca_f), recursive = TRUE, showWarnings = FALSE)
  fwrite(pca_dt, out_pca_f)

  message(sprintf("Done. Wrote %d cells x %d PCs to %s", nrow(pca_dt), n_dims, out_pca_f))
}
