suppressPackageStartupMessages({
  library("BPCells")
  library("Matrix")
  library("irlba")
  library("data.table")
  library("hdf5r")
})


# Ensure HDF5 has all datasets required by BPCells open_matrix_10x_hdf5().
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


.read_cell_ids <- function(cells_f) {
  cells <- fread(cells_f, select = 1L)[[1L]]
  if (!length(cells)) stop("Cell selection is empty: ", cells_f)
  if (anyNA(cells) || anyDuplicated(cells)) {
    stop("Cell selection must contain unique, non-missing IDs: ", cells_f)
  }
  as.character(cells)
}


.select_and_order_cells <- function(mat, cells_f = NULL) {
  if (is.null(cells_f) || !nzchar(cells_f)) return(mat)

  cells <- .read_cell_ids(cells_f)
  idx <- match(cells, colnames(mat))
  if (anyNA(idx)) {
    missing <- cells[is.na(idx)]
    stop(
      "Selected cell IDs are absent from the PCA matrix: ",
      paste(head(missing, 20L), collapse = ", "),
      if (length(missing) > 20L) " ..." else ""
    )
  }
  mat[, idx]
}


.get_library_sizes <- function(mat, coldata_f = NULL, exclude_mito = FALSE) {
  if (is.null(coldata_f) || !nzchar(coldata_f)) return(BPCells::colSums(mat))

  coldata <- fread(coldata_f)
  if (!"cell_id" %in% names(coldata)) stop("coldata has no cell_id column")
  if (anyDuplicated(coldata$cell_id)) stop("coldata contains duplicate cell_id values")

  idx <- match(colnames(mat), coldata$cell_id)
  if (anyNA(idx)) stop("Some PCA matrix cells are absent from coldata")

  required <- if (exclude_mito) c("total", "subsets_mito_sum") else "total"
  missing <- setdiff(required, names(coldata))
  if (length(missing)) stop("coldata is missing: ", paste(missing, collapse = ", "))

  lib_sizes <- coldata$total[idx]
  if (exclude_mito) lib_sizes <- lib_sizes - coldata$subsets_mito_sum[idx]
  if (anyNA(lib_sizes) || any(!is.finite(lib_sizes)) || any(lib_sizes <= 0)) {
    stop("PCA library sizes must be finite and greater than zero")
  }
  lib_sizes
}


.open_and_combine_10x <- function(counts_h5_fs) {
  counts_h5_fs <- as.character(counts_h5_fs)
  if (!length(counts_h5_fs)) stop("At least one counts HDF5 is required")

  mats <- lapply(counts_h5_fs, function(path) {
    .ensure_10x_hdf5_compat(path)
    open_matrix_10x_hdf5(path)
  })
  reference_features <- rownames(mats[[1L]])
  for (i in seq_along(mats)[-1L]) {
    if (!identical(reference_features, rownames(mats[[i]]))) {
      stop("All PCA matrices must contain identical features in identical order")
    }
  }
  if (length(mats) == 1L) mats[[1L]] else do.call(cbind, mats)
}


#' Compute PCA from one or more 10x HDF5 matrices using BPCells streaming.
#'
#' Normalization and scaling match the project integration path: library-size
#' normalization to 1e4, log1p, per-gene centring and unit-variance scaling,
#' followed by an upper cap of 10. The scaled matrix remains lazy.
#'
#' @param counts_h5_fs One or more compatible 10x HDF5 count matrices.
#' @param n_dims Number of principal components.
#' @param out_pca_f Output CSV.gz containing cell_id and pca_* columns.
#' @param coldata_f Optional cell metadata used for full-transcriptome library
#'   sizes. When absent, selected-matrix column sums are used (join behavior).
#' @param cells_f Optional one-column ordered cell selection file.
#' @param exclude_mito Subtract subsets_mito_sum from total library size.
#' @param scale_factor Library-size normalization target.
#' @param max_value Upper cap applied after gene scaling.
run_bpcells_pca <- function(
    counts_h5_fs, n_dims, out_pca_f, coldata_f = NULL, cells_f = NULL,
    exclude_mito = FALSE, scale_factor = 1e4, max_value = 10) {
  n_dims <- as.integer(n_dims)
  exclude_mito <- as.logical(exclude_mito)

  message("Opening count matrix/matrices with BPCells")
  mat <- .open_and_combine_10x(counts_h5_fs)
  mat <- .select_and_order_cells(mat, cells_f)
  message(sprintf("  Matrix dimensions: %d genes x %d cells", nrow(mat), ncol(mat)))
  if (n_dims < 1L || n_dims >= min(dim(mat))) {
    stop("n_dims must be positive and smaller than both matrix dimensions")
  }

  message("Normalizing: library size -> scale 1e4 -> log1p")
  lib_sizes <- .get_library_sizes(mat, coldata_f, exclude_mito)
  mat_norm <- multiply_cols(mat, scale_factor / lib_sizes) %>% log1p()

  # PCA makes many passes. Persist only the sparse, normalized representation;
  # centring, scaling and clipping remain queued to avoid a dense matrix.
  cache_dir <- tempfile("bpcells_pca_norm_")
  on.exit(unlink(cache_dir, recursive = TRUE, force = TRUE), add = TRUE)
  mat_norm <- write_matrix_dir(mat_norm, cache_dir, compress = FALSE)

  message("Computing per-gene mean and variance (streaming)")
  stats <- matrix_stats(mat_norm, row_stats = "variance")$row_stats
  gene_means <- stats["mean", ]
  gene_sds <- sqrt(stats["variance", ])
  gene_sds[!is.finite(gene_sds) | gene_sds == 0] <- 1

  message("Queuing gene centring, scaling and upper clipping")
  mat_scaled <- (mat_norm - gene_means) / gene_sds
  if (!is.null(max_value) && is.finite(max_value)) {
    mat_scaled <- min_scalar(mat_scaled, as.numeric(max_value))
  }

  message(sprintf("Running truncated SVD (n_dims = %d)", n_dims))
  set.seed(20230308L)
  svd_result <- irlba(t(mat_scaled), nv = n_dims, nu = n_dims)
  if (svd_result$iter <= 0L) stop("irlba did not perform any iterations")

  message("Writing PC scores")
  pc_scores <- sweep(svd_result$u, 2L, svd_result$d, "*")
  pca_colnames <- sprintf("pca_%d", seq_len(n_dims))
  pca_dt <- data.table(cell_id = colnames(mat))
  pca_dt[, (pca_colnames) := as.data.table(pc_scores)]

  dir.create(dirname(out_pca_f), recursive = TRUE, showWarnings = FALSE)
  fwrite(pca_dt, out_pca_f)
  message(sprintf("Done. Wrote %d cells x %d PCs to %s", nrow(pca_dt), n_dims, out_pca_f))
  invisible(pca_dt)
}


# Compatibility wrapper for the existing join rule.
run_join_pca <- function(counts_h5_f, n_dims, out_pca_f) {
  # Preserve the established join behavior, which did not cap scaled values.
  run_bpcells_pca(counts_h5_f, n_dims, out_pca_f, max_value = NULL)
}
