suppressPackageStartupMessages({
  library("BPCells")
  library("Matrix")
  library("irlba")
  library("data.table")
})


# A cells-by-genes low-rank ridge residual operator. The underlying expression
# matrix remains a BPCells IterableMatrix; only X and crossprod(X, Y) are held
# in memory. These methods are the public custom-matrix interface used by irlba.
setClass(
  "RidgeResidualMatrix",
  contains = "IterableMatrix",
  slots = c(base = "IterableMatrix", design = "matrix", inverse = "matrix", xty = "matrix")
)

setMethod("%*%", signature(x = "RidgeResidualMatrix", y = "numeric"),
  function(x, y) {
    base_product <- x@base %*% y
    base_product - x@design %*% x@inverse %*% crossprod(x@design, base_product)
  })

setMethod("%*%", signature(x = "numeric", y = "RidgeResidualMatrix"),
  function(x, y) {
    base_product <- x %*% y@base
    base_product - (x %*% y@design) %*% y@inverse %*% y@xty
  })


.cyclic_design <- function(tricycle_f, cell_ids, harmonics) {
  scores <- fread(tricycle_f)
  required <- c("cell_id", "sample_id", "tricycle_theta")
  missing <- setdiff(required, names(scores))
  if (length(missing)) stop("tricycle table is missing: ", paste(missing, collapse = ", "))
  if (anyDuplicated(scores$cell_id)) stop("tricycle table contains duplicate cell IDs")
  idx <- match(cell_ids, scores$cell_id)
  if (anyNA(idx)) stop("Some PCA cells are missing tricycle scores")
  scores <- scores[idx]
  center_groups <- paste0("sample:", scores$sample_id)
  unassigned <- is.na(scores$sample_id) | scores$sample_id == ""
  if (any(unassigned)) {
    if (!"tricycle_projection_group" %in% names(scores)) {
      stop("Unassigned tricycle cells require tricycle_projection_group")
    }
    center_groups[unassigned] <- scores$tricycle_projection_group[unassigned]
  }
  if (anyNA(center_groups) || any(center_groups == "")) {
    stop("tricycle centring groups contain missing values")
  }
  if (any(!is.finite(scores$tricycle_theta))) stop("tricycle_theta contains non-finite values")

  harmonics <- as.integer(harmonics)
  columns <- unlist(lapply(seq_len(harmonics), function(h) {
    list(sin(h * scores$tricycle_theta), cos(h * scores$tricycle_theta))
  }), recursive = FALSE)
  design <- do.call(cbind, columns)
  colnames(design) <- unlist(lapply(seq_len(harmonics), function(h) {
    c(sprintf("sin_%dtheta", h), sprintf("cos_%dtheta", h))
  }))
  groups <- split(seq_len(nrow(scores)), center_groups)
  for (rows in groups) {
    design[rows, ] <- sweep(design[rows, , drop = FALSE], 2L,
                            colMeans(design[rows, , drop = FALSE]), "-")
  }
  rms <- sqrt(colMeans(design^2))
  if (any(!is.finite(rms)) || any(rms < sqrt(.Machine$double.eps))) {
    stop("Cyclic design is degenerate after within-sample centring")
  }
  design <- sweep(design, 2L, rms, "/")
  max_sample_mean <- max(abs(unlist(lapply(groups, function(rows) {
    colMeans(design[rows, , drop = FALSE])
  }))))
  list(matrix = design, rms = rms, max_sample_mean = max_sample_mean)
}


.make_ridge_operator <- function(base, design, ridge_lambda) {
  gram <- crossprod(design)
  lambda_effective <- as.numeric(ridge_lambda) * sum(diag(gram)) / ncol(design)
  if (lambda_effective == 0 && qr(gram)$rank < ncol(gram)) {
    stop("The cyclic design is rank deficient; use a positive ridge_lambda")
  }
  inverse <- solve(gram + diag(lambda_effective, ncol(design)))
  xty <- do.call(rbind, lapply(seq_len(ncol(design)), function(j) {
    as.numeric(design[, j] %*% base)
  }))
  operator <- new(
    "RidgeResidualMatrix",
    base = base,
    design = design,
    inverse = inverse,
    xty = xty,
    dim = dim(base),
    transpose = FALSE,
    dimnames = dimnames(base)
  )
  list(
    matrix = operator,
    gram = gram,
    lambda_effective = lambda_effective,
    xty = xty,
    coefficients = inverse %*% xty
  )
}


.cell_set_fingerprint <- function(cell_ids) {
  path <- tempfile("pca_cells_", fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  writeLines(cell_ids, path, useBytes = TRUE)
  unname(tools::md5sum(path))
}


# Ensure HDF5 has all datasets required by BPCells open_matrix_10x_hdf5().
.ensure_10x_hdf5_compat <- function(h5_path) {
  h5 <- hdf5r::H5File$new(h5_path, mode = "r+")
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
    exclude_mito = FALSE, scale_factor = 1e4, max_value = 10,
    tricycle_f = NULL, harmonics = 2L, ridge_lambda = 0.1,
    out_regression_f = NULL, out_regression_genes_f = NULL) {
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

  pca_matrix <- t(mat_scaled)
  ridge <- NULL
  design_info <- NULL
  if (!is.null(tricycle_f) && nzchar(tricycle_f)) {
    message("Building within-sample cyclic design and lazy ridge operator")
    design_info <- .cyclic_design(tricycle_f, colnames(mat), harmonics)
    ridge <- .make_ridge_operator(pca_matrix, design_info$matrix, ridge_lambda)
    pca_matrix <- ridge$matrix
  }

  message(sprintf("Running truncated SVD (n_dims = %d)", n_dims))
  set.seed(20230308L)
  svd_result <- irlba(pca_matrix, nv = n_dims, nu = n_dims)
  if (svd_result$iter <= 0L) stop("irlba did not perform any iterations")

  message("Writing PC scores")
  pc_scores <- sweep(svd_result$u, 2L, svd_result$d, "*")
  pca_colnames <- sprintf("pca_%d", seq_len(n_dims))
  pca_dt <- data.table(cell_id = colnames(mat))
  pca_dt[, (pca_colnames) := as.data.table(pc_scores)]

  dir.create(dirname(out_pca_f), recursive = TRUE, showWarnings = FALSE)
  fwrite(pca_dt, out_pca_f)

  if (!is.null(out_regression_f) && nzchar(out_regression_f)) {
    if (is.null(ridge)) stop("out_regression_f requires tricycle regression")
    component_variance <- colSums(
      ridge$coefficients * (ridge$gram %*% ridge$coefficients)
    ) / max(1, nrow(design_info$matrix) - 1L)
    provenance <- data.table(
      n_cells = nrow(design_info$matrix),
      cell_set_md5 = .cell_set_fingerprint(colnames(mat)),
      harmonics = as.integer(harmonics),
      ridge_lambda = as.numeric(ridge_lambda),
      ridge_lambda_effective = ridge$lambda_effective,
      design_condition_number = kappa(ridge$gram),
      max_abs_within_sample_design_mean = design_info$max_sample_mean,
      design_rms = paste(signif(design_info$rms, 8L), collapse = ";"),
      coefficient_rms = sqrt(mean(ridge$coefficients^2)),
      cyclic_component_variance_mean = mean(component_variance),
      cyclic_component_variance_median = median(component_variance),
      cyclic_component_variance_max = max(component_variance)
    )
    dir.create(dirname(out_regression_f), recursive = TRUE, showWarnings = FALSE)
    fwrite(provenance, out_regression_f)
    if (!is.null(out_regression_genes_f) && nzchar(out_regression_genes_f)) {
      gene_provenance <- data.table(gene_id = rownames(mat))
      gene_provenance[, (rownames(ridge$coefficients)) := as.data.table(t(ridge$coefficients))]
      gene_provenance[, cyclic_component_variance := component_variance]
      dir.create(dirname(out_regression_genes_f), recursive = TRUE, showWarnings = FALSE)
      fwrite(gene_provenance, out_regression_genes_f)
    }
  }
  message(sprintf("Done. Wrote %d cells x %d PCs to %s", nrow(pca_dt), n_dims, out_pca_f))
  invisible(pca_dt)
}


# Compatibility wrapper for the existing join rule.
run_join_pca <- function(counts_h5_f, n_dims, out_pca_f) {
  # Preserve the established join behavior, which did not cap scaled values.
  run_bpcells_pca(counts_h5_f, n_dims, out_pca_f, max_value = NULL)
}
