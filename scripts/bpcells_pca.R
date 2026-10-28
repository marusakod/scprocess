suppressPackageStartupMessages({
  library("BPCells")
  library("Matrix")
  library("data.table")
})


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
#' @param n_threads Threads for native BPCells SVD matrix-vector products.
run_bpcells_pca <- function(
    counts_h5_fs, n_dims, out_pca_f, coldata_f = NULL, cells_f = NULL,
    exclude_mito = FALSE, scale_factor = 1e4, max_value = 10,
    tricycle_f = NULL, harmonics = 2L, ridge_lambda = 0.1,
    out_regression_f = NULL, out_regression_genes_f = NULL, n_threads = 1L,
    correction_centering = "within_sample") {
  if (!identical(correction_centering, "") ||
      (!is.null(tricycle_f) && nzchar(tricycle_f))) {
    correction_centering <- match.arg(correction_centering, c("within_sample", "global"))
  }
  if (length(n_threads) != 1L || !is.finite(n_threads) ||
      n_threads < 1 || n_threads != floor(n_threads)) {
    stop("n_threads must be a positive integer")
  }
  n_threads <- as.integer(n_threads)
  n_dims <- as.integer(n_dims)
  exclude_mito <- as.logical(exclude_mito)

  message("Opening count matrix/matrices with BPCells")
  mat <- .open_and_combine_10x(counts_h5_fs)
  mat <- .select_and_order_cells(mat, cells_f)
  message(sprintf("  Matrix dimensions: %d genes x %d cells", nrow(mat), ncol(mat)))
  if (min(dim(mat)) < 3L || n_dims < 1L || n_dims >= min(dim(mat))) {
    stop("PCA requires at least 3 cells and genes, and n_dims must be positive ",
         "and smaller than both matrix dimensions")
  }

  message("Normalizing: library size -> scale 1e4 -> log1p")
  lib_sizes <- .get_library_sizes(mat, coldata_f, exclude_mito)
  mat_norm <- log1p(multiply_cols(mat, scale_factor / lib_sizes))

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

  # BPCells SVD expects cells by genes here; U then contains cell scores.
  pca_matrix <- t(mat_scaled)
  ridge <- NULL
  design_info <- NULL
  if (!is.null(tricycle_f) && nzchar(tricycle_f)) {
    message("Building within-sample fit; correction centring: ", correction_centering)
    design_info <- .cyclic_design(tricycle_f, colnames(mat), harmonics, correction_centering)
    ridge <- .make_ridge_operator(
      pca_matrix, design_info$fit, design_info$correction, ridge_lambda
    )
    pca_matrix <- ridge$matrix
  }

  message(sprintf("Running threaded BPCells SVD (n_dims = %d, threads = %d)",
                  n_dims, n_threads))
  set.seed(20230308L)
  svd_result <- BPCells::svds(
    pca_matrix, k = n_dims, threads = n_threads, opts = list(tol = 1e-5)
  )
  .check_svd_convergence(svd_result, n_dims)

  message("Writing PC scores")
  pc_scores <- sweep(svd_result$u, 2L, svd_result$d, "*")
  pca_colnames <- sprintf("pca_%d", seq_len(n_dims))
  pca_dt <- data.table(cell_id = colnames(mat))
  pca_dt[, (pca_colnames) := as.data.table(pc_scores)]

  dir.create(dirname(out_pca_f), recursive = TRUE, showWarnings = FALSE)
  fwrite(pca_dt, out_pca_f)

  .write_cycle_regression_provenance(mat, ridge, design_info, harmonics,
    ridge_lambda, out_regression_f, out_regression_genes_f)
  message(sprintf("Done. Wrote %d cells x %d PCs to %s", nrow(pca_dt), n_dims, out_pca_f))
  invisible(pca_dt)
}


# TransformLinearResidual is the native low-rank subtraction used by regress_out.
# Its constructor is internal: keep the environment pinned and fail explicitly
# on unvalidated versions rather than silently changing the correction model.
.check_bpcells_residual_api <- function(version = utils::packageVersion("BPCells")) {
  if (as.character(version) != "0.3.1") {
    stop("Cell-cycle ridge correction requires the validated BPCells 0.3.1 ",
         "residual API; use envs/bpcells_pca.yaml (found ", version, ").")
  }
}


.write_cycle_regression_provenance <- function(mat, ridge, design_info, harmonics,
                                              ridge_lambda, out_regression_f,
                                              out_regression_genes_f) {
  if (!is.null(out_regression_f) && nzchar(out_regression_f)) {
    if (is.null(ridge)) stop("out_regression_f requires tricycle regression")
    component_variance <- colSums(
      ridge$coefficients * (ridge$correction_gram %*% ridge$coefficients)
    ) / max(1, nrow(design_info$fit) - 1L)
    group_mean_changes <- design_info$group_correction_means %*% ridge$coefficients
    provenance <- data.table(
      n_cells = nrow(design_info$fit),
      cell_set_md5 = .cell_set_fingerprint(colnames(mat)),
      harmonics = as.integer(harmonics),
      ridge_lambda = as.numeric(ridge_lambda),
      ridge_lambda_effective = ridge$lambda_effective,
      coefficient_fit_centring = "within_sample",
      correction_design_centring = design_info$correction_centering,
      design_condition_number = kappa(ridge$gram),
      max_abs_within_sample_design_mean = design_info$max_sample_mean,
      max_abs_global_correction_design_mean = design_info$max_global_correction_mean,
      design_rms = paste(signif(design_info$rms, 8L), collapse = ";"),
      coefficient_rms = sqrt(mean(ridge$coefficients^2)),
      n_regression_groups = nrow(design_info$group_correction_means),
      predicted_regression_group_mean_change_rms = sqrt(mean(group_mean_changes^2)),
      predicted_regression_group_mean_change_max = max(abs(group_mean_changes)),
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
}


.check_svd_convergence <- function(result, n_dims) {
  if (length(result$nconv) != 1L || !is.finite(result$nconv) ||
      result$nconv != n_dims || length(result$d) != n_dims ||
      any(!is.finite(result$d))) {
    stop("BPCells::svds did not converge to all ", n_dims, " requested PCs.")
  }
}


.cyclic_design <- function(tricycle_f, cell_ids, harmonics,
                           correction_centering = "within_sample") {
  correction_centering <- match.arg(correction_centering, c("within_sample", "global"))
  scores <- fread(tricycle_f)
  required <- c("cell_id", "sample_id", "tricycle_theta")
  missing <- setdiff(required, names(scores))
  if (length(missing)) stop("tricycle table is missing: ", paste(missing, collapse = ", "))
  if (anyDuplicated(scores$cell_id)) stop("tricycle table contains duplicate cell IDs")
  idx <- match(cell_ids, scores$cell_id)
  if (anyNA(idx)) stop("Some PCA cells are missing tricycle scores")
  scores <- scores[idx]

  # Centre within biological samples; unassigned cells use their run group.
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

  # Each harmonic contributes sine and cosine columns for the same cell order.
  harmonics <- as.integer(harmonics)
  columns <- unlist(lapply(seq_len(harmonics), function(h) {
    list(sin(h * scores$tricycle_theta), cos(h * scores$tricycle_theta))
  }), recursive = FALSE)
  raw_design <- do.call(cbind, columns)
  colnames(raw_design) <- unlist(lapply(seq_len(harmonics), function(h) {
    c(sprintf("sin_%dtheta", h), sprintf("cos_%dtheta", h))
  }))

  # Estimate shared slopes only from within-sample associations. This is the
  # fixed-effect fit without materializing a sample-indicator matrix.
  fit_design <- raw_design
  groups <- split(seq_len(nrow(scores)), center_groups)
  for (rows in groups) {
    fit_design[rows, ] <- sweep(
      fit_design[rows, , drop = FALSE], 2L,
      colMeans(fit_design[rows, , drop = FALSE]), "-"
    )
  }
  rms <- sqrt(colMeans(fit_design^2))
  if (any(!is.finite(rms)) || any(rms < sqrt(.Machine$double.eps))) {
    stop("Cyclic design is degenerate after within-sample centring")
  }
  fit_design <- sweep(fit_design, 2L, rms, "/")

  # By default preserve each sample mean and leave between-sample shifts to
  # integration. Global correction also subtracts predicted phase-composition
  # shifts between samples. Coefficient fitting is identical in both modes.
  correction_design <- fit_design
  if (correction_centering == "global") {
    correction_design <- sweep(raw_design, 2L, colMeans(raw_design), "-")
    correction_design <- sweep(correction_design, 2L, rms, "/")
  }
  max_sample_mean <- max(abs(unlist(lapply(groups, function(rows) {
    colMeans(fit_design[rows, , drop = FALSE])
  }))))
  group_correction_means <- do.call(rbind, lapply(groups, function(rows) {
    colMeans(correction_design[rows, , drop = FALSE])
  }))
  list(
    fit = fit_design,
    correction_centering = correction_centering,
    correction = correction_design,
    rms = rms,
    max_sample_mean = max_sample_mean,
    max_global_correction_mean = max(abs(colMeans(correction_design))),
    group_correction_means = group_correction_means
  )
}


.make_ridge_operator <- function(base, fit_design, correction_design, ridge_lambda) {
  .check_bpcells_residual_api()
  if (!is(fit_design, "matrix") || !is(correction_design, "matrix") ||
      !identical(dim(fit_design), dim(correction_design)) ||
      nrow(fit_design) != nrow(base) || ncol(fit_design) < 1L ||
      any(!is.finite(fit_design)) || any(!is.finite(correction_design))) {
    stop("Ridge designs must be finite matrices with matching cells and harmonics.")
  }
  if (length(ridge_lambda) != 1L || !is.finite(ridge_lambda) || ridge_lambda < 0) {
    stop("ridge_lambda must be finite and non-negative")
  }

  # Solve the small harmonic system; the cell-by-gene residual stays lazy.
  gram <- crossprod(fit_design)
  lambda_effective <- as.numeric(ridge_lambda) * sum(diag(gram)) / ncol(fit_design)
  if (lambda_effective == 0 && qr(gram)$rank < ncol(gram)) {
    stop("The cyclic design is rank deficient; use a positive ridge_lambda")
  }
  inverse <- solve(gram + diag(lambda_effective, ncol(fit_design)))
  xty <- do.call(rbind, lapply(seq_len(ncol(fit_design)), function(j) {
    as.numeric(fit_design[, j] %*% base)
  }))
  coefficients <- inverse %*% xty

  # For cells x genes, subtract correction_design %*% coefficients lazily.
  # Parameters follow physical storage order, not the logical transpose flag.
  row_params <- t(correction_design)
  col_params <- coefficients
  if (storage_order(base) == "row") {
    row_params <- coefficients
    col_params <- t(correction_design)
  }
  operator <- BPCells:::wrapMatrix(
    "TransformLinearResidual", convert_matrix_type(base, "double"),
    row_params = row_params,
    col_params = col_params,
    global_params = numeric(),
    vars_to_regress = if (is.null(colnames(fit_design))) character() else colnames(fit_design)
  )
  list(
    matrix = operator,
    gram = gram,
    correction_gram = crossprod(correction_design),
    lambda_effective = lambda_effective,
    xty = xty,
    coefficients = coefficients
  )
}


.cell_set_fingerprint <- function(cell_ids) {
  path <- tempfile("pca_cells_", fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  writeLines(cell_ids, path, useBytes = TRUE)
  unname(tools::md5sum(path))
}


# Add the ID dataset required by BPCells to project-generated count matrices.
.add_missing_10x_feature_ids <- function(h5_path) {
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
    .add_missing_10x_feature_ids(path)
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
