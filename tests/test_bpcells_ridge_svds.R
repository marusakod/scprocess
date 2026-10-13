# Run from the repository root in envs/bpcells_pca.yaml:
# Rscript --vanilla tests/test_bpcells_ridge_svds.R
source("scripts/bpcells_pca.R")

expect_close <- function(actual, expected, tolerance = 1e-8) {
  stopifnot(identical(dim(as.matrix(actual)), dim(as.matrix(expected))))
  error <- max(abs(as.matrix(actual) - as.matrix(expected)))
  if (!is.finite(error) || error > tolerance) {
    stop("Maximum absolute error ", error, " exceeds ", tolerance)
  }
  invisible(error)
}

expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

run_tests <- function() {
  .check_bpcells_residual_api()
  expect_error(.check_bpcells_residual_api("0.4.0"), "validated BPCells 0.3.1")
  .check_svd_convergence(list(nconv = 2L, d = c(2, 1)), 2L)
  for (result in list(list(), list(nconv = 1L, d = 2),
                      list(nconv = NA, d = c(2, 1)),
                      list(nconv = 2L, d = c(2, NaN)))) {
    expect_error(.check_svd_convergence(result, 2L), "did not converge")
  }
  for (threads in list(0, -1, 1.5, NA, Inf, c(1, 2))) {
    expect_error(run_bpcells_pca("unused", 2L, "unused", n_threads = threads),
                 "n_threads must be a positive integer")
  }

  tmp <- tempfile("scprocess-ridge-svds-test-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  set.seed(293)
  n <- 120L
  g <- 60L
  k <- 6L
  cell_ids <- paste0("cell", seq_len(n))
  genes <- paste0("gene", seq_len(g))
  samples <- rep(paste0("sample", 1:3), each = n / 3)
  theta <- rnorm(n, rep(c(0.4, 2, 4), each = n / 3), 0.8) %% (2 * pi)
  tricycle_f <- file.path(tmp, "tricycle.csv.gz")
  fwrite(data.table(cell_id = cell_ids, sample_id = samples,
                   tricycle_theta = theta), tricycle_f)
  design <- .cyclic_design(tricycle_f, cell_ids, 2L)
  stopifnot(design$correction_centering == "within_sample")
  expect_close(design$correction, design$fit)
  global_design <- .cyclic_design(tricycle_f, cell_ids, 2L, "global")
  expect_close(design$fit, global_design$fit)
  expect_error(.cyclic_design(tricycle_f, cell_ids, 2L, "invalid"), "arg")
  expect_error(run_bpcells_pca("unused", k, "unused", correction_centering = "invalid"), "arg")
  stopifnot(max(abs(global_design$group_correction_means)) > 0.1)
  stopifnot(max(abs(colMeans(design$correction))) < 1e-12)
  raw <- matrix(rpois(g * n, 2), g, dimnames = list(genes, cell_ids))
  raw[sample.int(length(raw), floor(length(raw) * 0.5))] <- 0
  cases <- 0L

  for (mode in c("within_sample", "global")) {
  design <- .cyclic_design(tricycle_f, cell_ids, 2L, mode)
  for (order in c("row", "col")) {
    counts <- as(Matrix(raw, sparse = TRUE), "IterableMatrix")
    counts <- convert_matrix_type(counts, "uint32_t")
    if (order == "col") counts <- transpose_storage_order(counts)
    counts <- write_matrix_dir(counts, file.path(tmp, paste0("counts-", mode, "-", order)))
    normalized <- log1p(multiply_cols(counts, 1e4 / BPCells::colSums(counts)))
    stats <- matrix_stats(normalized, row_stats = "variance")$row_stats
    scaled <- (normalized - stats["mean", ]) / sqrt(stats["variance", ])
    Y <- t(min_scalar(scaled, 10))
    stopifnot(storage_order(Y) == order)
    dense <- as.matrix(Y) # Only this tiny test materializes expression.

    for (penalty in c(0, 0.1, 1)) {
      ridge <- .make_ridge_operator(Y, design$fit, design$correction, penalty)
      stopifnot(is(ridge$matrix, "TransformLinearResidual"),
                identical(dimnames(ridge$matrix), dimnames(Y)))
      gram <- crossprod(design$fit)
      effective_penalty <- penalty * sum(diag(gram)) / ncol(design$fit)
      expected_beta <- solve(gram + diag(effective_penalty, ncol(gram)),
                             crossprod(design$fit, dense))
      corrected <- dense - design$correction %*% expected_beta
      expect_close(ridge$coefficients, expected_beta)
      expect_close(ridge$gram, gram)
      expect_close(ridge$correction_gram, crossprod(design$correction))
      expect_close(ridge$xty, crossprod(design$fit, dense))
      expect_close(ridge$lambda_effective, effective_penalty)
      expect_close(ridge$matrix, corrected)
      expect_close(colMeans(corrected), colMeans(dense))
      for (sample in unique(samples)) {
        rows <- which(samples == sample)
        mean_change <- colMeans(corrected[rows, ]) - colMeans(dense[rows, ])
        expected_change <- -colMeans(design$correction[rows, ]) %*% expected_beta
        expect_close(matrix(mean_change, nrow = 1), expected_change)
        if (mode == "within_sample") stopifnot(max(abs(mean_change)) < 1e-10)
      }
      # Add sample-specific technical offsets: within-sample fitting must
      # recover the same coefficients even when phase proportions differ.
      shifted <- dense + matrix(rep(rep(c(-3, 1, 5), each = n / 3), g), n)
      expect_close(crossprod(design$fit, shifted), crossprod(design$fit, dense))

      # The former R operator's forward and reverse multiplication formulas.
      v <- rnorm(g)
      u <- rnorm(n)
      inverse <- solve(gram + diag(effective_penalty, ncol(gram)))
      old_right <- dense %*% v - design$correction %*% inverse %*%
        crossprod(design$fit, dense %*% v)
      old_left <- u %*% dense - (u %*% design$correction) %*% inverse %*%
        crossprod(design$fit, dense)
      expect_close(ridge$matrix %*% v, old_right)
      expect_close(u %*% ridge$matrix, old_left)

      reference <- svd(corrected, nu = k, nv = k)
      for (threads in c(1L, 2L)) {
        fit <- BPCells::svds(ridge$matrix, k = k, threads = threads,
                             opts = list(tol = 1e-10))
        .check_svd_convergence(fit, k)
        expect_close(fit$d, reference$d[seq_len(k)], 1e-6)
        expect_close(tcrossprod(fit$u), tcrossprod(reference$u), 1e-6)
        expect_close(corrected %*% fit$v, sweep(fit$u, 2, fit$d, "*"), 1e-6)
        reverse <- BPCells::svds(t(ridge$matrix), k = k, threads = threads,
                                 opts = list(tol = 1e-10))
        .check_svd_convergence(reverse, k)
        expect_close(reverse$d, fit$d, 1e-6)
        expect_close(tcrossprod(reverse$v), tcrossprod(fit$u), 1e-6)
        cases <- cases + 1L
      }
    }
    singular <- cbind(design$fit[, 1], design$fit[, 1])
    expect_error(.make_ridge_operator(Y, singular, singular, 0), "rank deficient")
    stopifnot(is(.make_ridge_operator(Y, singular, singular, 0.1)$matrix,
                 "TransformLinearResidual"))
    expect_error(.make_ridge_operator(Y, design$fit, design$correction, -1),
                 "non-negative")
    expect_error(.make_ridge_operator(Y, design$fit[-1, ], design$correction, 0.1),
                 "matching cells")
  }
  }

  # End-to-end: multiple HDF5 inputs, selected/reordered cells, full-library
  # normalization, mitochondria exclusion, and regression diagnostic exports.
  first_h5 <- file.path(tmp, "counts-a.h5")
  second_h5 <- file.path(tmp, "counts-b.h5")
  write_matrix_10x_hdf5(convert_matrix_type(
    as(Matrix(raw[, 1:60], sparse = TRUE), "IterableMatrix"), "uint32_t"), first_h5)
  write_matrix_10x_hdf5(convert_matrix_type(
    as(Matrix(raw[, 61:120], sparse = TRUE), "IterableMatrix"), "uint32_t"), second_h5)
  coldata_f <- file.path(tmp, "coldata.csv.gz")
  fwrite(data.table(cell_id = cell_ids, total = colSums(raw) + 100,
                   subsets_mito_sum = 10), coldata_f)
  selected_ids <- rev(cell_ids[seq_len(90)])
  cells_f <- file.path(tmp, "cells.csv")
  fwrite(data.table(cell_id = selected_ids), cells_f)
  selected <- match(selected_ids, cell_ids)
  normalized <- log1p(sweep(raw[, selected], 2, colSums(raw)[selected] + 90, "/") * 1e4)
  reference_base <- t(pmin(sweep(sweep(normalized, 1, rowMeans(normalized), "-"),
                                 1, apply(normalized, 1, sd), "/"), 10))
  for (mode in c("none", "within_sample", "global")) {
    regression <- mode != "none"
    selected_design <- .cyclic_design(tricycle_f, selected_ids, 2L,
                                      if (regression) mode else "within_sample")
    pca_f <- file.path(tmp, paste0("pca-", mode, ".csv.gz"))
    diagnostic_f <- file.path(tmp, "regression.csv")
    genes_f <- file.path(tmp, "regression-genes.csv.gz")
    pcs <- run_bpcells_pca(c(first_h5, second_h5), k, pca_f,
      coldata_f = coldata_f, cells_f = cells_f, exclude_mito = TRUE,
      tricycle_f = if (regression) tricycle_f else "",
      harmonics = if (regression) 2L else NULL,
      ridge_lambda = if (regression) 0.1 else NULL,
      out_regression_f = if (regression) diagnostic_f else "",
      out_regression_genes_f = if (regression) genes_f else "",
      n_threads = 2L, correction_centering = if (regression) mode else "")
    stopifnot(identical(pcs$cell_id, selected_ids), file.exists(pca_f))
    corrected <- reference_base
    if (regression) {
      gram <- crossprod(selected_design$fit)
      lambda <- 0.1 * sum(diag(gram)) / ncol(gram)
      beta <- solve(gram + diag(lambda, ncol(gram)),
                    crossprod(selected_design$fit, reference_base))
      corrected <- reference_base - selected_design$correction %*% beta
      diagnostic <- fread(diagnostic_f)
      stopifnot(diagnostic$n_cells == length(selected_ids),
        diagnostic$coefficient_fit_centring == "within_sample",
        diagnostic$correction_design_centring == mode,
        diagnostic$cell_set_md5 == .cell_set_fingerprint(selected_ids))
      if (mode == "within_sample") {
        stopifnot(diagnostic$predicted_regression_group_mean_change_max < 1e-10)
      } else {
        stopifnot(diagnostic$predicted_regression_group_mean_change_max > 1e-5)
      }
      gene_diagnostic <- fread(genes_f)
      beta_columns <- colnames(selected_design$fit)
      # The production normalized cache can use float storage; the reference
      # above uses doubles throughout. Allow that roundoff in exported fits.
      expect_close(as.matrix(gene_diagnostic[, ..beta_columns]),
                   t(beta), 1e-6)
      stopifnot(identical(gene_diagnostic$gene_id, genes))
    }
    reference <- svd(corrected, nu = k, nv = k)
    pc_matrix <- as.matrix(pcs[, paste0("pca_", seq_len(k)), with = FALSE])
    expected_scores <- sweep(reference$u, 2, reference$d[seq_len(k)], "*")
    expect_close(tcrossprod(pc_matrix), tcrossprod(expected_scores), 1e-3)
  }
  message("PASS: ", cases, " native ridge/SVD configurations, both transpose ",
          "directions, both correction modes, invalid inputs and three end-to-end PCA runs")
}

run_tests()
