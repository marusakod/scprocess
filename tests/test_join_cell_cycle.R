# Run from the repository root in the workflow-managed BPCells environment.
source("scripts/bpcells_pca.R")
source("scripts/join_pca.R")

expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

run_tests <- function() {
  tmp <- tempfile("join-cycle-test-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  set.seed(718)
  n <- 120L
  g <- 60L
  k <- 6L
  ids <- paste0("cell", seq_len(n))
  raw <- matrix(rpois(g * n, 4), g,
                dimnames = list(paste0("gene", seq_len(g)), ids))
  counts <- as(Matrix(raw, sparse = TRUE), "IterableMatrix")
  counts_dir <- file.path(tmp, "counts")
  write_matrix_dir(convert_matrix_type(counts, "uint32_t"), counts_dir)
  files <- list.files(counts_dir, full.names = TRUE, recursive = TRUE)
  before <- tools::md5sum(files)
  before_mtime <- file.info(files)$mtime
  scores_f <- file.path(tmp, "scores.csv.gz")
  scores <- data.table(cell_id = ids,
    sample_id = rep(c("p1:s1", "p1:s2", "p2:s1"), each = 40L),
    tricycle_theta = runif(n, 0, 2 * pi))
  # Regression must align scores by ID rather than input row order.
  fwrite(scores[sample(.N)], scores_f)
  normalized <- log1p(sweep(raw, 2L, colSums(raw) / 1e4, "/"))
  scaled <- scale(t(normalized))
  for (mode in c("none", "within_sample", "global")) {
    expected <- scaled
    cycle_args <- list()
    if (mode != "none") {
      design <- .cyclic_design(scores_f, ids, 2L, mode)
      gram <- crossprod(design$fit)
      effective_lambda <- 0.1 * sum(diag(gram)) / ncol(gram)
      beta <- solve(gram + diag(effective_lambda, ncol(gram)),
                    crossprod(design$fit, scaled))
      expected <- scaled - design$correction %*% beta
      cycle_args <- list(tricycle_f = scores_f, correction_centering = mode,
        out_regression_f = file.path(tmp, paste0(mode, "-summary.csv.gz")),
        out_regression_genes_f = file.path(tmp, paste0(mode, "-genes.csv.gz")))
    }
    fit <- do.call(run_join_pca, c(list(counts_bp_dir = counts_dir, n_dims = k,
      n_threads = 2L, out_pca_f = file.path(tmp, paste0(mode, "-pca.csv.gz"))), cycle_args))
    stopifnot(identical(fit$cell_id, ids), ncol(fit) == k + 1L)
    actual <- as.matrix(fit[, -"cell_id"])
    reference <- svd(expected, nu = k, nv = k)
    target <- sweep(reference$u, 2L, reference$d[seq_len(k)], "*")
    stopifnot(max(abs(tcrossprod(actual) - tcrossprod(target))) < 1e-4)
    if (mode != "none") {
      summary <- fread(cycle_args$out_regression_f)
      gene_summary <- fread(cycle_args$out_regression_genes_f)
      stopifnot(summary$n_cells == n, summary$harmonics == 2L,
        summary$correction_design_centring == mode,
        summary$n_regression_groups == 3L, nrow(gene_summary) == g)
      actual_beta <- t(as.matrix(gene_summary[, c("sin_1theta", "cos_1theta",
                                                "sin_2theta", "cos_2theta")]))
      coefficient_error <- max(abs(actual_beta - beta))
      message("Dense-reference coefficient error: ", signif(coefficient_error, 5L))
      # BPCells' streaming normalized/scaled matrix may use float32 internally.
      stopifnot(coefficient_error < 1e-6)
      if (mode == "within_sample") {
        stopifnot(summary$predicted_regression_group_mean_change_max < 1e-10)
      }
    }
  }
  bad_f <- file.path(tmp, "missing.csv.gz")
  fwrite(scores[-1], bad_f)
  expect_error(run_join_pca(counts_dir, k, 1L, file.path(tmp, "bad.csv"),
                           tricycle_f = bad_f), "missing tricycle scores")
  scores[1, tricycle_theta := NA_real_]
  fwrite(scores, bad_f)
  expect_error(run_join_pca(counts_dir, k, 1L, file.path(tmp, "bad.csv"),
                           tricycle_f = bad_f), "non-finite")
  scores[1, tricycle_theta := 0.5]
  fwrite(rbind(scores, scores[1]), bad_f)
  expect_error(run_join_pca(counts_dir, k, 1L, file.path(tmp, "bad.csv"),
                           tricycle_f = bad_f), "duplicate cell IDs")
  stopifnot(identical(before, tools::md5sum(files)),
            identical(before_mtime, file.info(files)$mtime))
  message("PASS: baseline and both join correction modes match dense PCA; diagnostics and read-only inputs verified")
}
run_tests()
