# Run each tiny regression in its existing workflow environment:
# envs/rlibs_bpcells.yaml: Rscript tests/test_join_readonly_hdf5.R markers
# envs/bpcells_pca.yaml:   Rscript tests/test_join_readonly_hdf5.R pca
suppressPackageStartupMessages({
  library(BPCells)
  library(Matrix)
})
mode <- match.arg(commandArgs(trailingOnly = TRUE)[1], c("markers", "pca"))
if (mode == "markers") source("scripts/marker_genes_edger_bp.R")
if (mode == "pca") source("scripts/join_pca.R")

run_tests <- function() {
  tmp <- tempfile("join-readonly-test-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  raw <- matrix(c(1, 0, 2, 0, 3, 4, 5, 1, 0, 2, 4, 6), 3,
                dimnames = list(paste0("g", 1:3), paste0("c", 1:4)))
  if (mode == "markers") {
  paths <- file.path(tmp, c("a.h5ad", "b.h5ad"))
  for (i in 1:2) {
    invisible(write_matrix_anndata_hdf5(
      convert_matrix_type(Matrix(raw[, (2*i-1):(2*i)], sparse = TRUE), "uint32_t"), paths[i]))
    Sys.setFileTime(paths[i], as.POSIXct("2020-01-01", tz = "UTC"))
  }
  mtimes <- file.info(paths)$mtime
  hashes <- tools::md5sum(paths)
  integration <- file.path(tmp, "integration.csv")
  data.table::fwrite(data.table::data.table(cell_id = colnames(raw),
      sample_id = "s1", project_id = "p1", RNA_snn_res.0.5 = rep(c("cl1", "cl2"), 2)), integration)
  manifest <- file.path(tmp, "h5ads.yaml")
  yaml::write_yaml(as.list(setNames(paths, c("a", "b"))), manifest)
  # Exercise serial and forked reads, lazy evaluation, and cross-file aggregation.
  for (workers in c(1L, 2L)) {
    work <- file.path(tmp, paste0("tmp", workers))
    dir.create(work)
    cache <- build_pseudobulk_cache_bpcells(file.path(tmp, paste0("pb", workers)),
      integration, manifest, "0.5", "sample_id", 1L, workers, work)
    expected <- cbind(cl1 = raw[, 1] + raw[, 3], cl2 = raw[, 2] + raw[, 4])
    expected <- expected[, cache$coldata$cluster, drop = FALSE]
    stopifnot(isTRUE(all.equal(unname(as.matrix(cache$counts)), unname(expected))),
              identical(mtimes, file.info(paths)$mtime), identical(hashes, tools::md5sum(paths)))
  }
  } else {
  h5 <- file.path(tmp, "joint.h5")
  invisible(write_matrix_10x_hdf5(convert_matrix_type(Matrix(raw, sparse = TRUE), "uint32_t"), h5))
  Sys.setFileTime(h5, as.POSIXct("2020-01-01", tz = "UTC"))
  before_time <- file.info(h5)$mtime
  before_hash <- tools::md5sum(h5)
  native <- file.path(tmp, "joint.bpcells")
  write_join_bpcells_matrix(h5, native)
  stopifnot(identical(before_time, file.info(h5)$mtime), identical(before_hash, tools::md5sum(h5)),
            isTRUE(all.equal(unname(as.matrix(open_matrix_dir(native))), unname(raw))))
  # An incompatible input is rejected, never silently patched by a reader.
  fh <- hdf5r::H5File$new(h5, mode = "r+")
  fh[["matrix/features"]]$link_delete("id")
  fh$close_all()
  before_time <- file.info(h5)$mtime
  before_hash <- tools::md5sum(h5)
  error <- tryCatch(.ensure_10x_hdf5_compat(h5), error = identity)
  stopifnot(inherits(error, "error"), grepl("join_build_matrix", conditionMessage(error)),
            identical(before_time, file.info(h5)$mtime), identical(before_hash, tools::md5sum(h5)))
  }
  message("Read-only joined HDF5 regressions passed: ", mode)
}
run_tests()
