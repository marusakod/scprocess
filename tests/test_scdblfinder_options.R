#!/usr/bin/env Rscript

source('scripts/SampleQC.R')

captured_args = NULL
scDblFinder = function(sce, returnType, multiSampleMode, verbose, ...) {
  captured_args <<- list(...)
  n_real = ncol(sce)
  result = S4Vectors::DataFrame(
    type = c(rep('real', n_real), rep('doublet', 2L)),
    score = c(seq(0.1, 0.9, length.out = n_real), 0.8, 0.95),
    class = c(rep('singlet', max(0L, n_real - 1L)), 'doublet', NA, NA),
    row.names = c(colnames(sce), 'artificial_1', 'artificial_2'))
  S4Vectors::metadata(result)$scDblFinder.threshold = 0.7
  result
}

sce = SingleCellExperiment(assays = list(
  counts = Matrix::Matrix(matrix(1, nrow = 5L, ncol = 3L), sparse = TRUE)))
colnames(sce) = paste0('cell_', seq_len(ncol(sce)))
sce$sample_id = 'sample_a'
sce$detected = 5L

tmp_dir = tempfile('scdblfinder-options-')
dir.create(tmp_dir)
real_f = file.path(tmp_dir, 'real.csv.gz')
artificial_f = file.path(tmp_dir, 'artificial.csv.gz')
summary_f = file.path(tmp_dir, 'summary.csv.gz')

result = .run_scdblfinder(
  sce, 'sample_a', 'sample_id', 'none', real_f, artificial_f, summary_f,
  min_feats = 1L, min_cells = 1L, dbl_dbr_sd = 1, dbl_clusters = TRUE)
stopifnot(
  identical(captured_args$dbr.sd, 1),
  identical(captured_args$clusters, TRUE),
  nrow(result) == 3L,
  all(fread(artificial_f)$type == 'doublet'),
  fread(summary_f)$calling_threshold == 0.7)

invisible(.run_scdblfinder(
  sce, 'sample_a', 'sample_id', 'none', real_f, artificial_f, summary_f,
  min_feats = 1L, min_cells = 1L))
stopifnot(
  identical(captured_args$dbr.sd, 1),
  identical(captured_args$clusters, FALSE))

invisible(.run_scdblfinder(
  sce, 'sample_a', 'sample_id', 'none', real_f, artificial_f, summary_f,
  min_feats = 1L, min_cells = 1L, dbl_dbr_sd = NULL, dbl_clusters = NULL))
stopifnot(
  !('dbr.sd' %in% names(captured_args)),
  !('clusters' %in% names(captured_args)))

unlink(tmp_dir, recursive = TRUE)
