#!/usr/bin/env Rscript

source('scripts/doublet_id.R')

tmp_dir = tempfile('doublet-id-test-')
dir.create(tmp_dir)
real_f = file.path(tmp_dir, 'real.csv.gz')
artificial_f = file.path(tmp_dir, 'artificial.csv.gz')
summary_f = file.path(tmp_dir, 'summary.csv.gz')

real = data.table(
  sample_id = c(rep('sample_a', 4L), rep('sample_b', 5L)),
  cell_id = paste0('cell_', seq_len(9L)),
  type = 'real',
  score = c(0.01, 0.05, 0.7, 0.9, 0.01, 0.02, 0.03, 0.04, 0.8),
  scdbl_class = c(
    'singlet', 'singlet', 'doublet', 'doublet',
    'singlet', 'singlet', 'singlet', 'singlet', 'doublet'))
artificial = data.table(
  sample_id = rep(c('sample_a', 'sample_b'), each = 3L),
  cell_id = paste0('artificial_', seq_len(6L)),
  type = 'doublet',
  score = c(0.4, 0.7, 0.9, 0.3, 0.6, 0.95),
  scdbl_class = NA_character_)
summary = data.table(
  sample_id = c('sample_a', 'sample_b'),
  n_real = c(4L, 5L),
  n_artificial = c(3L, 3L),
  n_called_doublet = c(2L, 1L),
  pct_called_doublet = c(50, 20),
  calling_threshold = c(0.6, 0.7),
  dbl_dbr_sd = c(1, 1),
  dbl_clusters = c(FALSE, FALSE),
  scDblFinder_version = '1.24.0')
fwrite(real, real_f)
fwrite(artificial, artificial_f)
fwrite(summary, summary_f)

x = load_doublet_id_data(real_f, artificial_f, summary_f, 'sample_id')
overview_plot = plot_doublet_overview(x)
regression_dt = fit_doublet_median_regression(x$summary)
summary_table = doublet_summary_table(x)
x_pool = x
x_pool$run_var = 'pool_id'
pool_summary_table = doublet_summary_table(x_pool)
threshold_plot = plot_doublet_thresholds(x)
histogram_plot = plot_doublet_histograms(x, x$run_levels)
density_plot = plot_doublet_densities(x, x$run_levels)
run_labels = doublet_run_labels(x, x$run_levels)
stopifnot(
  nrow(x$real) == 9L,
  nrow(x$artificial) == 6L,
  identical(x$run_levels, c('sample_a', 'sample_b')),
  identical(x$run_var, 'sample_id'),
  identical(names(summary_table)[1L], 'sample_id'),
  identical(summary_table$sample_id, c('sample_a', 'sample_b')),
  identical(names(pool_summary_table)[1L], 'pool_id'),
  all.equal(x$summary$pct_score_gt_0_5, c(50, 20)) == TRUE,
  inherits(overview_plot, 'ggplot'),
  identical(overview_plot$layers[[1L]]$aes_params$colour, 'grey25'),
  identical(overview_plot$layers[[1L]]$aes_params$linetype, NULL),
  identical(overview_plot$layers[[2L]]$aes_params$colour, 'grey45'),
  identical(
    overview_plot$labels$caption,
    'Solid line shows the median quantile regression.'),
  nrow(regression_dt) == 100L,
  !'regression' %in% names(regression_dt),
  all(is.finite(regression_dt$fitted_pct)),
  identical(
    compact_count_labels(c(500, 3000, 12500)),
    c('500', '3k', '12.5k')),
  inherits(threshold_plot, 'ggplot'),
  identical(rlang::as_label(threshold_plot$mapping$x), 'calling_threshold'),
  identical(rlang::as_label(threshold_plot$mapping$y), 'pct_called_doublet'),
  identical(threshold_plot$layers[[1L]]$aes_params$colour, 'grey45'),
  is.null(threshold_plot$scales$get_scales('x')$limits),
  inherits(histogram_plot, 'ggplot'),
  inherits(density_plot, 'ggplot'),
  inherits(density_plot$layers[[2L]]$geom, 'GeomVline'),
  identical(density_plot$layers[[2L]]$aes_params$colour, 'grey45'),
  identical(
    unname(run_labels[c('sample_a', 'sample_b')]),
    c('sample_a (50.0%)', 'sample_b (20.0%)')))

unlink(tmp_dir, recursive = TRUE)
