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
x$score_summary = summarize_doublet_scores(x)
precision_plot = plot_doublet_beta_mean_precision(x)
set.seed(713)
beta_fit = fit_doublet_score_beta(rbeta(4000, 5, 15))
stopifnot(
  abs(beta_fit[['fitted_mean']] - 0.25) < 0.02,
  abs(beta_fit[['fitted_precision']] - 20) < 2,
  isTRUE(all.equal(unname(beta_fit[['fitted_variance']]),
    unname(beta_fit[['fitted_mean']] * (1 - beta_fit[['fitted_mean']]) /
      (beta_fit[['fitted_precision']] + 1)))),
  all(is.na(fit_doublet_score_beta(rep(0.5, 10)))),
  identical(rlang::as_label(precision_plot$mapping$y), 'fitted_precision'))
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
  isTRUE(all.equal(doublet_page_height(seq_len(4)), 3.4)),
  isTRUE(all.equal(doublet_page_height(seq_len(6)), 5.6)),
  isTRUE(all.equal(doublet_page_height(seq_len(20)), 10)))
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

# Exercise child chunks as well as the plot objects: page headings and device
# heights must survive knitting, including when the final page has fewer rows.
local({
  old_options = knitr::opts_chunk$get()
  on.exit(knitr::opts_chunk$restore(old_options))
  knitr::opts_chunk$set(fig.path = paste0(tmp_dir, '/'), dev = 'png')
  pages = split_doublet_runs(x, page_size = 1L)
  knitted = knitr::knit(text = c(
    '```{r test_pages, results="asis", echo=FALSE}',
    'render_doublet_pages(x, pages, plot_doublet_histograms, "test_hist", heading_level=5L)',
    '```'), envir = environment(), quiet = TRUE)
  stopifnot(
    grepl('##### Page 1', knitted, fixed = TRUE),
    grepl('##### Page 2', knitted, fixed = TRUE),
    length(list.files(tmp_dir, pattern = 'test_hist.*png$')) == 2L,
    !grepl('Error', knitted, fixed = TRUE))
})

# An older run's empty artificial table or missing threshold cannot supply the
# diagnostics required by current reports.
fwrite(artificial[0], artificial_f)
incomplete = tryCatch(load_doublet_id_data(real_f, artificial_f, summary_f, 'sample_id'), error = identity)
stopifnot(inherits(incomplete, 'error'), grepl('rerun QC', conditionMessage(incomplete)))
fwrite(artificial, artificial_f)
summary[, calling_threshold := NA_real_]
fwrite(summary, summary_f)
incomplete = tryCatch(load_doublet_id_data(real_f, artificial_f, summary_f, 'sample_id'), error = identity)
stopifnot(inherits(incomplete, 'error'), grepl('rerun QC', conditionMessage(incomplete)))

unlink(tmp_dir, recursive = TRUE)
