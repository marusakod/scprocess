suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

load_doublet_id_data <- function(real_f, artificial_f, summary_f, run_var) {
  real = fread(real_f)
  artificial = fread(artificial_f)
  summary = fread(summary_f)

  required_real = c(run_var, 'cell_id', 'type', 'score', 'scdbl_class')
  required_artificial = c(run_var, 'cell_id', 'type', 'score')
  required_summary = c(
    run_var, 'n_real', 'n_artificial', 'n_called_doublet',
    'pct_called_doublet', 'calling_threshold', 'dbl_dbr_sd',
    'dbl_clusters', 'scDblFinder_version')
  missing_real = setdiff(required_real, names(real))
  missing_artificial = setdiff(required_artificial, names(artificial))
  missing_summary = setdiff(required_summary, names(summary))
  if (length(missing_real))
    stop('Real-cell diagnostics lack: ', paste(missing_real, collapse = ', '))
  if (length(missing_artificial))
    stop('Artificial-doublet diagnostics lack: ',
      paste(missing_artificial, collapse = ', '))
  if (length(missing_summary))
    stop('Doublet summaries lack: ', paste(missing_summary, collapse = ', '))
  if (any(real$type != 'real') || any(artificial$type != 'doublet'))
    stop('Unexpected scDblFinder type values in diagnostic files')
  if (any(!is.finite(real$score)) || any(!is.finite(artificial$score)) ||
      any(real$score < 0 | real$score > 1) ||
      any(artificial$score < 0 | artificial$score > 1))
    stop('scDblFinder scores must be finite and between zero and one')
  if (anyDuplicated(summary[[run_var]]))
    stop('Doublet summary has more than one row per ', run_var)

  setnames(real, run_var, 'run_id')
  setnames(artificial, run_var, 'run_id')
  setnames(summary, run_var, 'run_id')
  high_score = real[, .(n_score_gt_0_5 = sum(score > 0.5)), by = run_id]
  summary = merge(summary, high_score, by = 'run_id', all.x = TRUE)
  summary[is.na(n_score_gt_0_5), n_score_gt_0_5 := 0L]
  summary[, pct_score_gt_0_5 := 100 * n_score_gt_0_5 / n_real]
  mean_score = real[, .(mean_score = mean(score)), by = run_id]
  summary = merge(summary, mean_score, by = 'run_id', all.x = TRUE)

  run_levels = summary[order(-pct_called_doublet, run_id), run_id]
  real[, run_id := factor(run_id, levels = run_levels)]
  artificial[, run_id := factor(run_id, levels = run_levels)]
  summary[, run_id := factor(run_id, levels = run_levels)]

  list(
    real = real,
    artificial = artificial,
    summary = summary,
    run_levels = run_levels,
    run_var = run_var)
}

doublet_summary_table <- function(x) {
  summary_table = x$summary[, .(
    run_id = as.character(run_id),
    real_cells = n_real,
    artificial_doublets = n_artificial,
    called_doublets = n_called_doublet,
    called_doublets_pct = pct_called_doublet,
    score_gt_0_5_pct = pct_score_gt_0_5,
    threshold = calling_threshold,
    mean_score)]
  setnames(summary_table, 'run_id', x$run_var)
  summary_table
}

plot_doublet_overview <- function(x) {
  regression_dt = fit_doublet_median_regression(x$summary)
  ggplot(x$summary, aes(x = n_real, y = pct_called_doublet)) +
    geom_line(
      data = regression_dt,
      aes(x = n_real, y = fitted_pct),
      inherit.aes = FALSE, colour = 'grey25', linewidth = 0.6) +
    geom_point(colour = 'grey45', alpha = 0.8, size = 2) +
    scale_x_continuous(
      breaks = pretty_breaks(), labels = compact_count_labels) +
    scale_y_continuous(
      limits = c(0, NA), breaks = pretty_breaks(), labels = label_percent(scale = 1),
      expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = 'Eligible real cells', y = 'Cells called as doublets',
      caption = 'Solid line shows the median quantile regression.') +
    theme_classic() +
    theme(plot.caption = element_text(colour = 'grey40', size = 8, hjust = 1))
}

fit_doublet_median_regression <- function(summary_dt) {
  if (nrow(summary_dt) < 2L || uniqueN(summary_dt$n_real) < 2L)
    stop('At least two distinct eligible-cell counts are required for quantile regression')
  prediction_dt = data.table(
    n_real = seq(min(summary_dt$n_real), max(summary_dt$n_real), length.out = 100L))
  fit = quantreg::rq(
    pct_called_doublet ~ n_real, tau = 0.5, data = summary_dt)
  prediction_dt[, fitted_pct := as.numeric(predict(fit, newdata = prediction_dt))]
  prediction_dt
}

plot_doublet_thresholds <- function(x) {
  ggplot(x$summary, aes(
      x = calling_threshold, y = pct_called_doublet)) +
    geom_point(colour = 'grey45', alpha = 0.8, size = 2) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(
      limits = c(0, NA), breaks = pretty_breaks(),
      labels = label_percent(scale = 1),
      expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = 'scDblFinder calling threshold',
      y = 'Cells called as doublets') +
    theme_classic()
}

compact_count_labels <- function(x) {
  labels = label_comma()(x)
  use_thousands = is.finite(x) & abs(x) >= 1000
  scaled = sprintf('%.1f', x[use_thousands] / 1000)
  scaled = sub('\\.0$', '', scaled)
  labels[use_thousands] = paste0(scaled, 'k')
  labels
}

plot_doublet_histograms <- function(x, runs, log_y = FALSE, binwidth = 0.02) {
  plot_dt = x$real[run_id %in% runs]
  threshold_dt = x$summary[run_id %in% runs]
  p = ggplot(plot_dt, aes(x = score)) +
    geom_histogram(
      binwidth = binwidth, boundary = 0, closed = 'left',
      fill = 'grey55', colour = 'white', linewidth = 0.15) +
    geom_vline(
      data = threshold_dt, aes(xintercept = calling_threshold),
      colour = '#B2182B', linetype = 'dashed', linewidth = 0.4) +
    facet_wrap(vars(run_id), ncol = 5, scales = 'free_y') +
    scale_x_continuous(limits = c(0, 1), breaks = pretty_breaks()) +
    labs(
      x = 'scDblFinder score', y = 'Real cells',
      caption = 'Dashed lines show the scDblFinder calling threshold.') +
    theme_classic() +
    theme(
      strip.background = element_rect(fill = 'white', colour = 'grey75'),
      strip.text = element_text(size = 8),
      plot.caption = element_text(colour = 'grey40', size = 8, hjust = 1))
  if (isTRUE(log_y))
    p + scale_y_log10(breaks = log_breaks(n = 4), labels = label_number())
  else
    p + scale_y_continuous(breaks = pretty_breaks(), labels = label_comma())
}

plot_doublet_densities <- function(x, runs, bandwidth = 0.01) {
  real = copy(x$real[run_id %in% runs])
  artificial = copy(x$artificial[run_id %in% runs])
  threshold_dt = x$summary[run_id %in% runs]
  real[, score_category := fifelse(
    scdbl_class == 'doublet', 'Called doublet', 'Called singlet')]
  artificial[, score_category := 'Artificial doublet']
  plot_dt = rbindlist(list(real, artificial), use.names = TRUE, fill = TRUE)
  category_levels = c('Called singlet', 'Called doublet', 'Artificial doublet')
  plot_dt[, score_category := factor(score_category, levels = category_levels)]
  counts = plot_dt[, .N, by = .(run_id, score_category)]
  plot_dt = plot_dt[counts[N >= 2L], on = .(run_id, score_category)]
  run_labels = doublet_run_labels(x, runs)

  ggplot(plot_dt, aes(x = score, colour = score_category)) +
    geom_density(linewidth = 0.55, bw = bandwidth, n = 1024) +
    geom_vline(
      data = threshold_dt, aes(xintercept = calling_threshold),
      inherit.aes = FALSE, colour = 'grey45', linetype = 'dashed',
      linewidth = 0.4) +
    facet_wrap(
      vars(run_id), ncol = 5, scales = 'free_y',
      labeller = labeller(run_id = as_labeller(run_labels))) +
    scale_colour_manual(values = c(
      'Called singlet' = '#4D4D4D',
      'Called doublet' = '#EF8A62',
      'Artificial doublet' = '#B2182B')) +
    scale_x_continuous(limits = c(0, 1), breaks = pretty_breaks()) +
    scale_y_continuous(
      trans = pseudo_log_trans(base = 10, sigma = 0.01),
      breaks = c(0, 0.1, 1, 10, 100, 1000), labels = label_number()) +
    labs(x = 'scDblFinder score', y = 'Density', colour = NULL) +
    theme_classic() +
    theme(
      strip.background = element_rect(fill = 'white', colour = 'grey75'),
      strip.text = element_text(size = 8),
      legend.position = 'bottom')
}

doublet_run_labels <- function(x, runs) {
  label_dt = x$summary[run_id %in% runs]
  setNames(
    sprintf('%s (%.1f%%)', as.character(label_dt$run_id),
      label_dt$pct_called_doublet),
    as.character(label_dt$run_id))
}

split_doublet_runs <- function(x, page_size = 20L) {
  split(x$run_levels, ceiling(seq_along(x$run_levels) / page_size))
}
