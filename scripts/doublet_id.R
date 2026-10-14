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

doublet_page_height <- function(runs) {
  1.2 + 2.2 * ceiling(length(runs) / 5)
}

render_doublet_pages <- function(x, pages, plot_fun, chunk_prefix,
                                heading_level = 4L, ...) {
  # Child chunks give each page its own device height, including a short last page.
  for (i in seq_along(pages)) {
    runs = pages[[i]]
    cat('\n', strrep('#', heading_level), ' Page ', i, '\n\n', sep = '')
    page_plot = plot_fun(x, runs, ...)
    child = c(
      sprintf(
        '```{r %s_%d, fig.width=14, fig.height=%s, echo=FALSE, cache=FALSE}',
        chunk_prefix, i, doublet_page_height(runs)),
      'print(page_plot)',
      '```')
    cat(knitr::knit_child(text = child, envir = environment(), quiet = TRUE),
      sep = '\n')
    cat('\n\n')
  }
}

fit_doublet_score_beta <- function(scores) {
  empty_fit = c(fitted_mean = NA_real_, fitted_variance = NA_real_)
  if (length(scores) < 3L || var(scores) == 0) return(empty_fit)
  n = length(scores)
  adjusted = (scores * (n - 1) + 0.5) / n
  mu = mean(adjusted)
  concentration = mu * (1 - mu) / var(adjusted) - 1
  if (!is.finite(concentration) || concentration <= 0) concentration = 2
  start = log(pmax(c(mu * concentration, (1 - mu) * concentration), 0.05))
  sum_log = sum(log(adjusted))
  sum_log_one_minus = sum(log1p(-adjusted))
  objective = function(log_parameters) {
    a = exp(log_parameters[1L])
    b = exp(log_parameters[2L])
    n * lbeta(a, b) - (a - 1) * sum_log - (b - 1) * sum_log_one_minus
  }
  fit = tryCatch(optim(start, objective, method = 'BFGS'), error = function(e) NULL)
  if (is.null(fit) || fit$convergence != 0L) return(empty_fit)
  parameters = exp(fit$par)
  total = sum(parameters)
  if (any(!is.finite(parameters)) || any(parameters <= 0)) return(empty_fit)
  c(fitted_mean = parameters[1L] / total,
    fitted_variance = prod(parameters) / (total^2 * (total + 1)))
}

summarize_doublet_scores <- function(x) {
  summary = x$real[, {
    entropy = numeric(.N)
    interior = score > 0 & score < 1
    entropy[interior] = -(score[interior] * log2(score[interior]) +
      (1 - score[interior]) * log2(1 - score[interior]))
    c(list(mean_binary_entropy = mean(entropy), n_cells = .N,
      n_high_scores = sum(score > 0.5)), as.list(fit_doublet_score_beta(score)))
  }, by = run_id]
  summary[, logit_high_score_prop := qlogis((n_high_scores + 0.5) / (n_cells + 1))]
  spread = mad(summary$logit_high_score_prop)
  cutoff = if (is.finite(spread) && spread > 0)
    median(summary$logit_high_score_prop) + 2 * spread else Inf
  summary[, high_score_outlier := logit_high_score_prop > cutoff]
  summary
}

plot_doublet_score_entropy <- function(x) {
  summary = copy(x$score_summary)
  summary[, run_id := factor(as.character(run_id),
    levels = as.character(summary[order(-mean_binary_entropy, run_id), run_id]))]
  ggplot(summary, aes(x = run_id, y = mean_binary_entropy)) +
    geom_hline(yintercept = median(summary$mean_binary_entropy),
      linewidth = 0.45, linetype = 'dashed', colour = 'grey55') +
    geom_point(size = 2.2, alpha = 0.85, colour = 'grey25') +
    scale_y_continuous(limits = c(0, NA), breaks = pretty_breaks(),
      expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = 'Mean binary entropy (bits)',
      caption = 'Samples are ordered by entropy; the dashed line is the across-sample median.') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.caption = element_text(colour = 'grey40', size = 8, hjust = 1))
}

plot_doublet_beta_mean_variance <- function(x) {
  fits = x$score_summary[is.finite(fitted_mean) & is.finite(fitted_variance)]
  outliers = fits[high_score_outlier == TRUE]
  omitted = nrow(x$score_summary) - nrow(fits)
  caption = 'Maximum-likelihood beta fits; red labels flag unusually high fractions of scores above 0.5.'
  if (omitted > 0) caption = paste(caption,
    sprintf('%d sample(s) omitted because a beta fit was not estimable.', omitted))
  ggplot(fits, aes(x = fitted_mean, y = fitted_variance)) +
    geom_point(size = 2.2, alpha = 0.8, colour = 'grey30') +
    geom_point(data = outliers, size = 2.4, colour = '#D7191C') +
    ggrepel::geom_label_repel(data = outliers, aes(label = run_id),
      seed = 1, size = 3, min.segment.length = 0,
      box.padding = 0.35, point.padding = 0.25,
      fill = alpha('white', 0.9), colour = 'grey20',
      label.size = 0.15, max.overlaps = Inf) +
    scale_x_continuous(breaks = pretty_breaks(), labels = label_number(accuracy = 0.01)) +
    scale_y_continuous(breaks = pretty_breaks(), labels = label_number(accuracy = 0.01)) +
    labs(x = 'Fitted beta mean', y = 'Fitted beta variance', caption = caption) +
    theme_classic() +
    theme(aspect.ratio = 1,
      plot.caption = element_text(colour = 'grey40', size = 8, hjust = 1))
}
