suppressPackageStartupMessages({
  library('anndataR')
})


read_xgboost_feature_selection <- function(model_card_f) {

  # Read the settings and final gene count from the fitted model, not current config.
  model_card = jsonlite::fromJSON(model_card_f)
  feature_selection = list(
    gain_threshold   = model_card$configuration$gain_threshold,
    n_selected_genes = model_card$n_genes)
  stopifnot(length(feature_selection$gain_threshold) == 1,
    length(feature_selection$n_selected_genes) == 1)
  return(feature_selection)
}


plot_train_test_split <- function(xgboost_dt) {

  # One row per cell; show each split using the same density limits and colour scale.
  density_limits = c(0.0001, 0.01)
  plot = ggplot(xgboost_dt) +
    aes(x = UMAP1, y = UMAP2) +
    facet_wrap(~split) +
    geom_bin2d(bins = 50,
      aes(fill = pmin(pmax(after_stat(density), density_limits[1]), density_limits[2]))) +
    scale_fill_distiller(palette = 'RdBu', trans = 'log10', limits = density_limits) +
    theme_classic() +
    theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank()) +
    labs(fill = 'density')
  return(plot)
}


plot_pred_true_umap <- function(predictions_dt, coarse = FALSE, label_colours = NULL) {

  # The caller supplies the desired cells; each input row contains both label calls.
  if (coarse) {
    reference_label_col = 'coarse_true'
    predicted_label_col = 'coarse_predicted'
  } else {
    reference_label_col = 'label'
    predicted_label_col = 'predicted_label'
  }
  label_levels = levels(predictions_dt[[reference_label_col]])

  # Reshape to two rows per cell: one reference label and one prediction.
  label_columns = c('cell_id', reference_label_col, predicted_label_col, 'UMAP1', 'UMAP2')
  plot_dt = predictions_dt[, label_columns, with = FALSE]
  plot_dt = melt(plot_dt,
    id.vars       = c('cell_id', 'UMAP1', 'UMAP2'),
    variable.name = 'label_type',
    value.name    = 'label')
  plot_dt[, label_type := fcase(label_type == reference_label_col,
    'True labels', default = 'Predicted labels')]
  plot_dt[, label_type := factor(label_type, levels = c('True labels', 'Predicted labels'))]
  if (!is.null(label_levels)) {
    plot_dt[, label := factor(label, levels = label_levels)]
  }

  # Preserve the supplied label ordering and named colour mapping in both panels.
  labels = levels(plot_dt$label)
  if (is.null(labels)) {
    labels = sort(unique(plot_dt$label))
  }
  label_cols = xgboost_colours(labels, label_colours)
  plot = ggplot(plot_dt) +
    aes(x = UMAP1, y = UMAP2, colour = label) +
    geom_point(size = 0.3) +
    facet_wrap(~label_type) +
    scale_colour_manual(values = label_cols) +
    theme_classic() +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
  return(plot)
}


get_metrics_dt <- function(predictions_dt, coarse = FALSE) {

  # Input is one row per cell; output is one row per reference label.
  if (coarse) {
    reference_label_col = 'coarse_true'
    predicted_label_col = 'coarse_predicted'
  } else {
    reference_label_col = 'label'
    predicted_label_col = 'predicted_label'
  }
  reference_labels = as.character(predictions_dt[[reference_label_col]])
  predicted_labels = as.character(predictions_dt[[predicted_label_col]])
  labels           = unique(reference_labels)

  # Calculate metrics for labels with reference cells in this subset.
  metrics_dt = rbindlist(lapply(labels, function(label) {
    tp        = sum(predicted_labels == label & reference_labels == label)
    fp        = sum(predicted_labels == label & reference_labels != label)
    fn        = sum(predicted_labels != label & reference_labels == label)
    precision = tp / max(tp + fp, 1)
    recall    = tp / max(tp + fn, 1)
    f1        = 2 * precision * recall / max(precision + recall, 1e-10)
    data.table(Class = label, N = sum(reference_labels == label),
      Precision = precision, Recall = recall, F1 = f1)
  }))
  return(metrics_dt)
}


plot_precision_recall <- function(metrics_dt, label_colours = NULL) {

  # Each point represents one reference label; point size shows its cell count.
  label_cols = xgboost_colours(as.character(metrics_dt$Class), label_colours)
  plot = ggplot(metrics_dt) +
    aes(x = Recall, y = Precision, colour = Class, size = N) +
    geom_hline(yintercept = 1, colour = 'grey70', linetype = 'dashed') +
    geom_vline(xintercept = 1, colour = 'grey70', linetype = 'dashed') +
    geom_point() +
    geom_label_repel(aes(label = Class), size = 3, max.overlaps = Inf, show.legend = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_colour_manual(values = label_cols, guide = 'none') +
    scale_size(range = c(1, 8)) +
    expand_limits(x = c(0, 1), y = c(0, 1)) +
    theme_classic() +
    labs(size = 'N cells')
  return(plot)
}


plot_gain_curve <- function(gain_dt, gain_threshold, n_selected_genes) {

  # One row per gene, already ordered by descending gain in the training output.
  # The final count includes the training minimum/maximum bounds on gene selection.
  stopifnot(gain_threshold > 0, gain_threshold <= 1,
    n_selected_genes >= 1, n_selected_genes <= nrow(gain_dt))
  plot_dt = copy(gain_dt)
  plot_dt[, gene_rank := .I]
  selection_colour = '#7A5195'

  plot = ggplot(plot_dt) +
    aes(x = gene_rank, y = cumulative_gain_frac) +
    geom_line(colour = 'grey30') +
    geom_hline(yintercept = gain_threshold, linetype = 'dashed', colour = 'grey60') +
    geom_vline(xintercept = n_selected_genes, linetype = 'dashed', colour = selection_colour) +
    annotate('text', x = n_selected_genes, y = 0.5,
      label = paste(n_selected_genes, 'genes'), hjust = -0.1, colour = selection_colour) +
    theme_classic() +
    labs(x = 'Gene rank', y = 'Cumulative gain fraction')
  return(plot)
}


plot_top_xgboost_genes <- function(selected_genes_dt, cpms_dt, label_order = NULL,
  pseudo_count = 10, ncol = 2, nrow = NULL, label_colours = NULL) {

  # selected_genes_dt has one row per gene; use the requested number of columns.
  if (is.null(nrow)) {
    nrow = max(5, ceiling(nrow(selected_genes_dt) / ncol))
  }

  # plot_dt has one row per gene, batch and label, with expression and cell counts.
  plot_dt = merge(cpms_dt, selected_genes_dt, by = 'gene_id')
  plot_dt[, symbol_lab := sprintf('%s (gain = %s)', symbol, signif(gain, 3))]

  # Collapse repeated batch/label rows to rank each gene facet by its gain.
  gene_gains_dt = plot_dt[, .(gain = gain[1]), by = symbol_lab][order(-gain)]
  plot_dt[, symbol_lab := factor(symbol_lab, levels = gene_gains_dt$symbol_lab)]

  # Keep the requested label order, or use alphabetical order when none is supplied.
  if (!is.null(label_order)) {
    labels = intersect(label_order, unique(plot_dt$cluster))
  } else {
    labels = sort(unique(plot_dt$cluster))
  }
  plot_dt[, cluster := factor(cluster, levels = labels)]
  label_cols = xgboost_colours(labels, label_colours)

  # logcpm uses natural logs; show CPM labels using the same supplied pseudocount.
  cpm_breaks = c(0, 10, 20, 50, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4)
  log_breaks = log(cpm_breaks + pseudo_count)
  cpm_labels = c('0', '10', '20', '50', '100', '200', '500',
    '1k', '2k', '5k', '10k', '20k', '50k')
  y_limits = log(c(0, 100) + pseudo_count)

  plot = ggplot(plot_dt) +
    aes(x = cluster, y = logcpm, fill = cluster, size = n_cells) +
    geom_quasirandom(colour = 'black', shape = 21) +
    scale_y_continuous(breaks = log_breaks, labels = cpm_labels) +
    expand_limits(y = y_limits) +
    scale_fill_manual(values = label_cols, guide = 'none') +
    scale_size(range = c(0, 6)) +
    expand_limits(size = 0) +
    facet_wrap(~symbol_lab, scales = 'free_y', ncol = ncol, nrow = nrow) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = 'label', y = 'counts per million', size = '# cells')
  return(plot)
}
