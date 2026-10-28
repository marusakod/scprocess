library(data.table)
source('scripts/utils.R')
source('scripts/train_xgboost.R')

# Per-label metrics retain the existing definitions for fine and coarse calls.
predictions_dt = data.table(
  cell_id = paste0('c', 1:4), UMAP1 = 1:4, UMAP2 = 4:1,
  label = factor(c('a', 'a', 'b', 'b'), levels = c('b', 'a')),
  predicted_label = c('a', 'b', 'b', 'b'),
  coarse_true = 'Neuron', coarse_predicted = 'Neuron')
metrics_dt = get_metrics_dt(predictions_dt)
stopifnot(identical(metrics_dt$Class, c('a', 'b')),
  identical(metrics_dt$N, c(2L, 2L)),
  isTRUE(all.equal(metrics_dt$Precision, c(1, 2 / 3))),
  isTRUE(all.equal(metrics_dt$Recall, c(.5, 1))),
  isTRUE(all.equal(metrics_dt$F1, c(2 / 3, .8))))
coarse_metrics_dt = get_metrics_dt(predictions_dt, coarse = TRUE)
stopifnot(coarse_metrics_dt$N == 4,
  coarse_metrics_dt$Precision == 1, coarse_metrics_dt$Recall == 1, coarse_metrics_dt$F1 == 1)

# UMAP reshaping keeps both calls per cell, the label order and the supplied colours.
label_colours = c(a = '#225588', b = '#BB4422')
umap = plot_pred_true_umap(predictions_dt, label_colours = label_colours)
stopifnot(nrow(umap$data) == 8,
  identical(levels(umap$data$label), c('b', 'a')),
  all(umap$data[, .N, by = cell_id]$N == 2L))
invisible(ggplot_build(umap))
invisible(ggplot_build(plot_precision_recall(metrics_dt, label_colours)))

# The plotted selection follows the fitted model, including exact and bounded cutoffs.
gain_dt = data.table(gene_id = paste0('g', 1:4), cumulative_gain_frac = c(.4, .8, .9, 1))
original_gain_dt = copy(gain_dt)
for (selection in list(c(.8, 2), c(.9, 4), c(.9, 2))) {
  gain_plot = plot_gain_curve(gain_dt, selection[1], selection[2])
  layers = ggplot_build(gain_plot)$data
  stopifnot(all(layers[[2]]$yintercept == selection[1]),
    all(layers[[3]]$xintercept == selection[2]),
    layers[[4]]$label == paste(selection[2], 'genes'),
    identical(gain_dt, original_gain_dt))
}
model_card_f = tempfile(fileext = '.json')
jsonlite::write_json(list(configuration = list(gain_threshold = .8), n_genes = 2),
  model_card_f, auto_unbox = TRUE)
feature_selection = read_xgboost_feature_selection(model_card_f)
stopifnot(feature_selection$gain_threshold == .8, feature_selection$n_selected_genes == 2)
unlink(model_card_f)

# Non-default columns and pseudocounts control the facet layout and axis limits.
selected_genes_dt = data.table(gene_id = paste0('g', 1:21), gain = 21:1, symbol = paste0('Gene', 1:21))
cpms_dt = CJ(gene_id = selected_genes_dt$gene_id, cluster = c('a', 'b'))
cpms_dt[, `:=`(logcpm = log(55), n_cells = 10)]
gene_plot = plot_top_xgboost_genes(selected_genes_dt, cpms_dt,
  pseudo_count = 5, ncol = 3, label_colours = label_colours)
stopifnot(gene_plot$facet$params$nrow == 7, gene_plot$facet$params$ncol == 3,
  gene_plot$scales$get_scales('y')$breaks[1] == log(5))
layers = ggplot_build(gene_plot)$data
stopifnot(isTRUE(all.equal(sort(unique(layers[[2]]$y)), log(c(5, 105)))))

cat('SUCCESS: label metrics/order, fitted gain cutoffs, facet layout and pseudocount scales passed\n')
