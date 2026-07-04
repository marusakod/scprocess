suppressPackageStartupMessages({
  library("anndataR")
})

plot_train_test_split <- function(xgboost_dt){

  ggplot(xgboost_dt) + aes(x = UMAP1, y = UMAP2) +
    facet_wrap(~split) +
    geom_bin2d( bins = 50, aes(fill = after_stat( density ) %>% pmin(0.01) %>% pmax(0.0001)) ) +
    scale_fill_distiller( palette = "RdBu", trans = "log10", limits = c(0.0001, 0.01) ) +
    theme_classic() +
    theme( aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank() ) +
    labs( fill = sprintf("density"))
}

plot_pred_true_umap <- function(val_dt, coarse = FALSE){

  if(coarse){
    true_lbl_var = 'coarse_true'
    pred_lbl_var = 'coarse_predicted'
  }else{
    true_lbl_var = 'label'
    pred_lbl_var = 'predicted_label'
  }

  # keep only validation data for plotting
  lbl_levels = levels(val_dt[[true_lbl_var]])
  val_dt_long = val_dt %>%
    .[ , c("cell_id", true_lbl_var, pred_lbl_var, "UMAP1", "UMAP2"), with = FALSE] %>%
    melt(id.vars = c('cell_id', 'UMAP1', 'UMAP2'), variable.name = 'label_type', value.name = 'label') %>%
    .[, label_type := fcase(label_type == true_lbl_var, 'True labels', default = 'Predicted labels')] %>%
    .[, label_type := factor(label_type, levels = c('True labels', 'Predicted labels'))]
  if (!is.null(lbl_levels)) {
    val_dt_long[, label := factor(label, levels = lbl_levels)]
  }

  lbls       = levels(val_dt_long$label)
  if (is.null(lbls)) lbls = val_dt_long[, label] %>% unique() %>% sort()
  lbl_cols   = nice_cols[ seq_along(lbls) ] %>% setNames(lbls)

  p = ggplot(val_dt_long) +
    aes(x = UMAP1, y = UMAP2, colour = label) +
    geom_point(size = 0.3) +
    facet_wrap(~label_type) +
    theme_classic() +
    scale_color_manual( values = lbl_cols) +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 3)))

}

get_metrics_dt <- function(val_dt){
  true_labels = as.character(val_dt$label)
  pred_labels = as.character(val_dt$predicted_label)
  classes     = unique(true_labels)
  metrics_dt  = rbindlist(lapply(classes, function(cls) {
    tp   = sum(pred_labels == cls & true_labels == cls)
    fp   = sum(pred_labels == cls & true_labels != cls)
    fn   = sum(pred_labels != cls & true_labels == cls)
    prec = tp / max(tp + fp, 1)
    rec  = tp / max(tp + fn, 1)
    f1   = 2 * prec * rec / max(prec + rec, 1e-10)
    data.table(Class = cls, N = sum(true_labels == cls),
      Precision = prec, Recall = rec, F1 = f1)
  }))

  return(metrics_dt)
}


plot_precision_recall <- function(metrics_dt) {
  lbl_cols = nice_cols[ seq_len(nrow(metrics_dt)) ] %>% setNames(metrics_dt$Class)
  ggplot(metrics_dt) +
    aes(x = Recall, y = Precision, colour = Class, size = N) +
    geom_hline(yintercept = 1, colour = "grey70", linetype = "dashed") +
    geom_vline(xintercept = 1, colour = "grey70", linetype = "dashed") +
    geom_point() +
    geom_label_repel(aes(label = Class), size = 3, max.overlaps = Inf, show.legend = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = lbl_cols, guide = "none") +
    scale_size(range = c(1, 8)) +
    expand_limits(x = c(0, 1), y = c(0, 1)) +
    theme_classic() +
    labs(size = "N cells")
}


plot_gain_curve <- function(gain_dt){
  gain_dt = gain_dt %>%
    .[, gene_rank := .I]

  p = ggplot(gain_dt) + aes(x = gene_rank, y = cumulative_gain_frac) +
  geom_line() +
  geom_hline(yintercept = 0.9, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = sum(gain_dt$cumulative_gain_frac <= 0.9) + 1,
    linetype = "dashed", colour = "blue") +
  theme_classic() +
  labs(x = "Gene rank", y = "Cumulative gain fraction") +
  annotate("text", x = sum(gain_dt$cumulative_gain_frac <= 0.9) + 1,
    y = 0.5, label = paste0(sum(gain_dt$cumulative_gain_frac <= 0.9) + 1, " genes"),
    hjust = -0.1, colour = "blue"
    )
}


plot_top_xgboost_genes <- function(sel_dt, cpms_dt, label_order = NULL,
  pseudo_count = 10, ncol = 2, nrow = NULL) {

  # infer number of rows
  if (is.null(nrow)) {
    nrow    = max(5, ceiling( nrow(sel_dt) / 2 ))
  }
  # get data; preserve gain order for facets
  plot_dt   = cpms_dt %>% merge( sel_dt, by = 'gene_id' ) %>%
   .[, symbol_lab := sprintf("%s (gain = %s)", symbol, signif(gain, 3)) ]
  lab_order = plot_dt[, .(gain = gain[1]), by = symbol_lab][order(-gain)]$symbol_lab
  plot_dt[, symbol_lab := factor(symbol_lab, levels = lab_order)]

  # order x-axis labels by frequency
  if (!is.null(label_order)) {
    cl_ls   = intersect(label_order, unique(plot_dt$cluster))
  } else {
    cl_ls   = plot_dt$cluster %>% unique %>% sort()
  }
  plot_dt[, cluster := factor(cluster, levels = cl_ls)]
  cl_cols   = nice_cols[ seq_along(cl_ls) ] %>% setNames(cl_ls)

  # plot!
  log_brks  = c(0, 10, 20, 50, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4) %>%
    `+`(pseudo_count) %>% log
  log_labs  = c('0', '10', '20', '50', '100', '200', '500',
    '1k', '2k', '5k', '10k', '20k', '50k')
  g = ggplot(plot_dt) +
    aes( x = cluster, y = logcpm, fill = cluster, size = n_cells ) +
    geom_quasirandom( colour = 'black', shape = 21 ) +
    scale_y_continuous( breaks = log_brks, labels = log_labs ) +
    expand_limits( y = (c(0, 100) + 10) %>% log ) +
    scale_fill_manual( values = cl_cols, guide = "none" ) +
    scale_size( range = c(0, 6) ) + expand_limits( size = 0 ) +
    facet_wrap( ~ symbol_lab, scales = 'free_y', ncol = ncol, nrow = nrow ) +
    theme_classic() +
    theme(
      axis.text.x       = element_text( angle = 90, hjust = 1, vjust = 0.5 )
    ) +
    labs(
      x = "cluster", y = "counts per million", size = "# cells"
    )

  return(g)
}
