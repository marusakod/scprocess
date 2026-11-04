
suppressPackageStartupMessages({
  library('magrittr')
  library('data.table')
  library('assertthat')
  library('scales')
  library('ggplot2')
  library('forcats')
})


plot_umap_labels <- function(umap_dt, clust_dt, name, label_order = NULL, label_colours = NULL,
  label_clusters = FALSE) {
  # join umap and clusters
  assert_that(
    all(c('UMAP1', 'UMAP2') %in% names(umap_dt)),
    'cell_id' %in% names(umap_dt),
    'cell_id' %in% names(clust_dt),
    'cluster' %in% names(clust_dt)
  )

  # define cluster name
  plot_dt     = merge(umap_dt, clust_dt, by = 'cell_id', all.x = TRUE) %>%
    .[, .(
      UMAP1   = rescale(UMAP1, to = c(0.05, 0.95)),
      UMAP2   = rescale(UMAP2, to = c(0.05, 0.95)),
      cluster
      )] %>% .[, cluster := factor(cluster) %>% fct_infreq ]

  if (!is.null(label_order))
    plot_dt[, cluster := factor(as.character(cluster), levels = unique(c(label_order, levels(cluster))))]

  # tweak if "ambiguous"
  if ("ambiguous" %in% levels(plot_dt$cluster) )
    plot_dt     = plot_dt[, cluster := cluster %>% fct_relevel("ambiguous", after = Inf) ]

  # define colours
  cl_lvls     = levels(plot_dt$cluster)
  cl_cols     = xgboost_colours(cl_lvls, label_colours)
  if ("ambiguous" %in% cl_lvls)
    cl_cols[[ "ambiguous" ]] = "grey"

  # make plot
  plot_dt     = plot_dt[ sample(.N, .N) ]
  g = ggplot(plot_dt) +
    aes( x = UMAP1, y = UMAP2, colour = cluster ) +
    geom_point(size = 0.1) +
    scale_colour_manual( values = cl_cols, guide = guide_legend(override.aes = list(size = 3)) ) +
    scale_x_continuous( breaks = pretty_breaks(), limits = c(0, 1) ) +
    scale_y_continuous( breaks = pretty_breaks(), limits = c(0, 1) ) +
    theme_bw() +
    theme( panel.grid = element_blank(), aspect.ratio = 1,
      axis.ticks = element_blank(), axis.text = element_blank() ) +
    labs( colour = name )

  # Print each cluster's name at its median position, so that clusters sharing
  # a repeated colour can still be told apart
  if (label_clusters) {
    centres_dt = plot_dt[, .(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), by = cluster]
    g = g + geom_text(data = centres_dt, aes(label = cluster), colour = 'black', size = 2.5)
  }

  return(g)
}

calc_labels_table <- function(guesses_dt, label_order = NULL) {
  ns_dt   = guesses_dt %>%
    .[, .(cell_id, unagg = predicted_label_naive, agg = predicted_label_agg)] %>%
    melt(id = "cell_id", variable.name = "label_method", value.name = "label") %>%
    .[, .N, by = .(label_method, label)] %>%
    dcast( label ~ label_method, value.var = "N", fill = 0) %>%
    .[ order(-agg, -unagg) ] %>% setcolorder(c("label", "agg", "unagg"))
  if (!is.null(label_order)) {
    ns_dt[, label := factor(label, levels = unique(c(label_order, as.character(label))))]
    setorder(ns_dt, label)
  }

  # Keep ambiguous calls after the named labels.
  if ("ambiguous" %in% ns_dt$label)
    ns_dt     = rbind( ns_dt[ label != "ambiguous" ], ns_dt[ label == "ambiguous" ] )

  # Use descriptive headings in the report table.
  ns_dt     = ns_dt %>% set_colnames(c("predicted\nlabel", "no. cells,\naggregated", "no. cells, not\naggregated"))

  return(ns_dt)
}

installed_label_style <- function(model_name) {
  registry_f = file.path(Sys.getenv("SCPROCESS_DATA_DIR"), "xgboost", "xgboost_models.csv")
  empty = list(order = NULL, colours = NULL, coarse_order = NULL,
    coarse_colours = NULL, mapping = character())
  if (!file.exists(registry_f)) return(empty)
  registry = fread(registry_f, na.strings = c("", "NA"))
  # model_name must differ from the registry's model column: inside data.table
  # brackets, a bare `model` would mean the column itself and match every row.
  entry = registry[model == model_name]
  # Models without specification tables keep the report's existing defaults.
  has_specs = nrow(entry) == 1L && "label_specs_f" %in% names(entry) &&
    !is.na(entry$label_specs_f)
  if (!has_specs) return(empty)

  coarse_specs_f = ""
  if ("coarse_label_specs_f" %in% names(entry) && !is.na(entry$coarse_label_specs_f)) {
    coarse_specs_f = entry$coarse_label_specs_f
  }
  read_xgboost_style(NULL, character(), entry$label_specs_f, coarse_specs_f)
}


# The tables shown in each labeller's tabs: the model's labels, then the
# grouped labels when the predictions have them. Each table has cell_id,
# predicted_label_naive and predicted_label_agg; the grouped table holds only
# the grouped columns under those names.
get_label_views <- function(guesses_dt, style) {
  views = list(Labels = list(
    labels_dt = guesses_dt[, .(cell_id, predicted_label_naive, predicted_label_agg)],
    order = style$order, colours = style$colours))
  if (all(c('coarse_predicted_label_naive', 'coarse_predicted_label_agg') %in% names(guesses_dt))) {
    views[['Grouped labels']] = list(
      labels_dt = guesses_dt[, .(cell_id, predicted_label_naive = coarse_predicted_label_naive,
        predicted_label_agg = coarse_predicted_label_agg)],
      order = style$coarse_order, colours = style$coarse_colours)
  }
  return(views)
}
